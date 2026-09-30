# -*- coding: utf-8 -*-
"""公共数值模块（项目编写；使用 NumPy、pandas 及 SciPy 线性代数）。

内容：烘房环境分段线性插值、一维柱坐标节点式有限体积网格（支持非均匀
网格与表面加密）、Thomas 追赶法、θ 格式时间推进单步（θ=0.5 为
Crank-Nicolson，θ=1 为后向 Euler；变系数时用预估-校正迭代）。

离散约定：节点 0=r_0<r_1<...<r_N=R，控制体以相邻节点中点为界（两端半
控制体），体积与面半径按单位长度并略去公因子 2π：
    V_0=rf_0^2/2,  V_i=(rf_i^2-rf_{i-1}^2)/2,  V_N=(R^2-rf_{N-1}^2)/2
半离散方程
    cap_i V_i du_i/dt = G_{i-1}(u_{i-1}-u_i)+G_i(u_{i+1}-u_i)+[i=N: B(amb-u_N)]+V_i S_i
其中面导通 G_i = rf_i*Γ_{i+1/2}/(r_{i+1}-r_i)，边界导通 B=R*β（β 为 h 或 h_m），
S 为可选体源（供制造解检验用）。内部通量两两抵消，格式天然离散守恒；
r=0 半控制体面积为零，无奇点问题。
"""
import numpy as np
import pandas as pd

ENV_XLSX = '附件/附件1.xlsx'
T_CONST, C_CONST = 50.0, 0.05   # 恒温干燥阶段设计值（附件1尾段均值约 49.99 / 0.0500）


def load_env(path=ENV_XLSX):
    """读附件1，返回环境温度 Ta(t)、水分浓度 Ca(t)。数据范围内分段线性插值，
    超出数据末端取恒温段设计值。"""
    df = pd.read_excel(path)
    t = df.iloc[:, 0].to_numpy(float)
    T = df.iloc[:, 1].to_numpy(float)
    C = df.iloc[:, 2].to_numpy(float)
    t_max = t[-1]

    def Ta(tt):
        tt = np.asarray(tt, float)
        return np.where(tt >= t_max, T_CONST, np.interp(tt, t, T))

    def Ca(tt):
        tt = np.asarray(tt, float)
        return np.where(tt >= t_max, C_CONST, np.interp(tt, t, C))

    return Ta, Ca


def make_graded_mesh(R, h_int, h_min, q, band, out_step):
    """生成含表面加密带的节点数组。

    [0, R-band] 均匀步长 h_int（要求整除，且 out_step 的整数倍都落在节点上）；
    (R-band, R] 步长按几何比 q 从内向外递减至约 h_min（表面边界层加密）。
    """
    ni = (R - band) / h_int
    assert abs(ni - round(ni)) < 1e-9, 'h_int 须整除 R-band'
    inner = np.linspace(0.0, R - band, int(round(ni)) + 1)
    m = int(np.ceil(np.log1p(band * (q - 1) / h_min) / np.log(q)))
    spac = h_min * q ** np.arange(m)[::-1]          # 靠内最大、靠表面最小
    spac *= band / spac.sum()                        # 精确铺满加密带
    nodes = np.concatenate([inner, R - band + np.cumsum(spac)])
    nodes[-1] = R
    return nodes


class CylGrid:
    """一维柱坐标节点式网格（可非均匀）及控制体几何量。"""

    def __init__(self, nodes):
        self.r = np.asarray(nodes, float)
        self.N = self.r.size - 1
        self.R = float(self.r[-1])
        self.dr = np.diff(self.r)                    # 节点间距，长度 N
        self.rf = 0.5 * (self.r[:-1] + self.r[1:])   # 面半径，长度 N
        V = np.empty(self.N + 1)
        V[0] = self.rf[0] ** 2 / 2
        V[1:-1] = (self.rf[1:] ** 2 - self.rf[:-1] ** 2) / 2
        V[-1] = (self.R ** 2 - self.rf[-1] ** 2) / 2
        self.V = V

    @classmethod
    def uniform(cls, R, N):
        return cls(np.linspace(0.0, R, N + 1))

    def out_index(self, r_targets):
        """输出半径（须恰为节点）的下标数组。"""
        idx = np.searchsorted(self.r, np.asarray(r_targets, float) - 1e-12)
        assert np.allclose(self.r[idx], r_targets, atol=1e-12), '输出点必须落在节点上'
        return idx


def thomas(sub, diag, sup, rhs):
    """追赶法解三对角方程组；sub[0]、sup[-1] 不参与。"""
    n = diag.size
    c = np.empty(n)
    d = np.empty(n)
    c[0] = sup[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - sub[i] * c[i - 1]
        c[i] = sup[i] / m
        d[i] = (rhs[i] - sub[i] * d[i - 1]) / m
    x = np.empty(n)
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


# 生产计算用 LAPACK 带状求解器加速（算法同追赶法；verify.py 中与自写
# thomas 交叉校验到机器精度）。scipy 不可用时自动回退。
try:
    from scipy.linalg import solve_banded as _solve_banded

    def tridiag_solve(sub, diag, sup, rhs):
        ab = np.zeros((3, diag.size))
        ab[0, 1:] = sup[:-1]
        ab[1] = diag
        ab[2, :-1] = sub[1:]
        return _solve_banded((1, 1), ab, rhs)
except ImportError:                                  # pragma: no cover
    tridiag_solve = thomas


def flux_semi(u, G, B, amb, V=None, S=None):
    """半离散右端 L(u)（含 Robin 边界项与可选体源）。"""
    L = np.zeros_like(u)
    q = G * (u[1:] - u[:-1])
    L[:-1] += q
    L[1:] -= q
    L[-1] += B * (amb - u[-1])
    if S is not None:
        L += V * S
    return L


def theta_step(u, dt, grid, gamma_face, cap, beta, amb0, amb1,
               theta=0.5, S0=None, S1=None):
    """θ 格式单步：M(u^{n+1}-u^n)/dt = θ L(u^{n+1};amb1,S1) + (1-θ) L(u^n;amb0,S0)。
    gamma_face: 面上的扩导系数；cap: 节点容量系数。返回 u^{n+1}。"""
    G = grid.rf * gamma_face / grid.dr
    B = grid.R * beta
    M = cap * grid.V / dt
    rhs = M * u + (1 - theta) * flux_semi(u, G, B, amb0, grid.V, S0)
    rhs[-1] += theta * B * amb1
    if S1 is not None:
        rhs += theta * grid.V * S1
    sub = np.zeros_like(u)
    sup = np.zeros_like(u)
    diag = M.copy()
    sub[1:] = -theta * G
    sup[:-1] = -theta * G
    diag[:-1] += theta * G
    diag[1:] += theta * G
    diag[-1] += theta * B
    return tridiag_solve(sub, diag, sup, rhs)


def theta_step_picard(u, dt, grid, gamma_nodes_of, cap, beta, amb0, amb1,
                      theta=0.5, S0=None, S1=None, tol=1e-6, itmax=20):
    """变系数 θ 格式的预估-校正（Picard）迭代。
    第 1 次以 u^n 计算物性作预估，随后以中点状态 (u^n+u^*)/2 的物性校正
    （θ=0.5 时保持二阶精度）。收敛判据：相邻两次迭代最大改变量 < tol。
    返回 (u^{n+1}, 迭代次数)。"""
    u_new = u.copy()
    for it in range(1, itmax + 1):
        g_node = gamma_nodes_of(0.5 * (u + u_new))
        g_face = 0.5 * (g_node[:-1] + g_node[1:])   # 面系数取相邻节点算术平均
        u_next = theta_step(u, dt, grid, g_face, cap, beta, amb0, amb1,
                            theta, S0, S1)
        delta = float(np.max(np.abs(u_next - u_new)))
        u_new = u_next
        if delta < tol:
            break
    else:
        raise RuntimeError(f'Picard iteration failed: delta={delta:.3e}, limit={itmax}')
    return u_new, it
