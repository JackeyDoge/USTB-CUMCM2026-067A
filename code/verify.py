# -*- coding: utf-8 -*-
"""问题1模型检验程序。

三项独立检验：
(1) 温度场半解析对照：常物性热方程用分离变量（Bessel 级数）+ Duhamel
    叠加原理处理时变环境温度（附件1分段线性），与数值解逐点对比；
(2) 传质求解器制造解（MMS）检验：为非线性扩散方程构造已知精确解，
    反推体源与等效环境函数，验证求解器二阶收敛；
(3) 自写 Thomas 追赶法与生产用 LAPACK 带状求解器交叉校验。

运行方式（项目根目录下）：python code/verify.py
结果追加写入 结果/p1_检验数据.md。
"""
import os
import numpy as np
import pandas as pd
from scipy import special, optimize

from common import CylGrid, thomas, tridiag_solve, theta_step_picard
from p1 import R0, RHO, CP, KH, H, HM, T0, C0, D_of_C

ALPHA = KH / (RHO * CP)
BI = H * R0 / KH


def eigenvalues(Bi, M):
    """特征方程 λ J1(λ) = Bi J0(λ) 的前 M 个正根（用 J0 零点分隔区间）。"""
    z = special.jn_zeros(0, M + 1)
    f = lambda lam: lam * special.j1(lam) - Bi * special.j0(lam)
    lams = []
    a = 1e-10
    for b in z[:M]:
        lams.append(optimize.brentq(f, a, b))
        a = b
    return np.array(lams)


def temp_semianalytic(r_eval, t_eval, M=200):
    """Bessel 级数 + Duhamel 叠加的温度半解析解（环境温度分段线性）。"""
    df = pd.read_excel('附件/附件1.xlsx')
    tt = df.iloc[:, 0].to_numpy(float)
    TT = df.iloc[:, 1].to_numpy(float)
    lam = eigenvalues(BI, M)
    mu = lam ** 2 * ALPHA / R0 ** 2
    Cn = 2 / lam * special.j1(lam) / (special.j0(lam) ** 2 + special.j1(lam) ** 2)
    J0r = special.j0(np.outer(np.asarray(r_eval) / R0, lam))     # (nr, M)
    out = np.empty((len(t_eval), len(r_eval)))
    for i, t in enumerate(t_eval):
        acc = np.full(len(r_eval), T0)
        for j in range(np.searchsorted(tt, t, side='left')):
            s0, s1 = tt[j], min(tt[j + 1], t)
            if s1 <= s0:
                continue
            slope = (TT[j + 1] - TT[j]) / (tt[j + 1] - tt[j])
            term = Cn * (np.exp(-mu * (t - s1)) - np.exp(-mu * (t - s0))) / mu
            acc = acc + slope * ((s1 - s0) - J0r @ term)
        out[i] = acc
    return out


def check_temperature():
    """生产数值解（p1_field.npz）与半解析解对比。"""
    dat = np.load(os.path.join('结果', 'p1_field.npz'))
    t_eval = np.arange(10, 1801, 10)
    r_eval = dat['r_cm'] / 100
    Tsa = temp_semianalytic(r_eval, t_eval)
    Tnum = dat['T'][t_eval - 1]
    err = np.abs(Tnum - Tsa)
    ij = np.unravel_index(np.argmax(err), err.shape)
    return dict(max=float(err.max()), mean=float(err.mean()),
                at=(int(t_eval[ij[0]]), float(r_eval[ij[1]] * 100)))


def mms_mass():
    """传质求解器制造解检验：C_m = 2.55 - 0.9 s(t) (r/R)^2, s=1-exp(-t/600)。"""
    R, Tend = R0, 1800.0

    def s_(t):
        return 1 - np.exp(-t / 600)

    def Cm(r, t):
        return 2.55 - 0.9 * s_(t) * (r / R) ** 2

    def S(r, t):
        sp = np.exp(-t / 600) / 600
        C = Cm(r, t)
        g = -1.8 * s_(t) * r / R ** 2                 # ∂C/∂r
        D = 7e-9 * np.exp(-0.89 / C)
        Dp = D * 0.89 / C ** 2                        # dD/dC
        return -0.9 * sp * (r / R) ** 2 + 1.8 * s_(t) / R ** 2 * (2 * D + r * Dp * g)

    def amb(t):
        Cs = Cm(R, t)
        return Cs + 7e-9 * np.exp(-0.89 / Cs) * (-1.8 * s_(t) / R) / HM

    errs = []
    for N, dt in ((40, 2.0), (80, 1.0), (160, 0.5)):
        grid = CylGrid.uniform(R, N)
        C = Cm(grid.r, 0.0)
        for n in range(int(round(Tend / dt))):
            th = 1.0 if n < 2 else 0.5
            t0n, t1n = n * dt, (n + 1) * dt
            C, _ = theta_step_picard(C, dt, grid, D_of_C, np.ones(grid.N + 1),
                                     HM, amb(t0n), amb(t1n), th,
                                     S(grid.r, t0n), S(grid.r, t1n))
        errs.append(float(np.max(np.abs(C - Cm(grid.r, Tend)))))
    p = [float(np.log2(a / b)) for a, b in zip(errs[:-1], errs[1:])]
    return dict(errs=errs, orders=p)


def check_solvers(n=200, trials=20, seed=1):
    """thomas 与 tridiag_solve 在随机对角占优三对角系统上的最大偏差。"""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(trials):
        sub = np.r_[0.0, -rng.random(n - 1)]
        sup = np.r_[-rng.random(n - 1), 0.0]
        diag = 2.0 + rng.random(n) - np.r_[sub[1:], 0] * 0 + np.abs(sub) + np.abs(sup)
        rhs = rng.standard_normal(n)
        worst = max(worst, float(np.max(np.abs(
            thomas(sub, diag, sup, rhs) - tridiag_solve(sub, diag, sup, rhs)))))
    return worst


def main():
    # 物理趋势断言：D 随 C 减小而减小（method.md 符号陷阱防线）
    assert D_of_C(np.array([0.15]))[0] < D_of_C(np.array([2.55]))[0]

    te = check_temperature()
    mm = mms_mass()
    sv = check_solvers()

    lines = [
        '',
        '## 问题1 独立检验（verify.py）',
        '',
        f'温度半解析对照（Bessel 级数 200 项 + Duhamel 叠加，t=10:10:1800 s × 21 个半径）：'
        f'最大偏差 {te["max"]:.2e} ℃（出现于 t={te["at"][0]} s, r={te["at"][1]:.1f} cm），'
        f'平均偏差 {te["mean"]:.2e} ℃。',
        '',
        f'传质求解器制造解检验（均匀网格 N=40/80/160，Δt=2/1/0.5 s）：'
        f'终端最大误差 ' + '、'.join(f'{e:.3e}' for e in mm['errs']) +
        f'；观测收敛阶 ' + '、'.join(f'{p:.2f}' for p in mm['orders']) + '（理论二阶）。',
        '',
        f'Thomas 追赶法 vs LAPACK 带状求解器（200 阶随机对角占优系统×20 次）：'
        f'最大偏差 {sv:.2e}。',
        '',
        'D(C) 单调性断言通过：D(0.15) < D(2.55)。',
    ]
    with open(os.path.join('结果', 'p1_检验数据.md'), 'a', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
