# -*- coding: utf-8 -*-
"""问题4：考虑径向收缩的动边界热质耦合模型。

模型（主模型 B，物质坐标观点）：药材失水收缩视为均匀（同源）收缩，
固体骨架上的物质点保持相对位置 ξ=r/R(t) 不变，水分依附骨架输运、
扩散相对骨架进行。对水分守恒方程取物质导数并作 Landau 变换后，
方程在 ξ∈[0,1] 固定域上与定域情形同构（无伪对流项）：
    R^2(t)·∂C/∂t|_ξ = (1/ξ)·∂/∂ξ( ξ·D(C,T)·∂C/∂ξ )
    R^2(t)·ρ(C)c_p(C)·∂T/∂t|_ξ = (1/ξ)·∂/∂ξ( ξ·k(C)·∂T/∂ξ )
表面（ξ=1）：ξD·∂C/∂ξ|_1 = R(t)·h_m·(C_a-C_s)，温度同理；中心对称。
R(t) 取附件2分段线性插值，物性取附录4。

对比模型 A：对固定空间扩散方程作坐标变换，保留 (ξR'/R)∂/∂ξ。
它与随骨架运动的模型B采用不同运动学假设，仅作结构对照。
这里的 rho(C) 用于有效热容量，并未同时求解局部干骨架密度守恒。

数值：ξ 域表面加密有限体积 + Crank-Nicolson（首两步后向 Euler）+
耦合预估-校正迭代；模型 A 的伪对流项按特征方向迎风、隐式并入三对角。

运行方式（项目根目录下）：python code/p4.py
输出：
  结果/result4.xlsx        每 60 s × 0.1 cm（r>R(t) 留空）+ 末列药材表面
  结果/p4_field.npz        采样场、半径、烘干时刻
  结果/p4_论文表格.md      论文表6
  结果/p4_检验数据.md      收敛性、模型A对比、附件2半径相容性检查
"""
import os
import shutil
import numpy as np
import pandas as pd
import openpyxl

from common import load_env, make_graded_mesh, CylGrid, theta_step, tridiag_solve

R0 = 0.02
H, HM = 25.0, 8e-7
T0, C0 = 28.0, 2.55
N_BE = 2
C_CRIT = 0.15
T_MAX4 = 130 * 3600.0
R_OUT = np.arange(21) * 1e-3            # 固定物理输出半径 0:0.1:2 cm

# ξ 域网格三级（由问题1-3物理网格按 R0 换算），表面加密带宽 0.025（=0.5 mm/R0）
XI_BAND = 0.025
XI_MESHES = [(6.25e-3, 2.5e-5, 1.3 ** 0.25),
             (3.125e-3, 1.25e-5, 1.3 ** 0.125),
             (1.5625e-3, 6.25e-6, 1.3 ** 0.0625)]
SCHED4_PROD = ((0.03125, 60.0), (0.25, 10800.0), (1.0, T_MAX4))
SCHED4_ALT = ((0.0625, 60.0), (0.5, 10800.0), (2.0, T_MAX4))


def rho_of(C):
    return 760.0 + 90.0 * C


def cp_of(C):
    return 1850.0 + 2150.0 * C / (C + 1.0)


def k_of(C):
    return 0.12 + 0.20 * C / (C + 1.0)


D_SCALE = 1.0                          # D 前置系数倍率（灵敏度分析用）


def set_params(fh=1.0, fhm=1.0, fD=1.0):
    """灵敏度分析入口：按倍率缩放 h、h_m 与 D 前置系数（恢复传 1,1,1）。"""
    global H, HM, D_SCALE
    H = 25.0 * fh
    HM = 8e-7 * fhm
    D_SCALE = fD


def D_of(C, T):
    """附录4扩散系数（T 为 ℃，式内换算 K；两个负号勿丢）。"""
    CC = np.maximum(C, 1e-9)
    return D_SCALE * 4.2e-4 * np.exp(-0.30 / CC) * np.exp(-3850.0 / (T + 273.15))


def load_radius(path=os.path.join('附件', '附件2.xlsx')):
    """附件2半径：分段线性插值函数（m），数据末端后取末值；并返回平台起始时刻。"""
    df = pd.read_excel(path)
    t = df.iloc[:, 0].to_numpy(float)
    R = df.iloc[:, 1].to_numpy(float) / 100.0
    i_pl = int(np.argmax(R <= R[-1] + 1e-8))

    def Rf(tt):
        return float(np.interp(tt, t, R))

    return Rf, float(t[i_pl]), float(R[-1])


def theta_step_adv(u, dt, grid, gamma_face, cap, beta, amb0, amb1, theta, w):
    """模型A专用：在 common.theta_step 基础上并入伪对流项 w·∂u/∂ξ（w≤0），
    特征方向迎风（取内侧差分），隐式 θ 加权，保持三对角结构。"""
    G = grid.rf * gamma_face / grid.dr
    B = grid.R * beta
    M = cap * grid.V / dt
    a = np.zeros_like(u)                    # 迎风对流系数 c_i（i≥1）
    a[1:] = cap[1:] * grid.V[1:] * w[1:] / grid.dr
    Ladv0 = np.zeros_like(u)
    Ladv0[1:] = a[1:] * (u[1:] - u[:-1])
    q = G * (u[1:] - u[:-1])
    L0 = np.zeros_like(u)
    L0[:-1] += q
    L0[1:] -= q
    L0[-1] += B * (amb0 - u[-1])
    rhs = M * u + (1 - theta) * (L0 + Ladv0)
    rhs[-1] += theta * B * amb1
    sub = np.zeros_like(u)
    sup = np.zeros_like(u)
    diag = M.copy()
    sub[1:] = -theta * G
    sup[:-1] = -theta * G
    diag[:-1] += theta * G
    diag[1:] += theta * G
    diag[-1] += theta * B
    diag[1:] += -theta * a[1:]              # w≤0 → -θa ≥ 0，保持对角占优
    sub[1:] += theta * a[1:]
    return tridiag_solve(sub, diag, sup, rhs)


def step_coupled4(T, C, dt, grid, th, Ta0, Ta1, Ca0, Ca1, Rm, w=None,
                  tol=1e-6, itmax=30):
    """ξ 域耦合单步（预估-校正）。Rm 为该步中点半径；w 为模型A的
    伪对流系数节点值（None 即模型B）。"""
    ones = np.ones(grid.N + 1)
    T_new, C_new = T.copy(), C.copy()
    for it in range(1, itmax + 1):
        T_mid = 0.5 * (T + T_new)
        C_mid = 0.5 * (C + C_new)
        capT = Rm * Rm * rho_of(C_mid) * cp_of(C_mid)
        capC = Rm * Rm * ones
        k_node = k_of(C_mid)
        D_node = D_of(C_mid, T_mid)
        kf = 0.5 * (k_node[:-1] + k_node[1:])
        Df = 0.5 * (D_node[:-1] + D_node[1:])
        bT, bC = Rm * H, Rm * HM
        if w is None:
            T_next = theta_step(T, dt, grid, kf, capT, bT, Ta0, Ta1, th)
            C_next = theta_step(C, dt, grid, Df, capC, bC, Ca0, Ca1, th)
        else:
            T_next = theta_step_adv(T, dt, grid, kf, capT, bT, Ta0, Ta1, th, w)
            C_next = theta_step_adv(C, dt, grid, Df, capC, bC, Ca0, Ca1, th, w)
        delta = max(float(np.max(np.abs(T_next - T_new))),
                    float(np.max(np.abs(C_next - C_new))))
        T_new, C_new = T_next, C_next
        if delta < tol:
            break
    else:
        raise RuntimeError(f'Moving-domain iteration failed: delta={delta:.3e}, limit={itmax}')
    return T_new, C_new, it, capT


def solve_p4(h_int, h_min, q, sched, model='B', c_crit=C_CRIT, t_max=T_MAX4):
    """求解至中心水分跌破 c_crit（延拓到其后首个 60 s 整点）。
    每 60 s 采样：固定物理半径网格上的水分（r>R(t) 为 NaN）、表面值、
    半径、平均含水率与几何加权量 W=R²ΣCV（不直接等于物理水质量）。"""
    xi = make_graded_mesh(1.0, h_int, h_min, q, XI_BAND, 1e-3)
    grid = CylGrid(xi)
    Ta, Ca = load_env()
    Rf, t_plateau, R_end = load_radius()

    T = np.full(grid.N + 1, T0)
    C = np.full(grid.N + 1, C0)
    nrow_max = int(round(t_max / 60))
    Cout = np.full((nrow_max, R_OUT.size), np.nan)
    Csurf = np.empty(nrow_max)
    Tsurf = np.empty(nrow_max)
    Rsamp = np.empty(nrow_max)
    Wsamp = np.empty(nrow_max)
    nrow = 0

    Vsum = float(np.sum(grid.V))
    mass_lhs = mass_rhs = 0.0
    picard_max = 0
    t_dry = None
    C_dry = None
    R_dry = None
    t_last_row = None
    center_is_max = True
    n_glob = 0
    t_now = 0.0
    done = False

    for dt, t_stop in sched:
        if done:
            break
        t_stop = min(t_stop, t_max)
        n_save = int(round(60.0 / dt))
        nstep = int(round((t_stop - t_now) / dt))
        tgrid = t_now + np.arange(nstep + 1) * dt
        Ta_all, Ca_all = Ta(tgrid), Ca(tgrid)
        for n in range(nstep):
            th = 1.0 if n_glob < N_BE else 0.5
            t0n, t1n = tgrid[n], tgrid[n + 1]
            R0_, R1_ = Rf(t0n), Rf(t1n)
            Rm = 0.5 * (R0_ + R1_)
            w = None
            if model == 'A':
                Rdot = (R1_ - R0_) / dt
                w = grid.r * Rdot / Rm
            Ca0, Ca1 = float(Ca_all[n]), float(Ca_all[n + 1])
            Cs_old = C[-1]
            C_old = C
            T, C, its, _ = step_coupled4(T, C, dt, grid, th,
                                         float(Ta_all[n]), float(Ta_all[n + 1]),
                                         Ca0, Ca1, Rm, w)
            picard_max = max(picard_max, its)
            mass_lhs += Rm * Rm * float(np.sum(grid.V * (C - C_old)))
            mass_rhs += dt * Rm * HM * (th * (Ca1 - C[-1]) + (1 - th) * (Ca0 - Cs_old))
            n_glob += 1
            t_cur = t1n
            if t_dry is None and C_old[0] >= c_crit > C[0]:
                frac = (C_old[0] - c_crit) / (C_old[0] - C[0])
                t_dry = t_cur - dt + frac * dt
                C_dry = C_old + frac * (C - C_old)
                R_dry = Rf(t_dry)
                t_last_row = np.ceil(t_dry / 60.0) * 60.0
            if C[0] < np.max(C) - 1e-12:
                center_is_max = False
            if (n + 1) % n_save == 0:
                Rc = Rf(t_cur)
                ok = R_OUT <= Rc + 1e-12
                Cout[nrow, ok] = np.interp(R_OUT[ok] / Rc, grid.r, C)
                Csurf[nrow] = C[-1]
                Tsurf[nrow] = T[-1]
                Rsamp[nrow] = Rc
                Wsamp[nrow] = Rc * Rc * float(np.sum(grid.V * C))
                nrow += 1
                if t_last_row is not None and t_cur >= t_last_row - 1e-9:
                    done = True
                    break
        t_now = t_stop

    mass_res = abs(mass_lhs - mass_rhs) / max(abs(mass_lhs), 1e-300)
    times = 60.0 * np.arange(1, nrow + 1)
    out = dict(times=times, r_cm=R_OUT * 100, C=Cout[:nrow],
               Csurf=Csurf[:nrow], Tsurf=Tsurf[:nrow], Rsamp=Rsamp[:nrow],
               W=Wsamp[:nrow], meanC=Wsamp[:nrow] / (Rsamp[:nrow] ** 2 * Vsum),
               t_dry=t_dry, R_dry=R_dry, mass_res=mass_res,
               picard_max=picard_max, n_nodes=grid.N + 1,
               center_is_max=center_is_max, t_plateau=t_plateau, R_end=R_end)
    if C_dry is not None:
        out['C_dry_xi'] = C_dry
        out['xi'] = grid.r
    return out


def write_result4(sol):
    os.makedirs('结果', exist_ok=True)
    dst = os.path.join('结果', 'result4.xlsx')
    shutil.copyfile(os.path.join('附件', '附件3', 'result4.xlsx'), dst)
    wb = openpyxl.load_workbook(dst)
    ws = wb['Sheet1']
    for j in range(21):
        ws.cell(1, 2 + j).value = round(0.1 * j, 1)
    ws.cell(1, 23).value = '药材表面'
    for i, tv in enumerate(sol['times']):
        ws.cell(2 + i, 1).value = int(round(tv))
        row = sol['C'][i]
        for j in range(21):
            v = row[j]
            if np.isfinite(v):
                cell = ws.cell(2 + i, 2 + j)
                cell.value = round(float(v), 4)
                cell.number_format = '0.0000'
        cell = ws.cell(2 + i, 23)
        cell.value = round(float(sol['Csurf'][i]), 4)
        cell.number_format = '0.0000'
    wb.save(dst)
    return dst


def md_table6(sol):
    """表6：每 6 h 及烘干结束行；列 0/0.5/1/1.5 cm（超出当前半径记 —）+ 表面。"""
    cols_cm = (0.0, 0.5, 1.0, 1.5)
    head = '| 时间/h | ' + ' | '.join(f'{c:g} cm' for c in cols_cm) + ' | 药材表面 |'
    lines = [head, '|' + '---|' * (len(cols_cm) + 2)]

    def fmt_row(label, cvals, csurf, Rc):
        cells = []
        for c_cm, v in zip(cols_cm, cvals):
            if c_cm / 100 <= Rc + 1e-12 and np.isfinite(v):
                cells.append(f'{v:.4f}')
            else:
                cells.append('—')
        return f'| {label} | ' + ' | '.join(cells) + f' | {csurf:.4f} |'

    hh = 6
    while hh * 3600 < sol['t_dry']:
        i = int(round(hh * 3600 / 60)) - 1
        vals = [sol['C'][i][int(round(c * 10))] for c in cols_cm]
        lines.append(fmt_row(f'{hh:d}', vals, sol['Csurf'][i], sol['Rsamp'][i]))
        hh += 6
    xi_t = np.array([0.0, 0.005, 0.01, 0.015]) / sol['R_dry']
    vals = [float(np.interp(x, sol['xi'], sol['C_dry_xi'])) if x <= 1 else np.nan
            for x in xi_t]
    lines.append(fmt_row(f'{sol["t_dry"]/3600:.2f}（烘干结束）', vals,
                         float(sol['C_dry_xi'][-1]), sol['R_dry']))
    return '\n'.join(lines)


def main():
    print('问题4：生产解（模型B，最细网格 + 三段步长）……')
    sol = solve_p4(*XI_MESHES[-1], SCHED4_PROD, model='B')
    print(f'  t_dry = {sol["t_dry"]/3600:.4f} h, R(t_dry) = {sol["R_dry"]*100:.3f} cm')
    print('问题4：网格收敛检验……')
    sol_c = solve_p4(*XI_MESHES[0], SCHED4_PROD, model='B')
    sol_m = solve_p4(*XI_MESHES[1], SCHED4_PROD, model='B')
    print(f'  粗 {sol_c["t_dry"]/3600:.4f} h，中 {sol_m["t_dry"]/3600:.4f} h')
    print('问题4：步长加倍检验……')
    sol_t = solve_p4(*XI_MESHES[-1], SCHED4_ALT, model='B')
    print(f'  步长加倍 {sol_t["t_dry"]/3600:.4f} h')
    print('问题4：对比模型A（含伪对流项，中网格）……')
    sol_a = solve_p4(*XI_MESHES[1], SCHED4_PROD, model='A')
    print(f'  模型A t_dry = {sol_a["t_dry"]/3600:.4f} h')

    os.makedirs('结果', exist_ok=True)
    np.savez(os.path.join('结果', 'p4_field.npz'),
             times=sol['times'], r_cm=sol['r_cm'], C=sol['C'],
             Csurf=sol['Csurf'], Rsamp=sol['Rsamp'], meanC=sol['meanC'],
             t_dry=sol['t_dry'], R_dry=sol['R_dry'])
    dst = write_result4(sol)
    tab6 = md_table6(sol)
    with open(os.path.join('结果', 'p4_论文表格.md'), 'w', encoding='utf-8') as f:
        f.write('## 表6 药材烘干过程的水分浓度（kg/kg）\n\n' + tab6 +
                f'\n\n烘干时长：t_dry = {sol["t_dry"]:.1f} s = '
                f'{sol["t_dry"]/3600:.4f} h ≈ {sol["t_dry"]/3600:.2f} h；'
                f'结束时半径 {sol["R_dry"]*100:.3f} cm\n')

    # 检验统计
    nr = min(sol['C'].shape[0], sol_c['C'].shape[0], sol_m['C'].shape[0],
             sol_t['C'].shape[0])

    def fdiff(a, b):
        d = np.abs(a['C'][:nr] - b['C'][:nr])
        return float(np.nanmax(d))

    dmesh1, dmesh2, ddt = fdiff(sol_c, sol_m), fdiff(sol_m, sol), fdiff(sol_t, sol)
    p_mesh = np.log2(dmesh1 / dmesh2) if dmesh2 > 0 else float('nan')
    # 输入半径与几何加权量的相容性描述，不是独立实验验证。
    nr2 = min(sol['W'].size, sol_a['W'].size)
    corrB = np.corrcoef(sol['Rsamp'][:nr2] ** 2, sol['W'][:nr2])[0, 1] ** 2
    corrA = np.corrcoef(sol_a['Rsamp'][:nr2] ** 2, sol_a['W'][:nr2])[0, 1] ** 2
    i_surf = int(np.argmax(sol['Csurf'] < C_CRIT))
    ver = [
        '## 问题4 数值检验记录', '',
        f'主模型B（物质坐标，无伪对流项）生产解：t_dry = {sol["t_dry"]:.1f} s = '
        f'{sol["t_dry"]/3600:.4f} h；结束时半径 {sol["R_dry"]*100:.4f} cm；'
        f'中心恒为最大值验证：{"通过" if sol["center_is_max"] else "未通过"}。', '',
        '| 检验 | t_dry/h | 与生产解之差/h | 采样场最大差 |', '|---|---|---|---|',
        f'| 粗网格({sol_c["n_nodes"]}节点) | {sol_c["t_dry"]/3600:.4f} | '
        f'{(sol_c["t_dry"]-sol["t_dry"])/3600:+.4f} | {dmesh1:.2e}(粗-中) |',
        f'| 中网格({sol_m["n_nodes"]}节点) | {sol_m["t_dry"]/3600:.4f} | '
        f'{(sol_m["t_dry"]-sol["t_dry"])/3600:+.4f} | {dmesh2:.2e}(中-细) |',
        f'| 步长加倍方案 | {sol_t["t_dry"]/3600:.4f} | '
        f'{(sol_t["t_dry"]-sol["t_dry"])/3600:+.4f} | {ddt:.2e} |', '',
        f'网格方向观测阶 p≈{p_mesh:.2f}；模型B水分方程离散收支残差 {sol["mass_res"]:.2e}；'
        f'Picard 最大迭代 {sol["picard_max"]}。', '',
        f'模型A（含伪对流项）对比：t_dry = {sol_a["t_dry"]/3600:.4f} h，'
        f'较模型B相差 {(sol_a["t_dry"]-sol["t_dry"])/3600:+.2f} h。'
        'A与B采用不同运动学假设；A是固定空间扩散场的坐标变换，'
        'B采用随骨架运动的有效扩散闭合。二者差异只作模型结构对照。', '',
        f'附件2相容性检查：半径平台起始 ≈ {sol["t_plateau"]/3600:.1f} h'
        f'（末值 {sol["R_end"]*100:.3f} cm），模型B烘干完成 {sol["t_dry"]/3600:.2f} h '
        f'早于最终平台，结束时半径已达 {sol["R_dry"]*100:.3f} cm'
        f'（与末值差 {abs(sol["R_dry"]-sol["R_end"])*100:.3f} cm）。'
        f'R(t)²与几何加权量W的相关平方：模型B {corrB:.4f}，模型A {corrA:.4f}。'
        'W不直接等于物理水质量，且R为模型输入，相关性不是独立准确性检验。', '',
        f'表面水分降至 {C_CRIT} 以下时刻 ≈ {sol["times"][i_surf]/3600:.2f} h。',
    ]
    with open(os.path.join('结果', 'p4_检验数据.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(ver) + '\n')

    print('=== 问题4 求解完成 ===')
    print(f'烘干时长（模型B）：{sol["t_dry"]:.1f} s = {sol["t_dry"]/3600:.4f} h '
          f'≈ {sol["t_dry"]/3600:.2f} h；结束半径 {sol["R_dry"]*100:.3f} cm')
    print(f'检验：粗/中/细 t_dry = {sol_c["t_dry"]/3600:.4f} / {sol_m["t_dry"]/3600:.4f}'
          f' / {sol["t_dry"]/3600:.4f} h；步长加倍 {sol_t["t_dry"]/3600:.4f} h；'
          f'模型A {sol_a["t_dry"]/3600:.4f} h')
    print(f'方程收支残差 {sol["mass_res"]:.1e}; 几何量相关平方: B {corrB:.4f} / A {corrA:.4f}; '
          f'结果文件 {dst}（{sol["C"].shape[0]} 行）')
    print('--- 表6 ---')
    print(tab6)


if __name__ == '__main__':
    main()
