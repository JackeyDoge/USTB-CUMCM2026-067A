# -*- coding: utf-8 -*-
"""问题2（及问题3共用求解器）：全烘干过程变物性热质耦合模型。

模型：一维柱坐标径向，ρ(C)、c_p(C)、k(C) 随水分变化，D(C,T) 呈
Arrhenius 型温度依赖（附录3），热-质双向耦合；边界条件同问题1，
环境温湿度 t≤14400 s 取附件1分段线性插值，之后取恒温段设计值
50 ℃ / 0.05 kg/kg（common.py）。

数值：有限体积（表面加密网格）+ Crank-Nicolson（首两步后向 Euler
启动）+ 温度/水分联立预估-校正迭代（每次迭代以中点状态更新全部
物性后依次解两方程，直至两场增量均小于阈值）。

运行方式（项目根目录下）：python code/p23.py
输出：
  结果/result2.xlsx          问题2完整结果（每 1 s × 每 0.1 cm，四位小数）
  结果/p2_field.npz          前 3 h 数值场（供绘图与问题3复用检查）
  结果/p2_论文表格.md        论文表3/表4数值
  结果/p2_检验数据.md        空间/时间收敛性、守恒与迭代记录
"""
import os
import shutil
import numpy as np
import openpyxl

from common import load_env, make_graded_mesh, CylGrid, theta_step

# ---- 题目与附录3给定参数（SI 单位；温度场以 ℃ 计算，仅 D 内换算 K） ----
R0 = 0.02
H, HM = 25.0, 8e-7
T0, C0 = 28.0, 2.55
BAND = 5e-4
N_BE = 2
R_OUT = np.arange(21) * 1e-3
T_END2 = 10800.0                     # 问题2输出时长 3 h

TAB_H = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)          # 论文表3/表4时间点 h
TAB_J = (0, 5, 10, 15, 20)                       # 论文表格距离下标

# 空间收敛研究的三级网格（固定步长方案）与时间收敛研究的三级步长方案
# 步长方案为分段表 ((dt_a, 换挡时刻), (dt_b, 终止时刻))：初始瞬态阶段用小步长
# 精细解析（表面节点早期变化剧烈），其后用常规步长。
MESHES = [(1.25e-4, 5e-7, 1.3 ** 0.25),
          (6.25e-5, 2.5e-7, 1.3 ** 0.125),
          (3.125e-5, 1.25e-7, 1.3 ** 0.0625)]
SCHEDS = [((0.125, 60.0), (1.0, T_END2)),
          ((0.0625, 60.0), (0.5, T_END2)),
          ((0.03125, 60.0), (0.25, T_END2))]
SCHED_PROD = SCHEDS[-1]


def rho_of(C):
    return 650.0 + 128.0 * C


def cp_of(C):
    return 1450.0 + 2736.0 * C / (C + 1.0)


def k_of(C):
    return 0.21 + 0.38 * C / (C + 1.0)


D_SCALE = 1.0                          # D 前置系数倍率（灵敏度分析用）


def set_params(fh=1.0, fhm=1.0, fD=1.0):
    """灵敏度分析入口：按倍率缩放 h、h_m 与 D 前置系数（恢复传 1,1,1）。"""
    global H, HM, D_SCALE
    H = 25.0 * fh
    HM = 8e-7 * fhm
    D_SCALE = fD


def D_of(C, T):
    """附录3扩散系数；C 下限保护避免 exp 溢出，T 为 ℃、式内换算 K。
    注意两个负号：D 随 C 减小而减小、随 T 升高而增大。"""
    CC = np.maximum(C, 1e-9)
    return D_SCALE * 2.4e-3 * np.exp(-0.45 / CC) * np.exp(-3850.0 / (T + 273.15))


def step_coupled(T, C, dt, grid, th, Ta0, Ta1, Ca0, Ca1, tol=1e-6, itmax=30):
    """耦合方程组的 θ 格式单步：预估-校正迭代，中点物性。
    返回 (T_new, C_new, 迭代次数, 末次迭代所用热容量 capT)。"""
    ones = np.ones(grid.N + 1)
    T_new, C_new = T.copy(), C.copy()
    for it in range(1, itmax + 1):
        T_mid = 0.5 * (T + T_new)
        C_mid = 0.5 * (C + C_new)
        capT = rho_of(C_mid) * cp_of(C_mid)
        k_node = k_of(C_mid)
        D_node = D_of(C_mid, T_mid)
        kf = 0.5 * (k_node[:-1] + k_node[1:])
        Df = 0.5 * (D_node[:-1] + D_node[1:])
        T_next = theta_step(T, dt, grid, kf, capT, H, Ta0, Ta1, th)
        C_next = theta_step(C, dt, grid, Df, ones, HM, Ca0, Ca1, th)
        delta = max(float(np.max(np.abs(T_next - T_new))),
                    float(np.max(np.abs(C_next - C_new))))
        T_new, C_new = T_next, C_next
        if delta < tol:
            break
    else:
        raise RuntimeError(f'Coupled iteration failed: delta={delta:.3e}, limit={itmax}')
    return T_new, C_new, it, capT


def solve_p2(h_int, h_min, q, sched, t_end=T_END2):
    """按分段步长方案 sched 求解 0~t_end，按每 1 s、0:0.1:2 cm 采样返回。
    后向 Euler 启动仅作用于全程最前 N_BE 步。"""
    grid = CylGrid(make_graded_mesh(R0, h_int, h_min, q, BAND, 1e-3))
    oi = grid.out_index(R_OUT)
    Ta, Ca = load_env()

    T = np.full(grid.N + 1, T0)
    C = np.full(grid.N + 1, C0)
    times = np.arange(1, int(t_end) + 1)
    Tout = np.empty((times.size, oi.size))
    Cout = np.empty((times.size, oi.size))

    Bt, Bm = grid.R * H, grid.R * HM
    mass0 = float(np.sum(C * grid.V))
    mass_in = 0.0
    ener_lhs = ener_rhs = 0.0
    picard_max = 0

    t_now = 0.0
    n_glob = 0
    for dt, t_stop in sched:
        t_stop = min(t_stop, t_end)
        n_save = int(round(1.0 / dt))
        assert abs(n_save * dt - 1.0) < 1e-12, 'dt 须整除 1 s'
        assert abs(round(t_now) - t_now) < 1e-9, '换挡时刻须为整秒'
        nstep = int(round((t_stop - t_now) / dt))
        tgrid = t_now + np.arange(nstep + 1) * dt
        Ta_all, Ca_all = Ta(tgrid), Ca(tgrid)
        for n in range(nstep):
            th = 1.0 if n_glob < N_BE else 0.5
            Ta0, Ta1 = float(Ta_all[n]), float(Ta_all[n + 1])
            Ca0, Ca1 = float(Ca_all[n]), float(Ca_all[n + 1])
            Ts_old, Cs_old = T[-1], C[-1]
            T_old = T
            T, C, its, capT = step_coupled(T, C, dt, grid, th,
                                           Ta0, Ta1, Ca0, Ca1)
            picard_max = max(picard_max, its)
            mass_in += dt * (th * Bm * (Ca1 - C[-1]) + (1 - th) * Bm * (Ca0 - Cs_old))
            ener_lhs += float(np.sum(capT * grid.V * (T - T_old)))
            ener_rhs += dt * (th * Bt * (Ta1 - T[-1]) + (1 - th) * Bt * (Ta0 - Ts_old))
            n_glob += 1
            if (n + 1) % n_save == 0:
                j = int(round(t_now + (n + 1) * dt)) - 1
                if 0 <= j < times.size:
                    Tout[j] = T[oi]
                    Cout[j] = C[oi]
        t_now = t_stop
        if t_now >= t_end:
            break

    mass_res = abs((float(np.sum(C * grid.V)) - mass0 - mass_in) / mass0)
    ener_res = abs(ener_lhs - ener_rhs) / max(abs(ener_lhs), 1e-300)
    return dict(times=times, r_cm=R_OUT * 100, T=Tout, C=Cout,
                mass_res=mass_res, ener_res=ener_res,
                picard_max=picard_max, n_nodes=grid.N + 1,
                T_final=T, C_final=C)


def write_result2(sol):
    os.makedirs('结果', exist_ok=True)
    dst = os.path.join('结果', 'result2.xlsx')
    shutil.copyfile(os.path.join('附件', '附件3', 'result2.xlsx'), dst)
    wb = openpyxl.load_workbook(dst)
    for name, data in (('温度', sol['T']), ('水分浓度', sol['C'])):
        ws = wb[name]
        for j in range(21):
            ws.cell(1, 2 + j).value = round(0.1 * j, 1)
        for i, tv in enumerate(sol['times']):
            ws.cell(2 + i, 1).value = int(tv)
            row = data[i]
            for j in range(21):
                cell = ws.cell(2 + i, 2 + j)
                cell.value = round(float(row[j]), 4)
                cell.number_format = '0.0000'
    wb.save(dst)
    return dst


def md_table(sol, field):
    head = '| 时间/h | 0 cm | 0.5 cm | 1 cm | 1.5 cm | 2 cm |'
    lines = [head, '|---|---|---|---|---|---|']
    for th_ in TAB_H:
        row = sol[field][int(round(th_ * 3600)) - 1]
        lines.append('| %.1f | %s |' % (th_, ' | '.join(f'{row[j]:.4f}' for j in TAB_J)))
    return '\n'.join(lines)


def diff_stats(a, b):
    dT = np.abs(a['T'] - b['T'])
    dC = np.abs(a['C'] - b['C'])
    it = np.ix_([int(round(h * 3600)) - 1 for h in TAB_H], list(TAB_J))
    return dict(T_all=dT.max(), C_all=dC.max(),
                T_tab=dT[it].max(), C_tab=dC[it].max())


def main():
    # 空间收敛（步长方案固定为生产方案）：粗、中网格 + 最细（即生产解）
    sols_s = [solve_p2(*m, SCHED_PROD) for m in MESHES]
    ds = [diff_stats(a, b) for a, b in zip(sols_s[:-1], sols_s[1:])]
    # 时间收敛（最细网格固定，两段步长方案整体逐级减半）
    sols_t = [solve_p2(*MESHES[-1], sc) for sc in SCHEDS[:-1]] + [sols_s[-1]]
    dt_ = [diff_stats(a, b) for a, b in zip(sols_t[:-1], sols_t[1:])]
    sol = sols_s[-1]                     # 生产解：最细网格 + 生产步长方案

    os.makedirs('结果', exist_ok=True)
    np.savez(os.path.join('结果', 'p2_field.npz'),
             times=sol['times'], r_cm=sol['r_cm'], T=sol['T'], C=sol['C'])
    dst = write_result2(sol)

    tab3, tab4 = md_table(sol, 'T'), md_table(sol, 'C')
    with open(os.path.join('结果', 'p2_论文表格.md'), 'w', encoding='utf-8') as f:
        f.write('## 表3 3小时内药材的温度（℃）\n\n' + tab3 +
                '\n\n## 表4 3小时内药材的水分浓度（kg/kg）\n\n' + tab4 + '\n')

    keys = (('T_all', '温度(全部输出点)'), ('C_all', '水分(全部输出点)'),
            ('T_tab', '温度(论文表格点)'), ('C_tab', '水分(论文表格点)'))

    def conv_block(title, pairs, labels):
        ln = [title, '| 相邻对 | ' + ' | '.join(n for _, n in keys) + ' |',
              '|' + '---|' * (len(keys) + 1)]
        for lab, d in zip(labels, pairs):
            ln.append(f'| {lab} | ' + ' | '.join(f'{d[k]:.3e}' for k, _ in keys) + ' |')
        if len(pairs) >= 2:
            p = {k: np.log2(pairs[-2][k] / pairs[-1][k]) for k, _ in keys}
            rich = {k: pairs[-1][k] / (2 ** p[k] - 1) if p[k] > 0 else float('nan')
                    for k, _ in keys}
            ln.append('| 观测阶 p | ' + ' | '.join(f'{p[k]:.2f}' for k, _ in keys) + ' |')
            ln.append('| 生产解误差估计 | ' + ' | '.join(f'{rich[k]:.2e}' for k, _ in keys) + ' |')
        return ln

    ver = ['## 问题2 数值检验记录', '',
           f'网格节点数（三级）：' + '、'.join(str(s['n_nodes']) for s in sols_s) +
           f'；生产解：最细网格 {sol["n_nodes"]} 节点，步长方案 '
           f'{SCHED_PROD[0][0]} s（前 {SCHED_PROD[0][1]:.0f} s）→ {SCHED_PROD[1][0]} s，'
           f'首 {N_BE} 步后向 Euler 启动。', '']
    ver += conv_block('### 空间收敛（步长方案固定，网格逐级减半）', ds,
                      ['粗-中', '中-细'])
    ver += ['']
    ver += conv_block('### 时间收敛（最细网格固定，两段步长整体逐级减半）', dt_,
                      ['S0-S1', 'S1-S2'])
    ver += ['',
            f'水分离散守恒残差（生产解）：{sol["mass_res"]:.2e}（相对值）。',
            f'温度方程逐步能量恒等式残差：{sol["ener_res"]:.2e}（相对值；'
            '变热容下按每步实际使用的 capT 累计）。',
            f'耦合预估-校正迭代最大次数：{sol["picard_max"]}（阈值 1e-6）。']
    with open(os.path.join('结果', 'p2_检验数据.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(ver) + '\n')

    print('=== 问题2 求解完成（生产解：%d 节点，步长 %g s→%g s）===' %
          (sol['n_nodes'], SCHED_PROD[0][0], SCHED_PROD[1][0]))
    print('\n'.join(conv_block('空间收敛：', ds, ['粗-中', '中-细'])))
    print('\n'.join(conv_block('时间收敛：', dt_, ['S0-S1', 'S1-S2'])))
    print(f'守恒/恒等式残差: 水分 {sol["mass_res"]:.1e}, 能量 {sol["ener_res"]:.1e}; '
          f'Picard 最大迭代 {sol["picard_max"]}')
    print(f'结果文件: {dst}')
    print('--- 表3 温度/℃ ---')
    print(tab3)
    print('--- 表4 水分浓度/(kg/kg) ---')
    print(tab4)


# ---------------------------------------------------------------------------
# 问题3：延拓至烘干结束（各处 C < 0.15，等价于中心 C(0,t) < 0.15）
# ---------------------------------------------------------------------------
C_CRIT = 0.15
T_MAX3 = 100 * 3600.0                 # 安全上限 100 h
SCHED3_PROD = ((0.03125, 60.0), (0.25, 10800.0), (1.0, T_MAX3))
SCHED3_ALT = ((0.0625, 60.0), (0.5, 10800.0), (2.0, T_MAX3))
TAB5_J = TAB_J


def solve_p3(h_int, h_min, q, sched, c_crit=C_CRIT, t_max=T_MAX3):
    """延拓求解至中心水分跌破 c_crit。每 60 s 采样水分/温度；
    跨越判据时对相邻两步状态线性插值，得烘干时刻 t_dry 与该时刻场。
    继续算到 ceil(t_dry/60)*60 以便结果文件末行覆盖烘干完成状态。"""
    grid = CylGrid(make_graded_mesh(R0, h_int, h_min, q, BAND, 1e-3))
    oi = grid.out_index(R_OUT)
    Ta, Ca = load_env()

    T = np.full(grid.N + 1, T0)
    C = np.full(grid.N + 1, C0)
    nrow_max = int(round(t_max / 60))
    Crows = np.empty((nrow_max, oi.size))
    Trows = np.empty((nrow_max, oi.size))
    nrow = 0

    Bm = grid.R * HM
    mass0 = float(np.sum(C * grid.V))
    mass_in = 0.0
    picard_max = 0
    t_dry = None
    C_dry = None                       # t_dry 时刻的全场水分（节点上）
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
        assert abs(n_save * dt - 60.0) < 1e-9, 'dt 须整除 60 s'
        nstep = int(round((t_stop - t_now) / dt))
        tgrid = t_now + np.arange(nstep + 1) * dt
        Ta_all, Ca_all = Ta(tgrid), Ca(tgrid)
        for n in range(nstep):
            th = 1.0 if n_glob < N_BE else 0.5
            Ca0, Ca1 = float(Ca_all[n]), float(Ca_all[n + 1])
            Cs_old = C[-1]
            C_old = C
            T, C, its, _ = step_coupled(T, C, dt, grid, th,
                                        float(Ta_all[n]), float(Ta_all[n + 1]),
                                        Ca0, Ca1)
            picard_max = max(picard_max, its)
            mass_in += dt * (th * Bm * (Ca1 - C[-1]) + (1 - th) * Bm * (Ca0 - Cs_old))
            n_glob += 1
            t_cur = t_now + (n + 1) * dt
            if t_dry is None and C_old[0] >= c_crit > C[0]:
                frac = (C_old[0] - c_crit) / (C_old[0] - C[0])
                t_dry = t_cur - dt + frac * dt
                C_dry = C_old + frac * (C - C_old)
                t_last_row = np.ceil(t_dry / 60.0) * 60.0
            if C[0] < np.max(C) - 1e-12:
                center_is_max = False
            if (n + 1) % n_save == 0:
                Crows[nrow] = C[oi]
                Trows[nrow] = T[oi]
                nrow += 1
                if t_last_row is not None and t_cur >= t_last_row - 1e-9:
                    done = True
                    break
        t_now = t_stop

    mass_res = abs((float(np.sum(C * grid.V)) - mass0 - mass_in) / mass0)
    times = 60.0 * np.arange(1, nrow + 1)
    return dict(times=times, r_cm=R_OUT * 100,
                C=Crows[:nrow], T=Trows[:nrow],
                t_dry=t_dry, C_dry_out=None if C_dry is None else C_dry[oi],
                mass_res=mass_res, picard_max=picard_max,
                n_nodes=grid.N + 1, center_is_max=center_is_max)


def write_result3(sol):
    os.makedirs('结果', exist_ok=True)
    dst = os.path.join('结果', 'result3.xlsx')
    shutil.copyfile(os.path.join('附件', '附件3', 'result3.xlsx'), dst)
    wb = openpyxl.load_workbook(dst)
    ws = wb['Sheet1']
    for j in range(21):
        ws.cell(1, 2 + j).value = round(0.1 * j, 1)
    for i, tv in enumerate(sol['times']):
        ws.cell(2 + i, 1).value = int(round(tv))
        row = sol['C'][i]
        for j in range(21):
            cell = ws.cell(2 + i, 2 + j)
            cell.value = round(float(row[j]), 4)
            cell.number_format = '0.0000'
    wb.save(dst)
    return dst


def md_table5(sol):
    head = '| 时间/h | 0 cm | 0.5 cm | 1 cm | 1.5 cm | 2 cm |'
    lines = [head, '|---|---|---|---|---|---|']
    hh = 6
    while hh * 3600 < sol['t_dry']:
        row = sol['C'][int(round(hh * 3600 / 60)) - 1]
        lines.append('| %d | %s |' % (hh, ' | '.join(f'{row[j]:.4f}' for j in TAB5_J)))
        hh += 6
    row = sol['C_dry_out']
    lines.append('| %.2f（烘干结束） | %s |' %
                 (sol['t_dry'] / 3600, ' | '.join(f'{row[j]:.4f}' for j in TAB5_J)))
    return '\n'.join(lines)


def main3():
    print('问题3：生产解（最细网格 + 三段步长）……')
    sol = solve_p3(*MESHES[-1], SCHED3_PROD)
    print(f'  t_dry = {sol["t_dry"]/3600:.4f} h')
    print('问题3：网格收敛检验（粗、中网格）……')
    sol_c = solve_p3(*MESHES[0], SCHED3_PROD)
    sol_m = solve_p3(*MESHES[1], SCHED3_PROD)
    print(f'  粗 {sol_c["t_dry"]/3600:.4f} h，中 {sol_m["t_dry"]/3600:.4f} h')
    print('问题3：步长收敛检验（最细网格，步长方案整体加倍）……')
    sol_t = solve_p3(*MESHES[-1], SCHED3_ALT)
    print(f'  加倍步长 {sol_t["t_dry"]/3600:.4f} h')

    os.makedirs('结果', exist_ok=True)
    np.savez(os.path.join('结果', 'p3_field.npz'),
             times=sol['times'], r_cm=sol['r_cm'], C=sol['C'], T=sol['T'],
             t_dry=sol['t_dry'])
    dst = write_result3(sol)
    tab5 = md_table5(sol)
    with open(os.path.join('结果', 'p3_论文表格.md'), 'w', encoding='utf-8') as f:
        f.write('## 表5 药材烘干过程的水分浓度（kg/kg）\n\n' + tab5 +
                f'\n\n烘干所需时间：t_dry = {sol["t_dry"]:.1f} s = '
                f'{sol["t_dry"]/3600:.4f} h ≈ {sol["t_dry"]/3600:.2f} h\n')

    nr = min(sol['C'].shape[0], sol_c['C'].shape[0], sol_m['C'].shape[0],
             sol_t['C'].shape[0])
    dmesh1 = float(np.max(np.abs(sol_c['C'][:nr] - sol_m['C'][:nr])))
    dmesh2 = float(np.max(np.abs(sol_m['C'][:nr] - sol['C'][:nr])))
    ddt = float(np.max(np.abs(sol_t['C'][:nr] - sol['C'][:nr])))
    p_mesh = np.log2(dmesh1 / dmesh2) if dmesh2 > 0 else float('nan')
    # 表面（最先达标）与中心的到达时刻，供论文干燥阶段划分
    i_surf = int(np.argmax(sol['C'][:, -1] < C_CRIT))
    ver = [
        '## 问题3 数值检验记录', '',
        f'烘干判据：各处 C<{C_CRIT}，等价于中心 C(0,t)<{C_CRIT}'
        f'（全程数值验证中心恒为最大值：{"通过" if sol["center_is_max"] else "未通过"}）。', '',
        f'生产解烘干时间 t_dry = {sol["t_dry"]:.1f} s = {sol["t_dry"]/3600:.4f} h。', '',
        '| 检验 | t_dry/h | 与生产解之差/h | 60s 采样场最大差 |', '|---|---|---|---|',
        f'| 粗网格({sol_c["n_nodes"]}节点) | {sol_c["t_dry"]/3600:.4f} | '
        f'{(sol_c["t_dry"]-sol["t_dry"])/3600:+.4f} | {dmesh1:.2e}(粗-中) |',
        f'| 中网格({sol_m["n_nodes"]}节点) | {sol_m["t_dry"]/3600:.4f} | '
        f'{(sol_m["t_dry"]-sol["t_dry"])/3600:+.4f} | {dmesh2:.2e}(中-细) |',
        f'| 步长加倍方案 | {sol_t["t_dry"]/3600:.4f} | '
        f'{(sol_t["t_dry"]-sol["t_dry"])/3600:+.4f} | {ddt:.2e} |', '',
        f'网格方向观测阶 p≈{p_mesh:.2f}；'
        f'生产解水分离散守恒残差 {sol["mass_res"]:.2e}；Picard 最大迭代 {sol["picard_max"]}。', '',
        f'表面（r=2 cm）最先降至 {C_CRIT} 以下的时刻 ≈ {sol["times"][i_surf]/3600:.2f} h，'
        f'其后内部扩散占主导，表面传质阻力仍保留在边界条件中。',
    ]
    with open(os.path.join('结果', 'p3_检验数据.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(ver) + '\n')

    print('=== 问题3 求解完成 ===')
    print(f'烘干所需时间：{sol["t_dry"]:.1f} s = {sol["t_dry"]/3600:.4f} h '
          f'≈ {sol["t_dry"]/3600:.2f} h')
    print(f'检验：粗/中/细网格 t_dry = {sol_c["t_dry"]/3600:.4f} / '
          f'{sol_m["t_dry"]/3600:.4f} / {sol["t_dry"]/3600:.4f} h；'
          f'步长加倍 {sol_t["t_dry"]/3600:.4f} h')
    print(f'守恒残差 {sol["mass_res"]:.1e}; Picard 最大迭代 {sol["picard_max"]}; '
          f'结果文件 {dst}（{sol["C"].shape[0]} 行）')
    print('--- 表5 水分浓度/(kg/kg) ---')
    print(tab5)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'p3':
        main3()
    else:
        main()
