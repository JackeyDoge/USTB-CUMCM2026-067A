# -*- coding: utf-8 -*-
"""问题1：预热平衡阶段常物性热质传递模型求解主程序。

模型：一维柱坐标径向热传导方程 + 非线性 Fick 扩散方程（D 依赖 C，二者解耦），
第三类边界条件，环境温湿度由附件1分段线性插值。
数值：有限体积（表面边界层加密网格）+ Crank-Nicolson（首两步后向 Euler
启动，抑制初始刚性振荡）+ 预估-校正迭代，见 common.py。

运行方式（项目根目录下）：python code/p1.py
输出：
  结果/result1.xlsx          完整结果（每 1 s × 每 0.1 cm，四位小数）
  结果/p1_field.npz          完整数值场（供绘图与检验复用）
  结果/p1_论文表格.md        论文表1/表2数值
  结果/p1_检验数据.md        网格收敛性与守恒残差记录
"""
import os
import shutil
import numpy as np
import openpyxl

from common import (load_env, make_graded_mesh, CylGrid,
                    theta_step, theta_step_picard)

# ---- 题目与附录2给定参数（SI 单位） ----
R0 = 0.02                      # 半径 m
RHO, CP, KH = 820.0, 2600.0, 0.36
H, HM = 25.0, 8e-7
T0, C0 = 28.0, 2.55
T_END = 1800.0
BAND = 5e-4                    # 表面加密带宽 m
N_BE = 2                       # Rannacher 启动：首两步用后向 Euler

R_OUT = np.arange(21) * 1e-3   # 输出半径 0:0.1:2 cm（m）
TAB_T = (100, 300, 600, 900, 1200, 1500, 1800)   # 论文表格时间点 s
TAB_J = (0, 5, 10, 15, 20)                       # 论文表格距离下标

# 网格/时间步收敛研究的五级参数：(h_int, h_min, q, dt)，逐级全部减半
LEVELS = [(5e-4, 2e-6, 1.3, 0.5),
          (2.5e-4, 1e-6, 1.3 ** 0.5, 0.25),
          (1.25e-4, 5e-7, 1.3 ** 0.25, 0.125),
          (6.25e-5, 2.5e-7, 1.3 ** 0.125, 0.0625),
          (3.125e-5, 1.25e-7, 1.3 ** 0.0625, 0.03125)]


def D_of_C(C):
    """附录2扩散系数经验公式（注意负号：D 随 C 减小而减小）。"""
    return 7e-9 * np.exp(-0.89 / np.maximum(C, 1e-9))


def solve(h_int, h_min, q, dt, t_end=T_END):
    """求解 0~t_end，按每 1 s、0:0.1:2 cm 采样返回。"""
    grid = CylGrid(make_graded_mesh(R0, h_int, h_min, q, BAND, 1e-3))
    oi = grid.out_index(R_OUT)
    Ta, Ca = load_env()
    n_save = int(round(1.0 / dt))
    assert abs(n_save * dt - 1.0) < 1e-12, 'dt 须整除 1 s'

    T = np.full(grid.N + 1, T0)
    C = np.full(grid.N + 1, C0)
    kf = np.full(grid.N, KH)
    capT = np.full(grid.N + 1, RHO * CP)
    capC = np.ones(grid.N + 1)

    nstep = int(round(t_end / dt))
    times = np.arange(1, int(t_end) + 1)
    Tout = np.empty((times.size, oi.size))
    Cout = np.empty((times.size, oi.size))

    Bt, Bm = grid.R * H, grid.R * HM
    mass0 = float(np.sum(C * grid.V))
    ener0 = float(np.sum(capT * T * grid.V))
    mass_in = ener_in = 0.0
    picard_max = 0

    for n in range(nstep):
        th = 1.0 if n < N_BE else 0.5
        t0n, t1n = n * dt, (n + 1) * dt
        Ta0, Ta1 = float(Ta(t0n)), float(Ta(t1n))
        Ca0, Ca1 = float(Ca(t0n)), float(Ca(t1n))
        Ts_old, Cs_old = T[-1], C[-1]
        T = theta_step(T, dt, grid, kf, capT, H, Ta0, Ta1, th)
        C, its = theta_step_picard(C, dt, grid, D_of_C, capC, HM,
                                   Ca0, Ca1, th)
        picard_max = max(picard_max, its)
        ener_in += dt * (th * Bt * (Ta1 - T[-1]) + (1 - th) * Bt * (Ta0 - Ts_old))
        mass_in += dt * (th * Bm * (Ca1 - C[-1]) + (1 - th) * Bm * (Ca0 - Cs_old))
        k = n + 1
        if k % n_save == 0:
            j = k // n_save - 1
            Tout[j] = T[oi]
            Cout[j] = C[oi]

    mass_res = abs((float(np.sum(C * grid.V)) - mass0 - mass_in) / mass0)
    ener_res = abs((float(np.sum(capT * T * grid.V)) - ener0 - ener_in) / ener0)
    return dict(times=times, r_cm=R_OUT * 100, T=Tout, C=Cout,
                mass_res=mass_res, ener_res=ener_res,
                picard_max=picard_max, n_nodes=grid.N + 1)


def write_result1(sol):
    os.makedirs('结果', exist_ok=True)
    dst = os.path.join('结果', 'result1.xlsx')
    shutil.copyfile(os.path.join('附件', '附件3', 'result1.xlsx'), dst)
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
    head = '| 时间/s | 0 cm | 0.5 cm | 1 cm | 1.5 cm | 2 cm |'
    lines = [head, '|---|---|---|---|---|---|']
    for tt in TAB_T:
        row = sol[field][tt - 1]
        lines.append('| %d | %s |' % (tt, ' | '.join(f'{row[j]:.4f}' for j in TAB_J)))
    return '\n'.join(lines)


def diff_stats(a, b):
    """两级解在输出点上的最大差：全域 / t≥60 s / 论文表格点。"""
    dT = np.abs(a['T'] - b['T'])
    dC = np.abs(a['C'] - b['C'])
    it = np.ix_([t - 1 for t in TAB_T], list(TAB_J))
    return dict(T_all=dT.max(), C_all=dC.max(),
                T_60=dT[59:].max(), C_60=dC[59:].max(),
                T_tab=dT[it].max(), C_tab=dC[it].max())


def main():
    sols = [solve(*lv) for lv in LEVELS]
    diffs = [diff_stats(a, b) for a, b in zip(sols[:-1], sols[1:])]
    sol = sols[-1]                       # 生产解：最细一级
    os.makedirs('结果', exist_ok=True)

    np.savez(os.path.join('结果', 'p1_field.npz'),
             times=sol['times'], r_cm=sol['r_cm'], T=sol['T'], C=sol['C'])
    dst = write_result1(sol)

    tab1, tab2 = md_table(sol, 'T'), md_table(sol, 'C')
    with open(os.path.join('结果', 'p1_论文表格.md'), 'w', encoding='utf-8') as f:
        f.write('## 表1 30分钟内药材的温度（℃）\n\n' + tab1 +
                '\n\n## 表2 30分钟内药材的水分浓度（kg/kg）\n\n' + tab2 + '\n')

    # 收敛表：逐对最大差与观测阶；末行给最细解的 Richardson 误差估计
    keys = (('T_all', '温度(全部输出点)'), ('C_all', '水分(全部输出点)'),
            ('T_tab', '温度(论文表格点)'), ('C_tab', '水分(论文表格点)'))
    conv_lines = ['| 相邻级 | ' + ' | '.join(n for _, n in keys) + ' |',
                  '|' + '---|' * (len(keys) + 1)]
    for i, d in enumerate(diffs):
        conv_lines.append(f'| L{i}-L{i+1} | ' +
                          ' | '.join(f'{d[k]:.3e}' for k, _ in keys) + ' |')
    orders = {k: np.log2(diffs[-2][k] / diffs[-1][k]) for k, _ in keys}
    rich = {k: diffs[-1][k] / (2 ** orders[k] - 1) for k, _ in keys}
    conv_lines.append('| 末两对观测阶 p | ' +
                      ' | '.join(f'{orders[k]:.2f}' for k, _ in keys) + ' |')
    conv_lines.append('| 最细解误差估计 | ' +
                      ' | '.join(f'{rich[k]:.2e}' for k, _ in keys) + ' |')

    ver = [
        '## 问题1 数值检验记录',
        '',
        f'网格：表面加密有限体积（加密带 {BAND*1e3:.1f} mm），五级细化节点数 '
        + '、'.join(str(s['n_nodes']) for s in sols) +
        f'；时间步长 ' + '、'.join(f'{lv[3]}' for lv in LEVELS) +
        f' s；首 {N_BE} 步后向 Euler 启动。生产解取最细一级。',
        '', *conv_lines, '',
        f'离散守恒残差（最细一级）：能量 {sol["ener_res"]:.2e}，水分 {sol["mass_res"]:.2e}（相对值）。',
        f'预估-校正迭代最大次数：{sol["picard_max"]}（收敛阈值 1e-6）。',
        '',
        f'特征数：Bi_热 = hR/k = {H*R0/KH:.4f}；'
        f'Bi_质(初值) = h_m R/D(C0) = {HM*R0/float(D_of_C(np.array([C0]))[0]):.4f}；'
        f'热扩散率 α = {KH/(RHO*CP):.4e} m²/s。',
    ]
    with open(os.path.join('结果', 'p1_检验数据.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(ver) + '\n')

    print('=== 问题1 求解完成（生产解=最细一级，节点 %d，Δt=%g s）===' %
          (sol['n_nodes'], LEVELS[-1][3]))
    print('\n'.join(conv_lines))
    print(f'守恒残差: 能量 {sol["ener_res"]:.1e}, 水分 {sol["mass_res"]:.1e}; '
          f'Picard 最大迭代 {sol["picard_max"]}')
    print(f'结果文件: {dst}')
    print('--- 表1 温度/℃ ---')
    print(tab1)
    print('--- 表2 水分浓度/(kg/kg) ---')
    print(tab2)


if __name__ == '__main__':
    main()
