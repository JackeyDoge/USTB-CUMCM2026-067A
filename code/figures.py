# -*- coding: utf-8 -*-
"""论文全套图件生成。

数据来源：结果/p1_field.npz、p2_field.npz、p3_field.npz、p4_field.npz、
sensitivity.npz，附件1/附件2 原始数据；网格收敛数字取自
结果/p1_检验数据.md（p1.py/verify.py 运行输出，此处按记录转录）。

设计约定（经调色板六项校验，Okabe-Ito 无障碍色系）：
  半径实体固定用色：r=0 蓝 #0072B2，r=1 绿 #009E73，r=1.5 橙 #E69F00，
  r=2(表面) 朱红 #D55E00；环境/参考线灰 #666；平均量天蓝 #56B4E9。
  数量场色带：温度 Oranges、水分 Blues（单色相、深=大）；不用彩虹色带；
  不使用双 y 轴（改用上下分栏）；网格线细淡置底。

运行方式（项目根目录下）：python code/figures.py
输出：paper/figures/*.pdf（矢量，LaTeX 引用）与 结果/图预览/*.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec

# ---- 全局样式 ----
plt.rcParams.update({
    'font.sans-serif': ['SimHei', 'Microsoft YaHei'],
    'axes.unicode_minus': False,
    'mathtext.fontset': 'dejavusans',
    'font.size': 10.5,
    'axes.linewidth': 0.8,
    'axes.grid': True,
    'grid.color': '#d9d9d9',
    'grid.linewidth': 0.5,
    'axes.axisbelow': True,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '#cccccc',
    'figure.dpi': 110,
})
BLUE, GREEN, ORANGE, VERMIL, SKY, GRAY = ('#0072B2', '#009E73', '#E69F00',
                                          '#D55E00', '#56B4E9', '#666666')
R_COLOR = {0.0: BLUE, 1.0: GREEN, 1.5: ORANGE, 2.0: VERMIL}
PDF_DIR = os.path.join('paper', 'figures')
PNG_DIR = os.path.join('结果', '图预览')


def savefig(fig, name):
    fig.savefig(os.path.join(PDF_DIR, name + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(PNG_DIR, name + '.png'), bbox_inches='tight', dpi=200)
    plt.close(fig)
    print('  已生成', name)


def fig_env():
    """图1 烘房环境条件（附件1）与两阶段划分。"""
    df = pd.read_excel(os.path.join('附件', '附件1.xlsx'))
    t = df.iloc[:, 0].to_numpy(float) / 3600
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 4.6), sharex=True)
    ext = np.array([4.0, 4.8])
    a1.plot(t, df.iloc[:, 1], color=VERMIL, lw=1.4, label='附件1实测')
    a1.plot(ext, [50, 50], color=VERMIL, lw=1.4, ls='--', label='恒温段设计值')
    a1.set_ylabel('烘房温度/℃')
    a1.legend(loc='lower right', fontsize=9)
    a2.plot(t, df.iloc[:, 2], color=BLUE, lw=1.4)
    a2.plot(ext, [0.05, 0.05], color=BLUE, lw=1.4, ls='--')
    a2.set_ylabel('烘房水分浓度/(kg/kg)')
    a2.set_xlabel('时间/h')
    for ax in (a1, a2):
        ax.axvline(4.0, color=GRAY, lw=0.9, ls=':')
        ax.set_xlim(0, 4.8)
    a1.text(1.9, 30.5, '预热平衡阶段', ha='center', fontsize=9.5, color=GRAY)
    a1.text(4.4, 44.5, '恒温干燥阶段', ha='center', fontsize=9.5, color=GRAY)
    savefig(fig, '图3')


def fig_p1_T():
    """图2 问题1温度场：云图 + 特征位置历程。"""
    d = np.load(os.path.join('结果', 'p1_field.npz'))
    t, r, T = d['times'], d['r_cm'], d['T']
    fig = plt.figure(figsize=(10.4, 3.7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1], wspace=0.28)
    ax = fig.add_subplot(gs[0])
    pc = ax.pcolormesh(t, r, T.T, cmap='Oranges', shading='gouraud')
    fig.colorbar(pc, ax=ax, label='温度/℃', pad=0.02)
    ax.set_xlabel('时间/s')
    ax.set_ylabel('到药材中心的距离/cm')
    ax.set_title('(a) 温度分布 T(r,t)', loc='left', fontsize=10.5)
    ax.grid(False)
    ax2 = fig.add_subplot(gs[1])
    env = pd.read_excel(os.path.join('附件', '附件1.xlsx'))
    ax2.plot(env.iloc[:31, 0], env.iloc[:31, 1], color=GRAY, lw=1.2, ls='--',
             label='烘房 $T_a$')
    for rv in (0.0, 1.0, 1.5, 2.0):
        j = int(round(rv * 10))
        ax2.plot(t, T[:, j], color=R_COLOR[rv], lw=1.8, label=f'r={rv:g} cm')
    ax2.set_xlabel('时间/s')
    ax2.set_ylabel('温度/℃')
    ax2.set_xlim(0, 1800)
    ax2.legend(fontsize=9, loc='upper left')
    ax2.set_title('(b) 特征位置温度历程', loc='left', fontsize=10.5)
    savefig(fig, 'fig02_p1_T')


def fig_p1_C():
    """图3 问题1水分场：云图 + 径向剖面（表面边界层）。"""
    d = np.load(os.path.join('结果', 'p1_field.npz'))
    t, r, C = d['times'], d['r_cm'], d['C']
    fig = plt.figure(figsize=(10.4, 3.7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1], wspace=0.28)
    ax = fig.add_subplot(gs[0])
    pc = ax.pcolormesh(t, r, C.T, cmap='Blues', shading='gouraud')
    fig.colorbar(pc, ax=ax, label='水分浓度/(kg/kg)', pad=0.02)
    ax.set_xlabel('时间/s')
    ax.set_ylabel('到药材中心的距离/cm')
    ax.set_title('(a) 水分浓度分布 C(r,t)', loc='left', fontsize=10.5)
    ax.grid(False)
    ax2 = fig.add_subplot(gs[1])
    times = (100, 300, 600, 900, 1200, 1500, 1800)
    shades = plt.cm.Blues(np.linspace(0.45, 1.0, len(times)))
    for tv, col in zip(times, shades):
        ax2.plot(r, C[tv - 1], color=col, lw=1.7)
        ax2.annotate(f'{tv}s', xy=(2.005, C[tv - 1, -1]), fontsize=8,
                     color='#333333', va='center')
    ax2.set_xlabel('到药材中心的距离/cm')
    ax2.set_ylabel('水分浓度/(kg/kg)')
    ax2.set_xlim(0, 2.25)
    ax2.set_title('(b) 各时刻径向剖面（表面薄边界层）', loc='left', fontsize=10.5)
    savefig(fig, 'fig03_p1_C')


def fig_convergence():
    """图4 数值检验：网格收敛（P1 场量、MMS 制造解）log-log 图。
    数字转录自 结果/p1_检验数据.md。"""
    lev = np.array([1, 2, 4, 8])          # 相邻级加密倍数（L0-L1 起）
    dT = np.array([3.149e-4, 7.879e-5, 1.970e-5, 4.924e-6])
    dC = np.array([3.938e-3, 1.004e-3, 2.519e-4, 6.304e-5])
    lev_m = np.array([1, 2])
    e_mms = np.array([3.328e-5, 8.346e-6, 2.088e-6])
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.loglog(lev, dT, 'o-', color=VERMIL, lw=1.8, ms=6, label='问题1 温度场相邻级差')
    ax.loglog(lev, dC, 's-', color=BLUE, lw=1.8, ms=6, label='问题1 水分场相邻级差')
    ax.loglog([1, 2, 4], e_mms, '^-', color=GREEN, lw=1.8, ms=6,
              label='制造解(MMS)误差')
    ref = 6e-3 * (1.0 / lev.astype(float)) ** 2
    ax.loglog(lev, ref, ls='--', color='#999999', lw=1.2, label='二阶参考斜率')
    ax.set_xlabel('网格/步长加密倍数')
    ax.set_ylabel('最大绝对差或误差')
    # 显式设置刻度标签（避免 \mathdefault 走中文字体缺失上标负号）
    from matplotlib.ticker import NullFormatter
    ax.set_xticks([1, 2, 4, 8], ['1', '2', '4', '8'])
    ax.set_yticks([1e-3, 1e-4, 1e-5],
                  [r'$10^{-3}$', r'$10^{-4}$', r'$10^{-5}$'])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.legend(fontsize=9)
    savefig(fig, 'fig04_convergence')


def fig_p2():
    """图5 问题2前3小时：温度与水分特征位置历程。"""
    d = np.load(os.path.join('结果', 'p2_field.npz'))
    th = d['times'] / 3600
    T, C = d['T'], d['C']
    env = pd.read_excel(os.path.join('附件', '附件1.xlsx'))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 3.7))
    a1.plot(env.iloc[:181, 0] / 3600, env.iloc[:181, 1], color=GRAY, lw=1.2,
            ls='--', label='烘房 $T_a$')
    for rv in (0.0, 1.0, 2.0):
        j = int(round(rv * 10))
        a1.plot(th, T[:, j], color=R_COLOR[rv], lw=1.8, label=f'r={rv:g} cm')
        a2.plot(th, C[:, j], color=R_COLOR[rv], lw=1.8, label=f'r={rv:g} cm')
    a1.set_xlabel('时间/h')
    a1.set_ylabel('温度/℃')
    a1.set_xlim(0, 3)
    a1.legend(fontsize=9, loc='lower right')
    a1.set_title('(a) 温度：约3 h 完成升温并趋于50 ℃', loc='left', fontsize=10.5)
    a2.set_xlabel('时间/h')
    a2.set_ylabel('水分浓度/(kg/kg)')
    a2.set_xlim(0, 3)
    a2.legend(fontsize=9, loc='lower left')
    a2.set_title('(b) 水分：由表及里进入整体下降', loc='left', fontsize=10.5)
    savefig(fig, 'fig05_p2')


def fig_p3():
    """图6 问题3干燥全程：干燥曲线 + 各时刻剖面。"""
    d = np.load(os.path.join('结果', 'p3_field.npz'))
    th = d['times'] / 3600
    C = d['C']
    t_dry = float(d['t_dry']) / 3600
    r = d['r_cm']
    r_m = r / 100
    meanC = np.trapezoid(C * r_m, r_m, axis=1) / np.trapezoid(r_m, r_m)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 3.7))
    a1.plot(th, C[:, 0], color=BLUE, lw=1.8, label='中心 r=0')
    a1.plot(th, meanC, color=SKY, lw=1.8, label='截面平均')
    a1.plot(th, C[:, -1], color=VERMIL, lw=1.8, label='表面 r=2 cm')
    a1.axhline(0.15, color=GRAY, lw=1.1, ls='--')
    a1.text(0.8, 0.165, '烘干判据 0.15', fontsize=9, color=GRAY, va='bottom')
    a1.axvline(t_dry, color=GRAY, lw=1.1, ls=':')
    a1.annotate(f'烘干结束\n{t_dry:.2f} h', xy=(t_dry, 1.15), fontsize=9,
                color='#333333', ha='right', xytext=(t_dry - 1.5, 1.15))
    a1.set_xlabel('时间/h')
    a1.set_ylabel('水分浓度/(kg/kg)')
    a1.set_xlim(0, 60)
    a1.legend(fontsize=9)
    a1.set_title('(a) 干燥曲线与烘干判据', loc='left', fontsize=10.5)
    times_h = (6, 12, 24, 36, 48, 57)
    shades = plt.cm.Blues(np.linspace(0.45, 1.0, len(times_h)))
    for tv, col in zip(times_h, shades):
        a2.plot(r, C[int(round(tv * 60)) - 1], color=col, lw=1.7,
                label=f'{tv} h')
    a2.axhline(0.15, color=GRAY, lw=1.1, ls='--')
    a2.set_xlabel('到药材中心的距离/cm')
    a2.set_ylabel('水分浓度/(kg/kg)')
    a2.set_xlim(0, 2.05)
    a2.legend(fontsize=8.5, ncols=2, title='时刻', title_fontsize=8.5)
    a2.set_title('(b) 各时刻径向剖面', loc='left', fontsize=10.5)
    savefig(fig, 'fig06_p3')


def fig_p4():
    """图7 问题4收缩烘干：动边界水分云图 + 半径/平均含水率演化。"""
    d = np.load(os.path.join('结果', 'p4_field.npz'))
    th = d['times'] / 3600
    C = d['C']
    R = d['Rsamp'] * 100
    t_dry = float(d['t_dry']) / 3600
    rad = pd.read_excel(os.path.join('附件', '附件2.xlsx'))
    fig = plt.figure(figsize=(10.4, 3.9))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.25, 1], wspace=0.26)
    ax = fig.add_subplot(gs[0])
    Cm = np.ma.masked_invalid(C)
    pc = ax.pcolormesh(th, d['r_cm'], Cm.T, cmap='Blues', shading='auto')
    ax.plot(th, R, color='#333333', lw=1.6)
    ax.text(20, 1.32, '药材表面 R(t)', fontsize=9, color='#333333')
    fig.colorbar(pc, ax=ax, label='水分浓度/(kg/kg)', pad=0.02)
    ax.axvline(t_dry, color=GRAY, lw=1.1, ls=':')
    ax.annotate(f'烘干结束 {t_dry:.2f} h', xy=(t_dry, 1.9), fontsize=9,
                ha='right', xytext=(t_dry - 1.2, 1.9), color='#333333')
    ax.set_xlabel('时间/h')
    ax.set_ylabel('到药材中心的距离/cm')
    ax.set_title('(a) 收缩域内水分分布与动边界', loc='left', fontsize=10.5)
    ax.grid(False)
    gs_r = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.12)
    b1 = fig.add_subplot(gs_r[0])
    b2 = fig.add_subplot(gs_r[1], sharex=b1)
    b1.plot(rad.iloc[::6, 0] / 3600, rad.iloc[::6, 1], 'o', ms=3.2,
            color='#333333', mfc='none', label='附件2实测')
    b1.plot(th, R, color=GREEN, lw=1.8, label='插值 R(t)')
    b1.axvline(t_dry, color=GRAY, lw=1.1, ls=':')
    b1.set_ylabel('半径/cm')
    b1.legend(fontsize=8.5, loc='upper right')
    b1.tick_params(labelbottom=False)
    b1.set_title('(b) 半径收缩与平均含水率', loc='left', fontsize=10.5)
    b2.plot(th, d['meanC'], color=SKY, lw=1.8)
    b2.axvline(t_dry, color=GRAY, lw=1.1, ls=':')
    b2.set_ylabel('平均含水率/(kg/kg)')
    b2.set_xlabel('时间/h')
    b2.set_xlim(0, 72)
    savefig(fig, 'fig07_p4')


def fig_compare():
    """问题3与问题4中心含水率对比，采用宋体和Times New Roman。"""
    from matplotlib.transforms import Bbox

    d3 = np.load(os.path.join('结果', 'p3_field.npz'))
    d4 = np.load(os.path.join('结果', 'p4_field.npz'))
    t3, t4 = float(d3['t_dry']) / 3600, float(d4['t_dry']) / 3600
    with plt.rc_context({'font.family': ['Times New Roman', 'SimSun'],
                         'pdf.fonttype': 42}):
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        ax.plot(d3['times'] / 3600, d3['C'][:, 0], color=BLUE, lw=1.8,
                label='问题3：尺寸固定')
        ax.plot(d4['times'] / 3600, d4['C'][:, 0], color=GREEN, lw=1.8,
                label='问题4：考虑收缩')
        ax.axhline(0.15, color=GRAY, lw=1.1, ls='--')
        ax.text(1.0, 0.17, '判据 0.15', fontsize=9, color=GRAY)
        for tv, col in ((t3, BLUE), (t4, GREEN)):
            ax.axvline(tv, color=col, lw=1.0, ls=':')
        ax.annotate(f'{t4:.2f} h', xy=(t4, 0.62), fontsize=9, ha='right',
                    xytext=(t4 - 1.0, 0.62), color=GREEN)
        ax.annotate(f'{t3:.2f} h', xy=(t3, 0.42), fontsize=9, ha='right',
                    xytext=(t3 - 1.0, 0.42), color=BLUE)
        ax.set_xlabel('时间(h)')
        ax.set_ylabel('中心水分浓度(kg/kg)')
        ax.set_xlim(0, 62)
        ax.legend(loc='upper center', fontsize=9.5, frameon=True,
                  facecolor='white', edgecolor='none', framealpha=1,
                  fancybox=False)
        # Keep the original PDF canvas and axis positions, in inches.
        bounds = Bbox.from_bounds(13.24375 / 72, -5.645 / 72,
                                 419.855507 / 72, 266.285 / 72)
        for folder, ext in ((PDF_DIR, '.pdf'), (PNG_DIR, '.png')):
            fig.savefig(os.path.join(folder, '图11' + ext),
                        bbox_inches=bounds, dpi=200)
        plt.close(fig)
        print('  已生成 图11')


def fig_sensitivity():
    """图9 灵敏度龙卷风图：±10%扰动下烘干时间的相对变化。"""
    d = np.load(os.path.join('结果', 'sensitivity.npz'), allow_pickle=True)
    cases = list(d['cases'])
    t3, t4 = d['t3'], d['t4']
    b3, b4 = t3[cases.index('基准')], t4[cases.index('基准')]

    def dp(tarr, base, name):
        return 100 * (tarr[cases.index(name)] - base) / base

    params = [('h', BLUE), ('hm', GREEN), ('D0', VERMIL)]
    labels = ['对流换热系数 h', '对流传质系数 $h_m$', '扩散系数前置系数 $D_0$']
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.4), sharey=True)
    for ax, tarr, base, ttl in ((axes[0], t3, b3, '(a) 问题3（尺寸固定）'),
                                (axes[1], t4, b4, '(b) 问题4（考虑收缩）')):
        y = np.arange(len(params))[::-1]
        for (pname, col), yy in zip(params, y):
            v_plus = dp(tarr, base, pname + '+10%')
            v_minus = dp(tarr, base, pname + '-10%')
            ax.barh(yy + 0.18, v_plus, height=0.32, color=col, alpha=0.95)
            ax.barh(yy - 0.18, v_minus, height=0.32, color=col, alpha=0.45)
            for v, yo in ((v_plus, yy + 0.18), (v_minus, yy - 0.18)):
                off = 0.35 if v >= 0 else -0.35
                ax.text(v + off, yo, f'{v:+.2f}%', va='center', fontsize=8,
                        ha='left' if v >= 0 else 'right', color='#333333')
        ax.axvline(0, color='#333333', lw=0.9)
        ax.set_yticks(y, labels)
        ax.set_xlim(-13, 13)
        ax.set_xlabel('烘干时间相对变化/%')
        ax.set_title(ttl, loc='left', fontsize=10.5)
    axes[1].text(6.2, 2.62, '深色: 参数+10%', fontsize=8.5, color='#333333')
    axes[1].text(6.2, 2.30, '浅色: 参数-10%', fontsize=8.5, color='#333333')
    savefig(fig, 'fig09_sensitivity')


def main():
    os.makedirs(PDF_DIR, exist_ok=True)
    os.makedirs(PNG_DIR, exist_ok=True)
    fig_env()
    fig_p1_T()
    fig_p1_C()
    fig_convergence()
    fig_p2()
    fig_p3()
    fig_p4()
    fig_compare()
    fig_sensitivity()
    print('全部图件已输出至', PDF_DIR, '与', PNG_DIR)


if __name__ == '__main__':
    main()
