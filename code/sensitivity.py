# -*- coding: utf-8 -*-
"""灵敏度分析：对流换热系数 h、对流传质系数 h_m、扩散系数前置系数 D0
各作 ±10% 单因素扰动，考察问题3（定域）与问题4（收缩，主模型B）
烘干时间 t_dry 的响应，计算归一化灵敏度系数
    S = (Δt_dry / t_dry,基准) / (Δp / p基准)   （即弹性，±10% 扰动下 S=Δ%/±10%）

网格采用粗一级（基准时长与细网格之差小于0.002 h，远小于主要参数扰动），
步长方案与生产解一致。

运行方式（项目根目录下）：python code/sensitivity.py
输出：结果/灵敏度分析.md、结果/sensitivity.npz
"""
import os
import time
import numpy as np

import p23
import p4

CASES = [('基准', 1.0, 1.0, 1.0),
         ('h+10%', 1.1, 1.0, 1.0), ('h-10%', 0.9, 1.0, 1.0),
         ('hm+10%', 1.0, 1.1, 1.0), ('hm-10%', 1.0, 0.9, 1.0),
         ('D0+10%', 1.0, 1.0, 1.1), ('D0-10%', 1.0, 1.0, 0.9)]

MESH3 = p23.MESHES[0]
MESH4 = p4.XI_MESHES[0]


def run_p3(fh, fhm, fD):
    p23.set_params(fh, fhm, fD)
    try:
        s = p23.solve_p3(*MESH3, p23.SCHED3_PROD)
    finally:
        p23.set_params(1.0, 1.0, 1.0)
    return s['t_dry'] / 3600.0


def run_p4(fh, fhm, fD):
    p4.set_params(fh, fhm, fD)
    try:
        s = p4.solve_p4(*MESH4, p4.SCHED4_PROD, model='B')
    finally:
        p4.set_params(1.0, 1.0, 1.0)
    return s['t_dry'] / 3600.0


def main():
    t3 = {}
    t4 = {}
    for name, fh, fhm, fD in CASES:
        tic = time.time()
        t3[name] = run_p3(fh, fhm, fD)
        t4[name] = run_p4(fh, fhm, fD)
        print(f'{name:8s}  P3 {t3[name]:8.4f} h   P4 {t4[name]:8.4f} h '
              f'  ({time.time()-tic:.0f} s)', flush=True)

    base3, base4 = t3['基准'], t4['基准']

    def row(name, dp):
        d3 = t3[name] - base3
        d4 = t4[name] - base4
        s3 = (d3 / base3) / dp
        s4 = (d4 / base4) / dp
        return (f'| {name} | {t3[name]:.4f} | {d3:+.4f} | {100*d3/base3:+.2f}% | '
                f'{s3:+.3f} | {t4[name]:.4f} | {d4:+.4f} | {100*d4/base4:+.2f}% | '
                f'{s4:+.3f} |')

    lines = [
        '## 灵敏度分析（h、h_m、D 前置系数 ±10% 单因素扰动）',
        '',
        f'响应量：烘干时间 t_dry。基准（粗网格）：问题3 {base3:.4f} h，'
        f'问题4 {base4:.4f} h（与生产解之差 <2e-3 h，见各问检验记录）。',
        '',
        '| 扰动 | P3 t_dry/h | Δ/h | Δ% | S | P4 t_dry/h | Δ/h | Δ% | S |',
        '|---|---|---|---|---|---|---|---|---|',
        row('h+10%', +0.1), row('h-10%', -0.1),
        row('hm+10%', +0.1), row('hm-10%', -0.1),
        row('D0+10%', +0.1), row('D0-10%', -0.1),
        '',
        'S 为归一化灵敏度系数（弹性）：S=-1 表示参数增 10% 则烘干时间缩短 10%。',
    ]
    os.makedirs('结果', exist_ok=True)
    with open(os.path.join('结果', '灵敏度分析.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    np.savez(os.path.join('结果', 'sensitivity.npz'),
             cases=np.array([c[0] for c in CASES]),
             t3=np.array([t3[c[0]] for c in CASES]),
             t4=np.array([t4[c[0]] for c in CASES]))
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
