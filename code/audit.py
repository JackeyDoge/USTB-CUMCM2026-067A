# -*- coding: utf-8 -*-
"""核对交付结果；可选重新求解，不覆盖既有结果工作簿和数值场。

python code/audit.py           # 工作簿、数值场及论文六张表一致性
python code/audit.py --solve   # 再运行四问生产网格及独立算法检验
"""
import argparse
import ast
import json
import os
import time
from pathlib import Path

import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '结果'


def workbook_check(problem):
    data = np.load(OUT / f'p{problem}_field.npz')
    fields = [('温度', 'T'), ('水分浓度', 'C')] if problem < 3 else [('Sheet1', 'C')]
    wb = openpyxl.load_workbook(OUT / f'result{problem}.xlsx', read_only=True)
    result = {'sheets': {}, 'checked_values': 0}
    for sheet, field in fields:
        expected = data[field]
        if problem == 4:
            expected = np.column_stack((expected, data['Csurf']))
        rows = wb[sheet].iter_rows()
        header = next(rows)
        assert np.allclose([c.value for c in header[1:22]], np.arange(21) / 10)
        if problem == 4:
            assert header[22].value == '药材表面'
        last_time = None
        n = 0
        for row in rows:
            if row[0].value is None:
                assert all(c.value is None for c in row)
                continue
            assert n < len(expected), f'extra data row in {sheet}'
            assert row[0].value == data['times'][n]
            last_time = row[0].value
            for cell, value in zip(row[1:], expected[n]):
                if np.isnan(value):
                    assert cell.value is None, (problem, cell.coordinate)
                else:
                    assert isinstance(cell.value, (int, float))
                    assert abs(cell.value - round(float(value), 4)) < 5e-13
                    assert cell.number_format == '0.0000', cell.coordinate
                    result['checked_values'] += 1
            assert all(c.value is None for c in row[1 + expected.shape[1]:])
            n += 1
        assert n == len(expected)
        result['sheets'][sheet] = {'rows': n, 'last_time_s': last_time}
    wb.close()
    return result


def paper_tables_check():
    """核对正文数字与原始求解表；末行小时标签允许提高显示精度。"""
    import re
    body = (ROOT / 'paper' / 'body.tex').read_text(encoding='utf-8')
    labels = ['p1T', 'p1C', 'p2T', 'p2C', 'p3', 'p4']
    expected_tables = []
    for name in ['p1', 'p2', 'p3', 'p4']:
        content = (OUT / f'{name}_论文表格.md').read_text(encoding='utf-8')
        for section in content.split('## ')[1:]:
            rows = []
            for line in section.splitlines():
                cells = [s.strip() for s in line.strip('|').split('|')]
                if line.startswith('| ') and re.match(r'^\d', cells[0]):
                    rows.append(cells[1:])
            if rows:
                expected_tables.append(rows)
    assert len(expected_tables) == len(labels)
    checks = {}
    for label, expected in zip(labels, expected_tables):
        part = body.split('\\label{tab:' + label + '}')[1].split('\\end{table}')[0]
        actual = []
        for line in part.splitlines():
            if re.match(r'^\d', line):
                cells = line.split('\\\\')[0].split('&')[1:]
                actual.append([s.strip().replace('\\textemdash', '—') for s in cells])
        assert actual == expected, label
        checks[label] = len(actual)
    return checks


def field_check(problem):
    """检查物理范围、固定物理位置域外掩码和未舍入终点判据。"""
    data = np.load(OUT / f'p{problem}_field.npz')
    checks = {}
    C = data['C']
    if problem == 4:
        outside = data['r_cm'][None, :] / 100 > data['Rsamp'][:, None] + 1e-12
        assert np.array_equal(np.isnan(C), outside)
        assert np.isfinite(data['Csurf']).all()
    else:
        assert np.isfinite(C).all()
    assert np.nanmin(C) >= 0 and np.nanmax(C) <= 2.55 + 1e-9
    checks['C_range'] = [float(np.nanmin(C)), float(np.nanmax(C))]
    if 'T' in data:
        assert np.isfinite(data['T']).all()
        assert data['T'].min() >= 28 - 1e-8 and data['T'].max() <= 51
        checks['T_range'] = [float(data['T'].min()), float(data['T'].max())]
    if problem >= 3:
        assert C[-1, 0] < 0.15
        assert np.all(C[:, 0] >= np.nanmax(C, axis=1) - 1e-12)
        checks['last_center_below_threshold'] = True
    return checks


def solve_again(problem):
    """以同一生产设置重新求解，与保存的未舍入场逐点比较。"""
    started = time.perf_counter()
    if problem == 1:
        import p1
        sol = p1.solve(*p1.LEVELS[-1])
    elif problem == 2:
        import p23
        sol = p23.solve_p2(*p23.MESHES[-1], p23.SCHED_PROD)
    elif problem == 3:
        import p23
        sol = p23.solve_p3(*p23.MESHES[-1], p23.SCHED3_PROD)
    else:
        import p4
        sol = p4.solve_p4(*p4.XI_MESHES[-1], p4.SCHED4_PROD)
    old = np.load(OUT / f'p{problem}_field.npz')
    diffs = {}
    for key in ['T', 'C', 'Csurf', 'meanC', 'Rsamp']:
        if key in old and key in sol:
            assert sol[key].shape == old[key].shape
            diffs[key] = float(np.nanmax(np.abs(sol[key] - old[key])))
            assert np.allclose(sol[key], old[key], atol=2e-9, rtol=0, equal_nan=True)
    record = {'max_abs_difference': diffs,
              'seconds': time.perf_counter() - started,
              'picard_max': sol['picard_max'],
              'mass_equation_residual': sol['mass_res']}
    if 't_dry' in sol:
        record['t_dry_h'] = sol['t_dry'] / 3600
        record['event_difference_s'] = sol['t_dry'] - float(old['t_dry'])
        record['center_is_max'] = sol['center_is_max']
        assert sol['center_is_max']
        assert abs(record['event_difference_s']) < 0.001
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solve', action='store_true')
    args = parser.parse_args()
    os.chdir(ROOT)
    target = OUT / 'audit_report.json'
    report = json.loads(target.read_text(encoding='utf-8')) if target.exists() else {}
    report['workbooks'], report['field_checks'] = {}, {}
    report.setdefault('recomputed', {})
    for source in (ROOT / 'code').glob('*.py'):
        ast.parse(source.read_text(encoding='utf-8'), filename=source.name)
    for p in range(1, 5):
        report['workbooks'][str(p)] = workbook_check(p)
        report['field_checks'][str(p)] = field_check(p)
    report['paper_tables'] = paper_tables_check()
    a, b = [np.load(OUT / f'p{i}_field.npz') for i in (3, 4)]
    report['drying'] = {
        'p3_h': float(a['t_dry']) / 3600,
        'p4_h': float(b['t_dry']) / 3600,
        'difference_h': float(a['t_dry'] - b['t_dry']) / 3600,
        'p3_last_center': float(a['C'][-1, 0]),
        'p4_last_center': float(b['C'][-1, 0])}

    def save():
        target.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')

    save()
    print('工作簿与六张论文表核对通过。', flush=True)
    if args.solve:
        import verify
        report['independent'] = {
            'temperature_semianalytic': verify.check_temperature(),
            'tridiagonal_difference': verify.check_solvers(),
            'manufactured_solution': verify.mms_mass()}
        save()
        for p in range(1, 5):
            print(f'重新求解问题{p}...', flush=True)
            report['recomputed'][str(p)] = solve_again(p)
            save()
            print(json.dumps(report['recomputed'][str(p)], ensure_ascii=False), flush=True)
    print('核对记录：结果/audit_report.json', flush=True)


if __name__ == '__main__':
    main()
