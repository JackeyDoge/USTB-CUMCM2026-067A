# 67A论文与支撑材料

当前论文入口为 `paper/67A论文.tex`，正文入口为 `paper/论文正文.tex`，附录入口为 `paper/论文附录.tex`，各节源文件位于 `paper/sec/`。请从当前入口编译，保留最新章节和图件。


## 本次归档（2026年9月30日）

- 论文：[最新版 PDF](paper/67A论文.pdf) · [TeX 入口](paper/67A论文.tex)。题名《基于热质耦合与收缩动边界的药材烘干过程预测》。
- PDF 共73页，原始编译时间为2026年9月13日19:20:09（北京时间），对应最终正文和附录；本次归档保留其原始字节，不重新计算或改写论文。
- 最新 PDF 的 SHA-256：`ef47999abd3d585ec8c5c93a4bd98331ab055429312cd9dc1f6b46e9a9b8d183`。
- 包含10个求解、核验、绘图和数据导出程序，4份结果工作簿、未舍入数值场、输入附件、图件及AI工具使用说明。
- [本次只读核查报告](verification/result_validation.json)：4份工作簿共645,916个数值、保存场、原稿六张表及10个程序语法检查通过。此核查不重新求解，也不代表当前分节论文的每张表均已逐表核验。
- [文件校验清单](MANIFEST.sha256.json)列出归档文件的大小和SHA-256，可用于核对下载内容；清单自身不递归列入清单。
- 未纳入旧交付PDF、旧压缩包、临时文件、第三方库/字体安装本体、优秀论文样本及未被正式论文采用的潜热独立试验。

## 文件与目录

- `paper/67A论文.tex`、`paper/67A论文.pdf`：当前论文源码入口与编译结果。
- `paper/ai-usage.tex`、`paper/AI工具使用详情.pdf`：AI使用详情源码与支撑材料PDF。
- `paper/figures/`：按当前论文图号命名的图件，以及绘图源文件和数据。图6由章节内TikZ直接绘制。
- `code/`：附录列明的10个Python程序。
- `结果/`：结果工作簿、未舍入数值场、检验和灵敏度分析记录。
- `附件/`：题目给定环境、半径记录及结果模板。
- `paper/body.tex`：保留供 `audit.py` 核对六张原始结果表的历史正文；它不是当前编译入口。旧版 main.tex、main.pdf 和 appendix.tex 未纳入本仓库。

仓库按项目相对路径组织，下载或克隆后在 `paper/` 目录编译。仓库包含当前论文、AI使用详情、源码、图件、输入附件、计算结果及本说明。

## 编译论文与AI使用详情

使用支持中文的完整TeX Live或兼容XeLaTeX环境，配置SimSun、SimHei、Consolas及Times New Roman等当前源码所用字体。在 `paper/`目录执行：

```text
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error 67A论文.tex
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error 67A论文.tex
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error -jobname=AI工具使用详情 ai-usage.tex
```

正文及附录中的相对路径依赖 `paper/`与 `code/`的同级目录关系。直接编译可使用已保存的图件，无需重新运行求解器或绘图程序。论文标题不随文件重命名改变。

## 输入与计算复现

Python 3.12。保持官方附件为 `附件/附件1.xlsx`、`附件/附件2.xlsx`、`附件/附件3/result1.xlsx`至`result4.xlsx`。从项目根目录运行：

```text
python -m pip install -r requirements.txt
python code/p1.py
python code/verify.py
python code/p23.py
python code/p23.py p3
python code/p4.py
python code/sensitivity.py
python code/figures.py
python code/audit.py --solve
```

求解脚本会更新 `结果/`下对应输出，保留官方模板。`audit.py --solve`复算并比较既有数据，不覆盖工作簿或采样场，检查报告写入 `结果/audit_report.json`；不带该选项时核对工作簿、原稿结果表与保存场，但仍会更新审核报告。`verify.py` 会追加检验记录；如需保留归档原始结果，请在副本中运行复现命令。网格、时间步和事件定位的具体设置以当前程序和论文附录为准。

## 当前图件的重绘

`code/figures.py`包含计算图件绘制及当前图11的样式设置；当前外部PNG不保证由运行该脚本自动复原。应保留包内已采用的图件。

重绘当前图3时，先从项目根目录运行 `python code/export_figure2_data.py`，再在 `paper/figures/`运行：

```text
xelatex -interaction=nonstopmode -halt-on-error -jobname=图3 fig01_env_style.tex
```

如先运行过 `figures.py`，随后执行上述步骤以恢复当前图3的字体、标签与图例。

当前图10(b)的辅助程序保留旧文件名 `code/export_figure8b_data.py`。运行后编译 `paper/figures/fig08b_style.tex`，按原比例输出为 `图10(b).png`，像素尺寸为2391×2298。具体步骤也见论文附录A。

## 结果与审核记录

固定尺寸与收缩工况的临界时长分别为57.4721 h、51.0905 h。问题一、二分别逐秒记录至1800 s、10800 s；问题三、四分别逐分钟记录至206940 s、183960 s，末个采样时刻晚于插值临界时刻。工作簿保留四位小数，达标判据使用未舍入结果；问题四域外位置留空，另列实际表面。

两问同时改变物性和几何，时长差不能全归于收缩；数值检验与独立复算也不等同于实测预测精度验证。

参赛队于2026年9月13日确认已完成全部代码和结果的人工审核及独立复算，结果可信。工具分工、AI参与范围和具体可查修改记录见 `paper/AI工具使用详情.pdf`。本次文件重命名与重新编译不重新计算四问数值结果。
