# -*- coding: utf-8 -*-
"""Export the original Figure 8(b) data for figures/fig08b_style.tex."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "paper" / "figures"
data = np.load(ROOT / "结果" / "p4_field.npz")
radius = pd.read_excel(ROOT / "附件" / "附件2.xlsx").iloc[::6]

# Keep the saved control-volume weighted mean, as in code/figures.py.
curves = np.column_stack(
    (data["times"] / 3600, data["Rsamp"] * 100, data["meanC"])
)
measurements = np.column_stack(
    (radius.iloc[:, 0].to_numpy(float) / 3600, radius.iloc[:, 1].to_numpy(float))
)
np.savetxt(OUTPUT / "fig08b_curves.dat", curves, fmt="%.17g",
           header="time radius mean", comments="")
np.savetxt(OUTPUT / "fig08b_measurements.dat", measurements, fmt="%.17g",
           header="time radius", comments="")
end_hours = float(data["t_dry"]) / 3600
(OUTPUT / "fig08b_endpoint.tex").write_text(
    "\\def\\DryingEndHour{" + format(end_hours, ".17g") + "}\n",
    encoding="utf-8",
)
print(f"Exported {len(curves)} curve samples and {len(measurements)} measurements.")
print(f"Drying endpoint: {end_hours:.12f} h")
