# -*- coding: utf-8 -*-
"""Export the original environment records for figures/fig01_env_style.tex."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "paper" / "figures"
records = pd.read_excel(ROOT / "附件" / "附件1.xlsx")
curves = records.iloc[:, :3].to_numpy(float)
curves[:, 0] /= 3600

# Straight segments between all recorded points retain the original linear
# interpolation. The constant extension is drawn separately in the TeX source.
np.savetxt(OUTPUT / "fig01_env_records.dat", curves, fmt="%.17g",
           header="time temperature concentration", comments="")
temperature = curves[:, 1]
concentration = curves[:, 2]
limits = {}
for name, values in (("Temperature", temperature), ("Concentration", concentration)):
    low, high = float(values.min()), float(values.max())
    margin = 0.05 * (high - low)
    limits[name + "Minimum"] = low - margin
    limits[name + "Maximum"] = high + margin
(OUTPUT / "fig01_env_limits.tex").write_text(
    "".join("\\def\\" + name + "{" + format(value, ".17g") + "}\n"
            for name, value in limits.items()),
    encoding="utf-8",
)
print(f"Exported {len(curves)} original environment records (0–4 h).")
