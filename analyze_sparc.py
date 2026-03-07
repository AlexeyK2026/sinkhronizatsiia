import numpy as np
import pandas as pd
from astroquery.vizier import Vizier

Pi = 1.2e-10  # м/с^2

# Забираем таблицы SPARC из VizieR (каталог J/AJ/152/157)
Vizier.ROW_LIMIT = -1

t1 = Vizier.get_catalogs("J/AJ/152/157/table1")[0].to_pandas()
t2 = Vizier.get_catalogs("J/AJ/152/157/table2")[0].to_pandas()

# Чистим и оставляем то, что нужно
t1 = t1[["Name","RHI","Vflat","Qual"]].copy()
t2 = t2[["Name","Rad","Vobs"]].copy()

# Берем галактики, где есть Vflat и RHI, и лучшее качество (Qual=1)
t1 = t1.replace({np.nan: None})
t1 = t1[(t1["Vflat"] > 0) & (t1["RHI"] > 0) & (t1["Qual"] == 1)].copy()

# Для "ядра": берем точку с минимальным Rad
t2 = t2[(t2["Rad"] > 0) & (t2["Vobs"] > 0)].copy()
core = (t2.sort_values(["Name","Rad"])
          .groupby("Name", as_index=False)
          .first()
          .rename(columns={"Rad":"Rad_core_kpc","Vobs":"C_core_kms"}))

df = t1.merge(core, on="Name", how="inner")

# SI-конверсии
kpc_to_m = 3.085677581e19
kms_to_ms = 1000.0

C_ms = df["Vflat"].to_numpy() * kms_to_ms
Ccore_ms = df["C_core_kms"].to_numpy() * kms_to_ms
P_m = df["RHI"].to_numpy() * kpc_to_m

# "Эффективный k", который требуется формуле для совпадения C
k_eff = C_ms / np.sqrt(Ccore_ms * Pi * P_m)

out = pd.DataFrame({
    "Name": df["Name"],
    "Vflat_km_s": df["Vflat"],
    "RHI_kpc": df["RHI"],
    "Rad_core_kpc": df["Rad_core_kpc"],
    "C_core_km_s": df["C_core_kms"],
    "k_eff": k_eff
})

out.to_csv("sparc_k_eff_check.csv", index=False)
print("Saved sparc_k_eff_check.csv with", len(out), "galaxies")
