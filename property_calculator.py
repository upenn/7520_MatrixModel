"""
property_calculator.py -- material properties and layer values for the thermal matrix model.

A materials table has one row per material, laid out like Table 1 of the ASHRAE Handbook--
Fundamentals, chapter 26 (SI): density, conductivity k (or, for products sold at one
thickness, a resistance R at a stated thickness), and specific heat. The built-in table
holds generic typical values; students add rows from the handbook or product data in their
own CSV (see TEMPLATE_COLUMNS), which is merged on top.

For a layer of a given area and thickness it gives the values the network needs
(conduction UA, capacitance C), split into segments for thick massive layers, plus formulas
ready to paste into the network editor.
"""
import difflib
import io
import math
import os
import urllib.request

import numpy as np
import pandas as pd

CALC_VERSION = "dev (not built)"

# Columns of a materials CSV (SI, as in ASHRAE Fundamentals ch. 26, Table 1).
# Give k, or R with the thickness it applies to (k is then thickness / R).
TEMPLATE_COLUMNS = ["material", "density_kg_m3", "k_W_mK", "R_m2K_W", "R_thickness_mm",
                    "cp_kJ_kgK", "source"]

# Generic typical values (same as the thermal matrix model's MATERIALS library).
GENERIC = [
    ("concrete", 2200, 1.75, None, None, 0.92),
    ("brick", 1920, 0.72, None, None, 0.84),
    ("adobe", 1600, 0.52, None, None, 0.90),
    ("rammed earth", 2000, 1.10, None, None, 0.90),
    ("stone", 2600, 2.30, None, None, 0.84),
    ("water", 1000, 0.60, None, None, 4.186),
    ("glass", 2500, 0.96, None, None, 0.84),
    ("gypsum board", 800, 0.16, None, None, 1.09),
    ("plaster", 1600, 0.70, None, None, 1.00),
    ("softwood", 510, 0.12, None, None, 1.38),
    ("hardwood", 720, 0.16, None, None, 1.26),
    ("CLT", 480, 0.13, None, None, 1.60),
    ("mineral wool", 30, 0.040, None, None, 0.84),
    ("EPS", 20, 0.035, None, None, 1.45),
    ("XPS", 35, 0.030, None, None, 1.45),
    ("straw bale", 110, 0.060, None, None, 2.00),
    ("steel", 7850, 50.0, None, None, 0.45),
    ("soil (moist)", 1800, 1.50, None, None, 1.30),
    ("air", 1.2, 0.026, None, None, 1.005),
]

# SI -> IP factors
K_IP = 0.577789        # W/m-K -> Btu/hr-ft-F
RHO_IP = 0.0624280     # kg/m3 -> lb/ft3
CP_IP = 1 / 4.1868     # kJ/kg-K -> Btu/lb-F
R_IP = 5.678263        # m2-K/W -> hr-ft2-F/Btu
M_PER_FT = 0.3048
UNITS = {
    "SI": dict(L="m", A="m²", k="W/m·K", rho="kg/m³", cp="J/kg·K", R="m²·K/W", U="W/m²·K",
               UA="W/K", C="J/K", Ca="kJ/m²·K"),
    "IP": dict(L="ft", A="ft²", k="Btu/hr·ft·°F", rho="lb/ft³", cp="Btu/lb·°F", R="hr·ft²·°F/Btu",
               U="Btu/hr·ft²·°F", UA="Btu/hr·°F", C="Btu/°F", Ca="Btu/ft²·°F"),
}
DAY = 86400.0


def generic_table():
    df = pd.DataFrame(GENERIC, columns=TEMPLATE_COLUMNS[:-1])
    df["source"] = "generic typical value"
    return df


def _clean(df):
    """Normalise a materials table: names, numbers, k from R where needed."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    missing = [c for c in ("material", "density_kg_m3", "cp_kJ_kgK") if c not in df]
    if missing:
        raise ValueError(f"materials table is missing column(s) {missing}; "
                         f"expected {', '.join(TEMPLATE_COLUMNS)}")
    for c in TEMPLATE_COLUMNS:
        if c not in df:
            df[c] = np.nan if c != "source" else ""
    df["material"] = df["material"].astype(str).str.strip()
    df = df[df["material"].ne("") & df["material"].ne("nan")]
    for c in ("density_kg_m3", "k_W_mK", "R_m2K_W", "R_thickness_mm", "cp_kJ_kgK"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    from_R = df["k_W_mK"].isna() & df["R_m2K_W"].gt(0) & df["R_thickness_mm"].gt(0)
    df.loc[from_R, "k_W_mK"] = df.loc[from_R, "R_thickness_mm"] / 1000 / df.loc[from_R, "R_m2K_W"]
    df["k_from_R"] = from_R
    df["source"] = df["source"].fillna("").astype(str)
    return df[TEMPLATE_COLUMNS + ["k_from_R"]].reset_index(drop=True)


def read_table(source):
    """Read a materials CSV from a path, URL, bytes or text."""
    if isinstance(source, bytes):
        source = source.decode("utf-8-sig")
    if isinstance(source, str) and source.startswith(("http://", "https://")):
        if "docs.google.com/spreadsheets" in source and "/edit" in source:   # Google Sheet link
            source = source.split("/edit")[0] + "/export?format=csv"
        with urllib.request.urlopen(source, timeout=60) as r:
            source = r.read().decode("utf-8-sig")
    if isinstance(source, str) and not os.path.exists(source):
        return _clean(pd.read_csv(io.StringIO(source)))
    return _clean(pd.read_csv(source, encoding="utf-8-sig"))


class Materials:
    """A materials table: generic values with the user's rows added (user rows win)."""

    def __init__(self, user=None):
        tables = [_clean(generic_table())]
        if user is not None:
            tables.append(user if "k_from_R" in user else _clean(user))
        df = pd.concat(tables, ignore_index=True)
        df["key"] = df["material"].str.lower()
        self.df = df.drop_duplicates("key", keep="last").drop(columns="key").reset_index(drop=True)
        self.n_user = 0 if user is None else len(user)

    def names(self):
        return list(self.df["material"])

    def find(self, name):
        """Row for a material name (case-insensitive; suggests close matches)."""
        key = str(name).strip().lower()
        hit = self.df[self.df["material"].str.lower() == key]
        if len(hit):
            return hit.iloc[0]
        part = self.df[self.df["material"].str.lower().str.contains(key, regex=False)]
        if len(part) == 1:
            return part.iloc[0]
        close = list(part["material"]) or difflib.get_close_matches(name, self.names(), n=8, cutoff=0.4)
        raise KeyError(f"'{name}' is not in the materials table." +
                       (f" Did you mean: {'; '.join(close[:12])}?" if close else ""))

    def search(self, text=""):
        t = str(text).strip().lower()
        df = self.df if not t else self.df[self.df["material"].str.lower().str.contains(t, regex=False)]
        return df

    def show(self, text="", units="SI"):
        df = self.search(text)
        out = pd.DataFrame({"material": df["material"]})
        k, rho, cp = props_in(df["k_W_mK"], df["density_kg_m3"], df["cp_kJ_kgK"], units)
        u = UNITS[units]
        out[f"density [{u['rho']}]"] = rho
        out[f"k [{u['k']}]"] = k
        out[f"c [{u['cp']}]"] = cp
        out["note"] = np.where(df["k_from_R"], "k from R at " + df["R_thickness_mm"].astype(str) + " mm", "")
        out["source"] = df["source"]
        return out.reset_index(drop=True)


def props_in(k_si, rho_si, cp_kj, units):
    """(k, rho, cp) in the network's units; cp in J/kg-K for SI (as the model's formulas use)."""
    if units == "SI":
        return k_si, rho_si, cp_kj * 1000.0
    return k_si * K_IP, rho_si * RHO_IP, cp_kj * CP_IP


def penetration_depth(k_si, rho_si, cp_kj, period_s=DAY):
    """Depth (m) at which a periodic temperature swing falls to 1/e: sqrt(alpha*P/pi)."""
    alpha = k_si / (rho_si * cp_kj * 1000.0)
    return math.sqrt(alpha * period_s / math.pi)


def _g(x, n=5):
    """Short number for formulas."""
    return f"{float(x):.{n}g}"


def layer(mats, material, area, thickness, units="SI", n_segments=0):
    """Values for one layer. area, thickness in m/m² (SI) or ft/ft² (IP).
    n_segments = 0 -> suggested from the daily penetration depth."""
    row = mats.find(material)
    if pd.isna(row.k_W_mK) or pd.isna(row.density_kg_m3) or pd.isna(row.cp_kJ_kgK):
        raise ValueError(f"'{row.material}' needs density, specific heat and k (or R and thickness).")
    k, rho, cp = props_in(row.k_W_mK, row.density_kg_m3, row.cp_kJ_kgK, units)
    L_m = thickness if units == "SI" else thickness * M_PER_FT
    delta = penetration_depth(row.k_W_mK, row.density_kg_m3, row.cp_kJ_kgK)
    # segments no thicker than half of delta; light layers (< 20 kJ/m2K) need only one
    light = row.density_kg_m3 * row.cp_kJ_kgK * L_m < 20.0
    suggested = 1 if light else max(1, math.ceil(2 * L_m / delta))
    n = int(n_segments) if n_segments and n_segments > 0 else suggested
    seg = thickness / n
    u = UNITS[units]
    R = thickness / k
    C = rho * cp * area * thickness
    rows = [
        ("thermal resistance R = L/k", R, u["R"]),
        ("U = k/L (no surface films)", 1 / R, u["U"]),
        ("whole layer: conduction UA = k·A/L", k * area / thickness, u["UA"]),
        ("whole layer: capacitance C = ρ·c·A·L", C, u["C"]),
        ("heat capacity per area ρ·c·L", rho * cp * thickness / (1000.0 if units == "SI" else 1.0), u["Ca"]),
        (f"per segment (L/{n}): conduction", k * area / seg, u["UA"]),
        ("per segment: capacitance", rho * cp * area * seg, u["C"]),
        ("surface → first segment centre (L/2n): conduction", k * area / (seg / 2), u["UA"]),
        ("daily penetration depth δ", delta if units == "SI" else delta / M_PER_FT, u["L"]),
    ]
    table = pd.DataFrame([(q, f"{v:,.4g}", un) for q, v, un in rows], columns=["quantity", "value", "units"])
    kf = f"k={_g(k)}"
    formulas = {
        "per segment conduction": f"conduction(A={_g(area)}, L={_g(seg)}, {kf})",
        "per segment capacitance": f"mass({_g(area)}, {_g(seg)}, rho={_g(rho)}, cp={_g(cp)})",
        "surface → first segment": f"conduction(A={_g(area)}, L={_g(seg / 2)}, {kf})",
        "whole layer conduction": f"conduction(A={_g(area)}, L={_g(thickness)}, {kf})",
        "whole layer capacitance": f"mass({_g(area)}, {_g(thickness)}, rho={_g(rho)}, cp={_g(cp)})",
    }
    info = dict(material=row.material, k=k, rho=rho, cp=cp, n=n, suggested=suggested,
                k_from_R=bool(row.k_from_R), source=row.source)
    return table, formulas, info


def assembly(mats, layers_text, area=1.0, units="SI", R_inside=None, R_outside=None):
    """Layers listed one per line as 'material, thickness' (inside to outside).
    Returns a per-layer table and totals. Surface resistances default to
    0.12 (inside) and 0.03 (outside) m²K/W, converted for IP."""
    R_in = (0.12 if units == "SI" else 0.12 * R_IP) if R_inside is None else R_inside
    R_out = (0.03 if units == "SI" else 0.03 * R_IP) if R_outside is None else R_outside
    rows = []
    for line in str(layers_text).splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        name, _, t = line.rpartition(",")
        if not name:
            raise ValueError(f"write each layer as 'material, thickness' (got '{line}')")
        t = float(t)
        _, _, info = layer(mats, name.strip(), area, t, units)
        R = t / info["k"]
        Ca = info["rho"] * info["cp"] * t
        rows.append(dict(material=info["material"], thickness=t, k=info["k"], R=R,
                         heat_capacity_per_area=Ca, suggested_segments=info["suggested"]))
    u = UNITS[units]
    df = pd.DataFrame(rows)
    R_total = R_in + R_out + (df["R"].sum() if len(df) else 0.0)
    totals = dict(R_surfaces=R_in + R_out, R_total=R_total, U=1 / R_total, UA=area / R_total,
                  C=(df["heat_capacity_per_area"].sum() * area) if len(df) else 0.0)
    df = df.rename(columns={"thickness": f"thickness [{u['L']}]", "k": f"k [{u['k']}]",
                            "R": f"R [{u['R']}]",
                            "heat_capacity_per_area": f"ρ·c·L [{u['C']} per {u['A']}]"})
    return df, totals
