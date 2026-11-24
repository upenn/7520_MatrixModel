"""
thermal_matrix.py -- a lumped-parameter (RC network) thermal "matrix model".

Python version of the spreadsheet method in Braham (1979), "A Brief Analysis of
the Performance of the CNRS Passive Solar House Collector: The Trombe Wall".

For every node i with capacitance C_i, an implicit (backward-difference) energy
balance over a time step dt gives

    (C_i/dt + sum_j UA_ij + UA_io) T_i(t) - sum_j UA_ij T_j(t)
        = (C_i/dt) T_i(t-dt) + UA_io T_o(t) + A_i * I(t) + Q_i(t)

or, for all nodes at once,     M1 . T(t) = M2(t)     ->     T(t) = M1^-1 . M2(t)

M1 depends only on the network (conductances and capacitances), so it is inverted
once; each hour M2 is rebuilt from last hour's temperatures and the weather.
Being fully implicit the scheme is unconditionally stable for any time step.

Units: a Network is either "IP" (F, Btu/hr-F, Btu/F, ft2, Btu/hr-ft2, Btu/hr;
time step in hours -- exactly as in the spreadsheet) or "SI" (C, W/K, J/K, m2,
W/m2, W; time step in seconds). Weather is stored in SI and converted on use.
"""
import copy
import io
import json
import math
import os
import re
import urllib.request
import uuid
import zipfile
from dataclasses import dataclass

# Stamped by build_notebook.py with the build date/time; recorded in saved networks and results.
MODEL_VERSION = "dev (not built)"

import numpy as np
import pandas as pd

__all__ = [
    "Network", "Weather", "Props", "Results", "MATERIALS", "MODEL_VERSION",
    "read_epw", "simulate", "plot_network", "edit_network",
    "add_mass_wall", "example_trombe_wall", "example_simple_room", "example_vent_room",
    "example_excel_weather",
]

# ---------------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------------
W_PER_BTUH = 0.29307107            # 1 Btu/hr = 0.293 W
WK_PER_BTUHF = 0.52752793          # 1 Btu/hr-F = 0.5275 W/K
JK_PER_BTUF = 1899.1005            # 1 Btu/F = 1899.1 J/K
M2_PER_FT2 = 0.09290304
WM2_PER_BTUHFT2 = 3.1545907        # 1 Btu/hr-ft2 = 3.1546 W/m2

UNIT_LABELS = {
    "IP": dict(T="°F", UA="Btu/hr·°F", C="Btu/°F", area="ft²", I="Btu/hr·ft²",
               Q="Btu/hr", E="kBtu", length="ft", volume="ft³", dt="hr"),
    "SI": dict(T="°C", UA="W/K", C="J/K", area="m²", I="W/m²",
               Q="W", E="kWh", length="m", volume="m³", dt="s"),
}


def c_to_f(t):
    return np.asarray(t) * 1.8 + 32.0


def f_to_c(t):
    return (np.asarray(t) - 32.0) / 1.8


# ---------------------------------------------------------------------------
# Materials and property helpers
# ---------------------------------------------------------------------------
# Typical values (SI): k [W/m-K], rho [kg/m3], cp [J/kg-K]. Approximate --
# students should replace with data for their own products where known.
MATERIALS = {
    "concrete":      dict(k=1.75,  rho=2200, cp=920),   # as in the Trombe spreadsheet
    "brick":         dict(k=0.72,  rho=1920, cp=840),
    "adobe":         dict(k=0.52,  rho=1600, cp=900),
    "rammed earth":  dict(k=1.10,  rho=2000, cp=900),
    "stone":         dict(k=2.30,  rho=2600, cp=840),
    "water":         dict(k=0.60,  rho=1000, cp=4186),
    "glass":         dict(k=0.96,  rho=2500, cp=840),
    "gypsum board":  dict(k=0.16,  rho=800,  cp=1090),
    "plaster":       dict(k=0.70,  rho=1600, cp=1000),
    "softwood":      dict(k=0.12,  rho=510,  cp=1380),
    "hardwood":      dict(k=0.16,  rho=720,  cp=1260),
    "CLT":           dict(k=0.13,  rho=480,  cp=1600),
    "mineral wool":  dict(k=0.040, rho=30,   cp=840),
    "EPS":           dict(k=0.035, rho=20,   cp=1450),
    "XPS":           dict(k=0.030, rho=35,   cp=1450),
    "straw bale":    dict(k=0.060, rho=110,  cp=2000),
    "steel":         dict(k=50.0,  rho=7850, cp=450),
    "soil (moist)":  dict(k=1.50,  rho=1800, cp=1300),
    "air":           dict(k=0.026, rho=1.2,  cp=1005),
}

SIGMA_SI = 5.670374e-8
DEFAULT_SOLAR_SOURCES = ["I_S", "I_E", "I_W", "I_N", "I_H"]


class Props:
    """Property calculators that return values in a network's unit system.

    Geometry is given in the network's units (ft / ft2 / ft3 for IP, m / m2 / m3
    for SI). Materials are looked up by name in MATERIALS (stored in SI) or given
    explicitly *in the network's units* (IP: k Btu/hr-ft-F, rho lb/ft3,
    cp Btu/lb-F, h and U Btu/hr-ft2-F; SI: W/m-K, kg/m3, J/kg-K, W/m2-K).
    """

    def __init__(self, units="SI"):
        units = units.upper()
        if units not in UNIT_LABELS:
            raise ValueError("units must be 'SI' or 'IP'")
        self.units = units

    # --- unit plumbing -----------------------------------------------------
    def _mat(self, material, k=None, rho=None, cp=None):
        """Return (k, rho, cp) in this unit system."""
        if material is not None:
            if material not in MATERIALS:
                raise KeyError(f"Unknown material '{material}'. Known: {', '.join(MATERIALS)}")
            m = MATERIALS[material]
            ks, rs, cs = m["k"], m["rho"], m["cp"]
            if self.units == "IP":
                ks, rs, cs = ks * 0.577789, rs * 0.0624280, cs / 4186.8
            k = ks if k is None else k
            rho = rs if rho is None else rho
            cp = cs if cp is None else cp
        return k, rho, cp

    # --- conductances (UA) -------------------------------------------------
    def conduction(self, A, L, k=None, material=None):
        """Conductance through a slab: k*A/L."""
        k, _, _ = self._mat(material, k=k)
        if k is None:
            raise ValueError("give k= or material=")
        return k * A / L

    def convection(self, A, h=None):
        """Surface film conductance h*A. Defaults: h = 3 W/m2K (0.53 Btu/hr-ft2-F),
        a typical still-air interior value; use ~15-25 W/m2K outdoors in wind."""
        if h is None:
            h = 3.0 if self.units == "SI" else 0.53
        return h * A

    def radiation(self, A, e1=0.9, e2=0.9, T_mean=None):
        """Linearised long-wave radiation between two parallel surfaces:
        4*sigma*Tm^3*A / (1/e1 + 1/e2 - 1). T_mean in network temperature units
        (default 20 C / 68 F)."""
        if self.units == "SI":
            Tk = (20.0 if T_mean is None else T_mean) + 273.15
            hr = 4 * SIGMA_SI * Tk ** 3
        else:
            Tr = (68.0 if T_mean is None else T_mean) + 459.67
            hr = 4 * 0.1714e-8 * Tr ** 3
        return hr * A / (1 / e1 + 1 / e2 - 1)

    def window(self, A, U):
        """U-value times area."""
        return U * A

    def ventilation(self, volume, ach):
        """Air exchange: rho*cp*volume*ACH/3600 (SI W/K) or 0.018*V*ACH (IP)."""
        if self.units == "SI":
            return 1.2 * 1005 * volume * ach / 3600.0
        return 0.018 * volume * ach

    def airflow(self, flow):
        """Air flow conductance: SI flow in m3/s (-> W/K), IP flow in cfm (-> Btu/hr-F)."""
        return 1206.0 * flow if self.units == "SI" else 1.08 * flow

    def series(self, *ua):
        """Conductances in series: 1 / sum(1/UA)."""
        return 1.0 / sum(1.0 / u for u in ua)

    def parallel(self, *ua):
        return float(sum(ua))

    # --- capacitances --------------------------------------------------------
    def mass(self, A, L, material=None, rho=None, cp=None):
        """Capacitance of a slab: rho*cp*A*L."""
        _, rho, cp = self._mat(material, rho=rho, cp=cp)
        if rho is None or cp is None:
            raise ValueError("give material= or rho= and cp=")
        return rho * cp * A * L

    def air_mass(self, volume, multiplier=1.0):
        """Capacitance of room air (x multiplier, often 3-10 to account for furniture)."""
        rc = 1.2 * 1005 if self.units == "SI" else 0.075 * 0.24
        return rc * volume * multiplier

    def aperture(self, A, tau=1.0, alpha=1.0, F=1.0):
        """Effective solar aperture A*tau*alpha*F: glazing area, glazing solar
        transmittance, solar absorptance of the receiving surface, and the fraction
        of the transmitted sun that reaches this node. Multiply by a solar source,
        e.g. aperture(6, 0.75, 0.9)*I_S."""
        return A * tau * alpha * F

    def namespace(self):
        """Names usable inside expressions typed into the network editor."""
        ns = {n: getattr(self, n) for n in (
            "conduction", "convection", "radiation", "window", "ventilation",
            "airflow", "series", "parallel", "mass", "air_mass", "aperture")}
        ns.update({k: getattr(math, k) for k in ("pi", "sqrt", "exp", "log", "sin", "cos")})
        ns["MATERIALS"] = MATERIALS
        return ns


# ---------------------------------------------------------------------------
# Weather
# ---------------------------------------------------------------------------
class Weather:
    """Hourly weather, stored in SI. Columns:
        T_out   dry-bulb [C]           T_ground  monthly ground temp [C]
        T_sky   sky temperature [C]    RH, wind  [%], [m/s]
        GHI, DNI, DHI  [W/m2]          I_<name>  irradiance on a surface [W/m2]
    """

    def __init__(self, df, lat=None, lon=None, tz=None, elevation=None, name=""):
        self.df = df
        self.lat, self.lon, self.tz, self.elevation, self.name = lat, lon, tz, elevation, name

    def __repr__(self):
        return f"<Weather '{self.name}': {len(self.df)} h, columns {list(self.df.columns)}>"

    # --- sources the network can refer to ----------------------------------
    @property
    def solar_sources(self):
        return [c for c in self.df.columns if c.startswith("I_")]

    @property
    def temperature_sources(self):
        return [c for c in self.df.columns if c.startswith("T_")]

    # --- solar geometry ----------------------------------------------------
    def sun_position(self):
        """Solar altitude and azimuth (degrees, azimuth clockwise from north) at
        the middle of each hour (EPW values are totals over the preceding hour)."""
        if self.lat is None:
            raise ValueError("this weather has no location; cannot compute sun position")
        t = self.df.index - pd.Timedelta(minutes=30)
        n = t.day_of_year.values
        hr = t.hour.values + t.minute.values / 60.0
        B = 2 * np.pi * (n - 1) / 365.0
        decl = (0.006918 - 0.399912 * np.cos(B) + 0.070257 * np.sin(B)
                - 0.006758 * np.cos(2 * B) + 0.000907 * np.sin(2 * B)
                - 0.002697 * np.cos(3 * B) + 0.00148 * np.sin(3 * B))
        eot = 229.18 * (0.000075 + 0.001868 * np.cos(B) - 0.032077 * np.sin(B)
                        - 0.014615 * np.cos(2 * B) - 0.040849 * np.sin(2 * B))
        solar_time = hr + (4 * (self.lon - 15 * self.tz) + eot) / 60.0
        w = np.radians(15 * (solar_time - 12))
        phi = np.radians(self.lat)
        east = -np.cos(decl) * np.sin(w)
        north = np.sin(decl) * np.cos(phi) - np.cos(decl) * np.cos(w) * np.sin(phi)
        up = np.sin(phi) * np.sin(decl) + np.cos(phi) * np.cos(decl) * np.cos(w)
        alt = np.degrees(np.arcsin(np.clip(up, -1, 1)))
        az = np.degrees(np.arctan2(east, north)) % 360
        return pd.DataFrame({"altitude": alt, "azimuth": az}, index=self.df.index)

    def add_surface(self, name, tilt=90.0, azimuth=180.0, albedo=0.2):
        """Add column I_<name>: total irradiance [W/m2] on a plane.
        tilt: 0 = horizontal, 90 = vertical.  azimuth: compass direction the
        surface faces (0 N, 90 E, 180 S, 270 W). Isotropic-sky model:
        beam*cos(theta) + DHI*(1+cos tilt)/2 + GHI*albedo*(1-cos tilt)/2."""
        need = {"GHI", "DNI", "DHI"}
        if not need <= set(self.df.columns):
            raise ValueError("this weather has no GHI/DNI/DHI columns; cannot add surfaces")
        sp = self.sun_position()
        a, z = np.radians(sp["altitude"].values), np.radians(sp["azimuth"].values)
        b, g = np.radians(tilt), np.radians(azimuth)
        sun = np.stack([np.cos(a) * np.sin(z), np.cos(a) * np.cos(z), np.sin(a)])
        nrm = np.array([np.sin(b) * np.sin(g), np.sin(b) * np.cos(g), np.cos(b)])
        cos_inc = np.clip(nrm @ sun, 0, None) * (np.sin(a) > 0)
        I = (self.df["DNI"].values * cos_inc
             + self.df["DHI"].values * (1 + np.cos(b)) / 2
             + self.df["GHI"].values * albedo * (1 - np.cos(b)) / 2)
        self.df["I_" + name] = I
        return self

    def add_standard_surfaces(self):
        for nm, az in (("S", 180), ("E", 90), ("W", 270), ("N", 0)):
            self.add_surface(nm, 90, az)
        self.add_surface("H", 0, 180)
        return self

    # --- conversion --------------------------------------------------------
    def in_units(self, units):
        """DataFrame converted to 'IP' or 'SI'."""
        df = self.df.copy()
        if units == "IP":
            for c in df.columns:
                if c.startswith("T_"):
                    df[c] = c_to_f(df[c])
                elif c.startswith("I_") or c in ("GHI", "DNI", "DHI"):
                    df[c] = df[c] / WM2_PER_BTUHFT2
        return df

    def plot(self, units="SI", start=None, end=None):
        import matplotlib.pyplot as plt
        df = self.in_units(units)
        if start or end:
            df = df.loc[_period_mask(df.index, start, end)]
        lab = UNIT_LABELS[units]
        fig, ax = plt.subplots(2, 1, figsize=(11, 5.5), sharex=True)
        ax[0].plot(df.index, df["T_out"], lw=0.8, color="tab:red", label="outdoor dry-bulb")
        if "T_ground" in df:
            ax[0].plot(df.index, df["T_ground"], lw=1.2, color="tab:brown", label="ground")
        ax[0].set_ylabel(lab["T"]); ax[0].legend(loc="upper left", fontsize=8)
        for c in self.solar_sources:
            ax[1].plot(df.index, df[c], lw=0.7, label=c)
        ax[1].set_ylabel(lab["I"]); ax[1].legend(loc="upper left", fontsize=8, ncol=6)
        ax[0].set_title(self.name)
        fig.tight_layout()
        return fig


def _period_mask(index, start=None, end=None):
    """Mask selecting 'MM-DD' to 'MM-DD' inclusive (ignores year)."""
    md = index.strftime("%m-%d")
    m = np.ones(len(index), bool)
    if start:
        m &= md >= start
    if end:
        m &= md <= end
    return m


def _read_text(source):
    """Read an EPW from a path, URL, or bytes; unzips .zip files."""
    if isinstance(source, (bytes, bytearray)):
        raw = bytes(source)
    elif re.match(r"^https?://", str(source)):
        req = urllib.request.Request(source, headers={"User-Agent": "thermal-matrix-model"})
        with urllib.request.urlopen(req, timeout=60) as r:
            raw = r.read()
    else:
        with open(source, "rb") as f:
            raw = f.read()
    if raw[:2] == b"PK":
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            names = [n for n in z.namelist() if n.lower().endswith(".epw")]
            if not names:
                raise ValueError("zip file contains no .epw")
            raw = z.read(names[0])
    return raw.decode("latin-1")


def read_epw(source, standard_surfaces=True):
    """Read an EnergyPlus weather file (path, URL, .zip or raw bytes)."""
    text = _read_text(source)
    lines = text.splitlines()
    loc = lines[0].split(",")
    name = ", ".join(s.strip() for s in loc[1:4] if s.strip())
    lat, lon, tz, elev = (float(loc[i]) for i in (6, 7, 8, 9))
    ground = None
    for ln in lines[:8]:
        if ln.upper().startswith("GROUND TEMPERATURES"):
            p = ln.split(",")
            try:
                if int(p[1]) > 0:
                    ground = [float(x) for x in p[6:18]]
            except (ValueError, IndexError):
                ground = None
    data = pd.read_csv(io.StringIO("\n".join(lines[8:])), header=None)
    n = len(data)
    year = 2001 if n < 8784 else 2004  # synthetic year: TMY files mix years
    idx = pd.Timestamp(year, 1, 1) + pd.to_timedelta(np.arange(1, n + 1), unit="h")
    df = pd.DataFrame(index=idx)
    df.index.name = "time (hour ending, local std.)"
    df["T_out"] = data[6].values.astype(float)
    df["RH"] = data[8].values.astype(float)
    ir = data[12].values.astype(float)
    df["T_sky"] = (np.clip(ir, 1, None) / SIGMA_SI) ** 0.25 - 273.15
    if ground:
        df["T_ground"] = np.array(ground)[(data[1].values.astype(int) - 1)]
    df["wind"] = data[21].values.astype(float)
    df["GHI"] = data[13].values.astype(float)
    df["DNI"] = data[14].values.astype(float)
    df["DHI"] = data[15].values.astype(float)
    w = Weather(df, lat, lon, tz, elev, name)
    if standard_surfaces:
        w.add_standard_surfaces()
    return w


# ---------------------------------------------------------------------------
# Network
# ---------------------------------------------------------------------------
class Network:
    """A thermal RC network.

    nodes:      {name, C, T0, solar, Q, heat_set, cool_set, x, y}
    boundaries: {name, source, x, y}  -- source is a weather column (T_out, T_ground,
                T_sky) or a constant temperature.
    links:      {a, b, UA, label}

    Any numeric field may be an expression string using Props helpers, e.g.
    "conduction(A=48, L=0.12, material='concrete')".

    `solar` is the absorbed solar gain, an expression in the weather's solar
    sources (I_S, I_E, I_W, I_N, I_H, ...) so a node can receive sun through
    windows in several directions, e.g. "aperture(6, 0.75, 0.9)*I_S + 3*0.6*I_E".
    (Older networks with solar_area/solar_source are converted on load.)
    """

    NODE_DEFAULTS = dict(C=0, T0=None, solar="", Q=0,
                         heat_set=None, cool_set=None, x=None, y=None)

    @staticmethod
    def _convert_legacy_solar(n):
        """solar_area + solar_source (older networks) -> solar expression."""
        area, src = n.pop("solar_area", None), n.pop("solar_source", None)
        if n.get("solar") in (None, "", 0) and src and area not in (None, "", 0, "0"):
            n["solar"] = f"{area}*{src}" if _is_number(area) else f"({area})*{src}"
        return n

    def __init__(self, units="SI", title=""):
        self.units = units.upper()
        self.title = title
        self.nodes, self.boundaries, self.links = [], [], []

    # --- building -----------------------------------------------------------
    @property
    def props(self):
        return Props(self.units)

    def _names(self):
        return [n["name"] for n in self.nodes] + [b["name"] for b in self.boundaries]

    def add_node(self, name, C=0, T0=None, solar="", Q=0,
                 heat_set=None, cool_set=None, x=None, y=None,
                 solar_area=None, solar_source=None):
        """solar: absorbed solar gain expression, e.g. "aperture(6, 0.75, 0.9)*I_S".
        solar_area/solar_source are accepted for older code and converted."""
        if name in self._names():
            raise ValueError(f"name '{name}' already used")
        n = dict(name=name, C=C, T0=T0, solar=solar, Q=Q, heat_set=heat_set,
                 cool_set=cool_set, x=x, y=y, solar_area=solar_area, solar_source=solar_source)
        self.nodes.append(self._convert_legacy_solar(n))
        return self

    def add_boundary(self, name, source="T_out", x=None, y=None):
        if name in self._names():
            raise ValueError(f"name '{name}' already used")
        self.boundaries.append(dict(name=name, source=source, x=x, y=y))
        return self

    def link(self, a, b, UA, label="", UA_open=None, sensor=None, open_above=None,
             outdoor_above=None, warmer_by=None, band=1.0):
        """Conductance UA between a and b.

        Controlled (switched) link, e.g. a temperature-responsive window: give
        UA_open, the conductance when open. The link opens when the `sensor` node
        is above `open_above` AND the outdoor air (T_out) is above `outdoor_above`
        (leave blank for no outdoor condition) AND the sensor is at least `warmer_by`
        degrees warmer than the outdoor air (blank = no indoor/outdoor test), so it
        doesn't let hot air in. It closes again when any of these drops `band` degrees
        below its set point (a dead band, like a real actuator, so it doesn't flutter
        open/shut every hour). UA is the closed value."""
        l = dict(a=a, b=b, UA=UA, label=label)
        if UA_open not in (None, ""):
            l.update(UA_open=UA_open, sensor=sensor or "", open_above=open_above,
                     outdoor_above=outdoor_above, warmer_by=warmer_by, band=band)
        self.links.append(l)
        return self

    @staticmethod
    def is_controlled(l):
        """True if a link has an 'open' conductance (switched ventilation)."""
        v = l.get("UA_open")
        return v is not None and str(v).strip() != ""

    def sensor_of(self, l):
        """Node whose temperature opens a controlled link (default: its first node end)."""
        s = str(l.get("sensor") or "").strip()
        if s:
            return s
        node_names = {n["name"] for n in self.nodes}
        return l["a"] if l["a"] in node_names else l["b"]

    def node(self, name):
        for n in self.nodes + self.boundaries:
            if n["name"] == name:
                return n
        raise KeyError(name)

    def set(self, name, **kw):
        """Change fields of a node or boundary, e.g. net.set('Room Air', heat_set=20)."""
        n = self.node(name)
        n.update(kw)
        if "solar_area" in kw or "solar_source" in kw:
            n["solar"] = ""
            self._convert_legacy_solar(n)
        return self

    def remove(self, name):
        self.nodes = [n for n in self.nodes if n["name"] != name]
        self.boundaries = [b for b in self.boundaries if b["name"] != name]
        self.links = [l for l in self.links if name not in (l["a"], l["b"])]
        return self

    # --- (de)serialisation --------------------------------------------------
    def to_dict(self):
        return dict(units=self.units, title=self.title, model_version=MODEL_VERSION,
                    nodes=copy.deepcopy(self.nodes),
                    boundaries=copy.deepcopy(self.boundaries), links=copy.deepcopy(self.links))

    @classmethod
    def from_dict(cls, d):
        net = cls(d.get("units", "SI"), d.get("title", ""))
        net._load(d)
        return net

    def _load(self, d):
        self.units = d.get("units", self.units).upper()
        self.title = d.get("title", self.title)
        self.nodes = [self._convert_legacy_solar({**self.NODE_DEFAULTS, **n}) for n in d.get("nodes", [])]
        self.boundaries = [{"source": "T_out", "x": None, "y": None, **b}
                           for b in d.get("boundaries", [])]
        self.links = [{"label": "", **l} for l in d.get("links", [])]
        self.saved_with = d.get("model_version")

    def to_json(self, path=None):
        s = json.dumps(self.to_dict(), indent=1)
        if path:
            with open(path, "w") as f:
                f.write(s)
        return s

    @classmethod
    def from_json(cls, s_or_path):
        if os.path.exists(str(s_or_path)):
            with open(s_or_path) as f:
                s_or_path = f.read()
        return cls.from_dict(json.loads(s_or_path))

    def copy(self):
        return Network.from_dict(self.to_dict())

    # --- evaluation ---------------------------------------------------------
    def value(self, v, what="value"):
        """Evaluate a number or expression string."""
        if v is None or (isinstance(v, str) and v.strip() == ""):
            return None
        if isinstance(v, (int, float, np.number)):
            return float(v)
        try:
            return float(eval(str(v), {"__builtins__": {}}, self.props.namespace()))
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"could not evaluate {what} = {v!r}: {e}") from None

    @staticmethod
    def has_solar(expr):
        return expr is not None and str(expr).strip() not in ("", "0", "0.0")

    def solar_gain(self, expr, sources, what="solar gain"):
        """Evaluate a solar gain expression. `sources` maps solar source names
        (I_S, I_E, ...) to irradiance values or arrays in the network's units."""
        if not self.has_solar(expr):
            return 0.0
        ns = self.props.namespace()
        ns.update(sources)
        try:
            return eval(str(expr), {"__builtins__": {}}, ns)
        except NameError as e:
            raise ValueError(f"could not evaluate {what} = {expr!r}: {e}. Solar sources "
                             f"available: {', '.join(sorted(sources))}") from None
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"could not evaluate {what} = {expr!r}: {e}") from None

    @staticmethod
    def solar_names(expr):
        """Solar source names (I_...) used in an expression."""
        return sorted(set(re.findall(r"\bI_\w+", str(expr or ""))))

    def evaluated(self):
        """Return (nodes_df, links_df) with every expression evaluated
        (the solar gain stays an expression; it depends on the weather)."""
        rows = []
        for n in self.nodes:
            rows.append(dict(name=n["name"], C=self.value(n["C"], f"C of {n['name']}") or 0.0,
                             T0=self.value(n["T0"]),
                             solar=str(n.get("solar") or "") if self.has_solar(n.get("solar")) else "",
                             Q=n["Q"] if isinstance(n["Q"], (list, tuple)) else (self.value(n["Q"]) or 0.0),
                             heat_set=self.value(n["heat_set"]), cool_set=self.value(n["cool_set"])))
        nodes = pd.DataFrame(rows, columns=["name", "C", "T0", "solar",
                                            "Q", "heat_set", "cool_set"]).set_index("name")
        lrows = [dict(a=l["a"], b=l["b"], UA=self.value(l["UA"], f"UA of {l['a']}–{l['b']}"),
                      label=l.get("label", "")) for l in self.links]
        cols = ["a", "b", "UA", "label"]
        if any(self.is_controlled(l) for l in self.links):
            cols += ["UA_open", "sensor", "open_above", "outdoor_above", "warmer_by", "band"]
            for r, l in zip(lrows, self.links):
                if self.is_controlled(l):
                    w = f"{l['a']}–{l['b']}"
                    band = self.value(l.get("band"), f"dead band of {w}")
                    r.update(UA_open=self.value(l["UA_open"], f"open UA of {w}"),
                             sensor=self.sensor_of(l),
                             open_above=self.value(l.get("open_above"), f"opening temperature of {w}"),
                             outdoor_above=self.value(l.get("outdoor_above"), f"outdoor temperature of {w}"),
                             warmer_by=self.value(l.get("warmer_by"), f"indoor–outdoor difference of {w}"),
                             band=1.0 if band is None else band)
        links = pd.DataFrame(lrows, columns=cols)
        return nodes, links

    def check(self, weather=None):
        """List of problems (empty list = OK)."""
        msgs = []
        names = self._names()
        if len(set(names)) != len(names):
            msgs.append("Duplicate names: " + ", ".join(sorted({x for x in names if names.count(x) > 1})))
        if not self.nodes:
            msgs.append("The network has no nodes.")
        try:
            nodes, links = self.evaluated()
        except ValueError as e:
            return msgs + [str(e)]
        bnames = {b["name"] for b in self.boundaries}
        for _, l in links.iterrows():
            for end in (l.a, l.b):
                if end not in names:
                    msgs.append(f"Link {l.a}–{l.b} refers to unknown '{end}'.")
            if l.UA is None or not np.isfinite(l.UA) or l.UA <= 0:
                msgs.append(f"Link {l.a}–{l.b} needs a positive UA (has {l.UA}).")
            if l.a == l.b:
                msgs.append(f"Link {l.a}–{l.b} connects a node to itself.")
            if l.a in bnames and l.b in bnames:
                msgs.append(f"Link {l.a}–{l.b} joins two boundaries; it has no effect.")
            if "UA_open" in links and pd.notna(l.UA_open):
                w = f"Controlled link {l.a}–{l.b}"
                if not np.isfinite(l.UA_open) or l.UA_open <= 0:
                    msgs.append(f"{w} needs a positive open UA (has {l.UA_open}).")
                if l.open_above is None or pd.isna(l.open_above):
                    msgs.append(f"{w} needs an opening temperature (open when the sensor is above …).")
                if l.sensor not in nodes.index:
                    msgs.append(f"{w}: sensor '{l.sensor}' is not a node.")
                if l.band < 0:
                    msgs.append(f"{w}: the dead band can't be negative.")
        for name, n in nodes.iterrows():
            if n.C < 0:
                msgs.append(f"Node '{name}' has negative capacitance.")
            if n.solar:
                avail = weather.solar_sources if weather is not None else DEFAULT_SOLAR_SOURCES
                used = self.solar_names(n.solar)
                if not used:
                    msgs.append(f"Node '{name}': solar gain '{n.solar}' doesn't use a solar source; "
                                f"multiply each aperture by one, e.g. aperture(6, 0.75, 0.9)*I_S.")
                else:
                    try:
                        self.solar_gain(n.solar, {c: 1.0 for c in set(avail) | set(used)},
                                        f"solar gain of {name}")
                    except ValueError as e:
                        msgs.append(str(e))
                    missing = [c for c in used if c not in avail]
                    if missing:
                        msgs.append(f"Node '{name}': solar source {', '.join(missing)} not in the weather "
                                    f"(available: {', '.join(avail)}; add other surfaces in step 1b).")
            if n.heat_set is not None and n.cool_set is not None and n.heat_set > n.cool_set:
                msgs.append(f"Node '{name}': heating set point is above cooling set point.")
        for b in self.boundaries:
            src = b["source"]
            if isinstance(src, str) and not _is_number(src):
                if weather is not None and src not in weather.df:
                    msgs.append(f"Boundary '{b['name']}': source '{src}' not in weather "
                                f"(available: {', '.join(weather.temperature_sources)} or a number).")
        # every group of connected nodes needs capacitance or a path to a boundary
        if not msgs:
            msgs += self._floating_groups(nodes, links, bnames)
        return msgs

    def _floating_groups(self, nodes, links, bnames):
        idx = {n: i for i, n in enumerate(nodes.index)}
        parent = list(range(len(idx)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        grounded = set()
        for _, l in links.iterrows():
            if l.a in idx and l.b in idx:
                parent[find(idx[l.a])] = find(idx[l.b])
        for _, l in links.iterrows():
            for x, y in ((l.a, l.b), (l.b, l.a)):
                if x in idx and y in bnames:
                    grounded.add(x)
        groups = {}
        for n, i in idx.items():
            groups.setdefault(find(i), []).append(n)
        out = []
        for members in groups.values():
            if not any(nodes.loc[m, "C"] > 0 for m in members) and not any(m in grounded for m in members):
                out.append("Nodes " + ", ".join(members) + " have no capacitance and no link to a "
                           "boundary, so their temperature is undefined (the matrix is singular).")
        return out

    # --- matrices -----------------------------------------------------------
    def matrices(self, dt=None):
        """Return dict with UA (node-node), UAb (node-boundary), C, M1 and labels.
        dt defaults to 1 hour."""
        if dt is None:
            dt = 1.0 if self.units == "IP" else 3600.0
        nodes, links = self.evaluated()
        names = list(nodes.index)
        bnames = [b["name"] for b in self.boundaries]
        ix, bx = {n: i for i, n in enumerate(names)}, {n: i for i, n in enumerate(bnames)}
        UA = np.zeros((len(names), len(names)))
        UAb = np.zeros((len(names), len(bnames)))
        for _, l in links.iterrows():
            a, b, u = l.a, l.b, l.UA
            if a in ix and b in ix:
                UA[ix[a], ix[b]] += u
                UA[ix[b], ix[a]] += u
            elif a in ix and b in bx:
                UAb[ix[a], bx[b]] += u
            elif b in ix and a in bx:
                UAb[ix[b], bx[a]] += u
        C = nodes["C"].values.astype(float)
        M1 = np.diag(UA.sum(1) + UAb.sum(1) + C / dt) - UA
        return dict(names=names, bnames=bnames, UA=UA, UAb=UAb, C=C, M1=M1, dt=dt, nodes=nodes,
                    links=links)

    def show(self, weather=None):
        """Print tables, check the network, and draw it."""
        from IPython.display import display, Markdown
        lab = UNIT_LABELS[self.units]
        msgs = self.check(weather)
        display(Markdown(f"### {self.title or 'Network'}  ({self.units} units)"))
        if msgs:
            display(Markdown("**Problems:**\n" + "\n".join(f"- ⚠️ {m}" for m in msgs)))
        else:
            display(Markdown("✅ Network is complete and solvable."))
        try:
            nodes, links = self.evaluated()
            nodes = nodes.rename(columns={"C": f"C [{lab['C']}]", "T0": f"T0 [{lab['T']}]",
                                          "solar": f"solar gain [{lab['Q']}]",
                                          "Q": f"Q [{lab['Q']}]"})
            display(nodes)
            display(links.rename(columns={"UA": f"UA [{lab['UA']}]"}))
            display(pd.DataFrame(self.boundaries)[["name", "source"]] if self.boundaries else
                    Markdown("_no boundaries_"))
        except ValueError as e:
            print(e)
        plot_network(self)

    def ua_table(self):
        """The UA_ij matrix with boundary columns, like the spreadsheet's Parameters sheet."""
        m = self.matrices()
        df = pd.DataFrame(m["UA"], index=m["names"], columns=m["names"])
        for j, b in enumerate(m["bnames"]):
            df[f"→ {b}"] = m["UAb"][:, j]
        df["C"] = m["C"]
        return df

    def __repr__(self):
        return (f"<Network '{self.title}' {self.units}: {len(self.nodes)} nodes, "
                f"{len(self.boundaries)} boundaries, {len(self.links)} links>")


def _is_number(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def add_mass_wall(net, prefix, A, thickness, n_segments, outside, inside,
                  material=None, k=None, rho=None, cp=None, T0=None, x0=None, y0=None, dx=None):
    """Discretise a homogeneous mass wall into n_segments capacitive nodes
    between two existing nodes `outside` and `inside` (usually surface nodes).
    Segment nodes sit at segment centres: outer links k*A/(L/2), inner k*A/L."""
    p = net.props
    L = thickness / n_segments
    ua_in = p.conduction(A, L, k=k, material=material)
    ua_half = p.conduction(A, L / 2, k=k, material=material)
    C = p.mass(A, L, material=material, rho=rho, cp=cp)
    names = [f"{prefix} {i + 1}" for i in range(n_segments)]
    for i, nm in enumerate(names):
        x = None if x0 is None else x0 + i * (dx or 80)
        net.add_node(nm, C=round(C, 6), T0=T0, x=x, y=y0)
    net.link(outside, names[0], round(ua_half, 6), "conduction")
    for a, b in zip(names[:-1], names[1:]):
        net.link(a, b, round(ua_in, 6), "conduction")
    net.link(names[-1], inside, round(ua_half, 6), "conduction")
    return names


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
@dataclass
class Results:
    T: pd.DataFrame            # node + boundary temperatures
    loads: pd.DataFrame        # heating (+) / cooling (-) per thermostat node
    solar: pd.DataFrame        # absorbed solar per node
    network: Network
    weather: pd.DataFrame      # weather in network units over the run
    dt: float
    vent: pd.DataFrame = None  # controlled links: fraction of each hour open (0 or 1)

    @property
    def units(self):
        return self.network.units

    def _vent_columns(self):
        """{link row: column of self.vent} for controlled links."""
        if self.vent is None or not len(self.vent.columns):
            return {}
        _, links = self.network.evaluated()
        rows = [i for i, l in links.iterrows() if "UA_open" in links and pd.notna(l.UA_open)]
        return dict(zip(rows, self.vent.columns))

    def flows(self):
        """Heat flow along each link, a -> b, in Q units (Btu/hr or W).
        Controlled links use their open or closed UA hour by hour."""
        _, links = self.network.evaluated()
        vc = self._vent_columns()
        out = {}
        for i, l in links.iterrows():
            ua = l.UA if i not in vc else l.UA + self.vent[vc[i]] * (l.UA_open - l.UA)
            out[f"{l.a} → {l.b}"] = ua * (self.T[l.a] - self.T[l.b])
        return pd.DataFrame(out, index=self.T.index)

    def vent_summary(self):
        """Hours open and heat carried by each controlled link."""
        lab = UNIT_LABELS[self.units]
        _, links = self.network.evaluated()
        vc = self._vent_columns()
        if not vc:
            return pd.DataFrame()
        fl = self.flows()
        hours = self.dt if self.units == "IP" else self.dt / 3600.0
        rows = []
        for i, col in vc.items():
            l = links.loc[i]
            v = self.vent[col]
            node_names = {n["name"] for n in self.network.nodes}
            inside = l.sensor if l.sensor in (l.a, l.b) else (l.a if l.a in node_names else l.b)
            q = fl[f"{l.a} → {l.b}"] * (1 if inside == l.a else -1)   # heat leaving the inside end
            rows.append({"link": col,
                         "opens when": f"{l.sensor} > {l.open_above:g}" +
                                       ("" if pd.isna(l.outdoor_above) else f" and T_out > {l.outdoor_above:g}") +
                                       ("" if pd.isna(l.warmer_by) else f" and {l.sensor} > T_out + {l.warmer_by:g}"),
                         "hours open": v.sum() * hours,
                         "% of time": 100 * v.mean(),
                         "times opened": self.vent.attrs.get("openings", {}).get(col, np.nan),
                         f"heat removed from {inside} while open [{lab['E']}]": self._energy(q * v),
                         f"… while closed [{lab['E']}]": self._energy(q * (1 - v))})
        return pd.DataFrame(rows).set_index("link").round(1)

    def plot_vent(self):
        """Sensor and outdoor temperature with the open periods shaded, per controlled link."""
        import matplotlib.pyplot as plt
        lab = UNIT_LABELS[self.units]
        _, links = self.network.evaluated()
        vc = self._vent_columns()
        if not vc:
            print("No controlled links in this network.")
            return None
        fig, axs = plt.subplots(len(vc), 1, figsize=(11, 3.2 * len(vc)), sharex=True, squeeze=False)
        for ax, (i, col) in zip(axs[:, 0], vc.items()):
            l = links.loc[i]
            v = self.vent[col]
            ax.fill_between(self.T.index, 0, 1, where=v > 0, step="pre", color="#7fb3d5", alpha=0.35,
                            transform=ax.get_xaxis_transform(), label="open")
            ax.plot(self.T.index, self.T[l.sensor], color="#c0392b", lw=1.2, label=f"{l.sensor}")
            ax.plot(self.T.index, self.weather["T_out"], "k--", lw=0.9, label="outdoor")
            ax.axhline(l.open_above, color="#c0392b", lw=0.8, ls=":", label=f"opens above {l.open_above:g}")
            if not pd.isna(l.outdoor_above):
                ax.axhline(l.outdoor_above, color="k", lw=0.8, ls=":", label=f"outdoor above {l.outdoor_above:g}")
            ax.set_ylabel(lab["T"]); ax.grid(alpha=0.3)
            ax.set_title(f"{col}: open {100 * v.mean():.0f}% of the time", fontsize=10)
            ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        return fig

    def _energy(self, q):
        """Sum of hourly power -> energy (kBtu or kWh)."""
        hours = self.dt if self.units == "IP" else self.dt / 3600.0
        return q.sum() * hours / 1000.0

    def summary(self):
        lab = UNIT_LABELS[self.units]
        nodes = [n["name"] for n in self.network.nodes]
        s = pd.DataFrame({f"min [{lab['T']}]": self.T[nodes].min(),
                          f"mean [{lab['T']}]": self.T[nodes].mean(),
                          f"max [{lab['T']}]": self.T[nodes].max(),
                          f"solar absorbed [{lab['E']}]": self._energy(self.solar[nodes])})
        for c in self.loads.columns:
            node = c.split(":", 1)[1]
            s.loc[node, f"heating [{lab['E']}]"] = self._energy(self.loads[c].clip(lower=0))
            s.loc[node, f"cooling [{lab['E']}]"] = abs(self._energy(self.loads[c].clip(upper=0)))
            s.loc[node, f"peak heating [{lab['Q']}]"] = self.loads[c].max()
        return s.round(2).astype(object).where(s.notna(), "")

    def plot(self, nodes=None, show_weather=True, title=None):
        import matplotlib.pyplot as plt
        lab = UNIT_LABELS[self.units]
        nodes = nodes or [n["name"] for n in self.network.nodes]
        nrow = 2 if show_weather else 1
        fig, axs = plt.subplots(nrow, 1, figsize=(11, 4 + 2 * (nrow - 1)), sharex=True,
                                gridspec_kw=dict(height_ratios=[3, 1][:nrow]))
        axs = np.atleast_1d(axs)
        ax = axs[0]
        cmap = plt.get_cmap("viridis", max(len(nodes), 2))
        for i, n in enumerate(nodes):
            ax.plot(self.T.index, self.T[n], lw=1.3, color=cmap(i), label=n)
        if show_weather:
            ax.plot(self.T.index, self.weather["T_out"], "k--", lw=0.9, label="outdoor")
        ax.set_ylabel(f"temperature [{lab['T']}]")
        ax.legend(fontsize=8, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5))
        ax.grid(alpha=0.3)
        ax.set_title(title or self.network.title)
        if show_weather:
            srcs = sorted({s for n in self.network.nodes for s in Network.solar_names(n.get("solar"))})
            for s in srcs:
                if s in self.weather:
                    axs[1].fill_between(self.T.index, self.weather[s], alpha=0.4, label=s)
            axs[1].set_ylabel(lab["I"])
            if srcs:
                axs[1].legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
            axs[1].grid(alpha=0.3)
        fig.tight_layout()
        return fig

    def plot_loads(self):
        import matplotlib.pyplot as plt
        lab = UNIT_LABELS[self.units]
        fig, ax = plt.subplots(figsize=(11, 3))
        for c in self.loads.columns:
            ax.plot(self.loads.index, self.loads[c], lw=1, label=c)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel(f"heating (+) / cooling (−) [{lab['Q']}]")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        fig.tight_layout()
        return fig

    def plot_average_day(self, nodes=None):
        """Mean daily profile (hour of day) for each node."""
        import matplotlib.pyplot as plt
        lab = UNIT_LABELS[self.units]
        nodes = nodes or [n["name"] for n in self.network.nodes]
        hod = (self.T.index - pd.Timedelta(minutes=30)).hour
        g = self.T.groupby(hod).mean()
        fig, ax = plt.subplots(figsize=(8, 4))
        cmap = plt.get_cmap("viridis", max(len(nodes), 2))
        for i, n in enumerate(nodes):
            ax.plot(g.index + 0.5, g[n], color=cmap(i), label=n)
        ax.plot(g.index + 0.5, self.weather["T_out"].groupby(hod).mean(), "k--", label="outdoor")
        ax.set_xlabel("hour of day"); ax.set_ylabel(lab["T"]); ax.set_xticks(range(0, 25, 3))
        ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5)); ax.grid(alpha=0.3)
        ax.set_title("average day")
        fig.tight_layout()
        return fig

    def to_csv(self, path):
        parts = [self.T, self.loads]
        if self.vent is not None and len(self.vent.columns):
            parts.append(self.vent.add_prefix("open: "))
        df = pd.concat(parts + [self.weather[["T_out"]].add_prefix("weather ")], axis=1)
        df.index.name = f"time (hour ending) - model version {MODEL_VERSION}"
        df.to_csv(path, encoding="utf-8")
        return path


def _day_indices(index, start, end):
    """Row positions for 'MM-DD'..'MM-DD' (inclusive), wrapping past Dec 31."""
    hours_start = (index - pd.Timedelta(minutes=30)).strftime("%m-%d")
    md = np.asarray(hours_start)
    pos = np.arange(len(index))
    if start is None:
        start = md[0]
    if end is None:
        end = md[-1]
    if start <= end:
        return pos[(md >= start) & (md <= end)]
    return np.concatenate([pos[md >= start], pos[md <= end]])


def simulate(net, weather, start=None, end=None, warmup_days=7, substeps=1, T_init=None):
    """Run the network through the weather.

    start, end : 'MM-DD' strings (inclusive), default the whole file. A period
                 with start > end wraps around the year end (e.g. '12-01'..'02-28').
    warmup_days: days simulated before `start` (results discarded) so the initial
                 temperatures stop mattering. The spreadsheet uses 0.
    substeps   : time steps per hour (weather held constant within the hour).
    T_init     : initial temperature for nodes without T0 (default: first outdoor temp).
    """
    msgs = net.check(weather)
    if msgs:
        raise ValueError("Network problems:\n  - " + "\n  - ".join(msgs))
    units = net.units
    dt_hour = 1.0 if units == "IP" else 3600.0
    dt = dt_hour / substeps
    m = net.matrices(dt)
    names, bnames, nodes = m["names"], m["bnames"], m["nodes"]
    M1, C = m["M1"], m["C"]
    Minv = np.linalg.inv(M1)

    wdf = weather.in_units(units)
    run = _day_indices(wdf.index, start, end)
    if len(run) == 0:
        raise ValueError(f"no weather between {start} and {end}")
    n_w = len(wdf)
    warm = np.arange(run[0] - int(warmup_days * 24), run[0]) % n_w if warmup_days else np.array([], int)
    seq = np.concatenate([warm, run]).astype(int)

    # boundary temperatures, solar, gains as arrays over the sequence
    Tb = np.zeros((len(seq), len(bnames)))
    for j, b in enumerate(net.boundaries):
        src = b["source"]
        Tb[:, j] = float(src) if _is_number(src) else wdf[src].values[seq]
    S = np.zeros((len(seq), len(names)))
    Q = np.zeros((len(seq), len(names)))
    hod = ((wdf.index[seq] - pd.Timedelta(minutes=30)).hour).values
    for i, (name, n) in enumerate(nodes.iterrows()):
        if n.solar:
            sources = {c: wdf[c].values[seq] for c in wdf.columns if c.startswith("I_")}
            S[:, i] = np.broadcast_to(net.solar_gain(n.solar, sources, f"solar gain of {name}"), len(seq))
        q = n.Q
        if isinstance(q, (list, tuple)):          # 24-hour profile
            Q[:, i] = np.asarray(q, float)[hod]
        elif isinstance(q, str) and q in wdf:     # weather column
            Q[:, i] = wdf[q].values[seq]
        else:
            Q[:, i] = float(q or 0.0)
    forcing = m["UAb"] @ Tb.T            # (n, steps)
    forcing = forcing.T + S + Q

    T0_default = T_init if T_init is not None else wdf["T_out"].values[seq[0]]
    T = np.array([T0_default if pd.isna(v) else v for v in nodes["T0"].values], float)

    heat = nodes["heat_set"].values
    cool = nodes["cool_set"].values
    tstat = [i for i in range(len(names)) if heat[i] is not None and not pd.isna(heat[i])
             or cool[i] is not None and not pd.isna(cool[i])]
    Cdt = C / dt

    # Controlled links (switched ventilation). M1 above uses the closed UA; opening
    # link c adds dUA = UA_open - UA to it. Each open/closed combination gets its own
    # M1, inverted once and cached. Any positive UA keeps M1 diagonally dominant, so
    # the implicit scheme stays stable whatever the switching; the decision uses the
    # temperatures at the end of the previous step, with a dead band against flutter.
    ctrl = []
    if "UA_open" in m["links"]:
        bx = {b: j for j, b in enumerate(bnames)}
        ix = {n: i for i, n in enumerate(names)}
        Tout = wdf["T_out"].values[seq]
        for li, l in m["links"].iterrows():
            if pd.isna(l.UA_open):
                continue
            d = l.UA_open - l.UA
            dM = np.zeros_like(M1)
            df_ = np.zeros(len(names))          # forcing per unit boundary temperature
            ends = [ix.get(l.a), ix.get(l.b)]
            if None not in ends:
                i, j = ends
                dM[i, i] += d; dM[j, j] += d; dM[i, j] -= d; dM[j, i] -= d
                jb = None
            else:
                i = ends[0] if ends[0] is not None else ends[1]
                jb = bx[l.b] if ends[0] is not None else bx[l.a]
                dM[i, i] += d
            oa = None if l.outdoor_above is None or pd.isna(l.outdoor_above) else float(l.outdoor_above)
            wb = None if l.warmer_by is None or pd.isna(l.warmer_by) else float(l.warmer_by)
            nm = f"{l.a} → {l.b}"
            if any(c["name"] == nm for c in ctrl):
                nm += f" ({l.label or li})"
            ctrl.append(dict(link=li, name=nm, dM=dM, node=i, jb=jb, d=d,
                             sensor=ix[l.sensor], on=float(l.open_above), out=oa, dT=wb,
                             band=float(l.band)))
    inv_cache = {(): Minv}
    state = [False] * len(ctrl)

    def inverse(st):
        key = tuple(c for c, s in enumerate(st) if s)
        if key not in inv_cache:
            inv_cache[key] = np.linalg.inv(M1 + sum(ctrl[c]["dM"] for c in key))
        return inv_cache[key]

    out_T = np.zeros((len(seq), len(names)))
    out_L = np.zeros((len(seq), len(names)))
    out_V = np.zeros((len(seq), len(ctrl)))
    out_N = np.zeros((len(seq), len(ctrl)))   # number of openings in each hour
    for k in range(len(seq)):
        for _ in range(substeps):
            f = forcing[k]
            if ctrl:
                for c, cl in enumerate(ctrl):
                    ts, to = T[cl["sensor"]], Tout[k]
                    if state[c]:   # open: close when any test falls below its set point - band
                        state[c] = not (ts < cl["on"] - cl["band"] or
                                        (cl["out"] is not None and to < cl["out"] - cl["band"]) or
                                        (cl["dT"] is not None and ts - to < cl["dT"] - cl["band"]))
                    else:          # closed: open when all tests pass
                        state[c] = (ts > cl["on"] and (cl["out"] is None or to > cl["out"]) and
                                    (cl["dT"] is None or ts - to > cl["dT"]))
                        out_N[k, c] += state[c]
                f = f.copy()
                for c, cl in enumerate(ctrl):
                    if state[c] and cl["jb"] is not None:
                        f[cl["node"]] += cl["d"] * Tb[k, cl["jb"]]
                    out_V[k, c] += state[c] / substeps
                Minv = inverse(state)
            Tn = Minv @ (Cdt * T + f)
            q = np.zeros(len(names))
            if tstat:
                Tn, q = _apply_thermostats(Minv, Tn, tstat, heat, cool)
            T = Tn
        out_T[k] = T
        out_L[k] = q if tstat else 0.0

    keep = slice(len(warm), None)
    idx = _continuous_index(wdf.index, run)
    Tdf = pd.DataFrame(out_T[keep], index=idx, columns=names)
    for j, b in enumerate(bnames):
        Tdf[b] = Tb[keep, j]
    loads = pd.DataFrame({f"load:{names[i]}": out_L[keep, i] for i in tstat}, index=idx)
    solar = pd.DataFrame(S[keep], index=idx, columns=names)
    wout = wdf.iloc[run].copy()
    wout.index = idx
    vent = pd.DataFrame({c["name"]: out_V[keep, n] for n, c in enumerate(ctrl)}, index=idx)
    vent.attrs["openings"] = {c["name"]: int(out_N[keep, n].sum()) for n, c in enumerate(ctrl)}
    return Results(Tdf, loads, solar, net.copy(), wout, dt_hour, vent)


def _apply_thermostats(Minv, Tfree, tstat, heat, cool):
    """Ideal, unlimited heating/cooling at thermostat nodes. Because the system
    is linear, T = Tfree + Minv[:, A] @ q_A; solve for q_A that holds each active
    node at its set point (active-set iteration)."""
    active = {}
    for _ in range(2 * len(tstat) + 2):
        q = np.zeros(len(Tfree))
        if active:
            A = list(active)
            G = Minv[np.ix_(A, A)]
            target = np.array([active[i] for i in A])
            qa = np.linalg.solve(G, target - Tfree[A])
            q[A] = qa
        T = Tfree + Minv @ q
        changed = False
        for i in tstat:
            h, c = heat[i], cool[i]
            h = None if h is None or pd.isna(h) else h
            c = None if c is None or pd.isna(c) else c
            if i in active:
                if (active[i] == h and q[i] < -1e-9) or (active[i] == c and q[i] > 1e-9):
                    del active[i]; changed = True   # wrong sign: release
            else:
                if h is not None and T[i] < h - 1e-9:
                    active[i] = h; changed = True
                elif c is not None and T[i] > c + 1e-9:
                    active[i] = c; changed = True
        if not changed:
            return T, q
    return T, q


def _continuous_index(index, run):
    """Monotonic timestamps for a run that may wrap past the end of the year."""
    t = index[run]
    jumps = np.where(np.diff(t.values) < np.timedelta64(0))[0]
    if len(jumps):
        off = np.zeros(len(t), "timedelta64[ns]")
        span = index[-1] - index[0] + pd.Timedelta(hours=1)
        for j in jumps:
            off[j + 1:] += np.timedelta64(span)
        t = pd.DatetimeIndex(t.values + off)
    return t


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def _auto_layout(net):
    """Give nodes without coordinates a simple left-to-right layout."""
    items = net.boundaries + net.nodes
    missing = [n for n in items if n.get("x") is None or n.get("y") is None]
    for i, n in enumerate(missing):
        n["x"] = 80 + (i % 9) * 95
        n["y"] = 90 + (i // 9) * 120


def plot_network(net, ax=None, values=True):
    """Static matplotlib drawing of the network (RC-circuit style)."""
    import matplotlib.pyplot as plt
    _auto_layout(net)
    try:
        nodes, links = net.evaluated()
    except ValueError:
        nodes, links = None, None
    pos = {n["name"]: (n["x"], -n["y"]) for n in net.nodes + net.boundaries}
    xs = [p[0] for p in pos.values()] or [0]
    ys = [p[1] for p in pos.values()] or [0]
    if ax is None:
        w = max(6, min(16, (max(xs) - min(xs)) / 70 + 2))
        h = max(3, min(10, (max(ys) - min(ys)) / 70 + 2))
        _, ax = plt.subplots(figsize=(w, h))
    for i, l in enumerate(net.links):
        if l["a"] in pos and l["b"] in pos:
            (x1, y1), (x2, y2) = pos[l["a"]], pos[l["b"]]
            ctl = Network.is_controlled(l)
            ax.plot([x1, x2], [y1, y2], color="#2f6aa3" if ctl else "#555", lw=1.2, zorder=1,
                    ls="--" if ctl else "-")
            if values:
                v = links.UA.iloc[i] if links is not None else l["UA"]
                txt = f"{v:,.4g}" if isinstance(v, (float, int)) else str(v)
                if ctl:
                    vo = links.UA_open.iloc[i] if links is not None else l["UA_open"]
                    txt += " ⇄ " + (f"{vo:,.4g}" if isinstance(vo, (float, int)) else str(vo))
                ax.text((x1 + x2) / 2, (y1 + y2) / 2, txt, fontsize=7, ha="center", va="center",
                        bbox=dict(fc="white", ec="none", pad=0.5), zorder=2)
    for n in net.nodes:
        x, y = pos[n["name"]]
        cval = nodes.loc[n["name"], "C"] if nodes is not None else 0
        fc = "#2b4c7e" if cval > 0 else "white"
        if cval > 0:   # capacitor to ground
            ax.plot([x, x], [y - 14, y - 26], color="#2b4c7e", lw=1, zorder=1)
            ax.plot([x - 7, x + 7], [y - 26, y - 26], color="#2b4c7e", lw=1.5)
            ax.plot([x - 7, x + 7], [y - 30, y - 30], color="#2b4c7e", lw=1.5)
        sol = Network.has_solar(n.get("solar"))
        ax.scatter([x], [y], s=260, c=fc, edgecolors="#e69500" if sol else "#2b4c7e",
                   linewidths=3 if sol else 1.5, zorder=3)
        ax.text(x, y + 16, n["name"], fontsize=8, ha="center", va="bottom", zorder=4)
    for b in net.boundaries:
        x, y = pos[b["name"]]
        ax.scatter([x], [y], s=300, marker="s", c="#d8e8f5", edgecolors="#3a7cb8", lw=1.5, zorder=3)
        ax.text(x, y + 16, b["name"], fontsize=8, ha="center", va="bottom", color="#3a7cb8")
        ax.text(x, y - 16, str(b["source"]), fontsize=7, ha="center", va="top", color="#3a7cb8")
    ax.set_aspect("equal")
    ax.margins(0.12)
    ax.axis("off")
    ax.set_title(f"{net.title}  ({net.units})", fontsize=10)
    return ax


def edit_network(net, weather=None, height=560, controls=False):
    """Interactive network editor (works in Google Colab; elsewhere it shows a
    JSON box you can copy back with Network.from_json).
    controls=True shows the switched-ventilation fields on links (they also
    appear for any link that already has an open UA)."""
    from IPython.display import HTML, display
    _auto_layout(net)
    uid = uuid.uuid4().hex[:8]
    cb = f"tmm.save_{uid}"
    solar_src = weather.solar_sources if weather is not None else DEFAULT_SOLAR_SOURCES
    lab = UNIT_LABELS

    def solar_info(expr):
        """Short summary of a solar gain expression for the editor panel."""
        if not Network.has_solar(expr):
            return None
        used = Network.solar_names(expr)
        try:
            aper = float(net.solar_gain(expr, {c: 1.0 for c in set(solar_src) | set(used)}))
            txt = f"total aperture {aper:,.4g} {lab[net.units]['area']}"
            if weather is not None and all(c in weather.df for c in used):
                wdf = weather.in_units(net.units)
                g = np.asarray(net.solar_gain(expr, {c: wdf[c].values for c in used}))
                txt += f" · peak gain {g.max():,.0f} {lab[net.units]['Q']}"
            return txt
        except (ValueError, TypeError) as e:
            return "error: " + str(e).split(": ", 1)[-1]
    temp_src = weather.temperature_sources if weather is not None else ["T_out", "T_ground", "T_sky"]

    def save(js):
        try:
            d = json.loads(js)
            net._load(d)
            msgs = net.check(weather)
            vals = {}
            for i, n in enumerate(net.nodes):
                for f in ("C", "T0", "Q", "heat_set", "cool_set"):
                    try:
                        v = net.value(n[f]) if not isinstance(n[f], list) else None
                    except ValueError as e:
                        v = "error: " + str(e).split(": ", 1)[-1]
                    vals[f"n{i}.{f}"] = v
                vals[f"n{i}.solar"] = solar_info(n.get("solar"))
            for i, l in enumerate(net.links):
                for f in ("UA", "UA_open", "open_above", "outdoor_above", "warmer_by", "band"):
                    if f != "UA" and not Network.is_controlled(l):
                        continue
                    try:
                        vals[f"l{i}.{f}"] = net.value(l.get(f))
                    except ValueError as e:
                        vals[f"l{i}.{f}"] = "error: " + str(e).split(": ", 1)[-1]
            res = dict(ok=not msgs, messages=msgs, values=vals)
        except Exception as e:  # noqa: BLE001
            res = dict(ok=False, messages=[f"{type(e).__name__}: {e}"], values={})
        from IPython.display import JSON
        return JSON(res)

    try:
        from google.colab import output  # type: ignore
        output.register_callback(cb, save)
    except ImportError:
        pass
    html = (_EDITOR_HTML
            .replace("__UID__", uid)
            .replace("__HEIGHT__", str(height))
            .replace("__CALLBACK__", cb)
            .replace("__NETWORK__", json.dumps(net.to_dict()))
            .replace("__SOLAR__", json.dumps(solar_src))
            .replace("__TEMPS__", json.dumps(temp_src))
            .replace("__VENT__", "true" if controls else "false")
            .replace("__UNITS__", json.dumps(UNIT_LABELS)))
    display(HTML(html))


_EDITOR_HTML = r"""
<style>
#tme-__UID__ { font-family: system-ui, -apple-system, Segoe UI, sans-serif; font-size: 13px; color: #1d2733; }
#tme-__UID__ .bar { display: flex; flex-wrap: wrap; gap: 6px; align-items: center; margin-bottom: 6px; }
#tme-__UID__ button { border: 1px solid #b8c4d0; background: #f6f8fa; border-radius: 6px; padding: 5px 10px; cursor: pointer; font-size: 13px; color: #1d2733; }
#tme-__UID__ button.on { background: #2b4c7e; color: #fff; border-color: #2b4c7e; }
#tme-__UID__ button:hover { border-color: #2b4c7e; }
#tme-__UID__ .wrap { display: flex; gap: 10px; align-items: flex-start; }
#tme-__UID__ svg { border: 1px solid #c9d3dd; border-radius: 8px; background: #fff; flex: 1 1 auto; min-width: 0; touch-action: none; }
#tme-__UID__ .panel { width: 290px; flex: 0 0 290px; border: 1px solid #c9d3dd; border-radius: 8px; padding: 10px; background: #fbfcfd; overflow-y: auto; max-height: __HEIGHT__px; box-sizing: border-box; }
#tme-__UID__ .panel h4 { margin: 0 0 8px; font-size: 14px; }
#tme-__UID__ .panel label { display: block; margin: 7px 0 2px; font-weight: 600; font-size: 12px; }
#tme-__UID__ .panel input, #tme-__UID__ .panel select, #tme-__UID__ textarea { width: 100%; box-sizing: border-box; padding: 4px 6px; border: 1px solid #b8c4d0; border-radius: 4px; font-size: 12px; font-family: ui-monospace, Consolas, monospace; background: #fff; color: #1d2733; }
#tme-__UID__ .val { font-size: 11px; color: #2b6e3f; min-height: 13px; font-family: ui-monospace, Consolas, monospace; }
#tme-__UID__ .val.err { color: #b3261e; }
#tme-__UID__ .hint { font-size: 11px; color: #5b6b7b; margin-top: 2px; }
#tme-__UID__ .status { margin-left: auto; font-size: 12px; color: #5b6b7b; }
#tme-__UID__ .msgs { margin-top: 6px; font-size: 12px; color: #b3261e; }
#tme-__UID__ .msgs.ok { color: #2b6e3f; }
#tme-__UID__ text { user-select: none; pointer-events: none; }
#tme-__UID__ .help { font-size: 12px; color: #5b6b7b; line-height: 1.5; }
</style>
<div id="tme-__UID__">
  <div class="bar">
    <button data-mode="select" class="on" title="Select and drag (V)">↖ Select / move</button>
    <button data-mode="node" title="Click the canvas to add a node (N)">● Node</button>
    <button data-mode="boundary" title="Click the canvas to add a boundary temperature (B)">■ Boundary</button>
    <button data-mode="link" title="Click two nodes to connect them (L)">— Link</button>
    <button data-act="delete" title="Delete selected (Del)">✕ Delete</button>
    <button data-act="json" title="Show / load the network as JSON">{ } JSON</button>
    <span class="status">…</span>
  </div>
  <div class="wrap">
    <svg viewBox="0 0 900 __HEIGHT__" style="width:100%;aspect-ratio:900/__HEIGHT__;height:auto"></svg>
    <div class="panel"></div>
  </div>
  <div class="msgs"></div>
  <div class="jsonbox" style="display:none;margin-top:6px">
    <textarea rows="8"></textarea>
    <div class="bar" style="margin-top:4px"><button data-act="loadjson">Load JSON into editor</button>
    <span class="hint">Outside Colab: copy this text and use <code>net = Network.from_json('''...''')</code>.</span></div>
  </div>
</div>
<script>
(function(){
const root = document.getElementById('tme-__UID__');
const svg = root.querySelector('svg'), panel = root.querySelector('.panel');
const statusEl = root.querySelector('.status'), msgsEl = root.querySelector('.msgs');
const jsonBox = root.querySelector('.jsonbox'), jsonText = jsonBox.querySelector('textarea');
const NS = 'http://www.w3.org/2000/svg';
const H = __HEIGHT__, W = 900, GRID = 10;
let net = __NETWORK__;
const SOLAR = __SOLAR__, TEMPS = __TEMPS__, UNITS = __UNITS__, VENT = __VENT__;
let mode = 'select', sel = null, linkFrom = null, drag = null, values = {}, saveTimer = null;
const inColab = !!(window.google && google.colab && google.colab.kernel);
net.nodes = net.nodes || []; net.boundaries = net.boundaries || []; net.links = net.links || [];

function U() { return UNITS[net.units] || UNITS.SI; }
function all() { return net.nodes.concat(net.boundaries); }
function byName(n) { return all().find(o => o.name === n); }
function isBoundary(o) { return net.boundaries.includes(o); }
function snap(v) { return Math.round(v / GRID) * GRID; }
function num(v) { const x = parseFloat(v); return isNaN(x) ? null : x; }
function shown(key, raw) {
  const v = values[key];
  if (typeof v === 'number') return Number(v.toPrecision(4)).toLocaleString();
  return (raw === null || raw === undefined) ? '' : String(raw);
}
function hasC(n, i) { const v = values['n' + i + '.C']; return (typeof v === 'number') ? v > 0 : (num(n.C) || 0) > 0 || (typeof n.C === 'string' && n.C.trim() !== '' && num(n.C) === null); }
function isCtrl(l) { return l.UA_open !== null && l.UA_open !== undefined && String(l.UA_open).trim() !== ''; }
function hasSolar(n, i) { const s = (n.solar === null || n.solar === undefined) ? '' : String(n.solar).trim(); return s !== '' && s !== '0'; }
function uniqueName(base) { let k = 1, nm = base + ' ' + k; while (byName(nm)) { k++; nm = base + ' ' + k; } return nm; }
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
function pt(evt) { const p = svg.createSVGPoint(); p.x = evt.clientX; p.y = evt.clientY; const q = p.matrixTransform(svg.getScreenCTM().inverse()); return { x: q.x, y: q.y }; }

function render() {
  svg.innerHTML = '';
  // grid
  const defs = el('defs', {}, svg);
  const pat = el('pattern', { id: 'g__UID__', width: 20, height: 20, patternUnits: 'userSpaceOnUse' }, defs);
  el('path', { d: 'M 20 0 L 0 0 0 20', fill: 'none', stroke: '#eef2f6', 'stroke-width': 1 }, pat);
  el('rect', { x: 0, y: 0, width: W, height: H, fill: 'url(#g__UID__)' }, svg);
  // links
  net.links.forEach((l, i) => {
    const a = byName(l.a), b = byName(l.b); if (!a || !b) return;
    const g = el('g', { style: 'cursor:pointer' }, svg);
    const selected = sel && sel.kind === 'link' && sel.i === i;
    el('line', { x1: a.x, y1: a.y, x2: b.x, y2: b.y, stroke: 'transparent', 'stroke-width': 14 }, g);
    // resistor zig-zag in the middle
    const dx = b.x - a.x, dy = b.y - a.y, L = Math.hypot(dx, dy) || 1, ux = dx / L, uy = dy / L, px = -uy, py = ux;
    const mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2, half = Math.min(18, L / 4);
    let d = `M ${a.x} ${a.y} L ${mx - ux * half} ${my - uy * half}`;
    for (let k = 1; k <= 6; k++) { const t = -half + (2 * half) * k / 6, s = (k === 6) ? 0 : (k % 2 ? 5 : -5); d += ` L ${mx + ux * t + px * s} ${my + uy * t + py * s}`; }
    d += ` L ${b.x} ${b.y}`;
    const ctl = isCtrl(l);
    el('path', { d, fill: 'none', stroke: selected ? '#e0662b' : (ctl ? '#2f6aa3' : '#4a5866'), 'stroke-width': selected ? 2.5 : 1.5, 'stroke-dasharray': ctl ? '6 3' : 'none' }, g);
    let txt = shown('l' + i + '.UA', l.UA) || '?';
    if (ctl) txt += ' ⇄ ' + (shown('l' + i + '.UA_open', l.UA_open) || '?');
    const tx = mx + px * 16, ty = my + py * 16;
    const t = el('text', { x: tx, y: ty + 4, 'text-anchor': 'middle', 'font-size': 11, fill: (values['l' + i + '.UA'] === undefined || typeof values['l' + i + '.UA'] === 'number') ? '#4a5866' : '#b3261e' }, g);
    t.textContent = txt;
    g.addEventListener('pointerdown', e => { e.stopPropagation(); e.preventDefault(); commitFocused(); if (mode === 'select' || mode === 'link') { sel = { kind: 'link', i }; render(); showPanel(); } });
  });
  // pending link
  if (linkFrom) { const a = byName(linkFrom); if (a) el('circle', { cx: a.x, cy: a.y, r: 22, fill: 'none', stroke: '#e0662b', 'stroke-dasharray': '4 3', 'stroke-width': 2 }, svg); }
  // boundaries
  net.boundaries.forEach(b => {
    const g = el('g', { style: 'cursor:move' }, svg);
    const selected = sel && sel.obj === b;
    el('rect', { x: b.x - 15, y: b.y - 15, width: 30, height: 30, rx: 3, fill: '#d8e8f5', stroke: selected ? '#e0662b' : '#3a7cb8', 'stroke-width': selected ? 3 : 1.5 }, g);
    for (let k = -10; k <= 10; k += 5) el('line', { x1: b.x + k - 4, y1: b.y + 10, x2: b.x + k + 4, y2: b.y - 10, stroke: '#9cc0e0', 'stroke-width': 1 }, g);
    const t1 = el('text', { x: b.x, y: b.y - 21, 'text-anchor': 'middle', 'font-size': 12, 'font-weight': 600, fill: '#2f6aa3' }, g); t1.textContent = b.name;
    const t2 = el('text', { x: b.x, y: b.y + 29, 'text-anchor': 'middle', 'font-size': 10, fill: '#3a7cb8' }, g); t2.textContent = String(b.source);
    hook(g, b);
  });
  // nodes
  net.nodes.forEach((n, i) => {
    const g = el('g', { style: 'cursor:move' }, svg);
    const selected = sel && sel.obj === n, cap = hasC(n, i), sol = hasSolar(n, i);
    if (cap) { // capacitor to ground
      el('line', { x1: n.x, y1: n.y + 14, x2: n.x, y2: n.y + 26, stroke: '#2b4c7e', 'stroke-width': 1.2 }, g);
      el('line', { x1: n.x - 8, y1: n.y + 26, x2: n.x + 8, y2: n.y + 26, stroke: '#2b4c7e', 'stroke-width': 2 }, g);
      el('line', { x1: n.x - 8, y1: n.y + 30, x2: n.x + 8, y2: n.y + 30, stroke: '#2b4c7e', 'stroke-width': 2 }, g);
    }
    if (sol) for (let k = 0; k < 8; k++) { const an = k * Math.PI / 4; el('line', { x1: n.x + Math.cos(an) * 16, y1: n.y + Math.sin(an) * 16, x2: n.x + Math.cos(an) * 22, y2: n.y + Math.sin(an) * 22, stroke: '#e69500', 'stroke-width': 2 }, g); }
    el('circle', { cx: n.x, cy: n.y, r: 12, fill: cap ? '#2b4c7e' : '#fff', stroke: selected ? '#e0662b' : '#2b4c7e', 'stroke-width': selected ? 3.5 : 1.8 }, g);
    if (n.heat_set !== null && n.heat_set !== undefined && n.heat_set !== '' || n.cool_set !== null && n.cool_set !== undefined && n.cool_set !== '') {
      const tt = el('text', { x: n.x + 15, y: n.y - 10, 'font-size': 11, fill: '#b3261e' }, g); tt.textContent = '🌡';
    }
    const t = el('text', { x: n.x, y: n.y - 18, 'text-anchor': 'middle', 'font-size': 12, 'font-weight': 600, fill: '#1d2733' }, g); t.textContent = n.name;
    hook(g, n);
  });
}

function hook(g, obj) {
  g.addEventListener('pointerdown', e => {
    e.stopPropagation(); e.preventDefault(); commitFocused();
    if (mode === 'link') {
      if (!linkFrom) { linkFrom = obj.name; sel = { kind: isBoundary(obj) ? 'boundary' : 'node', obj }; render(); showPanel(); return; }
      if (linkFrom !== obj.name) {
        const exists = net.links.findIndex(l => (l.a === linkFrom && l.b === obj.name) || (l.b === linkFrom && l.a === obj.name));
        if (exists >= 0) { sel = { kind: 'link', i: exists }; }
        else { net.links.push({ a: linkFrom, b: obj.name, UA: '', label: '' }); sel = { kind: 'link', i: net.links.length - 1 }; }
        linkFrom = null; render(); showPanel(true); schedule();
      }
      return;
    }
    sel = { kind: isBoundary(obj) ? 'boundary' : 'node', obj };
    const p = pt(e); drag = { obj, ox: obj.x - p.x, oy: obj.y - p.y, moved: false };
    svg.setPointerCapture(e.pointerId);
    render(); showPanel();
  });
}

svg.addEventListener('pointermove', e => {
  if (!drag) return;
  const p = pt(e); drag.obj.x = Math.max(20, Math.min(W - 20, snap(p.x + drag.ox))); drag.obj.y = Math.max(30, Math.min(H - 35, snap(p.y + drag.oy)));
  drag.moved = true; render();
});
svg.addEventListener('pointerup', e => { if (drag && drag.moved) schedule(); drag = null; });
svg.addEventListener('pointerdown', e => {
  e.preventDefault(); commitFocused();
  const p = pt(e);
  if (mode === 'node') {
    const n = { name: uniqueName('Node'), C: 0, T0: '', solar: '', Q: 0, heat_set: '', cool_set: '', x: snap(p.x), y: snap(p.y) };
    net.nodes.push(n); sel = { kind: 'node', obj: n }; render(); showPanel(true); schedule();
  } else if (mode === 'boundary') {
    const b = { name: uniqueName('Outdoor'), source: 'T_out', x: snap(p.x), y: snap(p.y) };
    net.boundaries.push(b); sel = { kind: 'boundary', obj: b }; render(); showPanel(true); schedule();
  } else { sel = null; linkFrom = null; render(); showPanel(); }
});

function field(label, key, value, opts) {
  opts = opts || {};
  const wrap = document.createElement('div');
  const lb = document.createElement('label'); lb.textContent = label; wrap.appendChild(lb);
  const inp = document.createElement(opts.multiline ? 'textarea' : 'input');
  if (opts.multiline) { inp.rows = 3; inp.style.resize = 'vertical'; }
  inp.value = (value === null || value === undefined) ? '' : value;
  if (opts.list) { const id = 'dl' + Math.random().toString(36).slice(2); inp.setAttribute('list', id); const dl = document.createElement('datalist'); dl.id = id; opts.list.forEach(s => { const o = document.createElement('option'); o.value = s; dl.appendChild(o); }); wrap.appendChild(dl); }
  if (opts.placeholder) inp.placeholder = opts.placeholder;
  let last = inp.value;
  const commit = () => { if (inp.value === last) return; last = inp.value; opts.set(inp.value); render(); schedule(); };
  inp.addEventListener('change', commit);
  inp.addEventListener('keydown', e => { if (e.key === 'Enter') { commit(); inp.blur(); } });
  inp._commit = commit;
  wrap.appendChild(inp);
  if (opts.vkey) { const v = document.createElement('div'); v.className = 'val'; v.dataset.vkey = opts.vkey; wrap.appendChild(v); }
  if (opts.hint) { const h = document.createElement('div'); h.className = 'hint'; h.innerHTML = opts.hint; wrap.appendChild(h); }
  panel.appendChild(wrap);
  return inp;
}
function coerce(v) { const s = String(v).trim(); if (s === '') return ''; const x = Number(s); return isNaN(x) ? s : x; }
function rename(obj, nm) {
  nm = nm.trim(); if (!nm || nm === obj.name) return;
  if (byName(nm)) { alert('Name "' + nm + '" is already used.'); return; }
  net.links.forEach(l => { if (l.a === obj.name) l.a = nm; if (l.b === obj.name) l.b = nm; });
  obj.name = nm;
}

function formulaHelp() {
  const u = U();
  const rows = [
    ['conduction(A, L, material="concrete")', 'k·A/L through a layer. <b>A</b> = face area the heat passes through; <b>L</b> = layer thickness in the direction of heat flow. Or give <code>k=</code> instead of a material.'],
    ['mass(A, L, "concrete")', 'ρ·c·A·L, heat stored in a slab of area <b>A</b> and thickness <b>L</b>. Or give <code>rho=</code>, <code>cp=</code>.'],
    ['air_mass(V, multiplier=5)', 'Heat stored in room air of volume <b>V</b>; the multiplier adds furniture and contents.'],
    ['convection(A, h)', 'h·A, surface film. <b>A</b> = surface area, <b>h</b> = film coefficient (default ' + (net.units === 'IP' ? '0.53 Btu/hr·ft²·°F' : '3 W/m²K') + ', still indoor air).'],
    ['radiation(A, e1, e2)', 'Long-wave (infrared) exchange between two facing surfaces. <b>A</b> = their area; <b>e1</b>, <b>e2</b> = their <i>infrared emissivities</i> (≈0.84–0.9 glass, masonry, paint; ≈0.1 low-e or metal). Optional <code>T_mean</code> (default ' + (net.units === 'IP' ? '68 °F' : '20 °C') + '). <i>Solar absorptance is not used here; it goes in a node\'s solar gain.</i>'],
    ['window(A, U)', 'U·A for glazing or any assembly with a known U-value.'],
    ['ventilation(V, ach)', 'Air exchange for room volume <b>V</b> at <b>ach</b> air changes per hour.'],
    ['airflow(flow)', 'Air-flow conductance; flow in ' + (net.units === 'IP' ? 'cfm' : 'm³/s') + '.'],
    ['series(ua1, ua2, …)', 'Conductances one after another: 1/Σ(1/UA). Use <code>+</code> for side-by-side paths.'],
    ['aperture(A, tau, alpha, F) * I_S', '<i>Solar gain only.</i> A·τ·α·F times a solar source: <b>A</b> = glazing area, <b>tau</b> = glazing solar transmittance, <b>alpha</b> = solar absorptance of the receiving surface, <b>F</b> = fraction of that window\'s sun reaching this node (default 1). Add one term per window: <code>aperture(6, 0.75, 0.9)*I_S + aperture(3, 0.75, 0.9)*I_E</code>.'],
  ];
  return '<div class="help">' + rows.map(r => '<p style="margin:6px 0"><code>' + r[0] + '</code><br>' + r[1] + '</p>').join('') + '</div>';
}

function showPanel(focus) {
  panel.innerHTML = '';
  const u = U();
  const h = document.createElement('h4'); panel.appendChild(h);
  let first = null;
  if (!sel) {
    h.textContent = 'Network';
    first = field('Title', null, net.title, { set: v => net.title = v });
    const lb = document.createElement('label'); lb.textContent = 'Units'; panel.appendChild(lb);
    const s = document.createElement('select'); ['IP', 'SI'].forEach(x => { const o = document.createElement('option'); o.value = x; o.textContent = x === 'IP' ? 'IP (°F, Btu/hr·°F, Btu/°F, ft²)' : 'SI (°C, W/K, J/K, m²)'; s.appendChild(o); });
    s.value = net.units; s.addEventListener('change', () => { net.units = s.value; schedule(); showPanel(); }); panel.appendChild(s);
    const p = document.createElement('div'); p.className = 'help';
    p.innerHTML = '<p><b>How to draw:</b> pick <b>● Node</b> and click the canvas to place thermal nodes (glazing, surfaces, wall layers, room air…). ' +
      'Pick <b>■ Boundary</b> for fixed temperatures (outdoor air, ground, sky, or a constant). Pick <b>— Link</b> and click two nodes to connect them with a conductance UA. ' +
      'Use <b>Select</b> to drag things and edit their values here.</p>' +
      '<p>Filled nodes have capacitance (drawn with a capacitor); hollow ones are massless. Orange rays = node absorbs sun.' + (VENT ? ' Dashed blue links are <b>controlled</b> (closed ⇄ open UA); select a link to set one up.' : '') + '</p>' +
      '<p><b>Any value may be a formula.</b> Lengths in ' + u.length + ', areas in ' + u.area + ', volumes in ' + u.volume + '.</p>' + formulaHelp() +
      '<p>Changes are sent to Python automatically.</p>';
    panel.appendChild(p);
  } else if (sel.kind === 'node') {
    const n = sel.obj, i = net.nodes.indexOf(n);
    h.textContent = 'Node';
    first = field('Name', null, n.name, { set: v => rename(n, v) });
    field('Capacitance C [' + u.C + ']', 'C', n.C, { set: v => n.C = coerce(v), vkey: 'n' + i + '.C', hint: '0 for a massless node (air gap, glass, surface). Formula: <code>mass(A, L, "concrete")</code> (A = area, L = thickness) or <code>air_mass(V, multiplier)</code>.' });
    field('Initial temperature T0 [' + u.T + ']', 'T0', n.T0, { set: v => n.T0 = coerce(v), vkey: 'n' + i + '.T0', placeholder: 'blank = first outdoor temp' });
    field('Solar gain [' + u.Q + ']', 'solar', n.solar, { set: v => n.solar = String(v).replace(/\s*\n\s*/g, ' ').trim(), vkey: 'n' + i + '.solar', placeholder: 'blank = no sun', multiline: true, hint: 'One term per window: <code>aperture(A, τ, α, F)*I_S</code>, added with <code>+</code>. A = glazing area, τ = glazing solar transmittance, α = <b>solar absorptance</b> of this surface (≈0.9 dark, 0.3–0.5 light), F = fraction of that window\'s sun landing on this node. E.g. <code>aperture(6, 0.75, 0.9)*I_S + aperture(3, 0.75, 0.9, 0.5)*I_E</code>. Sources: ' + SOLAR.join(', ') + '.' });
    field('Internal gain Q [' + u.Q + ']', 'Q', Array.isArray(n.Q) ? JSON.stringify(n.Q) : n.Q, { set: v => { const s = String(v).trim(); n.Q = s.startsWith('[') ? JSON.parse(s) : coerce(v); }, vkey: 'n' + i + '.Q', hint: 'Constant, or a list of 24 hourly values.' });
    field('Heating set point [' + u.T + ']', 'heat_set', n.heat_set, { set: v => n.heat_set = coerce(v), vkey: 'n' + i + '.heat_set', placeholder: 'blank = free-floating' });
    field('Cooling set point [' + u.T + ']', 'cool_set', n.cool_set, { set: v => n.cool_set = coerce(v), vkey: 'n' + i + '.cool_set', placeholder: 'blank = free-floating' });
  } else if (sel.kind === 'boundary') {
    const b = sel.obj;
    h.textContent = 'Boundary temperature';
    first = field('Name', null, b.name, { set: v => rename(b, v) });
    field('Source', 'source', b.source, { set: v => b.source = coerce(v), list: TEMPS.concat(['20', '68']), hint: 'Weather column (T_out outdoor air, T_ground, T_sky) or a constant temperature in ' + u.T + '.' });
  } else if (sel.kind === 'link') {
    const l = net.links[sel.i];
    h.textContent = 'Link  ' + l.a + ' — ' + l.b;
    first = field('Conductance UA [' + u.UA + ']', 'UA', l.UA, { set: v => l.UA = coerce(v), vkey: 'l' + sel.i + '.UA', hint: 'A number or a formula (see below).' });
    field('Label (optional)', 'label', l.label, { set: v => l.label = v });
    if (VENT || isCtrl(l)) ventFields(l, sel.i, u);
    const fh = document.createElement('div'); fh.innerHTML = formulaHelp(); panel.appendChild(fh);
  }
  updateValues();
  if (focus && first) setTimeout(() => { first.focus(); first.select(); }, 0);
}

function ventFields(l, i, u) {
  const hd = document.createElement('h4'); hd.style.marginTop = '14px'; hd.textContent = 'Controlled ventilation (optional)'; panel.appendChild(hd);
  const p = document.createElement('div'); p.className = 'help';
  p.innerHTML = 'For windows or vents that open by themselves. Give the <b>UA when open</b>; the UA above is then the <b>closed</b> value (e.g. infiltration). The link opens when the sensor node is above its set point, <b>and</b> the outdoor air is above its set point, <b>and</b> the sensor is warmer than outdoors by the given difference; it closes when any of them falls a dead band below.';
  panel.appendChild(p);
  const names = net.nodes.map(n => n.name);
  const defSensor = names.includes(l.a) ? l.a : l.b;
  field('UA when open [' + u.UA + ']', 'UA_open', l.UA_open, { set: v => { l.UA_open = coerce(v); if (isCtrl(l) && (l.band === undefined || l.band === null || l.band === '')) { l.band = 1; l.warmer_by = (net.units === 'IP') ? 3.6 : 2; } setTimeout(showPanel, 0); }, vkey: 'l' + i + '.UA_open', placeholder: 'blank = not controlled', hint: 'e.g. <code>ventilation(54, 8)</code> = 8 air changes per hour in a 54 ' + u.volume + ' room.' });
  if (!isCtrl(l)) return;
  field('Sensor node', 'sensor', l.sensor, { set: v => l.sensor = String(v).trim(), list: names, placeholder: 'blank = ' + defSensor, hint: 'The node whose temperature opens the window (usually the room air).' });
  field('Open when sensor is above [' + u.T + ']', 'open_above', l.open_above, { set: v => l.open_above = coerce(v), vkey: 'l' + i + '.open_above' });
  field('… and outdoor air is above [' + u.T + ']', 'outdoor_above', l.outdoor_above, { set: v => l.outdoor_above = coerce(v), vkey: 'l' + i + '.outdoor_above', placeholder: 'blank = any outdoor temperature' });
  field('… and sensor is warmer than outdoor by [' + u.T + ']', 'warmer_by', l.warmer_by, { set: v => l.warmer_by = coerce(v), vkey: 'l' + i + '.warmer_by', placeholder: 'blank = no indoor/outdoor test', hint: 'Keeps the window shut when it is hotter outside (e.g. 2 °C = 3.6 °F).' });
  field('Dead band [' + u.T + ']', 'band', (l.band === undefined || l.band === null) ? 1 : l.band, { set: v => l.band = coerce(v), vkey: 'l' + i + '.band', hint: 'Closes when any of the three tests falls this far below its set point, so it does not flutter open and shut.' });
}

function commitFocused() { const a = document.activeElement; if (a && a._commit && root.contains(a)) { a._commit(); a.blur(); } }
function updateValues() {
  panel.querySelectorAll('.val').forEach(d => {
    const v = values[d.dataset.vkey];
    if (v === undefined || v === null) { d.textContent = ''; d.className = 'val'; }
    else if (typeof v === 'number') { d.textContent = '= ' + Number(v.toPrecision(6)).toLocaleString(); d.className = 'val'; }
    else { d.textContent = String(v); d.className = String(v).startsWith('error') ? 'val err' : 'val'; }
  });
}

function schedule() { statusEl.textContent = 'editing…'; clearTimeout(saveTimer); saveTimer = setTimeout(save, 350); jsonText.value = JSON.stringify(net, null, 1); }
async function save() {
  if (!inColab) { statusEl.textContent = 'not in Colab — use { } JSON to copy'; return; }
  statusEl.textContent = 'saving…';
  try {
    const r = await google.colab.kernel.invokeFunction('__CALLBACK__', [JSON.stringify(net)], {});
    const d = r.data['application/json'];
    values = d.values || {};
    statusEl.textContent = '✓ saved to Python';
    if (d.messages && d.messages.length) { msgsEl.className = 'msgs'; msgsEl.innerHTML = '⚠ ' + d.messages.map(m => m.replace(/</g, '&lt;')).join('<br>⚠ '); }
    else { msgsEl.className = 'msgs ok'; msgsEl.textContent = '✓ network is complete and solvable'; }
    render(); updateValues();
  } catch (err) { statusEl.textContent = '⚠ could not reach Python (re-run this cell)'; }
}

root.querySelectorAll('button[data-mode]').forEach(b => b.addEventListener('click', () => {
  mode = b.dataset.mode; linkFrom = null;
  root.querySelectorAll('button[data-mode]').forEach(x => x.classList.toggle('on', x === b));
  svg.style.cursor = (mode === 'node' || mode === 'boundary') ? 'crosshair' : 'default';
  render();
}));
function del() {
  if (!sel) return;
  if (sel.kind === 'link') net.links.splice(sel.i, 1);
  else { const nm = sel.obj.name; net.nodes = net.nodes.filter(x => x !== sel.obj); net.boundaries = net.boundaries.filter(x => x !== sel.obj); net.links = net.links.filter(l => l.a !== nm && l.b !== nm); }
  sel = null; render(); showPanel(); schedule();
}
root.querySelector('button[data-act=delete]').addEventListener('click', del);
root.querySelector('button[data-act=json]').addEventListener('click', () => { jsonText.value = JSON.stringify(net, null, 1); jsonBox.style.display = jsonBox.style.display === 'none' ? 'block' : 'none'; });
root.querySelector('button[data-act=loadjson]').addEventListener('click', () => {
  try { const d = JSON.parse(jsonText.value); net = Object.assign({ units: net.units, title: '', nodes: [], boundaries: [], links: [] }, d); sel = null; render(); showPanel(); schedule(); }
  catch (e) { alert('Invalid JSON: ' + e.message); }
});
root.addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA' || e.target.tagName === 'SELECT') return;
  if (e.key === 'Delete' || e.key === 'Backspace') { del(); e.preventDefault(); }
  const m = { v: 'select', n: 'node', b: 'boundary', l: 'link' }[e.key.toLowerCase()];
  if (m) root.querySelector('button[data-mode=' + m + ']').click();
});
root.tabIndex = 0;
render(); showPanel(); save();
})();
</script>
"""


# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------
def example_trombe_wall():
    """The CNRS / Odeillo Trombe wall network from the spreadsheet (IP units).
    Solar source I_S: south vertical irradiance."""
    net = Network("IP", "Trombe wall (Braham 1979 / spreadsheet)")
    net.add_boundary("Outdoor", "T_out", x=50, y=180)
    y = 180
    for nm, x, yy in (("Glass 1", 140, y), ("Glass 2", 230, y), ("Air Channel", 320, 320),
                      ("Outer Surface", 400, y)):
        net.add_node(nm, C=0, T0=70, x=x, y=yy)
    for i in range(5):
        net.add_node(f"Wall {i + 1}", C=6150, T0=70, x=480 + i * 75, y=y)
    net.add_node("Inner Surface", C=0, T0=70, x=850, y=y)
    net.add_node("Room Air", C=7500, T0=70, x=650, y=420)
    net.set("Outer Surface", solar="480*I_S")   # the spreadsheet's AGF x S
    net.add_boundary("Outdoor (envelope)", "T_out", x=380, y=480)
    links = [("Outdoor", "Glass 1", 1080, "glass + film"),
             ("Glass 1", "Glass 2", 774, "air gap"),
             ("Glass 2", "Air Channel", 851, "convection"),
             ("Glass 2", "Outer Surface", 263, "radiation"),
             ("Air Channel", "Outer Surface", 851, "convection"),
             ("Air Channel", "Room Air", 748, "thermocirculation"),
             ("Outer Surface", "Wall 1", 2644, "conduction")]
    links += [(f"Wall {i}", f"Wall {i + 1}", 1322, "conduction") for i in range(1, 5)]
    links += [("Wall 5", "Inner Surface", 2644, "conduction"),
              ("Inner Surface", "Room Air", 851, "convection + radiation"),
              ("Room Air", "Outdoor (envelope)", 922, "envelope + infiltration")]
    for a, b, ua, lab in links:
        net.link(a, b, ua, lab)
    return net


def example_simple_room():
    """A small SI example: a room with south and west windows, a concrete floor
    slab on insulation over the ground, and lightweight walls -- built with formulas."""
    net = Network("SI", "Simple sunspace room")
    net.add_boundary("Outdoor", "T_out", x=80, y=120)
    net.add_boundary("Ground", "T_ground", x=620, y=470)
    net.add_node("Room Air", C="air_mass(4*5*2.7, multiplier=5)", T0=20, Q=[150] * 8 + [400] * 10 + [250] * 6,
                 x=300, y=180)
    net.add_node("Slab Surface", C=0, x=450, y=280,
                 solar="aperture(6, 0.75, 0.9, 0.8)*I_S + aperture(2, 0.75, 0.9, 0.8)*I_W")
    net.add_node("Slab", C="mass(20, 0.10, 'concrete')", x=450, y=370)
    net.link("Outdoor", "Room Air", "window(6 + 2, 1.4) + conduction(A=40, L=0.15, material='mineral wool')"
             " + ventilation(54, 0.5)", "window + walls + infiltration")
    net.link("Room Air", "Slab Surface", "convection(20, 3) + radiation(20, 0.9, 0.9)*0.5", "film")
    net.link("Slab Surface", "Slab", "conduction(A=20, L=0.05, material='concrete')", "half slab")
    net.link("Slab", "Ground", "series(conduction(A=20, L=0.05, material='concrete'), "
             "conduction(A=20, L=0.10, material='XPS'), conduction(A=20, L=0.5, material='soil (moist)'))",
             "half slab + insulation + soil")
    return net


def example_vent_room():
    """The simple SI room with temperature-responsive windows: the window link is
    infiltration only when closed and 15 air changes when open (windows wide open); it opens when the
    room air passes 24 °C, it is warmer than 16 °C outside, and the room is at least
    2 °C warmer than outside."""
    net = example_simple_room()
    net.title = "Simple room with operable windows"
    net.links[0].update(UA="window(6 + 2, 1.4) + conduction(A=40, L=0.15, material='mineral wool')",
                        label="window + walls")
    net.add_boundary("Outdoor (vent)", "T_out", x=80, y=300)
    net.link("Outdoor (vent)", "Room Air", "ventilation(54, 0.5)", "operable windows",
             UA_open="ventilation(54, 15)", sensor="Room Air", open_above=24, outdoor_above=16, warmer_by=2, band=1)
    return net


def example_excel_weather(T_F, S_btu):
    """Wrap spreadsheet-style hourly series (outdoor temp in F and solar on the
    collector in Btu/hr-ft2) as a Weather with columns T_out and I_S."""
    n = len(T_F)
    idx = pd.Timestamp(2001, 1, 1) + pd.to_timedelta(np.arange(1, n + 1), unit="h")
    df = pd.DataFrame({"T_out": f_to_c(np.asarray(T_F, float)),
                       "I_S": np.asarray(S_btu, float) * WM2_PER_BTUHFT2}, index=idx)
    return Weather(df, name="Spreadsheet weather (Trombe example)")
