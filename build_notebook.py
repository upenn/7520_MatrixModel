"""Build ThermalMatrixModel.ipynb (Colab) from thermal_matrix.py.

Run:  python build_notebook.py ["short note about what changed"]
      python build_notebook.py --vent ["note"]   -> ThermalMatrixModel_Vent.ipynb, the variant
                                                   with controlled (switched) ventilation links
Each build is stamped with its date and time; the stamp appears in the notebook title,
the setup cell's output, saved network .json files and results .csv files. The note
(if given) is added to CHANGELOG.md, which is also shown at the end of the notebook.
The engine source is embedded in the notebook's setup cell, so the notebook is
self-contained: students only need the .ipynb (and optionally their own .epw).
"""
import datetime
import json
import os
import re
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))


def excel_week():
    """First week of the spreadsheet's weather (outdoor F, solar Btu/hr-ft2)."""
    import openpyxl
    warnings.simplefilter("ignore")
    # MATRIX_XLSX can point at a copy, e.g. while the workbook is open in Excel
    path = os.environ.get("MATRIX_XLSX") or os.path.join(HERE, "MatrixModel_TrombeWall_25.xlsx")
    wb = openpyxl.load_workbook(path, data_only=True)
    w = wb["WeatherData"]
    T = [round(w.cell(r, 6).value, 4) for r in range(3, 3 + 168)]
    S = [round(w.cell(r, 7).value, 4) for r in range(3, 3 + 168)]
    return T, S


def md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src.strip("\n")}


def code(src, form=True):
    meta = {"cellView": "form"} if form else {}
    return {"cell_type": "code", "metadata": meta, "execution_count": None,
            "outputs": [], "source": src.strip("\n")}


CHANGELOG = os.path.join(HERE, "CHANGELOG.md")


def changelog_entries(version, note, vent=False, tag=None):
    """Existing CHANGELOG.md entries with this build's entry first."""
    entries = []
    if os.path.exists(CHANGELOG):
        entries = [ln for ln in open(CHANGELOG, encoding="utf-8").read().splitlines() if ln.startswith("- ")]
    if tag is None:
        tag = "[ventilation notebook] " if vent else ""
    return [f"- **{version}** - {tag}{note or '(rebuilt, no note)'}"] + entries


def write_changelog(entries):
    head = ("# Thermal Matrix Model: version history\n\n"
            "Newest first. Add a note with `python build_notebook.py \"what changed\"`.\n"
            "Entries marked [ventilation notebook] are builds of ThermalMatrixModel_Vent.ipynb "
            "(`python build_notebook.py --vent \"what changed\"`); both notebooks share the engine.\n"
            "Entries marked [property calculator] are builds of PropertyCalculator.ipynb "
            "(`python build_notebook.py --calc \"what changed\"`).\n\n")
    with open(CHANGELOG, "w", encoding="utf-8") as f:
        f.write(head + "\n".join(entries) + "\n")


DRIVE_HELPER = r'''

def drive_path(folder="ThermalMatrixModel"):
    """Folder in the student's Google Drive (mounted on first use); a local folder outside Colab."""
    if IN_COLAB:
        if not os.path.isdir("/content/drive/MyDrive"):
            from google.colab import drive
            drive.mount("/content/drive")
        path = os.path.join("/content/drive/MyDrive", folder.strip("/ ") or "ThermalMatrixModel")
    else:
        path = folder.strip("/ ") or "ThermalMatrixModel"
        print(f"(not in Colab: using the local folder '{path}' instead of Google Drive)")
    os.makedirs(path, exist_ok=True)
    return path
'''


def build(note="", vent=False):
    version = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    entries = changelog_entries(version, note, vent)
    engine = open(os.path.join(HERE, "thermal_matrix.py"), encoding="utf-8").read()
    engine, n = re.subn(r'^MODEL_VERSION = .*$', f'MODEL_VERSION = "{version}"', engine, count=1, flags=re.M)
    assert n == 1, "MODEL_VERSION line not found in thermal_matrix.py"
    T, S = excel_week()
    cells = []

    cells.append(md(f"""
# 🌡️ Thermal Matrix Model{" + controlled ventilation" if vent else ""}
**Version {version}**{(" - " + note) if note else ""}  (version history at the end)
""" + r"""
A lumped-parameter (resistance–capacitance) thermal network simulator for buildings,
adapted from the spreadsheet *Matrix Model* and W. Braham, *A Brief Analysis of the
Performance of the CNRS Passive Solar House Collector: The Trombe Wall* (1979).

**The method.** The building is reduced to a network of **nodes** (glass, surfaces, layers
of a wall, room air…) joined by **conductances** $UA_{ij}$. Each node may store heat
(capacitance $C_i$), absorb sun, receive internal gains, and connect to **boundary
temperatures** (outdoor air, ground, sky). An energy balance about each node, written as a
backward (implicit) difference over the time step $\Delta t$:

$$\Big(\tfrac{C_i}{\Delta t} + \sum_j UA_{ij} + UA_{io}\Big)\,T_i^{t} \;-\; \sum_j UA_{ij}\,T_j^{t}
\;=\; \tfrac{C_i}{\Delta t}\,T_i^{t-\Delta t} \;+\; UA_{io}\,T_o^{t} \;+\; A_i\,I^{t} \;+\; Q_i^{t}$$

For all nodes together this is $\;M_1\,T^{t} = M_2^{t}\;\Rightarrow\; T^{t} = M_1^{-1} M_2^{t}$.
$M_1$ depends only on the network, so it is inverted once; each hour $M_2$ is rebuilt from the
previous temperatures and the weather. The implicit scheme is stable for any time step.
""" + (r"""
**This version adds controlled ventilation.** A link can have two conductances, *closed* and
*open*, like a temperature-responsive window that opens when the room is warm and it isn't
too cold outside. Each open/closed combination has its own $M_1$, inverted once the first time
it occurs. See *Controlled ventilation* in step 2.
""" if vent else "") + r"""

### How to use this notebook
1. **Setup** – run the first cell once.
2. **Weather** – load an EnergyPlus `.epw` file (upload, URL, or the Philadelphia default).
3. **Network** – start from an example or a blank page, then **draw** your network: place
   nodes, boundaries and links, and type their properties (numbers or formulas such as
   `conduction(A=48, L=0.12, material="concrete")`).
4. **Simulate** – choose the dates.
5. **Results** – plots, tables, loads, and CSV download.

Run cells with ▶ (or *Shift+Enter*). Grey form fields on the right of each cell can be edited
without touching code. *Runtime ▸ Run all* runs the Trombe-wall example end to end.
"""))

    cells.append(code(
        "#@title 0 ▸ Setup: load the thermal-matrix engine (run once) { display-mode: \"form\" }\n"
        "# The full engine source follows. You don't need to read it, but you can.\n\n"
        + engine
        + "\n\n# --- first week of the spreadsheet's weather (outdoor °F, solar on collector Btu/hr·ft²)\n"
        + f"EXCEL_WEEK_T = {json.dumps(T)}\n"
        + f"EXCEL_WEEK_S = {json.dumps(S)}\n\n"
        + "import matplotlib.pyplot as plt\n"
        + "from IPython.display import display, Markdown\n"
        + "pd.set_option('display.width', 200)\n"
        + "try:\n"
        + "    from google.colab import files as _colab_files\n"
        + "    IN_COLAB = True\n"
        + "except ImportError:\n"
        + "    IN_COLAB = False\n"
        + DRIVE_HELPER
        + "print(f'✓ Thermal matrix engine loaded - model version {MODEL_VERSION}' + ('' if IN_COLAB else ' (not running in Colab: uploads and the drawing sync are disabled)'))\n"
    ))

    cells.append(md(r"""
## 1 ▸ Weather

Weather comes from an **EnergyPlus weather file** (`.epw`, or the `.zip` it usually comes in).
Free files for almost any location: <https://climate.onebuilding.org> (pick a *TMYx* file) or
<https://energyplus.net/weather>.

The file is converted into these hourly columns, which your network can refer to:

| column | meaning |
|---|---|
| `T_out` | outdoor dry-bulb temperature |
| `T_ground` | monthly ground temperature (from the file header, 0.5 m depth) |
| `T_sky` | effective sky temperature (from horizontal infrared) — for night-sky radiation |
| `I_S`, `I_E`, `I_W`, `I_N` | total solar on a **vertical** surface facing S, E, W, N |
| `I_H` | total solar on a **horizontal** surface |

Surface irradiance = beam·cos θ + isotropic diffuse + ground-reflected (albedo 0.2).
"""))

    cells.append(code(r'''
#@title 1a ▸ Load weather { display-mode: "form" }
weather_source = "Philadelphia (download TMYx)" #@param ["Philadelphia (download TMYx)", "Upload my own .epw or .zip", "Download from the URL below", "Spreadsheet week (Trombe example)"]
epw_url = "https://climate.onebuilding.org/WMO_Region_4_North_and_Central_America/USA_United_States_of_America/PA_Pennsylvania/USA_PA_Philadelphia.Intl.AP.724080_TMYx.zip" #@param {type:"string"}

PHILADELPHIA = ("https://climate.onebuilding.org/WMO_Region_4_North_and_Central_America/"
                "USA_United_States_of_America/PA_Pennsylvania/USA_PA_Philadelphia.Intl.AP.724080_TMYx.zip")
if weather_source.startswith("Philadelphia"):
    weather = read_epw(PHILADELPHIA)
elif weather_source.startswith("Upload"):
    if not IN_COLAB:
        raise RuntimeError("Uploading needs Colab. Locally use: weather = read_epw('path/to/file.epw')")
    up = _colab_files.upload()
    weather = read_epw(next(iter(up.values())))
elif weather_source.startswith("Download"):
    weather = read_epw(epw_url)
else:
    weather = example_excel_weather(EXCEL_WEEK_T, EXCEL_WEEK_S)

print(weather)
if weather.lat is not None:
    print(f"latitude {weather.lat}, longitude {weather.lon}, time zone {weather.tz:+g}, elevation {weather.elevation} m")
weather.plot(units="IP" if weather_source.startswith("Spreadsheet") else "SI");
plt.show()
'''))

    cells.append(code(r'''
#@title 1b ▸ (optional) Add solar on a surface with another tilt / orientation { display-mode: "form" }
#@markdown Creates a column `I_<name>` you can use as a node's solar source.
#@markdown Tilt: 0 = horizontal, 90 = vertical. Azimuth: the compass direction the surface faces (0 N, 90 E, 180 S, 270 W).
surface_name = "SE_roof" #@param {type:"string"}
tilt = 30 #@param {type:"slider", min:0, max:180, step:5}
azimuth = 135 #@param {type:"slider", min:0, max:355, step:5}
albedo = 0.2 #@param {type:"number"}
weather.add_surface(surface_name, tilt=tilt, azimuth=azimuth, albedo=albedo)
print("solar sources now available:", ", ".join(weather.solar_sources))
'''))

    cells.append(md(r"""
## 2 ▸ Describe the thermal network

A network has three kinds of things:

* **Nodes** (●) – a temperature you want to know. Give it a **capacitance C** if it stores heat
  (a wall layer, the room air and furniture); use C = 0 for massless nodes (a pane of glass, a
  surface, an air gap). A node can **absorb sun** through one or more windows: its
  *solar gain* is a formula with one term per window, each an aperture times a solar source
  (see below). A node can have **internal gains Q** (constant or 24 hourly values) and
  **heating / cooling set points** (an ideal thermostat that reports the load).
* **Boundaries** (■) – temperatures imposed from outside: `T_out`, `T_ground`, `T_sky`, or a
  constant number (e.g. a neighbour held at 20 °C).
* **Links** (—/\/\/—) – a conductance **UA** between two nodes or between a node and a boundary.

**Units.** *IP*: °F, Btu/hr·°F, Btu/°F, ft², Btu/hr — exactly as the spreadsheet.
*SI*: °C, W/K, J/K, m², W. Choose in the editor's *Network* panel (click empty canvas).

**Formulas** you can type in any value box. Geometry is in the network's units:
ft, ft², ft³ for IP; m, m², m³ for SI.

| formula | gives | arguments |
|---|---|---|
| `conduction(A, L, material="concrete")` or `conduction(A, L, k=1.75)` | k·A/L, conductance through a layer | **A** = face area the heat passes through (e.g. wall height × width); **L** = layer **thickness** in the direction of heat flow (not height or depth of the room); **k** = conductivity, if you don't name a material |
| `mass(A, L, "concrete")` or `mass(A, L, rho=…, cp=…)` | ρ·c·A·L, heat capacity of a slab | **A** = area, **L** = thickness (so A·L is its volume); **rho** = density, **cp** = specific heat |
| `air_mass(V, multiplier=5)` | heat capacity of room air | **V** = room volume; **multiplier** adds furniture and contents (typically 3–10) |
| `convection(A, h)` | h·A, surface film conductance | **A** = surface area; **h** = film coefficient (default 3 W/m²K = 0.53 Btu/hr·ft²·°F, still indoor air; ~15–25 W/m²K outdoors in wind) |
| `radiation(A, e1, e2, T_mean)` | 4σT³A / (1/e1 + 1/e2 − 1), long-wave radiation between two facing surfaces (the paper's equivalent radiation coefficient) | **A** = area of the facing surfaces; **e1, e2** = their **infrared emissivities** (≈0.84–0.9 for glass, masonry, paint; ≈0.1 or less for low-e coatings and bare metal); **T_mean** = typical temperature for linearising (default 20 °C / 68 °F) |
| `window(A, U)` | U·A | **A** = area; **U** = U-value of the assembly |
| `ventilation(V, ach)` | air-exchange conductance | **V** = room volume; **ach** = air changes per hour |
| `airflow(flow)` | air-flow conductance | **flow** in m³/s (SI) or cfm (IP) |
| `series(ua1, ua2, …)` | 1 / Σ(1/UA), conductances one after another | any conductances or formulas; add them with `+` for side-by-side paths |

**Solar gain (windows in several directions).** A node's *solar gain* box takes one term per
window, added together. Each term is an aperture times a **solar source**, the irradiance on a
surface facing that direction (`I_S`, `I_E`, `I_W`, `I_N` vertical, `I_H` horizontal, plus any
you add in step 1b):

| formula | gives | arguments |
|---|---|---|
| `aperture(A, tau, alpha, F) * I_S` | A·τ·α·F × irradiance, absorbed solar gain from one window | **A** = glazing area; **tau** = solar transmittance of the glazing; **alpha** = **solar absorptance** of the receiving surface (≈0.9 dark, 0.3–0.5 light); **F** = fraction of that window's sun that lands on this node (default 1) |

For example, a floor lit by a south window and an east window, with half the east sun landing
on it:
`aperture(6, 0.75, 0.9)*I_S + aperture(3, 0.75, 0.9, 0.5)*I_E`.
A plain number works too: `480*I_S` is the spreadsheet's AGF × S. The editor shows the total
aperture and the peak gain for the loaded weather.

**Emissivity is not solar absorptance.** `radiation()` is about infrared heat exchanged between
surfaces. Solar absorptance belongs only in the solar gain.

Materials: `concrete, brick, adobe, rammed earth, stone, water, glass, gypsum board, plaster,
softwood, hardwood, CLT, mineral wool, EPS, XPS, straw bale, steel, soil (moist), air`
(typical values — replace with product data when you have it).

**Property Calculator.** To look up materials (generic values, or your own table of ASHRAE
*Fundamentals* Table 1 values or product data) and work out each layer's conductance and
capacitance, split into segments for thick mass walls, use the companion notebook
[Property Calculator](https://colab.research.google.com/github/upenn/7520_MatrixModel/blob/main/PropertyCalculator.ipynb).
It prints formulas you can paste into the editor, e.g. `conduction(A=48, L=0.08, k=1.75)`.
"""))

    cells.append(code(r'''
#@title 2a ▸ Start a network { display-mode: "form" }
#@markdown ⚠️ Re-running this cell **replaces** the current network (save yours first with step 5, to Google Drive or your computer).
start_with = "Trombe wall (spreadsheet example, IP)" #@param ["Trombe wall (spreadsheet example, IP)", "Simple room with slab (SI, uses formulas)", "Blank network (IP)", "Blank network (SI)", "Open a saved network from my Google Drive", "Upload a saved network (.json)"]
title = "" #@param {type:"string"}
#@markdown For *Google Drive*: the folder in your Drive (as in step 5) and the file name; leave the name blank to list the files there.
drive_folder = "ThermalMatrixModel" #@param {type:"string"}
drive_file = "my_network.json" #@param {type:"string"}

if start_with.startswith("Trombe"):
    net = example_trombe_wall()
elif start_with.startswith("Simple"):
    net = example_simple_room()
elif start_with.startswith("Blank"):
    net = Network("IP" if "(IP)" in start_with else "SI")
    net.add_boundary("Outdoor", "T_out", x=80, y=200)
elif "Google Drive" in start_with:
    folder = drive_path(drive_folder)
    saved = sorted(f for f in os.listdir(folder) if f.endswith(".json"))
    if not drive_file or drive_file not in saved:
        raise FileNotFoundError(f"'{drive_file}' is not in {folder}. Networks there: {', '.join(saved) or 'none yet'}")
    net = Network.from_json(os.path.join(folder, drive_file))
    print("opened", os.path.join(folder, drive_file))
    if net.saved_with and net.saved_with != MODEL_VERSION:
        print(f"note: this network was saved with model version {net.saved_with}; you are running {MODEL_VERSION}")
else:
    if not IN_COLAB:
        raise RuntimeError("Uploading needs Colab. Locally use: net = Network.from_json('file.json')")
    up = _colab_files.upload()
    net = Network.from_json(next(iter(up.values())).decode("utf-8"))
    if net.saved_with and net.saved_with != MODEL_VERSION:
        print(f"note: this network was saved with model version {net.saved_with}; you are running {MODEL_VERSION}")
if title:
    net.title = title
print(net)
'''))

    cells.append(code(r'''
#@title 2b ▸ Draw / edit the network { display-mode: "form" }
#@markdown **● Node** / **■ Boundary**: click the canvas to place one. **— Link**: click two nodes.
#@markdown **Select**: drag to move, click to edit values in the side panel. Click empty canvas for title & units.
#@markdown Every change is sent to Python automatically (watch for *✓ saved to Python*).
#@markdown Keys: V select · N node · B boundary · L link · Del delete.
editor_height = 560 #@param {type:"integer"}
edit_network(net, weather, height=editor_height)
'''))

    cells.append(code(r'''
#@title 2c ▸ Check the network: tables, diagram and the UA matrix { display-mode: "form" }
show_matrix = True #@param {type:"boolean"}
net.show(weather)
plt.show()
if show_matrix and not net.check(weather):
    display(Markdown(f"**UA matrix** (node–node conductances, conductances to each boundary, and C) "
                     f"— the spreadsheet's *Parameters* sheet, in {UNIT_LABELS[net.units]['UA']}:"))
    display(net.ua_table().round(2))
'''))

    cells.append(md(r"""
<details><summary><b>Prefer typing to drawing?</b> Build or change the network in code (click to expand)</summary>

```python
net = Network("SI", "My studio")
net.add_boundary("Outdoor", "T_out")
net.add_boundary("Ground", "T_ground")
net.add_node("Room Air", C="air_mass(120, multiplier=5)", T0=20, heat_set=20, Q=[200]*8 + [600]*10 + [300]*6)
net.add_node("Floor Surface", solar="aperture(8, 0.7, 0.8)*I_S + aperture(3, 0.7, 0.8)*I_W")
net.add_node("Slab", C="mass(40, 0.15, 'concrete')")
net.link("Outdoor", "Room Air", "window(8 + 3, 1.2) + ventilation(120, 0.4)")
net.link("Room Air", "Floor Surface", "convection(40, 3)")
net.link("Floor Surface", "Slab", "conduction(A=40, L=0.075, material='concrete')")
net.link("Slab", "Ground", "series(conduction(A=40, L=0.075, material='concrete'), conduction(A=40, L=0.1, material='XPS'))")
net.set("Room Air", cool_set=26)      # change a value later
```
Then run **2b** to see and adjust it in the editor.
</details>
"""))

    cells.append(code("# your own code here (optional)\n", form=False))

    cells.append(md(r"""
## 3 ▸ Simulate
Dates are `MM-DD` and inclusive. A period like `12-01` → `02-28` wraps over the new year.
**Warm-up days** are simulated before the start and thrown away, so the guessed initial
temperatures don't distort the results (heavy mass walls may need 2–3 weeks).
The spreadsheet uses no warm-up and one step per hour.
"""))

    cells.append(code(r'''
#@title 3 ▸ Run the simulation { display-mode: "form" }
start_date = "01-01" #@param {type:"string"}
end_date = "01-07" #@param {type:"string"}
warmup_days = 7 #@param {type:"integer"}
steps_per_hour = 1 #@param {type:"integer"}

results = simulate(net, weather, start=start_date, end=end_date,
                   warmup_days=warmup_days, substeps=steps_per_hour)
_t = results.T.index - pd.Timedelta(minutes=30)
print(f"simulated {len(results.T)} hours, {_t[0]:%b %d} – {_t[-1]:%b %d} (+ {warmup_days} warm-up days), model version {MODEL_VERSION}")
display(results.summary())
'''))

    cells.append(md("## 4 ▸ Results"))

    cells.append(code(r'''
#@title 4a ▸ Temperatures over time { display-mode: "form" }
#@markdown Comma-separated node names to plot (blank = all nodes).
nodes_to_plot = "" #@param {type:"string"}
sel = [s.strip() for s in nodes_to_plot.split(",") if s.strip()] or None
results.plot(sel)
plt.show()
results.plot_average_day(sel)
plt.show()
'''))

    cells.append(code(r'''
#@title 4b ▸ Heating / cooling loads and heat flows { display-mode: "form" }
lab = UNIT_LABELS[net.units]
if len(results.loads.columns):
    results.plot_loads(); plt.show()
else:
    print("No thermostats: set a heating or cooling set point on a node (e.g. Room Air) to get loads.")
fl = results.flows()
display(Markdown(f"**Heat flow along each link** (positive = from first to second node), {lab['Q']}"))
display(pd.DataFrame({"mean": fl.mean(), "max": fl.max(), "min": fl.min()}).round(1))
'''))

    cells.append(code(r'''
#@title 4c ▸ Hourly table { display-mode: "form" }
display(results.T.round(2))
'''))

    cells.append(md("## 5 ▸ Save your work"))

    cells.append(code(r'''
#@title 5 ▸ Save your network (.json), results (.csv) and diagram (.png) { display-mode: "form" }
#@markdown **Google Drive** keeps your work between sessions: files go to *My Drive ▸ drive_folder*
#@markdown (Colab asks once for permission). Open a saved network again in step 2a.
file_stem = "my_network" #@param {type:"string"}
save_to = "Google Drive" #@param ["Google Drive", "Download to this computer", "Both"]
drive_folder = "ThermalMatrixModel" #@param {type:"string"}

folder = drive_path(drive_folder) if save_to != "Download to this computer" else "."
base = os.path.join(folder, file_stem)
outs = [base + ".json", base + "_results.csv", base + "_network.png"]
net.to_json(outs[0])
results.to_csv(outs[1])
fig = plt.figure(figsize=(12, 6)); plot_network(net, ax=fig.gca()); fig.text(0.99, 0.01, f"model version {MODEL_VERSION}", ha="right", fontsize=7, color="gray"); fig.savefig(outs[2], dpi=150, bbox_inches="tight"); plt.close(fig)
print("saved:\n  " + "\n  ".join(outs))
if IN_COLAB and save_to != "Google Drive":
    for f in outs:
        _colab_files.download(f)
'''))

    cells.append(md(r"""
---
### Notes and limits of the method
* All heat transfer is **linear** (constant UA). Radiation and natural convection are
  linearised about a typical temperature, as in the 1979 paper. Things that switch — vent
  dampers closing at night, movable insulation, shading — are not modelled; a common
  work-around is to run day and night cases separately.
* Weather is hourly and held constant within each hour. `steps_per_hour > 1` refines the time
  step; it matters only for light nodes with small time constants.
* Thermostats are **ideal** (unlimited capacity, exact set point). The load reported is the
  heat that must be added (+) or removed (−) at that node.
* Solar uses an isotropic sky; it doesn't account for shading by overhangs or the glazing's
  angle-dependent transmittance — fold those into each window's aperture (τ or F).
* With the *Spreadsheet week* weather and the *Trombe wall* example (no warm-up, 1 step/hour),
  this notebook reproduces the spreadsheet's temperatures exactly.
"""))

    cells.append(md("### Version history\n" + "\n".join(entries)))

    if vent:
        cells = vent_variant(cells)

    nb = {"nbformat": 4, "nbformat_minor": 0,
          "metadata": {"colab": {"provenance": [], "toc_visible": True},
                       "kernelspec": {"name": "python3", "display_name": "Python 3"},
                       "language_info": {"name": "python"}},
          "cells": cells}
    for c in nb["cells"]:   # notebook format wants source as list of lines
        c["source"] = c["source"].splitlines(keepends=True)
    out = os.path.join(HERE, "ThermalMatrixModel_Vent.ipynb" if vent else "ThermalMatrixModel.ipynb")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
    write_changelog(entries)
    print("wrote", out, "- version", version)


VENT_MD = r"""
### Controlled ventilation (temperature-responsive windows)
Select a link in the editor and fill in **UA when open**. The link's ordinary UA becomes the
**closed** value (infiltration through the shut window), and it is drawn as a dashed blue
line labelled *closed ⇄ open*.

| field | meaning |
|---|---|
| **UA when open** | conductance with the windows open, e.g. `ventilation(54, 15)` (15 air changes per hour in a 54 m³ room) or `airflow(0.3)` |
| **sensor node** | the node whose temperature opens it (blank = the link's node end, usually *Room Air*) |
| **open when sensor is above** | indoor set point |
| **… and outdoor air is above** | outdoor set point (blank = no outdoor condition); keeps the windows shut in cold weather |
| **… and sensor is warmer than outdoor by** | indoor–outdoor difference (blank = no test); keeps the windows shut when it is hotter outside, so they don't let heat in. New controlled links start at 2 °C (3.6 °F) |
| **dead band** | the link closes again when any of the three tests falls this far below its set point (e.g. with a 2 °C difference and a 1 °C band, it closes once the room is less than 1 °C warmer than outside) |

The link **opens when all three tests pass**: the sensor is above its set point, the outdoor air
is above its set point, and the sensor is warmer than the outdoor air by the given difference.
Step 4d shows the heat each link removed (negative = it brought heat in).

Put the controlled ventilation on its **own link**, separate from the walls and glazing; giving
it its own outdoor boundary (as in the example) keeps the drawing clear.

**How it is solved, and why it stays stable.** Each step decides open/closed from the
temperatures at the end of the previous step, then uses the $M_1$ for that combination (with
*n* controlled links there are at most $2^n$ of them, each inverted once and kept). Switching
can't make the solution blow up: with any positive UA, $M_1$ stays diagonally dominant and
the backward difference is stable for any time step. What switching *can* do is **cycle**. Room air
is light and can swing several degrees in an hour, so with hourly steps a window near its set
point may flip open and shut every hour (a saw-tooth in the room temperature). Use
**steps_per_hour = 4** or more (the default here): the window can then switch every 15 minutes,
and each hour's result shows the fraction of the hour it was open, which behaves like a
window that opens part way. The dead band also reduces cycling.

In code:
```python
net.link("Outdoor (vent)", "Room Air", "ventilation(54, 0.5)", "operable windows",
         UA_open="ventilation(54, 15)", sensor="Room Air", open_above=24, outdoor_above=16, warmer_by=2, band=1)
```
"""

VENT_RESULTS = r'''
#@title 4d ▸ Controlled ventilation: when was it open? { display-mode: "form" }
vs = results.vent_summary()
if len(vs):
    display(Markdown("**Controlled links**: hours open, number of openings, and the heat each removed "
                     "from the room while open and while closed (negative = it brought heat in)."))
    display(vs)
    results.plot_vent(); plt.show()
else:
    print("No controlled links: select a link in 2b and give it a 'UA when open'.")
'''


def vent_variant(cells):
    """Turn the standard notebook's cells into the controlled-ventilation notebook."""
    def edit(start, *pairs):
        for c in cells:
            if c["source"].lstrip().startswith(start):
                for a, b in pairs:
                    assert a in c["source"], (start, a)
                    c["source"] = c["source"].replace(a, b)
                return
        raise AssertionError(f"cell {start!r} not found")

    def insert_after(start, cell):
        i = next(k for k, c in enumerate(cells) if c["source"].lstrip().startswith(start))
        cells.insert(i + 1, cell)

    edit("# 🌡️", ("*Runtime ▸ Run all* runs the Trombe-wall example end to end.",
                 "*Runtime ▸ Run all* runs the operable-window example for July end to end."))
    edit("#@title 2a", (
        'start_with = "Trombe wall (spreadsheet example, IP)" #@param ["Trombe wall',
        'start_with = "Room with operable windows (SI, controlled ventilation)" #@param '
        '["Room with operable windows (SI, controlled ventilation)", "Trombe wall'),
        ('if start_with.startswith("Trombe"):',
         'if start_with.startswith("Room with operable"):\n    net = example_vent_room()\n'
         'elif start_with.startswith("Trombe"):'))
    edit("#@title 2b", (
        "#@markdown Keys: V select · N node · B boundary · L link · Del delete.",
        "#@markdown Keys: V select · N node · B boundary · L link · Del delete.\n"
        "#@markdown **Controlled ventilation**: select a link and fill in *UA when open* (dashed blue = controlled)."),
        ("edit_network(net, weather, height=editor_height)",
         "edit_network(net, weather, height=editor_height, controls=True)"))
    insert_after("## 2 ▸", md(VENT_MD))
    edit("## 3 ▸", ("The spreadsheet uses no warm-up and one step per hour.",
                    "The spreadsheet uses no warm-up and one step per hour. With controlled ventilation use "
                    "4 or more steps per hour (see *Controlled ventilation* above). The default dates are July; "
                    "the *Spreadsheet week* weather only covers 01-01 to 01-07."))
    edit("#@title 3", ('start_date = "01-01"', 'start_date = "07-01"'),
         ('end_date = "01-07"', 'end_date = "07-31"'),
         ("steps_per_hour = 1 #@param", "steps_per_hour = 4 #@param"))
    insert_after("#@title 4b", code(VENT_RESULTS))
    edit("---\n### Notes", (
        "Things that switch — vent\n  dampers closing at night, movable insulation, shading — are not modelled; a common\n"
        "  work-around is to run day and night cases separately.",
        "Links can switch between two\n  values (controlled ventilation, above); other things that change — movable insulation,\n"
        "  shading — are not modelled."))
    return cells


CALC_INTRO = r"""
Material properties and layer values for the **Thermal Matrix Model**. Look up a material, then
get the conductance UA and capacitance C of a layer, split into segments for thick, massive
layers, with formulas you can paste straight into the network editor.

**Your own materials table.** The built-in table holds generic typical values. To use
**ASHRAE Handbook—Fundamentals, chapter 26, Table 1** (*Thermal Properties of Common Building and
Insulating Materials*) or manufacturers' data, put the rows you need in a CSV or Google Sheet
with these columns (SI, as in Table 1):

| column | meaning |
|---|---|
| `material` | your name for it, e.g. *hollow clay block 25 cm* |
| `density_kg_m3` | density ρ, kg/m³ |
| `k_W_mK` | conductivity k, W/(m·K) |
| `R_m2K_W` and `R_thickness_mm` | for products listed by **resistance at a given thickness** (boards, blocks), give R and that thickness instead of k; k is worked out as thickness / R |
| `cp_kJ_kgK` | specific heat, kJ/(kg·K) |
| `source` | where the numbers came from (e.g. *ASHRAE Fund. 2025 ch. 26 Table 1*) |

Step 1 can save a template with these columns to your Google Drive; open it with Google Sheets,
add rows, and load it back. Your rows are added to the generic ones; a row with the same name
replaces the generic value. ASHRAE's handbook is copyrighted, so keep your table in your own
Drive (or the course site) rather than posting it publicly. Penn users can read the handbook
online through the library.
"""

CALC_CELLS = [
    ("md", "## 1 ▸ Materials table"),
    ("code", r'''
#@title 1 ▸ Load the materials table { display-mode: "form" }
#@markdown *Google Drive*: the file `csv_name` in *My Drive ▸ drive_folder*. *Link*: a Google Sheet shared as "anyone with the link", or a CSV URL.
materials_from = "Generic values only" #@param ["Generic values only", "My CSV in Google Drive", "Upload a CSV", "Google Sheet or CSV link"]
drive_folder = "ThermalMatrixModel" #@param {type:"string"}
csv_name = "materials.csv" #@param {type:"string"}
link = "" #@param {type:"string"}
#@markdown Tick to save a template (generic rows + the right columns) as `csv_name` in your Drive folder.
save_template = False #@param {type:"boolean"}

if save_template:
    path = os.path.join(drive_path(drive_folder), csv_name)
    if os.path.exists(path):
        print(f"{path} already exists - not overwritten (pick another csv_name)")
    else:
        generic_table()[TEMPLATE_COLUMNS].to_csv(path, index=False)
        print("template saved to", path, "- open it with Google Sheets, add rows, then load it here")

user = None
if materials_from.startswith("My CSV"):
    user = read_table(os.path.join(drive_path(drive_folder), csv_name))
elif materials_from.startswith("Upload"):
    if not IN_COLAB:
        raise RuntimeError("Uploading needs Colab. Locally use: user = read_table('file.csv')")
    up = _colab_files.upload()
    user = read_table(next(iter(up.values())))
elif materials_from.startswith("Google Sheet"):
    user = read_table(link)
mats = Materials(user)
print(f"{len(mats.df)} materials" + (f" ({mats.n_user} rows read from your table, added to the generic values)" if mats.n_user else " (generic values)"))
'''),
    ("code", r'''
#@title 2 ▸ Look up materials { display-mode: "form" }
#@markdown Part of a name (blank = all), shown in SI or IP units.
search = "" #@param {type:"string"}
units = "SI" #@param ["SI", "IP"]
display(mats.show(search, units).round(4))
'''),
    ("md", r"""
## 2 ▸ One layer
Geometry in the network's units: m and m² for SI, ft and ft² for IP. **area** = face area the
heat passes through; **thickness** = layer thickness in the direction of heat flow.

**How many segments?** A massive layer is modelled as several capacitive nodes (like the
spreadsheet's *Wall 1–5*). A daily temperature swing dies away over the **penetration depth**
δ = √(α·P/π) (α = k/ρc, P = 1 day); segments no thicker than δ/2 follow it well, so that is the
suggestion when `n_segments` is 0. Light or insulating layers need only one.
"""),
    ("code", r'''
#@title 3 ▸ Layer values { display-mode: "form" }
material = "concrete" #@param {type:"string"}
units = "SI" #@param ["SI", "IP"]
area = 48 #@param {type:"number"}
thickness = 0.4 #@param {type:"number"}
n_segments = 0 #@param {type:"integer"}

table, formulas, info = layer(mats, material, area, thickness, units, n_segments)
u = UNITS[units]
print(f"{info['material']}: k = {info['k']:.4g} {u['k']}, ρ = {info['rho']:.4g} {u['rho']}, "
      f"c = {info['cp']:.4g} {u['cp']}" + ("  (k from R)" if info["k_from_R"] else "") +
      (f"  [{info['source']}]" if info["source"] else ""))
print(f"{info['n']} segment(s)" + ("" if info["n"] == info["suggested"] else f"  (suggested: {info['suggested']})"))
display(table.round(3))
display(Markdown("**Formulas to paste into the network editor** (" + units + " network):\n\n" +
                 "\n".join(f"- {k}: `{v}`" for k, v in formulas.items())))
'''),
    ("md", r"""
## 3 ▸ An assembly of layers
List the layers from **inside to outside**, one per line as `material, thickness`. The table
gives each layer's R and heat capacity, and the totals with surface films (defaults: inside
0.12, outside 0.03 m²·K/W; set to change them). Use it for the steady-state U-value of a wall
and to see where its mass sits.
"""),
    ("code", r'''
#@title 4 ▸ Assembly { display-mode: "form" }
units = "SI" #@param ["SI", "IP"]
area = 1 #@param {type:"number"}
R_inside = 0 #@param {type:"number"}
R_outside = 0 #@param {type:"number"}
#@markdown Edit the layers in the code (*Show code*, or double-click this cell):
layers = """
gypsum board, 0.0125
mineral wool, 0.15
brick, 0.10
"""
df, tot = assembly(mats, layers, area, units, R_inside or None, R_outside or None)
u = UNITS[units]
display(df.round(4))
print(f"R total (with surface films {tot['R_surfaces']:.3g}) = {tot['R_total']:.3g} {u['R']}   "
      f"U = {tot['U']:.3g} {u['U']}   UA = {tot['UA']:.4g} {u['UA']}   C = {tot['C']:.4g} {u['C']} for {area:g} {u['A']}")
'''),
]


def build_calculator(note=""):
    version = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    entries = changelog_entries(version, note, tag="[property calculator] ")
    src = open(os.path.join(HERE, "property_calculator.py"), encoding="utf-8").read()
    src, n = re.subn(r'^CALC_VERSION = .*$', f'CALC_VERSION = "{version}"', src, count=1, flags=re.M)
    assert n == 1
    cells = [md(f"# 🧱 Property Calculator\n**Version {version}**" + (f" - {note}" if note else "") + "\n" + CALC_INTRO)]
    cells.append(code(
        "#@title 0 ▸ Setup (run once) { display-mode: \"form\" }\n" + src
        + "\nfrom IPython.display import display, Markdown\n"
        + "try:\n    from google.colab import files as _colab_files\n    IN_COLAB = True\n"
        + "except ImportError:\n    IN_COLAB = False\n" + DRIVE_HELPER
        + "print(f'✓ property calculator loaded - version {CALC_VERSION}')\n"))
    for kind, text in CALC_CELLS:
        cells.append(md(text) if kind == "md" else code(text, form=not text.lstrip().startswith("#@title 4")))
    cells.append(md("### Version history\n" + "\n".join(e for e in entries if "[property calculator]" in e)))
    nb = {"nbformat": 4, "nbformat_minor": 0,
          "metadata": {"colab": {"provenance": []},
                       "kernelspec": {"name": "python3", "display_name": "Python 3"},
                       "language_info": {"name": "python"}},
          "cells": cells}
    for c in nb["cells"]:
        c["source"] = c["source"].splitlines(keepends=True)
    out = os.path.join(HERE, "PropertyCalculator.ipynb")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
    write_changelog(entries)
    print("wrote", out, "- version", version)


if __name__ == "__main__":
    args = sys.argv[1:]
    note = " ".join(a for a in args if not a.startswith("--")).strip()
    if "--calc" in args:
        build_calculator(note)
    else:
        build(note, vent="--vent" in args)
