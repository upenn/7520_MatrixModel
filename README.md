# Thermal Matrix Model (Python / Colab)

Python version of the spreadsheet `MatrixModel_TrombeWall_25.xlsx` (method from
`TrombeWall_1979.pdf`), with a Colab notebook for students.

## Open in Colab
| notebook | |
|---|---|
| **Thermal Matrix Model** (standard) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/upenn/7520_MatrixModel/blob/main/ThermalMatrixModel.ipynb) |
| **Thermal Matrix Model + controlled ventilation** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/upenn/7520_MatrixModel/blob/main/ThermalMatrixModel_Vent.ipynb) |
| **Property Calculator** (materials, layer values, assemblies) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/upenn/7520_MatrixModel/blob/main/PropertyCalculator.ipynb) |

The links always open the latest version on GitHub. Step 5 saves your network (`.json`), results
(`.csv`) and diagram (`.png`) to *My Drive ▸ ThermalMatrixModel* in your Google Drive (or downloads
them); step 2a opens a saved network again. To keep edits to the notebook itself, use
*File ▸ Save a copy in Drive*.

## Files

| file | what it is |
|---|---|
| `ThermalMatrixModel.ipynb` | **The notebook to give students.** Self-contained: upload it to Colab (File ▸ Upload notebook) or put it on Drive / GitHub. |
| `ThermalMatrixModel_Vent.ipynb` | **Variant with controlled ventilation**: a link can switch between a closed and an open UA (temperature-responsive windows). Same engine; opens on the operable-window example in July. |
| `PropertyCalculator.ipynb` | **Stand-alone property calculator**: generic materials plus the student's own table (CSV or Google Sheet laid out like ASHRAE *Fundamentals* ch. 26 Table 1); layer UA and C split into segments, paste-ready formulas, assembly U-values. Built from `property_calculator.py` with `python build_notebook.py --calc "note"`. |
| `thermal_matrix.py` | The engine (source of truth). Network, EPW reader, solar on surfaces, property helpers, implicit solver, thermostats, plots, drawing editor. |
| `build_notebook.py` | Rebuilds the notebook from `thermal_matrix.py` (run after editing the engine): `python build_notebook.py "what changed"` (add `--vent` for the ventilation notebook). Stamps the build date/time as the model version. |
| `CHANGELOG.md` | Version history, one line per build (also shown at the end of the notebook). |
| `tests/test_against_excel.py` | Checks the Python model reproduces the spreadsheet's hourly temperatures (max difference ~1e-13 °F). |

## Student workflow in the notebook
1. **Setup** cell (engine is embedded, collapsed).
2. **Weather**: Philadelphia TMYx is downloaded by default; or upload an `.epw`/`.zip`, paste a URL
   (climate.onebuilding.org), or use the spreadsheet's week. Adds solar on S/E/W/N vertical and
   horizontal surfaces; any other tilt/azimuth can be added.
3. **Network**: start from the Trombe-wall example, a simple SI room, a blank page or a saved `.json`,
   then **draw it** in the editor (nodes, boundaries, links; values can be formulas like
   `conduction(A=48, L=0.12, material="concrete")`). Edits sync to Python live (Colab only).
   A check cell shows tables, a diagram and the UA matrix (like the *Parameters* sheet).
4. **Simulate** any date range (wraps over New Year) with warm-up days and optional sub-hourly steps.
5. **Results**: temperatures, average day, heating/cooling loads for nodes with set points,
   link heat flows, CSV + JSON + PNG download.

## Beyond the spreadsheet
Several boundary temperatures (outdoor, ground, sky, constant), internal gains with 24-hour profiles,
ideal heating/cooling thermostats (exact, by superposition on M1⁻¹), a materials library, and
SI or IP units per network.

## Controlled ventilation (ThermalMatrixModel_Vent.ipynb)
A link with a **UA when open** switches between its closed UA and the open one. It opens when a
sensor node (usually Room Air) is above a set point **and** the outdoor air is above an outdoor
set point **and** the sensor is warmer than outdoors by `warmer_by` (default 2 °C), and closes
when any of these falls a dead band below. Each step decides from the previous
step's temperatures; each open/closed combination has its own M1, inverted once and cached.
Positive UAs keep M1 diagonally dominant, so switching cannot destabilise the implicit scheme,
but light nodes can cycle open/shut every hour: use 4+ steps per hour (the notebook's default),
which turns cycling into a fraction of each hour open. Results add `vent`, `vent_summary()` and
`plot_vent()`; `flows()` uses the open or closed UA hour by hour. Networks without controlled
links behave exactly as before (the spreadsheet test still matches).

## Versions
Every build is stamped with its date and time (e.g. `2026-11-19 19:31`). The stamp appears in the
notebook title, the setup cell's output, the simulation summary, saved network `.json` files,
the header of results `.csv` files and the corner of the network `.png`. When a student uploads
a network saved with a different version, the notebook says so. Builds and their notes are listed
in `CHANGELOG.md`. If the spreadsheet is open in Excel, point `MATRIX_XLSX` at a copy of it.
