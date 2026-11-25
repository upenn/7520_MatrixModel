# Thermal Matrix Model: version history

Newest first. Add a note with `python build_notebook.py "what changed"`.
Entries marked [ventilation notebook] are builds of ThermalMatrixModel_Vent.ipynb (`python build_notebook.py --vent "what changed"`); both notebooks share the engine.
Entries marked [property calculator] are builds of PropertyCalculator.ipynb (`python build_notebook.py --calc "what changed"`).

- **2026-11-24 19:15** - [property calculator] Read the materials table straight from a Google Sheet in your Drive by name (no need to save as CSV); template instructions updated
- **2026-11-24 18:53** - [property calculator] First version: materials table (generic values plus your own CSV or Google Sheet in the layout of ASHRAE Fundamentals ch. 26 Table 1), layer values with suggested segments and paste-ready formulas, and assembly U-values
- **2026-11-24 18:52** - [ventilation notebook] Property calculator moved to its own notebook (PropertyCalculator.ipynb); link to it in step 2
- **2026-11-24 18:52** - Property calculator moved to its own notebook (PropertyCalculator.ipynb); link to it in step 2
- **2026-11-24 18:11** - [ventilation notebook] Save networks, results and diagrams to Google Drive (default) or download them; open a saved network from Google Drive in step 2a
- **2026-11-24 18:11** - Save networks, results and diagrams to Google Drive (default) or download them; open a saved network from Google Drive in step 2a
- **2026-11-23 16:47** - [ventilation notebook] Indoor-outdoor difference test: a controlled link opens only when the sensor is warmer than outdoors by warmer_by (default 2 C / 3.6 F), so windows stay shut when it is hotter outside
- **2026-11-23 16:43** - [ventilation notebook] First version: links can switch between a closed and an open UA (temperature-responsive windows), opening when a sensor node and the outdoor air are above set points, with a dead band; defaults to 4 steps per hour
- **2026-11-21 15:31** - Solar gain per node is now a formula with one term per window, e.g. aperture(6, 0.75, 0.9)*I_S + aperture(3, 0.75, 0.9)*I_E; older networks convert automatically
- **2026-11-19 19:31** - Version stamp (date/time) on the notebook, saved networks, results and diagrams; version history
- **2026-11-19** - Clearer formula help in the editor and notebook (L = thickness, e1/e2 = infrared emissivities, solar absorptance goes in the solar aperture); removed the "insert layered mass wall" cell.
- **2026-11-19** - First version: Python engine, EPW weather, drawing editor, Trombe wall and simple room examples; matches the spreadsheet exactly.
