"""Check thermal_matrix.py reproduces MatrixModel_TrombeWall_25.xlsx exactly."""
import os, sys, warnings
import numpy as np, openpyxl
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from thermal_matrix import *

warnings.simplefilter("ignore")
# MATRIX_XLSX can point at a copy, e.g. while the workbook is open in Excel
_path = os.environ.get("MATRIX_XLSX") or os.path.join(HERE, "..", "MatrixModel_TrombeWall_25.xlsx")
wb = openpyxl.load_workbook(_path, data_only=True)
w = wb["WeatherData"]
T = [w.cell(r, 6).value for r in range(3, 8763)]
S = [w.cell(r, 7).value for r in range(3, 8763)]
calc = wb["Calculations"]
excel = np.array([[calc.cell(r, c).value for c in range(3, 3 + 168)] for r in range(42, 53)]).T

wx = example_excel_weather(T, S)
net = example_trombe_wall()
res = simulate(net, wx, start="01-01", end="01-07", warmup_days=0)
names = [n["name"] for n in net.nodes]
py = res.T[names].values
err = np.abs(py - excel).max()
print("hours compared:", len(excel), " max |difference| (F):", err)
assert err < 1e-6, err
print("OK: Python matches the spreadsheet")
