# -*- coding: utf-8 -*-
"""
NACA 0012 airfoil mesh study: Gmsh boundary-layer mesh + SU2 RANS (Spalart-Allmaras)

What it does, in order
  1. installs Gmsh and downloads SU2 (Linux build) if they are missing
  2. estimates the first cell height from the target y+
  3. builds the mesh with Gmsh and checks it
  4. writes the SU2 settings file and runs SU2
  5. reports lift and drag, checks convergence, and adds the run to results_log.csv

How to use (Google Colab or any Linux machine)
  - Upload this file, then in a notebook cell run:      %run -i naca0012_study.py
  - Change settings without editing this file, in a cell of its own:
        OVERRIDES = {"LC_AIRFOIL": 0.0025, "LOW_MACH_CORR": "NO"}
        %run -i naca0012_study.py
    Names that are not in DEFAULTS below are rejected, so a typo cannot pass silently.
  - Build and check only the mesh (about a minute, no solver):
        OVERRIDES = {"MESH_ONLY": True}
  - The run's results are also left in the notebook variable LAST_RESULT.
"""

import csv
import datetime
import importlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import urllib.request
import zipfile

# =============================================================================
# 0. SETTINGS (change them with OVERRIDES, or edit the values here)
# =============================================================================
DEFAULTS = {
    "CASE": "naca0012_bl",        # name used for the mesh file

    # Flow conditions: NASA TMR NACA 0012 validation case
    "MACH": 0.15,
    "REYNOLDS": 6.0e6,            # based on chord = 1
    "TEMPERATURE": 300.0,         # kelvin
    "AOA": 0.0,                   # angle of attack in degrees (TMR also uses 10 and 15)

    # Mesh
    "YPLUS": 1.0,                 # target y+ of the first cell (flat-plate estimate)
    "GROWTH_RATIO": 1.2,          # growth between boundary-layer cells
    "BL_THICKNESS": 0.04,         # total thickness of the layered block, in chords
    "LC_AIRFOIL": 0.005,          # cell size along the airfoil surface, in chords
    "LC_FAR": 20.0,               # cell size at the outer circle
    "R_FAR": 100.0,               # outer circle radius in chords (NASA TMR grids reach ~500)
    "N_POINTS": 100,              # points defining each side of the airfoil outline
    "WAKE_LC": 0.005,             # cell size at the trailing edge (None switches this off)
    "WAKE_GROWTH": 0.1,           # that cell size grows by this much per unit distance

    # Solver
    "ITERATIONS": 3000,
    "CFL": 1.0,
    "MUSCL_TURB": "YES",                    # second-order turbulence advection
    "FLOW_LIMITER": "VENKATAKRISHNAN_WANG",
    "VENKAT_COEFF": 0.05,
    # Low-Mach correction of the Roe scheme. On NASA's 449x129 structured grid "NO" matches
    # the reference (+1.4%) and "YES" adds about 13% drag. On the coarser Gmsh meshes "YES"
    # is much closer to the reference; on the finest Gmsh mesh both agree within about 1.5%.
    "LOW_MACH_CORR": "YES",

    # NASA TMR reference values at alpha = 0 (CFL3D, Spalart-Allmaras). They are only used
    # when AOA is 0. Set both to None to switch the comparison off.
    "CL_REF": 0.0,
    "CD_REF": 0.00819,

    # Control
    "MESH_ONLY": False,           # True = build and check the mesh, then stop
    "GMSH_VERSION": "4.15.2",     # versions the published results were produced with
    "SU2_VERSION": "v8.5.0",      # use "latest" to take the newest SU2 release
}

YES_NO = ("YES", "NO")


def apply_overrides(defaults, overrides):
    """Return the settings with the overrides applied. Unknown names are an error."""
    unknown = sorted(k for k in overrides if k not in defaults)
    if unknown:
        raise ValueError("Unknown setting name(s): %s.\nValid names are: %s"
                         % (", ".join(unknown), ", ".join(sorted(defaults))))
    merged = dict(defaults)
    merged.update(overrides)
    return merged


def validate_settings(s):
    """Raise ValueError with a clear message if any setting is unusable."""
    problems = []
    for key in ("LOW_MACH_CORR", "MUSCL_TURB"):
        if s[key] not in YES_NO:
            problems.append("%s must be \"YES\" or \"NO\" (got %r)" % (key, s[key]))
    for key in ("MACH", "REYNOLDS", "TEMPERATURE", "YPLUS", "BL_THICKNESS", "LC_AIRFOIL",
                "LC_FAR", "R_FAR", "CFL", "VENKAT_COEFF", "WAKE_GROWTH", "ITERATIONS"):
        value = s[key]
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not value > 0:
            problems.append("%s must be a positive number (got %r)" % (key, value))
    if not isinstance(s["GROWTH_RATIO"], (int, float)) or not s["GROWTH_RATIO"] > 1:
        problems.append("GROWTH_RATIO must be a number above 1 (got %r)" % (s["GROWTH_RATIO"],))
    if not isinstance(s["N_POINTS"], int) or s["N_POINTS"] < 20:
        problems.append("N_POINTS must be a whole number of at least 20 (got %r)"
                        % (s["N_POINTS"],))
    if s["WAKE_LC"] is not None and not (isinstance(s["WAKE_LC"], (int, float))
                                         and s["WAKE_LC"] > 0):
        problems.append("WAKE_LC must be a positive number or None (got %r)" % (s["WAKE_LC"],))
    if isinstance(s["R_FAR"], (int, float)) and s["R_FAR"] <= 5:
        problems.append("R_FAR must be well above the airfoil size (at least 5 chords)")
    if not isinstance(s["CASE"], str) or not s["CASE"].strip():
        problems.append("CASE must be a non-empty text")
    if problems:
        raise ValueError("Invalid settings:\n  " + "\n  ".join(problems))


# Overrides come from the notebook. When this file is run with `%run -i` or exec(), the
# notebook's own variables are visible here; when it is imported (tests) there are none.
OVERRIDES_USED = dict(globals().get("OVERRIDES", {}) or {})
SETTINGS = apply_overrides(DEFAULTS, OVERRIDES_USED)
validate_settings(SETTINGS)


# =============================================================================
# 1. SMALL HELPERS (no Gmsh or SU2 needed, covered by tests/test_study.py)
# =============================================================================
def reference_values(s):
    """NASA reference CL and CD, or (None, None) when they do not apply."""
    if s["AOA"] != 0:
        return None, None
    return s["CL_REF"], s["CD_REF"]


def percent_difference(value, reference):
    return 100.0 * (value - reference) / reference


def first_cell_height(reynolds, yplus):
    """Flat-plate estimate of the first cell height (in chords) for a target y+."""
    cf = 0.026 * reynolds ** (-1.0 / 7.0)
    y1 = yplus * math.sqrt(2.0) / (reynolds * math.sqrt(cf))
    return cf, y1


def expected_layers(y1, thickness, ratio):
    """Number of boundary-layer cells needed to reach the block thickness."""
    return math.log(1.0 + thickness * (ratio - 1.0) / y1) / math.log(ratio)


def render_config(s, mesh_file):
    """Text of the SU2 settings file for the given settings."""
    mach, aoa, temp = s["MACH"], s["AOA"], s["TEMPERATURE"]
    reynolds = "%.1f" % s["REYNOLDS"]
    cfl, iters = s["CFL"], int(s["ITERATIONS"])
    low_mach, muscl_turb = s["LOW_MACH_CORR"], s["MUSCL_TURB"]
    limiter, venkat = s["FLOW_LIMITER"], s["VENKAT_COEFF"]
    return f"""SOLVER= RANS
KIND_TURB_MODEL= SA
MATH_PROBLEM= DIRECT
MACH_NUMBER= {mach}
AOA= {aoa}
FREESTREAM_TEMPERATURE= {temp}
REYNOLDS_NUMBER= {reynolds}
REYNOLDS_LENGTH= 1.0
REF_ORIGIN_MOMENT_X= 0.25
REF_ORIGIN_MOMENT_Y= 0.0
REF_ORIGIN_MOMENT_Z= 0.0
REF_LENGTH= 1.0
REF_AREA= 1.0
REF_DIMENSIONALIZATION= FREESTREAM_PRESS_EQ_ONE
FLUID_MODEL= STANDARD_AIR
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716E-5
MU_T_REF= 273.15
SUTHERLAND_CONSTANT= 110.4
MARKER_HEATFLUX= ( airfoil, 0.0 )
MARKER_FAR= ( farfield )
MARKER_PLOTTING= ( airfoil )
MARKER_MONITORING= ( airfoil )
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CFL_NUMBER= {cfl}
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.5, 1.1, 0.5, 20.0 )
ITER= {iters}
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-6
LINEAR_SOLVER_ITER= 10
CONV_NUM_METHOD_FLOW= ROE
LOW_MACH_CORR= {low_mach}
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= {limiter}
VENKAT_LIMITER_COEFF= {venkat}
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= {muscl_turb}
SLOPE_LIMITER_TURB= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_TURB= EULER_IMPLICIT
CONV_RESIDUAL_MINVAL= -10
MESH_FILENAME= {mesh_file}
MESH_FORMAT= SU2
SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, LIFT, DRAG)
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF)
VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)
OUTPUT_FILES= (RESTART, PARAVIEW, SURFACE_CSV)
"""


def read_mesh_summary(path):
    """Read node/element counts and boundary-marker sizes from an SU2 mesh file."""
    info = {"npoin": None, "nelem": None, "markers": {}}
    tag = None
    with open(path) as f:
        for line in f:
            if line.startswith("NPOIN="):
                info["npoin"] = int(line.split("=")[1].split()[0])
            elif line.startswith("NELEM="):
                info["nelem"] = int(line.split("=")[1].split()[0])
            elif line.startswith("MARKER_TAG="):
                tag = line.split("=")[1].strip()
            elif line.startswith("MARKER_ELEMS=") and tag is not None:
                info["markers"][tag] = int(line.split("=")[1].split()[0])
                tag = None
    return info


def read_history(path):
    """Read SU2's history.csv into {column name: list of numbers}."""
    with open(path, newline="") as f:
        reader = csv.reader(f, skipinitialspace=True)
        header = [name.strip().strip('"').strip() for name in next(reader)]
        columns = {name: [] for name in header}
        for row in reader:
            if len(row) < len(header):
                continue
            for name, value in zip(header, row):
                try:
                    columns[name].append(float(value))
                except ValueError:
                    columns[name].append(float("nan"))
    return columns


def convergence_report(cd_series, final_residual, min_residual=-6.0,
                       max_rel_cd_change=0.005, window=100):
    """Judge whether a run settled. Returns {"converged", "problems", "cd_rel_change"}."""
    tail = list(cd_series)[-window:]
    last = tail[-1]
    rel_change = (max(tail) - min(tail)) / abs(last) if last != 0 else float("inf")
    problems = []
    if final_residual is not None and final_residual > min_residual:
        problems.append("density residual is only %.2f (want below %.1f)"
                        % (final_residual, min_residual))
    if rel_change > max_rel_cd_change:
        problems.append("CD still varied by %.2f%% over the last %d iterations (want below %.1f%%)"
                        % (100 * rel_change, len(tail), 100 * max_rel_cd_change))
    return {"converged": not problems, "problems": problems, "cd_rel_change": rel_change}


def merge_into_log(path, row):
    """Add one run to the results table, keeping old rows even if columns changed."""
    old_rows = []
    if os.path.exists(path):
        with open(path, newline="") as f:
            old_rows = [{k: v for k, v in r.items() if k is not None}
                        for r in csv.DictReader(f)]
    all_rows = old_rows + [row]
    fields = []
    for r in all_rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(all_rows)
    return all_rows


# =============================================================================
# 2. TOOLS: Gmsh and SU2
# =============================================================================
def ensure_gmsh(version):
    """Import gmsh, installing it (and its one system library) if needed."""
    try:
        import gmsh
    except (ImportError, OSError):
        print("installing gmsh %s ..." % version)
        done = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "gmsh==" + version])
        if done.returncode != 0:
            print("could not install gmsh==%s; trying the newest version instead" % version)
            subprocess.run([sys.executable, "-m", "pip", "install", "-q", "gmsh"], check=True)
        subprocess.run("apt-get install -y -q libglu1-mesa", shell=True, check=False)
        importlib.invalidate_caches()
        import gmsh
    if gmsh.__version__ != version:
        print("note: Gmsh %s is installed; the published results used %s"
              % (gmsh.__version__, version))
    return gmsh


def ensure_su2(version, su2_dir="su2"):
    """Return (path to SU2_CFD, version tag), downloading the Linux release if needed."""
    su2_dir = os.path.abspath(su2_dir)
    su2_bin = os.path.join(su2_dir, "bin", "SU2_CFD")
    version_file = os.path.join(su2_dir, "VERSION.txt")
    if os.path.exists(su2_bin):
        tag = "unknown (installed by an earlier run)"
        if os.path.exists(version_file):
            tag = open(version_file).read().strip()
        return su2_bin, tag
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "AMD64"):
        raise RuntimeError(
            "This script downloads SU2's Linux build, and this machine is %s/%s. Run it in "
            "Google Colab or on Linux, or install SU2 yourself and place SU2_CFD at %s"
            % (platform.system(), platform.machine(), su2_bin))
    api = "https://api.github.com/repos/su2code/SU2/releases/"
    api += "latest" if version == "latest" else "tags/" + version
    data = json.load(urllib.request.urlopen(api))
    urls = [a["browser_download_url"] for a in data["assets"] if "linux64-omp" in a["name"]]
    if not urls:
        raise RuntimeError("no linux64-omp download found in SU2 release %s" % data["tag_name"])
    print("downloading SU2", data["tag_name"], "...")
    urllib.request.urlretrieve(urls[0], "su2.zip")
    with zipfile.ZipFile("su2.zip") as z:
        z.extractall(su2_dir)
    inner = os.path.join(su2_dir, "linux64-omp.zip")      # the release holds a second zip
    if os.path.exists(inner) and not os.path.exists(su2_bin):
        with zipfile.ZipFile(inner) as z:
            z.extractall(su2_dir)
    if not os.path.exists(su2_bin):
        raise RuntimeError("SU2_CFD not found after unpacking; look inside the %s folder" % su2_dir)
    bin_dir = os.path.join(su2_dir, "bin")
    for name in os.listdir(bin_dir):
        os.chmod(os.path.join(bin_dir, name), 0o755)
    with open(version_file, "w") as f:
        f.write(data["tag_name"])
    return su2_bin, data["tag_name"]


# =============================================================================
# 3. MESH with Gmsh
# =============================================================================
def build_mesh(gmsh, s, y1, mesh_file):
    """Build the airfoil mesh, write it to mesh_file, and return counts and plot data."""
    import numpy as np

    n_points = int(s["N_POINTS"])
    lc_airfoil, lc_far, r_far = s["LC_AIRFOIL"], s["LC_FAR"], s["R_FAR"]

    if gmsh.isInitialized():
        gmsh.finalize()
    gmsh.initialize()
    try:
        gmsh.model.add(s["CASE"])
        geo = gmsh.model.geo
        fld = gmsh.model.mesh.field

        # airfoil outline (sharp trailing edge version of the NACA 0012 formula)
        t = 0.12
        x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n_points)))
        yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                      + 0.2843 * x**3 - 0.1036 * x**4)
        yt[0] = 0.0
        yt[-1] = 0.0

        p = [geo.addPoint(x[i], yt[i], 0, lc_airfoil) for i in range(n_points)]
        q = [geo.addPoint(x[i], -yt[i], 0, lc_airfoil) for i in range(1, n_points - 1)]
        upper = geo.addSpline(p)
        lower = geo.addSpline([p[-1]] + q[::-1] + [p[0]])

        # outer circle
        c = geo.addPoint(0.5, 0, 0, lc_far)
        f1 = geo.addPoint(0.5 + r_far, 0, 0, lc_far)
        f2 = geo.addPoint(0.5, r_far, 0, lc_far)
        f3 = geo.addPoint(0.5 - r_far, 0, 0, lc_far)
        f4 = geo.addPoint(0.5, -r_far, 0, lc_far)
        a1 = geo.addCircleArc(f1, c, f2)
        a2 = geo.addCircleArc(f2, c, f3)
        a3 = geo.addCircleArc(f3, c, f4)
        a4 = geo.addCircleArc(f4, c, f1)

        outer = geo.addCurveLoop([a1, a2, a3, a4])
        inner = geo.addCurveLoop([upper, lower])
        surf = geo.addPlaneSurface([outer, inner])
        geo.synchronize()

        # thin stretched layers next to the airfoil
        bl = fld.add("BoundaryLayer")
        fld.setNumbers(bl, "CurvesList", [upper, lower])
        fld.setNumbers(bl, "PointsList", [p[-1]])        # sharp trailing edge
        fld.setNumber(bl, "Size", y1)
        fld.setNumber(bl, "Ratio", s["GROWTH_RATIO"])
        fld.setNumber(bl, "Thickness", s["BL_THICKNESS"])
        fld.setNumber(bl, "Quads", 1)
        fld.setAsBoundaryLayer(bl)

        # optional: small cells around the trailing edge (x = 1, y = 0), growing with distance
        if s["WAKE_LC"] is not None:
            expr = "%g + %g*Sqrt((x-1)*(x-1) + y*y)" % (s["WAKE_LC"], s["WAKE_GROWTH"])
            wk = fld.add("MathEval")
            fld.setString(wk, "F", expr)
            fld.setAsBackgroundMesh(wk)

        # labels the solver reads
        gmsh.model.addPhysicalGroup(1, [upper, lower], name="airfoil")
        gmsh.model.addPhysicalGroup(1, [a1, a2, a3, a4], name="farfield")
        gmsh.model.addPhysicalGroup(2, [surf], name="fluid")

        gmsh.model.mesh.generate(2)

        # collect the mesh for counting and plotting
        tags, coords, _ = gmsh.model.mesh.getNodes()
        pts = coords.reshape(-1, 3)
        index = {int(tg): i for i, tg in enumerate(tags)}
        etypes, _, enodes = gmsh.model.mesh.getElements(2)
        pieces = []
        n_tri = 0
        n_quad = 0
        for et, en in zip(etypes, enodes):
            if et == 2:                                   # triangles
                a = np.array([index[int(v)] for v in en]).reshape(-1, 3)
                pieces.append(a)
                n_tri = len(a)
            elif et == 3:                                 # quads: split into 2 triangles
                a = np.array([index[int(v)] for v in en]).reshape(-1, 4)
                pieces.append(a[:, [0, 1, 2]])
                pieces.append(a[:, [0, 2, 3]])
                n_quad = len(a)
        tri = np.vstack(pieces)

        gmsh.write(mesh_file)
    finally:
        gmsh.finalize()
    return {"nodes": len(tags), "triangles": n_tri, "quads": n_quad, "pts": pts, "tri": tri}


def plot_mesh(pts, tri, path="mesh.png"):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    views = [("Whole airfoil", (-0.1, 1.1), (-0.15, 0.15)),
             ("Nose close-up", (-0.03, 0.08), (-0.05, 0.05)),
             ("Tail close-up", (0.9, 1.05), (-0.04, 0.04))]
    for ax, (title, xl, yl) in zip(axs, views):
        ax.triplot(pts[:, 0], pts[:, 1], tri, linewidth=0.3)
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.set_title(title)
    plt.savefig(path, dpi=120)
    plt.show()


def plot_history(columns, path="convergence.png"):
    import matplotlib.pyplot as plt
    rho_name = next((k for k in columns if k.startswith("rms[Rho]")), None)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4))
    if rho_name:
        axs[0].plot(columns[rho_name])
    axs[0].set_title("Residual rms[Rho]")
    axs[1].plot(columns["CL"])
    axs[1].set_title("CL (lift)")
    axs[2].plot(columns["CD"])
    axs[2].set_title("CD (drag)")
    for ax in axs:
        ax.set_xlabel("iteration")
    plt.savefig(path, dpi=120)
    plt.show()


# =============================================================================
# 4. SU2
# =============================================================================
STALE_OUTPUTS = ("history.csv", "surface.csv", "vol_solution.vtu", "restart.dat", "su2_run.log")


def run_su2(su2_bin, cfg="rans.cfg"):
    """Run SU2 and return the run time in minutes. Raises if the run did not finish cleanly."""
    for name in STALE_OUTPUTS:                      # never read results of an earlier run
        if os.path.exists(name):
            os.remove(name)
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(su2_bin) + ":" + env.get("PATH", "")
    start = time.time()
    proc = subprocess.run([su2_bin, cfg], env=env, capture_output=True, text=True)
    minutes = (time.time() - start) / 60.0
    with open("su2_run.log", "w") as f:
        f.write(proc.stdout + "\n" + proc.stderr)
    print("\n".join(proc.stdout.splitlines()[-25:]))
    if "diverged" in proc.stdout:
        raise RuntimeError("SU2 diverged: see su2_run.log (try a smaller CFL)")
    if proc.returncode != 0 or "Exit Success" not in proc.stdout:
        raise RuntimeError("SU2 did not finish cleanly: see su2_run.log")
    return minutes


# =============================================================================
# 5. THE WHOLE RUN
# =============================================================================
def main():
    s = SETTINGS
    cl_ref, cd_ref = reference_values(s)
    if OVERRIDES_USED:
        print("overrides applied:", OVERRIDES_USED)
    if s["AOA"] != 0:
        print("note: the NASA reference values are for AOA = 0, so no comparison is made")

    gmsh = ensure_gmsh(s["GMSH_VERSION"])
    su2_bin, su2_tag = None, "not needed (MESH_ONLY)"
    if not s["MESH_ONLY"]:
        su2_bin, su2_tag = ensure_su2(s["SU2_VERSION"])
    print("tools ready (Gmsh %s, SU2 %s)" % (gmsh.__version__, su2_tag))

    cf, y1 = first_cell_height(s["REYNOLDS"], s["YPLUS"])
    layers = expected_layers(y1, s["BL_THICKNESS"], s["GROWTH_RATIO"])
    print("skin friction estimate Cf = %.5f" % cf)
    print("first cell height (chords) = %.3e" % y1)
    print("boundary-layer cells expected: %d" % round(layers))

    mesh_file = s["CASE"] + ".su2"
    mesh = build_mesh(gmsh, s, y1, mesh_file)
    print("mesh: nodes = %d, triangles = %d, quads = %d"
          % (mesh["nodes"], mesh["triangles"], mesh["quads"]))
    print("mesh file saved:", mesh_file)

    summary = read_mesh_summary(mesh_file)
    for name in ("airfoil", "farfield"):
        if summary["markers"].get(name, 0) <= 0:
            raise RuntimeError("the mesh file has no boundary elements named %r; found %s"
                               % (name, summary["markers"]))
    print("boundary elements: airfoil = %d, farfield = %d"
          % (summary["markers"]["airfoil"], summary["markers"]["farfield"]))
    if summary["npoin"] is not None and summary["npoin"] != mesh["nodes"]:
        print("note: the mesh file lists %d points, Gmsh counted %d nodes"
              % (summary["npoin"], mesh["nodes"]))
    plot_mesh(mesh["pts"], mesh["tri"])

    result = {"nodes": mesh["nodes"], "triangles": mesh["triangles"], "quads": mesh["quads"],
              "y1": y1, "settings": dict(s)}
    if s["MESH_ONLY"]:
        print("MESH_ONLY is on: stopping before the solver")
        return result

    with open("rans.cfg", "w") as f:
        f.write(render_config(s, mesh_file))
    print("settings file written: rans.cfg")
    print("running SU2 (this can take several minutes) ...")
    minutes = run_su2(su2_bin)

    columns = read_history("history.csv")
    for name in ("CL", "CD"):
        if name not in columns:
            raise RuntimeError("history.csv has no %s column; columns found: %s"
                               % (name, list(columns)))
    rho_name = next((k for k in columns if k.startswith("rms[Rho]")), None)
    cl, cd = columns["CL"][-1], columns["CD"][-1]
    residual = columns[rho_name][-1] if rho_name else None
    last = slice(-min(100, len(columns["CD"])), None)
    cl_var = max(columns["CL"][last]) - min(columns["CL"][last])
    cd_var = max(columns["CD"][last]) - min(columns["CD"][last])
    report = convergence_report(columns["CD"], residual)

    print("\niterations run:", len(columns["CD"]))
    print("final CL = %.6f   (variation over last 100 iterations: %.2e)" % (cl, cl_var))
    print("final CD = %.6f   (variation over last 100 iterations: %.2e)" % (cd, cd_var))
    if residual is not None:
        print("final density residual (log10) = %.2f" % residual)
    if cd_ref is not None:
        print("CD vs reference: %.6f vs %.6f  (difference %.2f%%)"
              % (cd, cd_ref, percent_difference(cd, cd_ref)))
    if cl_ref is not None:
        print("CL vs reference: %.6f vs %.6f  (difference %.6f)" % (cl, cl_ref, cl - cl_ref))
    if report["converged"]:
        print("convergence check: OK")
    else:
        print("convergence check: WARNING - " + "; ".join(report["problems"]))
        print("  treat this result with caution")

    row = {
        "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "AOA": s["AOA"], "y1": "%.3e" % y1, "growth": s["GROWTH_RATIO"],
        "bl_thickness": s["BL_THICKNESS"], "lc_airfoil": s["LC_AIRFOIL"],
        "lc_far": s["LC_FAR"], "R_far": s["R_FAR"], "nodes": mesh["nodes"],
        "muscl_turb": s["MUSCL_TURB"], "limiter": s["FLOW_LIMITER"],
        "low_mach": s["LOW_MACH_CORR"], "cfl": s["CFL"],
        "wake_lc": "" if s["WAKE_LC"] is None else s["WAKE_LC"],
        "yplus_target": s["YPLUS"], "iterations": len(columns["CD"]),
        "CL": "%.6f" % cl, "CD": "%.6f" % cd,
        "log10_rms_rho": "" if residual is None else "%.2f" % residual,
        "converged": "YES" if report["converged"] else "NO",
        "su2_minutes": "%.1f" % minutes,
        "gmsh_version": gmsh.__version__, "su2_version": su2_tag,
    }
    merge_into_log("results_log.csv", row)
    print("run added to results_log.csv")
    plot_history(columns)

    result.update({"CL": cl, "CD": cd, "residual": residual, "converged": report["converged"],
                   "iterations": len(columns["CD"]), "su2_minutes": minutes})
    return result


if __name__ == "__main__":
    LAST_RESULT = main()
