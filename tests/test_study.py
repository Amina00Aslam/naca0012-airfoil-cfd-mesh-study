"""
Tests for naca0012_study.py that do not need Gmsh or SU2.

Run from the repository folder:   python tests/test_study.py
(They also work with pytest if you have it installed.)

The Gmsh parts are checked with a stand-in object that records the calls, so these tests
catch typos and wrong variable names, but only a real run proves the mesh itself.
"""
import csv
import importlib.util
import math
import os
import sys
import tempfile
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "naca0012_study.py")
spec = importlib.util.spec_from_file_location("study", SCRIPT)
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)          # imports only; main() is not run


def expect_error(exc_type, func, *args, text=None):
    try:
        func(*args)
    except exc_type as e:
        if text is not None:
            assert text in str(e), "message %r does not contain %r" % (str(e), text)
        return
    raise AssertionError("expected %s" % exc_type.__name__)


# ---------------------------------------------------------------- settings
def test_overrides_reject_typos():
    expect_error(ValueError, study.apply_overrides, study.DEFAULTS, {"LC_AIRFOL": 0.01},
                 text="LC_AIRFOL")
    merged = study.apply_overrides(study.DEFAULTS, {"LC_AIRFOIL": 0.01})
    assert merged["LC_AIRFOIL"] == 0.01 and study.DEFAULTS["LC_AIRFOIL"] == 0.005


def test_validate_settings():
    study.validate_settings(dict(study.DEFAULTS))
    for bad in ({"LOW_MACH_CORR": "yes"}, {"LC_AIRFOIL": -1}, {"GROWTH_RATIO": 1.0},
                {"N_POINTS": 5}, {"WAKE_LC": 0}, {"R_FAR": 2}, {"ITERATIONS": "many"},
                {"CASE": ""}):
        s = dict(study.DEFAULTS)
        s.update(bad)
        expect_error(ValueError, study.validate_settings, s)
    s = dict(study.DEFAULTS)
    s["WAKE_LC"] = None                      # switching the tail refinement off is allowed
    study.validate_settings(s)


def test_reference_values_only_at_zero_aoa():
    assert study.reference_values(dict(study.DEFAULTS)) == (0.0, 0.00819)
    s = dict(study.DEFAULTS)
    s["AOA"] = 10.0
    assert study.reference_values(s) == (None, None)


# ---------------------------------------------------------------- numbers
def test_first_cell_height_matches_the_runs():
    cf, y1 = study.first_cell_height(6.0e6, 1.0)
    assert abs(cf - 0.0027968297) < 1e-9
    assert abs(y1 - 4.456878e-6) < 1e-11
    assert round(study.expected_layers(y1, 0.04, 1.2)) == 41
    assert round(study.expected_layers(y1, 0.04, 1.1)) == 71        # as seen in earlier runs
    assert abs(study.percent_difference(0.010521, 0.00819) - 28.5) < 0.1


# ---------------------------------------------------------------- settings file
KNOWN_GOOD_CONFIG = """SOLVER= RANS
KIND_TURB_MODEL= SA
MATH_PROBLEM= DIRECT
MACH_NUMBER= 0.15
AOA= 0.0
FREESTREAM_TEMPERATURE= 300.0
REYNOLDS_NUMBER= 6000000.0
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
CFL_NUMBER= 1.0
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.5, 1.1, 0.5, 20.0 )
ITER= 3000
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-6
LINEAR_SOLVER_ITER= 10
CONV_NUM_METHOD_FLOW= ROE
LOW_MACH_CORR= YES
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG
VENKAT_LIMITER_COEFF= 0.05
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= YES
SLOPE_LIMITER_TURB= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_TURB= EULER_IMPLICIT
CONV_RESIDUAL_MINVAL= -10
MESH_FILENAME= naca0012_bl.su2
MESH_FORMAT= SU2
SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, LIFT, DRAG)
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF)
VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)
OUTPUT_FILES= (RESTART, PARAVIEW, SURFACE_CSV)
"""


def test_config_equals_the_file_used_for_the_published_runs():
    # this text is the settings file the final runs were made with
    assert study.render_config(dict(study.DEFAULTS), "naca0012_bl.su2") == KNOWN_GOOD_CONFIG


def test_config_changes_follow_the_settings():
    s = study.apply_overrides(study.DEFAULTS, {"LOW_MACH_CORR": "NO", "AOA": 2.0,
                                               "ITERATIONS": 500, "CFL": 2.5})
    text = study.render_config(s, "other.su2")
    for line in ("LOW_MACH_CORR= NO", "AOA= 2.0", "ITER= 500", "CFL_NUMBER= 2.5",
                 "MESH_FILENAME= other.su2"):
        assert line in text.splitlines(), line
    assert "{" not in text and "}" not in text


# ---------------------------------------------------------------- reading files
def test_read_mesh_summary():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "m.su2")
        open(path, "w").write(
            "NDIME= 2\nNELEM= 1900\n5 1 2 3 0\nNPOIN= 1185 1185\n0.1 0.2 0\n"
            "NMARK= 2\nMARKER_TAG= airfoil\nMARKER_ELEMS= 408\n3 1 2\n"
            "MARKER_TAG= farfield\nMARKER_ELEMS= 60\n3 5 6\n")
        info = study.read_mesh_summary(path)
    assert info == {"npoin": 1185, "nelem": 1900, "markers": {"airfoil": 408, "farfield": 60}}


def write_history(path, n=300, cd_end=0.0099, residual_end=-8.0):
    with open(path, "w") as f:
        f.write('"Time_Iter", "Outer_Iter", "Inner_Iter", "rms[Rho]", "CD", "CL"\n')
        for i in range(n):
            frac = i / (n - 1)
            f.write("0, 0, %d, %.4f, %.6f, %.6f\n"
                    % (i, -2 + (residual_end + 2) * frac,
                       cd_end + 0.01 * math.exp(-20 * frac),      # settles to cd_end
                       -0.002 * (1 - math.exp(-20 * frac))))


def test_read_history_handles_quotes_and_spaces():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "history.csv")
        write_history(path)
        cols = study.read_history(path)
    assert list(cols)[:4] == ["Time_Iter", "Outer_Iter", "Inner_Iter", "rms[Rho]"]
    assert len(cols["CD"]) == 300 and abs(cols["CD"][-1] - 0.0099) < 1e-12


def test_convergence_report():
    good = study.convergence_report([0.0099] * 150 + [0.00990, 0.009902], -8.0)
    assert good["converged"] and good["problems"] == []
    stalled = study.convergence_report([0.0099] * 200, -4.5)
    assert not stalled["converged"] and "residual" in stalled["problems"][0]
    drifting = study.convergence_report([0.0099 + 1e-4 * i / 100 for i in range(200)], -8.0)
    assert not drifting["converged"] and "CD still varied" in drifting["problems"][0]


def test_log_keeps_old_rows_when_columns_change():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "results_log.csv")
        open(path, "w").write("a,b\n1,2\n3,4,5\n")          # includes a row with an extra field
        rows = study.merge_into_log(path, {"a": "9", "b": "8", "new_column": "x"})
        text = open(path).read().splitlines()
    assert len(rows) == 3
    assert text[0] == "a,b,new_column" and text[3] == "9,8,x" and text[1] == "1,2,"


# ---------------------------------------------------------------- Gmsh calls (stand-in)
def make_fake_gmsh(write_marker_file=True):
    g = mock.MagicMock()
    g.__version__ = "4.15.2"
    g.isInitialized.return_value = False
    g.model.mesh.getNodes.return_value = (np.arange(1, 7), np.arange(18, dtype=float), None)
    g.model.mesh.getElements.return_value = (
        [2, 3], [None, None], [np.array([1, 2, 3]), np.array([1, 2, 3, 4])])
    if write_marker_file:
        def fake_write(path):
            open(path, "w").write(
                "NDIME= 2\nNELEM= 2\nNPOIN= 6 6\nNMARK= 2\n"
                "MARKER_TAG= airfoil\nMARKER_ELEMS= 10\nMARKER_TAG= farfield\nMARKER_ELEMS= 4\n")
        g.write.side_effect = fake_write
    return g


def test_build_mesh_uses_the_settings():
    g = make_fake_gmsh(write_marker_file=False)
    s = dict(study.DEFAULTS)
    info = study.build_mesh(g, s, 4.457e-6, "unused.su2")
    assert (info["nodes"], info["triangles"], info["quads"]) == (6, 1, 1)
    assert info["tri"].shape == (3, 3)                       # 1 triangle + 1 quad split in two
    names = [c.kwargs.get("name") for c in g.model.addPhysicalGroup.call_args_list]
    assert names == ["airfoil", "farfield", "fluid"]
    fld = g.model.mesh.field
    set_numbers = [c.args[1:] for c in fld.setNumber.call_args_list]
    assert ("Size", 4.457e-6) in set_numbers and ("Ratio", 1.2) in set_numbers
    assert ("Thickness", 0.04) in set_numbers and ("Quads", 1) in set_numbers
    assert fld.setString.call_args.args[1:] == ("F", "0.005 + 0.1*Sqrt((x-1)*(x-1) + y*y)")
    g.finalize.assert_called()                                # Gmsh is always closed again
    # the airfoil points use the surface cell size, the outer circle uses the far cell size
    sizes = [c.args[3] for c in g.model.geo.addPoint.call_args_list]
    assert sizes.count(0.005) == 100 + 98 and sizes.count(20.0) == 5


def test_build_mesh_without_tail_refinement_and_after_an_error():
    g = make_fake_gmsh(write_marker_file=False)
    s = dict(study.DEFAULTS)
    s["WAKE_LC"] = None
    study.build_mesh(g, s, 4.457e-6, "unused.su2")
    g.model.mesh.field.setString.assert_not_called()
    g2 = make_fake_gmsh(write_marker_file=False)
    g2.model.mesh.generate.side_effect = RuntimeError("meshing failed")
    expect_error(RuntimeError, study.build_mesh, g2, dict(study.DEFAULTS), 4.457e-6, "x.su2")
    g2.finalize.assert_called()                               # still closed after the failure


# ---------------------------------------------------------------- the whole run, with fakes
def run_main_with_fakes(overrides, residual_end=-8.0, **fake_kwargs):
    """Run main() in a temporary folder with Gmsh and SU2 replaced by stand-ins."""
    old_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            fake_gmsh = make_fake_gmsh()

            def fake_run_su2(su2_bin, cfg="rans.cfg"):
                assert os.path.exists(cfg)
                write_history("history.csv", residual_end=residual_end)
                return 0.5

            settings = study.apply_overrides(study.DEFAULTS, overrides)
            with mock.patch.object(study, "SETTINGS", settings), \
                 mock.patch.object(study, "OVERRIDES_USED", overrides), \
                 mock.patch.object(study, "ensure_gmsh", return_value=fake_gmsh), \
                 mock.patch.object(study, "ensure_su2", return_value=("/x/SU2_CFD", "v8.5.0")), \
                 mock.patch.object(study, "run_su2", side_effect=fake_run_su2), \
                 mock.patch("matplotlib.pyplot.show"):
                result = study.main()
            files = sorted(os.listdir(d))
            log = list(csv.DictReader(open("results_log.csv"))) if "results_log.csv" in files else []
            cfg = open("rans.cfg").read() if "rans.cfg" in files else None
        finally:
            os.chdir(old_dir)
    return result, files, log, cfg


def test_main_full_run():
    result, files, log, cfg = run_main_with_fakes({"LC_AIRFOIL": 0.0025})
    assert {"mesh.png", "convergence.png", "rans.cfg", "results_log.csv",
            "naca0012_bl.su2"} <= set(files)
    assert result["converged"] and abs(result["CD"] - 0.0099) < 1e-9
    assert len(log) == 1 and log[0]["lc_airfoil"] == "0.0025" and log[0]["converged"] == "YES"
    assert log[0]["gmsh_version"] == "4.15.2" and log[0]["su2_version"] == "v8.5.0"
    assert log[0]["low_mach"] == "YES" and log[0]["nodes"] == "6"


def test_main_flags_a_run_that_did_not_settle():
    result, files, log, cfg = run_main_with_fakes({}, residual_end=-4.5)
    assert not result["converged"] and log[0]["converged"] == "NO"


def test_main_mesh_only_stops_before_the_solver():
    result, files, log, cfg = run_main_with_fakes({"MESH_ONLY": True})
    assert "rans.cfg" not in files and "results_log.csv" not in files
    assert "mesh.png" in files and result["nodes"] == 6 and "CD" not in result


def test_main_rejects_a_mesh_without_the_wall():
    old_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            g = make_fake_gmsh(write_marker_file=False)
            g.write.side_effect = lambda p: open(p, "w").write(
                "NPOIN= 6\nMARKER_TAG= farfield\nMARKER_ELEMS= 4\n")
            with mock.patch.object(study, "ensure_gmsh", return_value=g), \
                 mock.patch.object(study, "ensure_su2", return_value=("/x", "v")), \
                 mock.patch("matplotlib.pyplot.show"):
                expect_error(RuntimeError, study.main, text="airfoil")
        finally:
            os.chdir(old_dir)


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    failed = 0
    for name, func in tests:
        try:
            func()
            print("PASS", name)
        except Exception as e:                  # report every failure, not just the first
            failed += 1
            print("FAIL", name, "->", type(e).__name__, e)
    print("\n%d tests, %d failed" % (len(tests), failed))
    sys.exit(1 if failed else 0)
