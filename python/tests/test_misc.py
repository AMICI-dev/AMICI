"""Miscellaneous AMICI Python interface tests"""

import os
import subprocess
from tempfile import TemporaryDirectory

import amici
import libsbml
import pytest
import sympy as sp


def test_parameter_scaling_from_int_vector():
    """Ensure we can generate a ParameterScaling vector from Python"""

    scale_vector = amici.parameterScalingFromIntVector(
        [
            amici.ParameterScaling_log10,
            amici.ParameterScaling_ln,
            amici.ParameterScaling_none
        ])

    assert scale_vector[0] == amici.ParameterScaling_log10
    assert scale_vector[1] == amici.ParameterScaling_ln
    assert scale_vector[2] == amici.ParameterScaling_none


def test_sbml2amici_no_observables():
    """Test model generation works for model without observables"""

    # test model
    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits('mole')
    c1 = model.createCompartment()
    c1.setId('C1')
    model.addCompartment(c1)
    s1 = model.createSpecies()
    s1.setId('S1')
    s1.setCompartment('C1')
    model.addSpecies(s1)

    sbml_importer = amici.sbml_import.SbmlImporter(sbml_source=model,
                                                   from_file=False)
    tmpdir = TemporaryDirectory()
    sbml_importer.sbml2amici(modelName="test",
                             output_dir=tmpdir.name,
                             observables=None)


def test_hill_function_dwdx():
    """Kinetic laws with Hill functions, may lead to NaNs in the Jacobian
    if involved states are zero if not properly arranged symbolically.
    Test that what we are applying the right sympy simplification."""

    w = sp.sympify('Pow(x1, p1) / (Pow(x1, p1) + a)')
    dwdx = w.diff(sp.Symbol('x1'))

    # Verify that without simplification we fail
    with pytest.raises(ZeroDivisionError):
        with sp.evaluate(False):
            res = dwdx.subs({'x1': 0.0})
        _ = str(res)

    # Test that powsimp does the job
    dwdx = sp.powsimp(dwdx)
    with sp.evaluate(False):
        res = dwdx.subs({'x1': 0.0})
    _ = str(res)


def test_cmake_compilation(sbml_example_presimulation_module):
    """Check that CMake build succeeds for one of the models generated during
    Python tests"""

    source_dir = os.path.dirname(sbml_example_presimulation_module.__path__[0])

    cmd = f"set -e; cd {source_dir}; mkdir -p build; cd build; "\
          "cmake ..; make"
    subprocess.run(cmd, shell=True, capture_output=True, check=True)
