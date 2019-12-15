#!/usr/bin/env python3

"""Miscellaneous AMICI Python interface tests"""

import amici
import sys
import copy
import os
import unittest
import sympy as sp
import libsbml
from tempfile import TemporaryDirectory

from amici.petab_import import constant_species_to_parameters


class TestAmiciMisc(unittest.TestCase):
    """TestCase class various AMICI Python interface functions"""

    def setUp(self):
        self.default_path = copy.copy(sys.path)
        self.resetdir = os.getcwd()

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        self.test_parameterScalingFromIntVector()

    def test_parameterScalingFromIntVector(self):
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

    def test_csc_matrix(self):
        """Test sparse CSC matrix creation"""
        matrix = sp.Matrix([[1, 0], [2, 3]])
        symbolColPtrs, symbolRowVals, sparseList, symbolList, sparseMatrix = \
            amici.ode_export.csc_matrix(matrix, 'a')

        assert symbolColPtrs == [0, 2, 3]
        assert symbolRowVals == [0, 1, 1]
        assert sparseList == sp.Matrix([[1], [2], [3]])
        assert symbolList == ['a0', 'a1', 'a2']
        assert str(sparseMatrix) == 'Matrix([[a0, 0], [a1, a2]])'

    def test_csc_matrix_empty(self):
        """Test sparse CSC matrix creation for empty matrix"""
        matrix = sp.Matrix()
        symbolColPtrs, symbolRowVals, sparseList, symbolList, sparseMatrix = \
            amici.ode_export.csc_matrix(matrix, 'a')
        print(symbolColPtrs, symbolRowVals, sparseList, symbolList, sparseMatrix)
        assert symbolColPtrs == []
        assert symbolRowVals == []
        assert sparseList == sp.Matrix(0, 0, [])
        assert symbolList == []
        assert str(sparseMatrix) == 'Matrix(0, 0, [])'

    def test_csc_matrix_vector(self):
        """Test sparse CSC matrix creation from matrix slice"""
        matrix = sp.Matrix([[1, 0], [2, 3]])
        symbolColPtrs, symbolRowVals, sparseList, symbolList, sparseMatrix = \
            amici.ode_export.csc_matrix(matrix[:, 0], 'a')

        assert symbolColPtrs == [0, 2]
        assert symbolRowVals == [0, 1]
        assert sparseList == sp.Matrix([[1], [2]])
        assert symbolList == ['a0', 'a1']
        assert str(sparseMatrix) == 'Matrix([[a0], [a1]])'

        '''Test continuation of numbering of symbols'''
        symbolColPtrs, symbolRowVals, sparseList, symbolList, sparseMatrix = \
            amici.ode_export.csc_matrix(matrix[:, 1], 'a',
                                        base_index=len(symbolList))

        assert symbolColPtrs == [0, 1]
        assert symbolRowVals == [1]
        assert sparseList == sp.Matrix([[3]])
        assert symbolList == ['a2']
        assert str(sparseMatrix) == 'Matrix([[0], [a2]])'

    def test_sbml2amici_no_observables(self):
        """Test model generation works for model without observables"""

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

    def test_constant_species_to_parameters(self):
        """test conversion from species to constant parameters"""

        document = libsbml.SBMLDocument(3, 1)
        model = document.createModel()
        model.setTimeUnits("second")
        model.setExtentUnits("mole")
        model.setSubstanceUnits('mole')

        s = model.createSpecies()
        s.setId('x1')
        s.setConstant(True)
        s.setInitialConcentration(1.0)

        # Ensure species is replaced in reactions
        r = model.createReaction()
        r.setId('r1')
        # Add multiple instances to ensure all are removed
        for (coeff, name) in [(1, 'x1'), (1, 'x1')]:
            species_ref = r.createReactant()
            species_ref.setSpecies(name)
            species_ref.setStoichiometry(coeff)
        for (coeff, name) in [(1, 'x1'), (1, 'x1')]:
            species_ref = r.createProduct()
            species_ref.setSpecies(name)
            species_ref.setStoichiometry(coeff)
        for name in ['x1', 'x1']:
            species_ref = r.createModifier()
            species_ref.setSpecies(name)

        constant_species_to_parameters(model)

        assert len(list(model.getListOfParameters())) == 1
        assert len(list(model.getListOfSpecies())) == 0
        assert len(list(r.getListOfReactants())) == 0
        assert len(list(r.getListOfProducts())) == 0
        assert len(list(r.getListOfModifiers())) == 0

    def test_hill_function_dwdx(self):
        """Kinetic laws with Hill functions, may lead to NaNs in the Jacobian
        if involved states are zero if not properly arranged symbolically.
        Test that what we are applying the right sympy simplification."""

        w = sp.sympify('Pow(x1, p1) / (Pow(x1, p1) + a)')
        dwdx = w.diff(sp.Symbol('x1'))

        # Verify that without simplification we fail
        with self.assertRaises(ZeroDivisionError):
            with sp.evaluate(False):
                res = dwdx.subs({'x1': 0.0})
            _ = str(res)

        # Test that powsimp does the job
        dwdx = sp.powsimp(dwdx)
        with sp.evaluate(False):
            res = dwdx.subs({'x1': 0.0})
        _ = str(res)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciMisc())
    unittest.main()
