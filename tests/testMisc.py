#!/usr/bin/env python3

"""Miscellaneous AMICI Python interface tests"""

import amici
import sys
import copy
import os
import unittest
import sympy as sp

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


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciMisc())
    unittest.main()
