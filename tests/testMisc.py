#!/usr/bin/env python3

"""Miscellaneous AMICI Python interface tests"""

import amici
import sys
import copy
import os
import unittest


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


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciMisc())
    unittest.main()
