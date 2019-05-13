#!/usr/bin/env python3

import amici
import os
import unittest


class TestTextODEImport(unittest.TestCase):
    """
    TestCase class for testing ODE import from a generic text file and simulation from AMICI python interface
    """

    def test_text_ode_import(self):

        ode_file = os.path.join('test_ODE_import', 'ode_input1.txt')
        sbml_file = os.path.join('test_ODE_import', 'sbml_output.xml')
        expected_result_file = os.path.join('test_ODE_import', 'sbml_output.xml')

        amici.import_from_txt(ode_file, sbml_file)

        with open(sbml_file, 'r') as f_in:
            sbml_contents = f_in.read()

        with open(expected_result_file, 'r') as f_in:
            expected_sbml_contents = f_in.read()

        self.assertEqual(expected_sbml_contents, sbml_contents)

        os.remove(sbml_file)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestTextODEImport())
    unittest.main()
