import unittest as ut
import numpy as np
import numpy.testing as npt
from DOSclass import (DensityOfStates, H_TO_EV)

class TestDOSClass(ut.TestCase):

    def setUp(self):
        self.pdos = DensityOfStates.from_cp2k_pdos_output("./testDOS/el-k1-1.pdos")
        self.dos_1 = DensityOfStates.from_cp2k_output("./testDOS/outfile_1.out", "1")
        self.dos_2 = DensityOfStates.from_cp2k_output("./testDOS/outfile_2.out", "1", sigma=0.001, stepsize=1)
        return super().setUp()

    def test_bins(self) -> None:
        self.assertRaises(ValueError, DensityOfStates.calculate_bins, "1", 2)
        self.assertRaises(ValueError, DensityOfStates.calculate_bins, 1, "2")
        self.assertRaises(ValueError, DensityOfStates.calculate_bins, 2, 1)
        self.assertRaises(ValueError, DensityOfStates.calculate_bins, 1, 2, "1")
        self.assertRaises(ValueError, DensityOfStates.calculate_bins, 1, 2, -1)
        self.assertTrue(isinstance(DensityOfStates.calculate_bins(1,10,1), int))
        self.assertEqual(DensityOfStates.calculate_bins(1,10,1), 12)

    def test_pdos_values(self) -> None:
        # el-k1-1.pdos has exactly 1540 molecular orbital on which atomic orbital's coefficients are evaluated.
        self.assertEqual(len(self.pdos.parameter_vector), 1540)
        # testing the first and the last element of the energy values.
        self.assertAlmostEqual(self.pdos.parameter_vector[0], -0.806592*H_TO_EV)
        self.assertAlmostEqual(self.pdos.parameter_vector[-1], 0.144732*H_TO_EV)
        # for each MO, there is a specific weight given by the squared sum of its AO coefficients.
        self.assertEqual(len(self.pdos.weights), 1540)
        self.assertAlmostEqual(self.pdos.weights[0], 0.000052932754740)
        self.assertAlmostEqual(self.pdos.weights[-1], 0.0489148595823616)
        # no values for the Fermi energy is evaluated.
        self.assertEqual(len(self.pdos.fermi_energies), 0)

    def test_out_values(self) -> None:
        # outfile_1 has 1800 KS and 2 Fermi energy
        self.assertEqual(len(self.dos_1.parameter_vector), 1800)
        self.assertEqual(len(self.dos_1.fermi_energies), 2)
        self.assertEqual(self.dos_1.fermi_energies[0], 1.135131)
        self.assertEqual(self.dos_1.fermi_energies[-1], 1.130812)
        self.assertEqual(self.dos_1.parameter_vector[0], (-0.71032018*H_TO_EV))
        self.assertEqual(self.dos_1.parameter_vector[-1], (0.24501299*H_TO_EV))
        self.assertEqual(self.dos_1.norm_factor, 1800)
        self.assertEqual(self.dos_1.stepsize, DensityOfStates.DOS_STEP)
        self.assertEqual(self.dos_1.sigma, DensityOfStates.DOS_SIGMA)
        self.assertEqual(self.dos_1.parameter_minimum, np.min(self.dos_1.parameter_vector))
        self.assertEqual(self.dos_1.parameter_maximum, np.max(self.dos_1.parameter_vector))
        self.assertEqual(len(self.dos_1.weights), 1800)
        npt.assert_equal(self.dos_1.weights, np.ones(1800))
        # outfile_2 has 108000 KS and 12 Fermi energies. Sigma and stepsize are different than the default values.
        self.assertEqual(len(self.dos_2.parameter_vector), 10800)
        self.assertEqual(len(self.dos_2.fermi_energies), 12)
        npt.assert_allclose(self.dos_2.fermi_energies, 1.135131)
        self.assertEqual(self.dos_2.parameter_vector[0], (-0.71032018*H_TO_EV))
        self.assertEqual(self.dos_2.parameter_vector[-1], (0.24381002*H_TO_EV))
        self.assertEqual(self.dos_2.norm_factor, 10800)
        self.assertEqual(self.dos_2.norm_factor, 10800)
        self.assertEqual(self.dos_2.stepsize, 1)
        self.assertEqual(self.dos_2.sigma, 0.001)
        self.assertEqual(self.dos_2.parameter_minimum, np.min(self.dos_2.parameter_vector))
        self.assertEqual(self.dos_2.parameter_maximum, np.max(self.dos_2.parameter_vector))
        self.assertEqual(len(self.dos_2.weights), 10800)
        npt.assert_equal(self.dos_2.weights, np.ones(10800))
    
    def tearDown(self):
        return super().tearDown()

if __name__ == "__main__":
    ut.main(verbosity=3)