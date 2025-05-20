import unittest
import math
from cnd_calculator import (
    calculate_Rd,
    calculate_geometric_mean,
    calculate_V_values,
    calculate_CND_indices,
    calculate_CND_r_squared,
    calculate_cnd
)

class TestCNDCalculator(unittest.TestCase):

    def test_calculate_Rd(self):
        self.assertAlmostEqual(calculate_Rd({'N': 3.0, 'P': 0.5, 'K': 2.0}), 100.0 - (3.0 + 0.5 + 2.0), places=7)
        self.assertAlmostEqual(calculate_Rd({}), 100.0, places=7)
        self.assertAlmostEqual(calculate_Rd({'N': 5.0}), 95.0, places=7)
        self.assertAlmostEqual(calculate_Rd({'N': 100.0}), 0.0, places=7)
        self.assertAlmostEqual(calculate_Rd({'N': 110.0}), -10.0, places=7)


    def test_calculate_geometric_mean(self):
        # Typical case
        data1 = {'N': 3.0, 'P': 0.5, 'K': 2.0, 'Rd': 94.5}
        expected1 = (3.0 * 0.5 * 2.0 * 94.5)**(1/4)
        self.assertAlmostEqual(calculate_geometric_mean(data1), expected1, places=7)

        # Single component
        data2 = {'N': 10.0}
        self.assertAlmostEqual(calculate_geometric_mean(data2), 10.0, places=7)

        # Includes a zero
        data3 = {'N': 3.0, 'P': 0, 'K': 2.0, 'Rd': 95.0}
        self.assertAlmostEqual(calculate_geometric_mean(data3), 0.0, places=7)
        
        # Empty input - current implementation returns 0.0
        data4 = {}
        self.assertAlmostEqual(calculate_geometric_mean(data4), 0.0, places=7)


    def test_calculate_V_values(self):
        # Typical case
        data1 = {'N': 3.0, 'P': 0.5, 'K': 2.0, 'Rd': 94.5}
        g1 = (3.0 * 0.5 * 2.0 * 94.5)**(1/4)
        result1 = calculate_V_values(data1, g1)
        self.assertAlmostEqual(result1['N'], math.log(3.0 / g1), places=7)
        self.assertAlmostEqual(result1['P'], math.log(0.5 / g1), places=7)
        self.assertAlmostEqual(result1['K'], math.log(2.0 / g1), places=7)
        self.assertAlmostEqual(result1['Rd'], math.log(94.5 / g1), places=7)

        # G is 0
        data2 = {'N': 3.0, 'P': 0.5}
        result2 = calculate_V_values(data2, 0)
        self.assertTrue(math.isnan(result2['N']))
        self.assertTrue(math.isnan(result2['P']))

        # Concentration is 0 or negative
        data3 = {'N': 0, 'P': -0.5, 'K': 2.0}
        g3 = 10.0 # Assuming a non-zero G for this test part
        result3 = calculate_V_values(data3, g3)
        self.assertTrue(math.isnan(result3['N'])) # ln(0/G)
        self.assertTrue(math.isnan(result3['P'])) # ln(-0.5/G)
        self.assertAlmostEqual(result3['K'], math.log(2.0/g3), places=7)

        # G = 1
        data4 = {'N': 3.0, 'P': 0.5}
        g4 = 1.0
        result4 = calculate_V_values(data4, g4)
        self.assertAlmostEqual(result4['N'], math.log(3.0), places=7)
        self.assertAlmostEqual(result4['P'], math.log(0.5), places=7)
        
        # Concentration is non-positive, G non-zero
        data5 = {'N': 0.0, 'P': 2.0}
        g5 = calculate_geometric_mean(data5) # This will be 0.0
        result5 = calculate_V_values(data5, g5) # G will be 0, so all NaN
        self.assertTrue(math.isnan(result5['N']))
        self.assertTrue(math.isnan(result5['P']))

        data6 = {'N': -1.0, 'P': 2.0}
        # g6 = calculate_geometric_mean(data6) # This would be complex or error
        # Instead, let's test calculate_V_values directly with a placeholder G
        # as geometric mean of negative numbers is problematic for this model.
        # The function itself guards against log(negative or zero).
        result6 = calculate_V_values(data6, 10.0)
        self.assertTrue(math.isnan(result6['N']))
        self.assertAlmostEqual(result6['P'], math.log(2.0/10.0))


    def test_calculate_CND_indices(self):
        V_values = {'N': 1.0, 'P': -0.5, 'K': 0.0, 'X': 2.0, 'Y': math.nan}
        V_star = {'N': 0.8, 'P': -0.6, 'K': 0.1, 'Z': 3.0} # X missing, Z not in V_values
        SD_star = {'N': 0.2, 'P': 0.1, 'K': 0.0, 'Z': 0.5} # K has SD_star = 0

        result = calculate_CND_indices(V_values, V_star, SD_star)

        self.assertAlmostEqual(result['N'], (1.0 - 0.8) / 0.2, places=7)
        self.assertAlmostEqual(result['P'], (-0.5 - (-0.6)) / 0.1, places=7)
        self.assertTrue(math.isnan(result['K'])) # SD_star is 0
        self.assertTrue(math.isnan(result['X'])) # X missing in V_star/SD_star
        self.assertTrue(math.isnan(result['Y'])) # V_value was nan

        # Nutrient in V_star/SD_star but not in V_values (should not appear in result)
        self.assertNotIn('Z', result)
        
        # Empty V_values
        self.assertEqual(calculate_CND_indices({}, V_star, SD_star), {})
        
        # V_values has item, but not in V_star
        self.assertTrue(math.isnan(calculate_CND_indices({'A':1}, V_star, SD_star)['A']))


    def test_calculate_CND_r_squared(self):
        # Typical case
        cnd_indices1 = {'N': 1.0, 'P': -2.0, 'K': 0.5}
        expected1 = 1.0**2 + (-2.0)**2 + 0.5**2
        self.assertAlmostEqual(calculate_CND_r_squared(cnd_indices1), expected1, places=7)

        # One index is math.nan
        cnd_indices2 = {'N': 1.0, 'P': math.nan, 'K': 0.5}
        expected2 = 1.0**2 + 0.5**2
        self.assertAlmostEqual(calculate_CND_r_squared(cnd_indices2), expected2, places=7)

        # All indices are math.nan
        cnd_indices3 = {'N': math.nan, 'P': math.nan}
        self.assertAlmostEqual(calculate_CND_r_squared(cnd_indices3), 0.0, places=7)

        # Empty CND_indices
        cnd_indices4 = {}
        self.assertAlmostEqual(calculate_CND_r_squared(cnd_indices4), 0.0, places=7)


    def test_calculate_cnd_main_function(self):
        # Comprehensive test case (example values)
        sample_conc = {'N': 2.8, 'P': 0.25, 'K': 2.5}
        # Rd = 100 - (2.8 + 0.25 + 2.5) = 100 - 5.55 = 94.45
        # nutrient_concentrations_with_Rd = {'N': 2.8, 'P': 0.25, 'K': 2.5, 'Rd': 94.45}
        # G = (2.8 * 0.25 * 2.5 * 94.45)^(1/4) = (165.2875)^(1/4) approx 3.5859
        # V_N = ln(2.8 / G) = ln(2.8 / 3.5859) = ln(0.78078) approx -0.2475
        # V_P = ln(0.25 / G) = ln(0.25 / 3.5859) = ln(0.06971) approx -2.6637
        # V_K = ln(2.5 / G) = ln(2.5 / 3.5859) = ln(0.69714) approx -0.3608
        # V_Rd = ln(94.45 / G) = ln(94.45 / 3.5859) = ln(26.338) approx 3.2710

        # Manually calculated for this example (approximations)
        # Using more precise G: (2.8 * 0.25 * 2.5 * 94.45)**(1/4) = 3.585900199
        # V_N  = math.log(2.8 / 3.585900199)   # -0.247511
        # V_P  = math.log(0.25 / 3.585900199)  # -2.663697
        # V_K  = math.log(2.5 / 3.585900199)   # -0.360787
        # V_Rd = math.log(94.45 / 3.585900199) #  3.271995

        high_yield_means = {'N': -0.20, 'P': -2.70, 'K': -0.30, 'Rd': 3.25}
        high_yield_sds = {'N': 0.1, 'P': 0.15, 'K': 0.12, 'Rd': 0.10}

        # I_N = (-0.247511 - (-0.20)) / 0.1 = -0.047511 / 0.1 = -0.47511
        # I_P = (-2.663697 - (-2.70)) / 0.15 = 0.036303 / 0.15 = 0.24202
        # I_K = (-0.360787 - (-0.30)) / 0.12 = -0.060787 / 0.12 = -0.506558
        # I_Rd = (3.271995 - 3.25) / 0.10 = 0.021995 / 0.10 = 0.21995

        # r_squared = (-0.47511)^2 + (0.24202)^2 + (-0.506558)^2 + (0.21995)^2
        #           = 0.2257295 + 0.0585737 + 0.256601 + 0.048378
        #           = 0.5892822

        cnd_indices, r_squared = calculate_cnd(sample_conc, high_yield_means, high_yield_sds)

        self.assertAlmostEqual(cnd_indices['N'], -0.47511, places=5)
        self.assertAlmostEqual(cnd_indices['P'], 0.24202, places=5)
        self.assertAlmostEqual(cnd_indices['K'], -0.506558, places=5) # Corrected expected sign based on V_K and V_star_K
        self.assertAlmostEqual(cnd_indices['Rd'], 0.21995, places=5)
        self.assertAlmostEqual(r_squared, 0.58928, places=5)


        # Test edge case: empty sample_nutrient_concentrations
        # Rd = 100.0
        # G = (100.0)**(1/1) = 100.0
        # V_Rd = ln(100.0 / 100.0) = 0
        # I_Rd = (0 - V_star_Rd) / SD_star_Rd = (0 - 3.25) / 0.10 = -32.5
        # r_squared = (-32.5)^2 = 1056.25
        # Other indices will be nan if their V* and SD* are not defined for non-existing sample values.
        # The current CND_indices function will only compute for keys in V_values.
        # So, only Rd will be in cnd_indices.
        
        empty_sample_conc = {}
        cnd_indices_empty, r_squared_empty = calculate_cnd(empty_sample_conc, high_yield_means, high_yield_sds)
        
        self.assertIn('Rd', cnd_indices_empty)
        self.assertNotIn('N', cnd_indices_empty) # N, P, K are not in the sample
        self.assertNotIn('P', cnd_indices_empty)
        self.assertNotIn('K', cnd_indices_empty)

        expected_I_Rd_empty = (0.0 - high_yield_means['Rd']) / high_yield_sds['Rd']
        self.assertAlmostEqual(cnd_indices_empty['Rd'], expected_I_Rd_empty, places=7)
        self.assertAlmostEqual(r_squared_empty, expected_I_Rd_empty**2, places=7)

        # Test case where a nutrient is in sample but not in reference data
        sample_conc_missing_ref = {'N': 2.8, 'P': 0.25, 'K': 2.5, 'Xtra': 1.0} # Xtra not in high_yield_means/sds
        # Rd = 100 - (2.8+0.25+2.5+1.0) = 100 - 6.55 = 93.45
        cnd_indices_missing, r_squared_missing = calculate_cnd(sample_conc_missing_ref, high_yield_means, high_yield_sds)
        self.assertTrue(math.isnan(cnd_indices_missing['Xtra']))
        
        # Check that r_squared calculation ignores this nan
        # Recalculate expected r_squared without Xtra
        # G_missing = (2.8 * 0.25 * 2.5 * 1.0 * 93.45)**(1/5) approx 3.4019
        # V_N_m  = math.log(2.8 / G_missing) approx -0.1959
        # V_P_m  = math.log(0.25 / G_missing) approx -2.6097
        # V_K_m  = math.log(2.5 / G_missing) approx -0.3085
        # V_Rd_m = math.log(93.45 / G_missing) approx 3.3134
        # V_Xtra_m = math.log(1.0 / G_missing) approx -1.2243

        # I_N_m = (-0.1959 - (-0.20)) / 0.1 = 0.041
        # I_P_m = (-2.6097 - (-2.70)) / 0.15 = 0.602
        # I_K_m = (-0.3085 - (-0.30)) / 0.12 = -0.0708
        # I_Rd_m = (3.3134 - 3.25) / 0.10 = 0.634
        # I_Xtra_m = nan
        # r_sq_m_expected = (0.041)**2 + (0.602)**2 + (-0.0708)**2 + (0.634)**2
        #                 = 0.001681 + 0.362404 + 0.00501264 + 0.401956 = 0.77105364
        
        # Using actual calculated values from the function for N, P, K, Rd
        r_sq_m_expected = 0.0
        for nutrient in ['N', 'P', 'K', 'Rd']:
            if not math.isnan(cnd_indices_missing[nutrient]):
                 r_sq_m_expected += cnd_indices_missing[nutrient]**2
        self.assertAlmostEqual(r_squared_missing, r_sq_m_expected, places=5)


if __name__ == '__main__':
    unittest.main()
