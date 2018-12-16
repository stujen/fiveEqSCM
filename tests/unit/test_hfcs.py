import numpy as np

from U_FaIR.concentrations import calculate_hfc_conc

def test_hfc_impulse_response():
    time = np.array([0, 1, 2, 3])
    input_emissions = np.array([10, 0, 0, 0])
    expected = 10 * np.exp(-time)

    result = calculate_hfc_conc(input_emissions, time, lifetime=1.0)

    # assert result == expected
    np.testing.assert_allclose(result, expected)

# test under constant emissions
# test where pulse isn't in year zero
