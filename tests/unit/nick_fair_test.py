import numpy as np

from U_FaIR.fiveFAIR import *

from U_FaIR.testFAIR import *

def test_hfc_impulse_response():
    time = np.array([0, 1, 2, 3])
    input_emissions = np.array([10, 0, 0, 0])
    expected = 10 * np.exp(-time)

    result = calculate_hfc_conc(input_emissions, time, lifetime=1.0)

    # assert result == expected
    np.testing.assert_allclose(result, expected)


def test_oxfair_impulse_response_to_carbon():
    time = np.array([0, 1, 2, 3])
    input_emissions = np.array([[10, 0, 0, 0],[0,0,0,0],[0,0,0,0]])

    Cexpected,RFexpected,Texpected = oxfair_test(input_emissions, time, lifetime=1.0)
    Cresult,RFresult,Tresult = oxfair(input_emissions, time, lifetime=1.0)

    # assert result == expected
    np.testing.assert_allclose(Cresult, Cexpected)

# test under constant emissions
# test where pulse isn't in year zero

