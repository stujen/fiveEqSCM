import numpy as np

from U_FaIR.fiveFAIR import *

from temp_test_fair.testFAIR import *


def test_oxfair_impulse_response_to_carbon():
    time = np.array([0, 1, 2, 3])
    input_emissions = np.array([[10, 0, 0, 0],[0,0,0,0],[0,0,0,0]])

    Cexpected,RFexpected,Texpected = oxfair_test(input_emissions)
    Cresult,RFresult,Tresult = oxfair(input_emissions)

    # assert result == expected
    np.testing.assert_allclose(Cresult, Cexpected)

# test under constant emissions
# test where pulse isn't in year zero

