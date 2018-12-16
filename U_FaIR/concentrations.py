import numpy as np


def calculate_hfc_conc(emissions, time, lifetime):
    return emissions[0] * np.exp(-time)
