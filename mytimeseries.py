import numpy as np

import pdb

def power_spectral_density(var):
    """
    calculate power spectral density
    """

    n = len(var)
    Y = np.fft.fft(var)/n  # fft computing and normalization
    Y = Y[range(int(n/2))]
    k = np.arange(n)
    Ts = 30*60  # sampling interval (sec)
    Fs = 1./Ts
    T = n/Fs
    frq = k/T  # two sides frequency range
    frq = frq[range(int(n/2))]  # one side frequency range

    return frq, abs(Y)
