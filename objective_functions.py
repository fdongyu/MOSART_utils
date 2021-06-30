import numpy as np
from datetime import datetime

import pdb

def convert2date(intime):
    """
    convert datetime to date
    """
    intime = [itime.date() for itime in intime]
    return np.asarray(intime)

def intersection(lst1, lst2):
    """
    find the matches values in two lists
    """
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def getOverlap(time1, var1, time2, var2):
    """
    return the matches variables with the same length
    time1: discharge
    time2: conc
    need to convert conc to date first
    """
    #time1_date = convert2date(time1)
    #time2_date = convert2date(time2)

    #time1_date = time1_date.tolist()
    #time2_date = time2_date.tolist()

    time1 = time1.tolist()
    time2 = time2.tolist()

    both = intersection(time1, time2)

    ind1 = [time1.index(x) for x in both]
    ind2 = [time2.index(x) for x in both]

    return np.asarray(time1)[ind1], var1[ind1], np.asarray(time2)[ind2], var2[ind2]

def R2(time1, var1, time2, var2):
    """
    calculate the R2 metric for the simulation
    """
    time1, var1, time2, var2 = getOverlap(time1, var1, time2, var2)
    #pdb.set_trace()
    return (np.corrcoef(var1, var2)[0,1])**2

def nashsutcliffe(t1, evaluation, t2, simulation):
    """
    https://github.com/thouska/spotpy/blob/master/spotpy/objectivefunctions.py
    Nash-Sutcliffe model efficinecy
        .. math::
         NSE = 1-\\frac{\\sum_{i=1}^{N}(e_{i}-s_{i})^2}{\\sum_{i=1}^{N}(e_{i}-\\bar{e})^2} 
    :evaluation: Observed data to compared with simulation data.
    :type: list
    :simulation: simulation data to compared with evaluation data
    :type: list
    :return: Nash-Sutcliff model efficiency
    :rtype: float
    """
    t1, evaluation, t2, simulation = getOverlap(t1, evaluation, t2, simulation)

    if len(evaluation) == len(simulation):
        s, e = np.array(simulation), np.array(evaluation)
        # s,e=simulation,evaluation
        mean_observed = np.nanmean(e)
        # compute numerator and denominator
        numerator = np.nansum((e - s) ** 2)
        denominator = np.nansum((e - mean_observed)**2)
        # compute coefficient
        return 1 - (numerator / denominator)

    else:
        logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan

