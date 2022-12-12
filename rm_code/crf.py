import math
import numpy as np


def resolution(tR1, tR2, W1, W2):
    """Return the resolution of 2 peaks, given tR and W for both peaks."""
    resolution = ((2*abs(tR2-tR1))/(W1+W2))
    return(resolution)


def sigmoid(x, a, b):
    """Return a sigmoidal transformation of x."""
    sigmoid = 1/(1 + np.exp(-a*x + b))
    return(sigmoid)


def sort_peaks(retention_times, peak_widths):
    """
    Sort peaks based on retention time
    and return sorted retention time list and peak width list.
    """
    number_of_peaks = len(retention_times)

    # Create a list of tuples, one for each peak (rt, W)
    peak_tuple_list = []
    for i in range(number_of_peaks):
        peak_tuple = (retention_times[i], peak_widths[i])
        peak_tuple_list.append(peak_tuple)
    # Sort according to first element
    peak_tuples_sorted = sorted(peak_tuple_list, key=lambda x: x[0])

    retention_times = []
    peak_widths = []
    for i in range(number_of_peaks):
        retention_times.append(peak_tuples_sorted[i][0])
        peak_widths.append(peak_tuples_sorted[i][1])

    return(retention_times, peak_widths)


def crf(retention_times, peak_widths, phi_list):
    """
    Return CRF score for a chromatogram characterized by a list of retention
    times and a corresponding list of peak widths.
    """

    N = len(retention_times)

    # Parameters sigmoidal transformations
    b0 = 3.93
    b1 = 3.66
    b2 = -0.0406
    b3 = -4.646

    resolutions = []
    sigmoid_resolutions = []

    # Sort retention times and peak widths
    retention_times, peak_widths = sort_peaks(retention_times, peak_widths)

    prod_S = 1

    # Loop over all neighboring peak pairs and get S. Multiply together.
    for i in range(N - 1):
        tR1 = retention_times[i]
        tR2 = retention_times[i+1]
        W1 = peak_widths[i]
        W2 = peak_widths[i+1]

        R = resolution(tR1, tR2, W1, W2)
        S = sigmoid(R, b0, b1)
        prod_S = prod_S * S

        sigmoid_resolutions.append(S)
        resolutions.append(R)

    # Create f and g
    f = prod_S ** (1/(N-1))

    # Get T
    tR_last = retention_times[-1]
    W_last = peak_widths[-1]
    T = tR_last + 0.5*(W_last)

    g = sigmoid(T, b2, b3)

    score = f * g

    # Optional penalty for gradient segments with negative slope
    # Comment out to remove penalty:
    if(sorted(phi_list) != phi_list):
        return(0.8 * score)

    return(score)

def capped_sum_of_resolutions(retention_times, peak_widths, phi_list=None, max_time=60, min_res=1, max_res=1.5):
    """
    Resolution equation as defined in eq. 5 and 6 of
     https://chemrxiv.org/engage/chemrxiv/article-details/62e2a383e7fc8f9e388caabc
    Uses symmetric resolution equation.

    :param retention_times: ndarray containing retention_times
    :param peak_widths: ndarray containing peak widths
    :param phi_list: list of phi points
    :param max_time: int maximum allowed time
    :param min_res: float minimum required resolution
    :param max_res: float maximum required resolution
    :return: float score
    """

    N = len(retention_times)

    resolutions = np.zeros(len(retention_times))

    mask = retention_times < max_time
    for i in range(N - 1):
        # check if tR1 is later then max_time, if yes we can stop
        tR1 = retention_times[i]
        tR2 = retention_times[i+1]
        W1 = peak_widths[i]
        W2 = peak_widths[i+1]

        resolutions[i] = resolution(tR1, tR2, W1, W2)
        if resolutions[i] < min_res:
            resolutions[i]=0
        if min_res < resolutions[i] < max_res:
            resolutions[i] = resolutions[i]/max_res
        if resolutions[i] > max_res:
            resolutions[i]=1

    # zero out scores of peaks that have eluted after max_time
    resolutions = resolutions*mask
    #print(resolutions, 'resolutions')

    score = resolutions.sum()
    # Optional penalty for gradient segments with negative slope
    # Comment out to remove penalty:
    if phi_list is not None:
        if(sorted(phi_list) != phi_list):
            return(0.8 * score)

    return(score)

def tyteca_eq_11(retention_times, peak_widths, max_time=60, min_time=2, max_res=1.6, prefacs=[1,1,1]):
    """
    Implements the CRF defined in Tyteca 2014, 10.1016/J.CHROMA.2014.08.014, Category II-A, equation 11.

    :param retention_times: ndarray containing retention_times
    :param peak_widths: ndarray containing peak widths
    :param max_time: int maximum allowed time
    :param min_time: int minimum allowed time
    :param max_res: float maximum required resolution
    :param prefacs: prefactors that dictate importance of each term.
    :return: float CRF score
    """
    N = len(retention_times)
    nobs_term = N**prefacs[0]

    resolutions = np.zeros(N)
    for i in range(N - 1):
        # check if tR1 is later then max_time, if yes we can stop
        tR1 = retention_times[i]
        tR2 = retention_times[i+1]
        W1 = peak_widths[i]
        W2 = peak_widths[i+1]

        resolutions[i] = resolution(tR1, tR2, W1, W2)
        if resolutions[i] > max_res:
            resolutions[i] = max_res

    res_term = resolutions.sum()
    max_time_term = prefacs[1] + np.abs(max_time - np.max(retention_times))
    min_time_term = prefacs[2] * (min_time-np.min(retention_times))

    return nobs_term + res_term - max_time_term + min_time_term

def tyteca_eq_24(retention_times, peak_widths, max_res=1.6):
    """
    Implements the CRF defined in Tyteca 2014, 10.1016/J.CHROMA.2014.08.014, Category I-B, equation 24.
    :param retention_times: ndarray containing retention_times
    :param peak_widths: ndarray containing peak widths
    :param max_res: float maximum required resolution
    :return: float CRF score
    """
    N = len(retention_times)

    resolutions = np.zeros(N)
    for i in range(N - 1):
        # check if tR1 is later then max_time, if yes we can stop
        tR1 = retention_times[i]
        tR2 = retention_times[i+1]
        W1 = peak_widths[i]
        W2 = peak_widths[i+1]

        resolutions[i] = resolution(tR1, tR2, W1, W2)
        if resolutions[i] > max_res:
            resolutions[i] = max_res
        res_term = np.sum(resolutions)

        return N + (res_term / (max_res*(N-1)))
