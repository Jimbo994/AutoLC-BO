import numpy as np
import math
import scipy.stats as ss
import heapq


"""
This file contains implementations of a number of chromatographic response
functions from the literature on this subject.
"""


def resolution(tR1, tR2, W1, W2):
    """Return resolution of two peaks."""
    average_peak_width = (W1 + W2)/2
    res = abs((tR2 - tR1)) / average_peak_width
    return(res)


def critical_resolution(tR_array, W_array):
    """
    Return the critical resolution (resolution of the
    least well-separated peak pair) of a collection of peaks.
    """

    resolutions = []

    # Calculate resolutions between all peak combinations
    for i in range(len(tR_array)-1):
        tR1 = tR_array[i]
        W1 = W_array[i]

        for j in range(i+1, len(tR_array)):
            tR2 = tR_array[j]
            W2 = W_array[j]

            res = resolution(tR1, tR2, W1, W2)
            print("res", res)
            resolutions.append(res)

    critical_res = min(resolutions)
    #return(critical_res)
    last_tR = max(tR_array)
    print("critical res", critical_res)
    return(min(critical_res, 2.5))


def order_tR_W(tR_array, W_array):
    """Order retention times and peak widths according to retention time."""

    zipped = list(zip(tR_array, W_array))
    sorted_zip = sorted(zipped, key=lambda x: x[0])
    unzipped = [list(t) for t in zip(*sorted_zip)]

    tR_list = unzipped[0]
    W_list = unzipped[1]

    tR_array = np.array(tR_list, dtype="float128")
    W_array = np.array(W_list, dtype="float128")

    return(tR_array, W_array)


def sum_of_adjacent_resolutions(tR_array, W_array):
    """Return the sum of resolutions of adjacent peaks."""

    # Define threshold for time penalty as a constant.
    threshold = 20

    adjacent_resolutions = []
    tR_array, W_array = order_tR_W(tR_array, W_array)

    # Loop over all peaks except the last one
    for i in range(len(tR_array) - 1):
        tR1 = tR_array[i]
        W1 = W_array[i]

        if(tR1 < threshold):

            tR2 = tR_array[i+1]
            W2 = W_array[i+1]

            res = resolution(tR1, tR2, W1, W2)
            adjacent_resolutions.append(res)

    return(sum(adjacent_resolutions))


def product_of_resolutions(tR_array, W_array):
    """Return the product of resolution of a collection of peaks."""

    prod = 1

    # Time PENALTY
    pen_threshold = 30

    # Compare all peaks
    for i in range(len(tR_array)-1):
        tR1 = tR_array[i]
        W1 = W_array[i]

        for j in range(i+1, len(tR_array)):
            tR2 = tR_array[j]
            W2 = W_array[j]
            # Time PENALTY
            if(tR1 < pen_threshold and tR2 < pen_threshold):
                res = resolution(tR1, tR2, W1, W2)
                prod = prod * res
            else:
                prod = prod * 0.001
    print("prod found: ", prod)
    return(prod)


def sum_of_resolutions(tR_array, W_array):
    """Return the sum of resolutions of all peaks."""

    sum = 0

    # Define time penalty
    pen_threshold = 30

    # Sum of resolutions between all peak combinations
    for i in range(len(tR_array)-1):
        tR1 = tR_array[i]
        W1 = W_array[i]

        for j in range(i+1, len(tR_array)):
            tR2 = tR_array[j]
            W2 = W_array[j]

            # Time penalty
            if(tR1 < pen_threshold and tR2 < pen_threshold):
                res = resolution(tR1, tR2, W1, W2)
                sum = sum + res

    return(sum)


def hao(tR_list, W_list, phi_list):
    # Penalty voor constraint geen negatieve slopes.
    if(phi_list != sorted(phi_list)):
        return(-100000)

    t_s = 50
    r_min = critical_resolution(tR_list, W_list)
    tR_max = max(tR_list)
    cf = round(r_min, 2) * 10**5 + round(t_s - tR_max, 2) * 10
    return(cf)


################################################################################
################################################################################
################################# PEAK PURITY ##################################
################################################################################
################################################################################


#https://stats.stackexchange.com/questions/311592/how-to-find-the-point-where-two-normal-distributions-intersect
def gaussian_intersection(m1, m2, std1, std2):
    '''
    Parameters:
        m1, m2: Means of Gaussians
        std1, std2: Standard deviations of Gaussians

    Returns:
        Points of intersection of 2 Gaussian curves
    '''
    a = 1/(2*std1**2) - 1/(2*std2**2)
    b = m2/(std2**2) - m1/(std1**2)
    c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
    return np.roots([a,b,c])


def get_non_overlapping_area_on_interval(interval_points, mus, sigmas):
    point1 = interval_points[0]
    point2 = interval_points[1]

    # Pick a random point on the interval (maybe in the middle)
    if(point1 == -math.inf):
        mid_point = point2 - 0.001
    elif(point2 == math.inf):
        mid_point = point1 + 0.001
    else:
        mid_point = point1 + (point2 - point1)/2

    curve_list = []
    # Evaluate all functions at that point
    for i in range(len(mus)):
        mu = mus[i]
        sigma = sigmas[i]
        y = ss.norm.pdf(mid_point, mu, sigma)
        curve_list.append(y)

    # Find the area underneath the highest and subtract area underneath
    # second highest

    # Integrate highest over interval
    max_value = max(curve_list)
    max_index = curve_list.index(max_value)
    mu_h = mus[max_index]
    sigma_h = sigmas[max_index]

    # If the retention time of the highest peak is larger than 20, we dont want to count it
    # this is to encourage faster retention times.
    delta_t = 20

    if(mu_h > delta_t):
        return(0)
    else:
        integral1 = ss.norm.cdf(point2, mu_h, sigma_h) - ss.norm.cdf(point1, mu_h, sigma_h)
        #print(integral1)

        # Integrate second highest over interval
        two_largest = heapq.nlargest(2, curve_list)
        second_max_value = min(two_largest)
        second_max_index = curve_list.index(second_max_value)
        mu_h2 = mus[second_max_index]
        sigma_h2 = sigmas[second_max_index]

        integral2 = ss.norm.cdf(point2, mu_h2, sigma_h2) - ss.norm.cdf(point1, mu_h2, sigma_h2)
        #print(integral2)

        # Return non-overlapping area between interval
        non_overlapping_area = integral1 - integral2
        #return((1/mu_h)*non_overlapping_area)
        return(non_overlapping_area)


def peak_purity(tR_list, W_list):
    # Cast to array
    tR_list = np.array(tR_list)
    W_list = np.array(W_list)

    mus = tR_list
    sigmas = W_list/4

    number_of_curves = len(mus)

    # For every curve, find the intersections with all other curves
    # Yields a list of intersections from left to right
    intersections_list = []
    i_range = range(number_of_curves - 1)
    for i in i_range:
        mu1 = mus[i]
        sigma1 = sigmas[i]

        if(mu1 < 30):

            range1 = range(i + 1, number_of_curves)

            for j in range1:
                mu2 = mus[j]
                sigma2 = sigmas[j]
                intersections = gaussian_intersection(mu1, mu2, sigma1, sigma2)
                intersections_list.extend(intersections)

    intersections_list.sort()
    # Add -inf and inf
    intersections_list.insert(0, -math.inf)
    intersections_list.append(math.inf)

    total_non_overlapping_area = 0

    # For each interval this yields, determine which curve is the highest
    # and which curve is the second highest in this interval
    for i in range(len(intersections_list) - 1):
        interval = [intersections_list[i], intersections_list[i + 1]]
        area = get_non_overlapping_area_on_interval(interval, mus, sigmas)
        total_non_overlapping_area = total_non_overlapping_area + area

    return(total_non_overlapping_area)
