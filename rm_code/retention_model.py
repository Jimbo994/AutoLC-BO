import math
import rm_code.peak_width as peak_width
from rm_code.crf import sort_peaks
import numpy as np

def get_k(k_0, S, phi):
    """Return retention factor k for a given phi."""
    k = k_0 * math.exp(-S * phi)
    return(k)


def get_B_list(phi_list, t_list):
    """Return list of slopes; 1 for each gradient segment."""
    B_list = []
    # Loop over all gradient segments except the last one (that is parallel to t-axis).
    for i in range(len(phi_list) - 1):
        # Calculate the slope of each gradient segment.
        B = (phi_list[i + 1] - phi_list[i])/(t_list[i + 1] - t_list[i])
        B_list.append(B)

    return(B_list)


def isocratic_retention_time(t_0, k):
    """Return retention time for a given analyte under isocratic elution."""
    t_R = t_0 * (1 + k)
    return(t_R)


def retention_time_multisegment_gradient(k_0, S, t_0, t_D, t_init, phi_list, t_list, N):
    """
    Return retention time for a given analyte under multisegment gradient elution.
    :k_0: k0 value for the compound.
    :S: S value for the compound.
    :t_0: Column dead time (defined in globals.py).
    :t_D: Dwell time (defined in globals.py).
    :t_init: Length in time of initial isocratic segment.
    :phi_list: List of phi values; one for each turning point in the gradient profile.
    :t_list: List of t values; one for each turning point in the gradient profile.
    :N: Column efficieny (defined in globals.py).
    :return: Retention time and peak width for the compound under the specified
             gradient profile.
    """
    # k_init is the retention factor before the gradient starts.
    phi_init = phi_list[0]
    k_init = get_k(k_0, S, phi_init)

    fractional_migr_before_gradient = (t_D + t_init)/(t_0 * k_init)

    # CASE 1: compound elutes during first isocratic segment.
    if(fractional_migr_before_gradient >= 1):

        # Calculate retention time
        # The analyte is eluted before the gradient reaches it.
        t_R = isocratic_retention_time(t_0, k_init)

        # Calculate peak width
        # Peak width is calculated using isocratic formula
        W = peak_width.get_peak_width_isocratic(N, t_0, k_init)
        return(t_R, W)

    # Otherwise, the analyte is eluted during the gradient program.

    # Get a list of slopes, one for each gradient segment except the final isocratic segment.
    B_list = get_B_list(phi_list, t_list)

    # Loop over slopes to calculate all a_i's.
    as_list = []
    a_ws_list = []
    for i, B in enumerate(B_list):
        phi_i = phi_list[i]
        phi_i_1 = phi_list[i+1]
        t_i = t_list[i]
        t_i_1 = t_list[i+1]

        k_phi_i = get_k(k_0, S, phi_i)
        k_phi_i_1 = get_k(k_0, S, phi_i_1)


        if(B == 0):
            # Retention time addend if B = 0
            a_i = (t_i_1 - t_i)/(t_0 * k_phi_i)
            # Peak width addend if B = 0
            a_i_w = ((t_i_1 - t_i) * (1 + k_phi_i)**2)/(k_phi_i**3)

        else:
            # Retention time addend if B != 0
            a_i = (1/(t_0 * S * B)) * ((1/k_phi_i_1) - (1/k_phi_i))

            # Peak width addend if B != 0
            a_i_w = (1/(S * B)) * ( (1/3)*((1/k_phi_i_1**3) - (1/k_phi_i**3)) +
                                    (1/k_phi_i_1**2 - 1/k_phi_i**2) +
                                    (1/k_phi_i_1 - 1/k_phi_i) )


        as_list.append(a_i)
        a_ws_list.append(a_i_w)

    cum_sum = 0
    cum_sum_w = 0
    elution_before_last_segment = False
    # Loop over a list until the sum reaches 1 and remember the last segment.
    for i, a in enumerate(as_list):
        if((cum_sum + a + fractional_migr_before_gradient) < 1):
            cum_sum = cum_sum + a
            cum_sum_w = cum_sum_w + a_ws_list[i]
        else:
            elution_segment = i
            elution_before_last_segment = True
            # At this point we want to break out of the loop
            break

    # CASE 2 (2 subcases) Compound is eluted during one of the other gradient segments.
    if(elution_before_last_segment == True):
        # Analyte elutes during another gradient segment.
        # Calculate retention time
        t_n = t_list[elution_segment]
        phi_n = phi_list[elution_segment]
        k_phi_n = get_k(k_0, S, phi_n)
        B_n = B_list[elution_segment]

        # Elution segment is isocratic.
        if(B_n == 0):
            # Slope of gradient segment during which analyte elutes is 0
            t_R = t_D + t_init + t_0 + t_n + t_0*k_phi_n * (1 - fractional_migr_before_gradient - cum_sum)

            # Calculate peak width
            cum_sum_w = cum_sum_w + ((t_R - t_0 - t_D - t_init - t_n) * (1 + k_phi_n)**2) / k_phi_n**3
            k_phi_n = get_k(k_0, S, phi_n)
            G = peak_width.get_G(cum_sum_w, k_phi_n, k_init, t_D, t_0, t_init)
            W = peak_width.get_peak_width_gradient(N, t_0, k_phi_n, G)
            return(t_R, W)

        # Elution segment has a gradient.
        else:
            # Slope of gradient segment during which analyte elutes is not 0
            # Calculate retention time
            t_R = t_0 + t_D + t_init + t_n + (1/(S*B_n)) * np.log(1 + t_0*S*B_n*k_phi_n*(1 - fractional_migr_before_gradient - cum_sum))

            # Calculate peak width
            phi_at_elution = phi_n + B_n * (t_R - t_0 - t_D - t_init - t_n)
            k_at_elution = get_k(k_0, S, phi_at_elution)

            cum_sum_w = cum_sum_w + (1/(S * B_n)) * (
                                    (1/3)* ((1/k_at_elution**3) - (1/k_phi_n**3)) +
                                    (1/k_at_elution**2 - 1/k_phi_n**2) +
                                    (1/k_at_elution - 1/k_phi_n)
                                    )

            G = peak_width.get_G(cum_sum_w, k_at_elution, k_init, t_D, t_0, t_init)
            #G = peak_width.poppe(k_0, S, B_n, t_0)
            W = peak_width.get_peak_width_gradient(N, t_0, k_at_elution, G)
            return(t_R, W)

    # CASE 3: Compound is eluted during last isocratic segment.
    else:
        # Analyte elutes during last isocratic segment.
        phi_n = phi_list[-1]
        k_phi_n = get_k(k_0, S, phi_n)
        t_n = t_list[-1]
        t_R = t_D + t_init + t_0 + t_n + t_0*k_phi_n * (1 - fractional_migr_before_gradient - cum_sum)

        # Add last part to peak peak_width
        cum_sum_w = cum_sum_w + ((t_R - t_0 - t_D - t_init - t_n) * (1 + k_phi_n)**2) / (k_phi_n**3)
        G = peak_width.get_G(cum_sum_w, k_phi_n, k_init, t_D, t_0, t_init)
        W = peak_width.get_peak_width_gradient(N, t_0, k_phi_n, G)
        return(t_R, W)


def compute_chromatogram(k0_list, S_list, t_0, t_D, t_init, phi_list, t_list, N):
    tR_list = []
    W_list = []

    # Calculate retention times and peak widths
    for i in range(len(k0_list)):
        k_0 = k0_list[i]
        S = S_list[i]
        tR, W = retention_time_multisegment_gradient(k_0, S, t_0, t_D, t_init, phi_list, t_list, N)
        tR_list.append(tR)
        W_list.append(W)

        # We need to do sorting so in the list are neighbours.
        tR_list, W_list = sort_peaks(tR_list, W_list)
    return np.array(tR_list), np.array(W_list)