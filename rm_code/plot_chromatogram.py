import sys
import rm_code.crf as crf
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


def plot_normal(x_range, mu=0, sigma=1, plot_individual_peaks=True):
    """Plot the normal distribution function for a given x range."""
    gaussian = norm.pdf(x_range, mu, scale=sigma)
    if(plot_individual_peaks):
        plt.plot(x_range, gaussian, linewidth=1)
    return(gaussian)


def plot_chromatogram(t_R_list, pw_list, phi_list, t_list, t_D, t_0, t_init, plot_individual_peaks=False):
    """
    Plot chromatogram and gradient profile given a list of retention times
    and peak widths and a gradient profile.

    :t_R_list: List of retention times.
    :pw_list: List of peak widths
    :t_0: Column dead time (defined in globals.py).
    :t_D: Dwell time (defined in globals.py).
    :t_init: Length in time of initial isocratic segment.
    :phi_list: List of phi values; one for each turning point in the gradient profile.
    :t_list: List of t values; one for each turning point in the gradient profile.
    :score: CRF score of the chromatogram.
    """
    #print(t_R_list)
    t_R_list, pw_list = crf.sort_peaks(t_R_list, pw_list)
    time_before_gradient_reaches_detector = t_D + t_init + t_0 # +t_0 voor graph van einde kolom
    end_time_profile = t_list[-1] + time_before_gradient_reaches_detector
    end_time_chromatogram = t_R_list[-1] + 2*pw_list[-1]

    # The last peak elutes after the last turning point of the gradient profile
    if(end_time_chromatogram > end_time_profile):
        case = 1
    # The last peak elutes before the last turning point of the gradient profile
    else:
        case = 2

    """Plot chromatogram."""
    fig, ax = plt.subplots()

    # The standard deviation is the pw/4
    sigmas = [pw/4 for pw in pw_list]

    if(case == 1):
        x = np.linspace(0, end_time_chromatogram, 100000)
    else:
        x = np.linspace(0, end_time_profile , 100000)

    #ax.set_xlabel(score)
    sum_of_ys = np.zeros(100000)

    for i, mean in enumerate(t_R_list):
        sigma = sigmas[i]
        y = plot_normal(x, mean, sigma, plot_individual_peaks)
        sum_of_ys = np.add(sum_of_ys, y)

    ax.set_xlim(xmin=0)
    ax.set_xlabel('t (min)')
    # Plot the chromatogram
    ax.plot(x, sum_of_ys, linewidth=1)


    """Plot gradient profile."""
    #phi_list = [phi + 0.01 for phi in phi_list]
    ax2=ax.twinx()
    t_list = [t + time_before_gradient_reaches_detector for t in t_list]

    # Plot dotted red line for time_before_gradient_reaches_detector
    ax2.plot([0, time_before_gradient_reaches_detector], [phi_list[0], phi_list[0]], linestyle="-", linewidth=0.8, color="red", alpha=0.7)

    # Plot the rest of the gradient profile
    ax2.plot(t_list, phi_list, marker=".", markersize=4, linestyle="-", linewidth=0.8, color="black", alpha=0.7)


    if(case == 1):
        ax2.plot([t_list[-1], end_time_chromatogram], [phi_list[-1], phi_list[-1]], marker=".", markersize=4, linestyle="-", linewidth=0.8, color="gray", alpha=0.7)
        ax2.set_xlim(xmin=0, xmax=end_time_chromatogram)
    else:
        ax2.set_xlim(xmin=0, xmax=end_time_profile)

    ax2.set_ylim(ymin=0, ymax=1.1)
    ax2.set_ylabel(r'$\phi$')

    #plt.savefig('fig.png', dpi = 300)
    plt.show()
