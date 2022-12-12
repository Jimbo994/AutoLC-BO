import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths
import matplotlib.pyplot as plt
from scipy.stats import norm

def create_signal(tR, W, tlim):
    """
    Generate signal from

    :param tR: List of retention times
    :param W: List of respective peak widths
    :param tlim: Length of signal
    :return: Signal with Gaussian peaks
    """
    sigmas = np.array(W)/4
    # resolution of the signal might be excessive, not sure
    x = np.linspace(0, tlim, 100000)
    signal = np.zeros(100000)
    for i in range(len(tR)):
        signal += norm.pdf(x, loc=tR[i], scale=sigmas[i])
    return x, signal

def detect_peaks(x, y, height_thresh=0, plot=False):
    """
    peak detection on a 1D signal.

    Parameters
    ----------
    :param x: ndarray Measurement timepoints
    :param y: ndarray Measurement intensities
    :param height_thresh: int Required height of peaks (see find_peaks function from Scipy
    :param plot: Bool If true, will show a plot of the peaks and widths on signal

    Returns
    ----------
    :peaks: ndarray timepoints of detected peaks in 'x' that satisfy all given conditions
    :widths: ndarray widths of detected peaks

    Note: Code relies heavily on Scipy functions, look into documentation for more information.
     Height threshold might need some tuning.
    """

    peaks, heights = find_peaks(y, height=height_thresh)  # , width=[0, width_thresh], threshold=0)
    prominences, left_bases, right_bases = peak_prominences(y, peaks)

    # Set peak bases correctly, Scipy seems to get this wrong
    for i in range(1, len(peaks)):
        if heights['peak_heights'][i - 1] > heights['peak_heights'][i]:
            right_bases[i - 1] = left_bases[i]
        elif heights['peak_heights'][i - 1] < heights['peak_heights'][i]:
            left_bases[i] = right_bases[i - 1]

    # Set peak prominences to peak heights
    offset = heights['peak_heights']

    # Find widths at respective heights
    # Calculate widths at x[peaks] - offset * rel_height at height 0.135 (4 sigma)
    widths, h_eval, left_ips, right_ips = peak_widths(
        y, peaks,
        rel_height=0.865,  # This is 4 Sigma
        prominence_data=(offset, left_bases, right_bases)
    )

    dt = x[-1] / len(x)

    if plot == True:
        # With time-axis
        xvals = x

        # Plot utility
        fig, ax = plt.subplots(1, 1)
        dt = x[-1] / len(x)

        ax.plot(x, y, linewidth=1)
        ax.set_ylabel('Absorption intensity')
        ax.set_xlabel('time (mins)')

        ax.plot(xvals[peaks], y[peaks], ".")
        ax.hlines(h_eval, left_ips * dt + xvals[0], xvals[peaks], color='C1', label='i_a_i')
        ax.hlines(h_eval, xvals[peaks], right_ips * dt + xvals[0], color='C2', label='i_b_i')
        ax.vlines(xvals[peaks], 0, heights['peak_heights'], alpha=0.5, color='black', linewidth=1)

        plt.show()
    return x[peaks], widths * dt


# To test the code
if __name__ == "__main__":
    from chromatogram_given_profile import plot_chromatogram_given_gradient_profile
    from crf import crf, capped_sum_of_resolutions

    # Well separated sample, essentially tRs = peaks, widths_pred = widths
    tRs = np.array([11, 14, 20, 34,  44, 50, 65])
    widths_pred = np.array([1, 2, 2.5, 3, 4, 3, 2])
    x, y = create_signal(tRs, widths_pred, 70)

    peaks, widths = detect_peaks(x,y,height_thresh=0, plot=True)

    #print(tRs, widths_pred)
    #print(peaks, widths)

    print(capped_sum_of_resolutions(tRs, widths_pred, [0, 0.5]))
    print(capped_sum_of_resolutions(peaks, widths, [0, 0.5]))

    # Poorly separated sample, essentially tRs != peaks, widths_pred != widths
    tRs = np.array([11, 14, 20, 34,  44, 50, 65])
    widths_pred = np.array([1, 2, 2.5, 3, 4, 3, 2])*2
    x, y = create_signal(tRs, widths_pred, 70)

    peaks, widths = detect_peaks(x,y,height_thresh=0, plot=True)

    # print(tRs, widths_pred)
    # print(peaks, widths)

    print(capped_sum_of_resolutions(tRs, widths_pred, [0, 0.5]))
    print(capped_sum_of_resolutions(peaks, widths, [0, 0.5]))

