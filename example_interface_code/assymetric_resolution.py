import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths
import matplotlib.pyplot as plt

def gen_res_24(x, y, iteration, pars, bound=1.5, height_thresh=1000, width_thresh=0, plot=False):
    """"
    TODO: fill docstring with reference, etc
    """
    peaks, heights = find_peaks(y, height=height_thresh, width=[0, width_thresh], threshold=0)

    peaks_time, heights_time = find_peaks(y, height=height_thresh*2.5, width=[0, width_thresh], threshold=0)
    prominences, left_bases, right_bases = peak_prominences(y, peaks)

    # Set bases correctly
    for i in range(1, len(peaks)):

        if heights['peak_heights'][i - 1] > heights['peak_heights'][i]:
            right_bases[i - 1] = left_bases[i]
        elif heights['peak_heights'][i - 1] < heights['peak_heights'][i]:
            left_bases[i] = right_bases[i - 1]

    # set prominences to peak heights
    offset = heights['peak_heights']

    # find widths at respective heights
    # all_widths, all_h_eval, all_left_ips, all_right_ips = [], [], [], []
    # Calculate widths at x[peaks] - offset * rel_height at height 0.135
    widths, h_eval, left_ips, right_ips = peak_widths(
        y, peaks,
        rel_height=0.865,
        prominence_data=(offset, left_bases, right_bases)
    )

    a = peaks - left_ips
    b = right_ips - peaks

    N = 16 * ((peaks / (a + b)) ** 2)

    A = b / a

    resolutions_graph = np.ones(len(peaks))
    residues = 0
    h = heights['peak_heights']
    all_res = []

    # plot utility
    fig, ax = plt.subplots(1, 1)
    dt = x.values[-1] / len(x.values)


    #With time-axis
    xvals = x.values
    ax.plot(x, y, linewidth=1)
    ax.set_ylabel('Absorption intensity')
    ax.set_xlabel('time (mins)')

    ax.plot(xvals[peaks], y[peaks], ".")
    ax.plot(xvals[np.max(peaks_time)], y[np.max(peaks_time)], ".", color='red')
    ax.plot(xvals[np.max(peaks)], y[np.max(peaks)], ".", color='pink')
    ax.hlines(h_eval, left_ips*dt + xvals[0], xvals[peaks], color='C1', label='i_a_i')
    ax.hlines(h_eval, xvals[peaks], right_ips*dt + xvals[0], color='C2', label='i_b_i')
    ax.vlines(xvals[peaks], 0, heights['peak_heights'], alpha=0.5, color='black', linewidth=1)

    # with no time
    # ax.plot(peaks, y[peaks], "x")
    # ax.hlines(h_eval, left_ips, peaks, color='C1', label='i_a_i')
    # ax.hlines(h_eval, peaks, right_ips, color='C2', label='i_b_i')
    # ax.vlines(peaks, 0, heights['peak_heights'], alpha=0.5, color='black', linewidth=1)

    # # NOTE!! Only works for certain segments now (only triple) TODO: add more functionality
    # time_matrix = [0, 0.25, pars[3]+0.25, pars[4]+0.25, pars[5]+0.25, 20.25]
    # phi_matrix = [0,0, pars[0], pars[1], pars[2], 100]
    num_pars = len(pars)
    num_phi_pars = int(len(pars)/2)
    time_matrix = np.concatenate((np.concatenate((np.array([0, 0.25]), pars[num_phi_pars:].flatten()+0.25)), np.array([20.25])))
    phi_matrix = np.concatenate((np.concatenate((np.array([0, 0]), pars[:num_phi_pars].flatten())), np.array([100])))

    ax2 = ax.twinx()
    ax2.plot(time_matrix, phi_matrix, color='black', alpha=0.8, linewidth=1)
    ax2.scatter(time_matrix[2:num_phi_pars+2], phi_matrix[2:num_phi_pars+2], color='blue')
    ax2.set_ylabel('mobile phase composition B')
    ax2.set_ylim(-2,100)

    plt.savefig('past_runs/measurement' +str(iteration) + '.png', dpi=300)
    plt.savefig('past_runs/measurement' +str(iteration) + '.pdf', dpi=300)

    if plot == True:
        plt.show()
    else:
        plt.close()

    #plt.show()

    for i in range(1, len(peaks)):

        if (h[i] / h[i - 1]) <= np.exp(-2):
            sqrt_term = 0
        else:
            sqrt_term = np.sqrt(1 + 0.5 * np.log(h[i] / h[i - 1]))

        # Eq 24.
        i_Rs_ji = ((peaks[i] - peaks[i - 1]) * (1 + A[i - 1]) * (1 + A[i]) * np.sqrt(N[i - 1] * N[i])) / (
                    ((4 * A[i - 1] * peaks[i - 1]) * (1 + A[i]) * np.sqrt(N[i])) + 4 * peaks[i] * (
                        1 + A[i - 1]) * np.sqrt(N[i - 1]) * sqrt_term)
        j_Rs_ji = ((peaks[i] - peaks[i - 1]) * (1 + A[i - 1]) * (1 + A[i]) * np.sqrt(N[i - 1] * N[i])) / (
                    ((4 * A[i - 1] * peaks[i - 1]) * (1 + A[i]) * np.sqrt(N[i]) * sqrt_term) + 4 * peaks[i] * (
                        1 + A[i - 1]) * np.sqrt(N[i - 1]))

        if i_Rs_ji < j_Rs_ji:
            all_res.append(i_Rs_ji)
            res = i_Rs_ji
        else:
            all_res.append(j_Rs_ji)
            res = j_Rs_ji
        if i_Rs_ji >= bound:
            resolutions_graph[i - 1] = 1
        elif 1 < res < bound:
            resolutions_graph[i - 1] = 0
            residues += res / 1.5
        elif res < 1:
            resolutions_graph[i - 1] = 0
            residues += res / 1.5
    print('tot peaks', len(peaks))
    con_components = np.sum(resolutions_graph)
    total = con_components + residues
    print('con components', con_components, 'residues', total, 'peaks', len(peaks))
    return con_components, total, len(peaks), all_res, - xvals[np.max(peaks)]
