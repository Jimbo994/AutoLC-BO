import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
from botorch.sampling import SobolQMCNormalSampler
from botorch.utils.multi_objective import is_non_dominated, Hypervolume
import os

import shutil

from MaternSortKernel import MaternSortKernel

from gpytorch.kernels import ScaleKernel, MaternKernel
from gpytorch.priors import GammaPrior
import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from botorch.utils import standardize
from botorch.utils.transforms import normalize, unnormalize
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import ExpectedImprovement

from botorch.optim.optimize import optimize_acqf, optimize_acqf_list
from botorch.utils.multi_objective.box_decompositions.non_dominated import NondominatedPartitioning
from botorch.acquisition.multi_objective.monte_carlo import qExpectedHypervolumeImprovement

from baseline_als import baseline_als
from Slack_utils import send_message_to_slack, send_file_to_slack
from assymetric_resolution import gen_res_24
from infer_reference_point import infer_reference_point

def BO_loop(iteration):
    # Load in data
    scores = pd.read_csv('past_runs/scores.txt', header=None, dtype=float)
    pars = pd.read_csv('past_runs/pars.txt', delimiter=',', header=None, dtype=float)
    pars_dim = len(pars.values[0])

    # These need to be lower and upper bounds of the parameters  this needs some insight into the sample
    bounds = torch.stack([torch.tensor([0., 0., 0., 0., 0., 0.]), torch.tensor([100., 100., 100, 19., 19., 19.])])
    bounds_norm = torch.stack([torch.zeros(6), torch.ones(6)])
    dim = len(bounds_norm[0])
    X = torch.from_numpy(pars.values)
    train_X = normalize(X, bounds)
    Y = torch.from_numpy(scores.values)
    train_Y = standardize(Y)

    # Initialize Gaussian Process and MarginalLogLikelihood
    gp = SingleTaskGP(train_X, train_Y)

    # custom GP
    gp.covar_module = ScaleKernel(
        MaternSortKernel(
            nu=2.5,
            ard_num_dims=dim,
            batch_shape=gp._aug_batch_shape,
            lengthscale_prior=GammaPrior(3.0, 6.0),
        ),
        batch_shape=gp._aug_batch_shape,
        outputscale_prior=GammaPrior(2.0, 0.15),
    )

    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)

    # Fit parameters
    fit_gpytorch_model(mll, max_retries=40)

    # Saving kernel parameters to file, for error checking
    df = pd.DataFrame(
        data=[gp.covar_module.outputscale.detach().numpy(), gp.likelihood.noise_covar.noise.detach().numpy(),
              gp.covar_module.base_kernel.lengthscale.detach().numpy()],
        index=['outputscale', 'noise', 'lengthscale']
    )
    if iteration == 0:
        df.to_csv('past_runs/kernel_parameters.csv', mode='w', sep=',', encoding='utf-8', header=True, index=False)
    else:
        df.to_csv('past_runs/kernel_parameters.csv', mode='a', sep=',', encoding='utf-8', header=False, index=False)

    # # Initialize Acquisition functions
    EI = ExpectedImprovement(gp, train_Y.max())
    # Optimize acq. fun.
    candidate, acq_value = optimize_acqf(
        EI, bounds=bounds_norm, q=1, num_restarts=20, raw_samples=512,
    )
    candidate[0, :int(dim/2)], _ = torch.sort(candidate[0, :int(dim/2)])
    candidate[0, int(dim / 2):], _ = torch.sort(candidate[0, int(dim / 2):])

    pars = unnormalize(candidate, bounds).numpy()
    timepoints = pars[0, int(dim/2):]
    for i in range(len(timepoints)-1):
        if np.allclose(timepoints[i], timepoints[i+1], rtol=0, atol=1e-2):
            timepoints[i+1] += np.random.uniform(low=0.01, high=1)
    pars[0, int(dim/2):] = timepoints
    return pars

def BO_loop_MO(iteration):
    # Load in data
    scores = pd.read_csv('past_runs/scores.txt', header=None, dtype=float)
    pars = pd.read_csv('past_runs/pars.txt', delimiter=',', header=None, dtype=float)
    pars_dim = len(pars.values[0])

    bounds = torch.stack([torch.tensor([0., 0., 0., 0., 0., 0.]), torch.tensor([100., 100., 100., 19., 19., 19.])])    # phi_init (t_init fixed at 0.25) , phi_1 , dt_1 (Then actual time should be dt_1 + 0.25), phi_final ,dt_2 (Then actual time should be dt_2 + dt_1 + 0.25
    bounds_norm = torch.stack([torch.zeros(6), torch.ones(6)])
    dim = len(bounds_norm[0])

    X = torch.from_numpy(pars.values)
    train_X = normalize(X, bounds)
    Y = torch.from_numpy(scores.values)
    train_Y = standardize(Y)

    # Initialize Gaussian Process and MarginalLogLikelihood
    gp = SingleTaskGP(train_X, train_Y)

    # custom GP
    gp.covar_module = ScaleKernel(
        MaternSortKernel(
            nu=2.5,
            ard_num_dims=dim,
            batch_shape=gp._aug_batch_shape,
            lengthscale_prior=GammaPrior(3.0, 6.0),
        ),
        batch_shape=gp._aug_batch_shape,
        outputscale_prior=GammaPrior(2.0, 0.15),
    )

    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)

    # Fit parameters
    fit_gpytorch_model(mll, max_retries=40)

    # Saving kernel parameters to file, for error checking
    df = pd.DataFrame(
        data=[gp.covar_module.outputscale.detach().numpy(), gp.likelihood.noise_covar.noise.detach().numpy(),
              gp.covar_module.base_kernel.lengthscale.detach().numpy()],
        index=['outputscale', 'noise', 'lengthscale']
    )
    if iteration == 0:
        df.to_csv('past_runs/kernel_parameters.csv', mode='w', sep=',', encoding='utf-8', header=True, index=False)
    else:
        df.to_csv('past_runs/kernel_parameters.csv', mode='a', sep=',', encoding='utf-8', header=False, index=False)

    # pick reference point to be 10% worse than the previous measurement (see qEHVI paper from Balandat)
    # using the infere_reference_point function
    pareto_mask = is_non_dominated(Y)
    pareto_y = Y[pareto_mask]
    r_point = infer_reference_point(pareto_y)
    hv = Hypervolume(ref_point=r_point)
    volume = hv.compute(pareto_y)

    open("past_runs/ref_points_and_hv.txt", "a").write(
        str(r_point.numpy()[0]) + ',' + str(r_point.numpy()[1]) + ',' +str(volume) + '\n')

    qehvi_sampler = SobolQMCNormalSampler(num_samples=128)

    # partition non-dominated space into disjoint rectangles
    partitioning = NondominatedPartitioning(ref_point=r_point, Y=train_Y)
    acq_func = qExpectedHypervolumeImprovement(
        model=gp,
        ref_point=r_point,  # use known reference point
        partitioning=partitioning,
        sampler=qehvi_sampler,
    )
    # optimize
    candidate, _ = optimize_acqf(
        acq_function=acq_func,
        bounds=bounds_norm,
        q=1,
        num_restarts=20,
        raw_samples=512,  # used for intialization heuristic
        options={"batch_limit": 5, "maxiter": 200, "nonnegative": True},
        sequential=True,
    )
    candidate[0, :int(dim/2)], _ = torch.sort(candidate[0, :int(dim/2)])
    candidate[0, int(dim/2):], _ = torch.sort(candidate[0, int(dim / 2):])
    pars = unnormalize(candidate, bounds).numpy()
    timepoints = pars[0, int(dim/2):]
    for i in range(len(timepoints)-1):
        if np.allclose(timepoints[i], timepoints[i+1], rtol=0, atol=1e-2) or timepoints[i] > timepoints[i+1]:
            timepoints[i+1] += np.random.uniform(low=0.01, high=1)
    pars[0, int(dim/2):] = timepoints
    return pars

print(os.getcwd())
os.chdir('D:\Jim_Boelrijk\AutoLC\AutoLC-BO')
channel = 'autolc-bo'
mode = 'QUADRUPLE' #Can also be double
obj = 'SO'

iterations = 36
iteration = 0

if mode == 'QUADRUPLE':
    initial_parameters = np.loadtxt('past_runs/initial_parameters_4step.txt')

initial_iterations = len(initial_parameters)

try:
    while iteration < iterations:
        whosturn = np.loadtxt('WhosTurn.txt', dtype=str, encoding='utf-16')
        print('iteration is:', iteration, 'turn is', whosturn)

        if whosturn == 'LC':
            print('LC machine is running, current iteration:', iteration)
            time.sleep(30)
        elif whosturn=='busy':
            print('Measurement running')
            time.sleep(30)
        elif whosturn=='finished':
            print('We are finished')
            break
        elif whosturn == 'BO':
            # Check if we are still in scanning experiments
            print('turn is BO')
                # Check if very first run
            if iteration == 0:
                # write first parameters to file
                print('loading first initial experiment')
                if mode == 'QUADRUPLE':
                    np.savetxt('parameters_4step.txt', initial_parameters[iteration])
                    os.system('D:\Jim_Boelrijk\AutoLC\Bridge-BO-4step\Debug\Bridge-BO-4step.exe')
                # Set turn to LC
                turn = np.array(['LC'])
                np.savetxt('WhosTurn.txt', turn, fmt="%s", encoding='utf-16')
                iteration += 1
            else:
                print('storing  previous runs')
                # appending previous parameters
                if mode == 'QUADRUPLE':
                    par = pd.read_csv('parameters_4step.txt', header=None)

                # pars = pd.read_csv('parameters.txt', header=None)
                par.T.to_csv('past_runs/pars.txt', mode='a', header=None, index=False)

                ####### CODE FOR TOTAL SPECTRUM #########
                filename = 'MegaMix.csv'

                # Inspiration
                f_measurement = pd.read_csv(filename, delimiter=',', encoding='utf-16le', dtype=float, index_col=0)
                print(f_measurement.head())

                plt.plot(np.sum(f_measurement.values, axis=1), linewidth=1, label='measurement')
                plt.legend()
                plt.savefig('past_runs/uncorrected' +str(iteration) + '.pdf', dpi=300)
                plt.savefig('past_runs/uncorrected' +str(iteration) + '.png', dpi=300)

                plt.close()

                f_times = pd.read_csv(filename, delimiter=',', encoding='utf-16le', header=None, dtype=float, skiprows=1,
                                      usecols=[0])

                ### CODE to take out last 3.1 minutes of measurement (in that 3.1 minutes, the gradient is reset),
                # not necessary when dealing with blank
                cutoff = 3.1
                t_length = len(f_times.values)
                dt = f_times.values[-1]/t_length
                idx_tocut = int(np.round(cutoff/dt)[0])
                f_wavelengths = pd.read_csv(filename, delimiter=',', encoding='utf-16le', header=None, dtype=float,
                                            usecols=[i for i in range(1, 902)])
                wavelengths = f_wavelengths.iloc[0].values
                f_intensities = pd.read_csv(filename, delimiter=',', encoding='utf-16le', header=None, dtype=float, skiprows=1,
                                            usecols=[i for i in range(1, 902)])  #columns of wavelengths
                intensities = f_intensities.values

                # sum intensities at all wavelengths and cut off the cutoff.
                total_spectrum = np.sum(intensities, axis=1)[:t_length-idx_tocut]

                # Code for baseline correction, not used at the moment because we measure a blank
                # compute baseline
                baseline = baseline_als(total_spectrum, 100000, 0.001)
                correct_tot_spec = total_spectrum - baseline

                plt.plot(baseline, label='baseline', linewidth=0.7)
                plt.plot(total_spectrum, linewidth=0.7)
                plt.legend()
                plt.savefig('past_runs/corrected' +str(iteration) + '.pdf', dpi=300)
                plt.savefig('past_runs/corrected' +str(iteration) + '.png', dpi=300)

                plt.close()

                # gen_res_24 also does plotting for you!
                if obj == 'SO':
                    con_comps, total, num_peaks, _, last_peak = gen_res_24(f_times[:t_length-idx_tocut], correct_tot_spec,  iteration, par.values, height_thresh=1000, width_thresh=100) # I think we need to lower this
                if obj == 'MO':
                    con_comps, total, num_peaks, _, last_peak = gen_res_24(f_times[:t_length-idx_tocut], correct_tot_spec,  iteration, par.values, height_thresh=1000, width_thresh=100, plot=False) # I think we need to lower it slightly

                # Save files
                shutil.copyfile(filename, 'past_runs/measurement' +str(iteration) + '.csv')
                pd.DataFrame(total_spectrum).to_csv('past_runs/tot_spec' + str(iteration) + '.txt', header=False)
                pd.DataFrame(correct_tot_spec).to_csv('past_runs/corrected_tot_spec' + str(iteration) + '.txt', header=False)
                f_times.to_csv('past_runs/time_steps' + str(iteration) + '.txt', header=False)
                ###### END CODE FOR TOTAL SPECTRUM #########

                if mode == 'QUADRUPLE' and obj != 'MO':
                    open("past_runs/all_scores.txt", "a").write(
                        str(con_comps) + ',' + str(total) + ',' + str(num_peaks) + ',' + str(last_peak[0]) + '\n')
                    open("past_runs/scores.txt", "a").write(str(total) + '\n')
                if mode == 'QUADRUPLE' and obj == 'MO':
                    open("past_runs/all_scores.txt", "a").write(
                        str(con_comps) + ',' + str(total) + ',' + str(num_peaks) + ',' + str(last_peak[0]) + '\n')
                    open("past_runs/scores.txt", "a").write(str(total) + ',' + str(last_peak[0]) + '\n')

                # Sending information to Slack
                if mode == 'QUADRUPLE':
                    send_message_to_slack('Using parameters: phi1 ' + str(par.values[0][0]) + ' phi2: ' + str(
                        par.values[1][0]) + ' phi3 ' + str(par.values[2][0]) + ' t1 ' + str(par.values[3][0]) + ' t2 ' + str(par.values[4][0]) + ' t3 ' + str(par.values[5][0]), channel)

                send_message_to_slack('Number of separated peaks: ' + str(con_comps) + ' Sum of resolutions: ' + str(total) + ' Number of peaks detected: ' + str(num_peaks), channel)
                send_file_to_slack('past_runs/measurement' +str(iteration) + '.png', channel)
                send_file_to_slack('past_runs/uncorrected' + str(iteration) + '.png', channel)
                send_file_to_slack('past_runs/corrected' + str(iteration) + '.png', channel)

                if iteration < initial_iterations:
                    print('load next initial experiment')
                    if mode == 'QUADRUPLE':
                        np.savetxt('parameters_4step.txt', initial_parameters[iteration])
                        os.system('D:\Jim_Boelrijk\AutoLC\Bridge-BO-4step\Debug\Bridge-BO-4step.exe')
                    iteration += 1
                    turn = np.array(['LC'])
                    np.savetxt('WhosTurn.txt', turn, fmt="%s", encoding='utf-16')
                else:
                    print('execute BO algorithm')
                    if mode == 'QUADRUPLE' and obj != 'MO':
                        new_pars = BO_loop(iteration)
                        print('new_pars', new_pars, new_pars[0][0], type(new_pars))
                        np.savetxt('parameters_4step.txt', new_pars, delimiter='\n')
                        os.system('D:\Jim_Boelrijk\AutoLC\Bridge-BO-4step\Debug\Bridge-BO-4step.exe')
                    if mode == 'QUADRUPLE' and obj == 'MO':
                        new_pars = BO_loop_MO(iteration)
                        np.savetxt('parameters_4step.txt', new_pars, delimiter='\n')
                        os.system('D:\Jim_Boelrijk\AutoLC\Bridge-BO-4step\Debug\Bridge-BO-4step.exe')
                    iteration += 1
                    turn = np.array(['LC'])
                    np.savetxt('WhosTurn.txt', turn, fmt="%s", encoding='utf-16')
    turn = np.array(['finished'])
    np.savetxt('WhosTurn.txt', turn, fmt="%s", encoding='utf-16')
    print(turn)
    print('done')
except Exception as exception:
    print(exception)
    send_message_to_slack('Error so going into standby method! Errormessage: ' + str(exception), channel)
    # os.system('D:\Jim_Boelrijk\AutoLC\BridgeBO-standbye\Debug\BridgeBO-standbye.exe')
    exit()