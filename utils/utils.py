import random, os
import numpy as np
import torch

from botorch.models import SingleTaskGP
from gpytorch.mlls import ExactMarginalLogLikelihood

from botorch.utils.sampling import draw_sobol_samples
from botorch.cross_validation import gen_loo_cv_folds
from botorch.cross_validation import batch_cross_validation
from botorch.utils.transforms import normalize, unnormalize

from botorch.utils.multi_objective.box_decompositions.dominated import DominatedPartitioning


import matplotlib.pyplot as plt

def seed_everything(seed: int):
    # set all random seeds.
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = True

def generate_initial_sorted_data(bounds, n=6):
    # generate initial data using sobol sampling, and sort parameters to exclude negative gradients
    random_pars = draw_sobol_samples(bounds=bounds,n=n, q=1).squeeze(1)
    dim = len(bounds[0])
    random_pars_sorted = np.concatenate((np.sort(random_pars[:, :int(dim/2)]), np.sort(random_pars[:, int(dim/2):])), axis=1)
    return random_pars_sorted


def bo_to_rm(pars, fixed_phi_pars, fixed_time_pars):
    """
    convenience function to generate friendly format for the retention modeling code.
    The retention modeling code was created seperately from the bo code.
    """
    shape = pars.shape
    if shape[1] == 2:
        dim_half = 1
    if shape[1] == 3:
        dim_half = 2
    else:
        dim_half = int(shape[1] / 2)

    if len(fixed_phi_pars) > 1:
        phi_pars = np.concatenate((np.repeat(fixed_phi_pars[0], shape[0]).reshape(shape[0], -1), pars[:, :dim_half]),
                                  axis=1)
        phi_pars = np.concatenate((phi_pars, np.repeat(fixed_phi_pars[1], shape[0]).reshape(shape[0], -1)), axis=1)
    else:
        phi_pars = np.concatenate((pars[:, :dim_half], np.repeat(fixed_phi_pars[0], shape[0]).reshape(shape[0], -1)),
                                  axis=1)

    time_pars = np.concatenate((np.repeat(fixed_time_pars[0], shape[0]).reshape(shape[0], -1), pars[:, dim_half:]),
                               axis=1)
    time_pars = np.concatenate((time_pars, np.repeat(fixed_time_pars[1], shape[0]).reshape(shape[0], -1)), axis=1)
    return phi_pars, time_pars


def best_so_far(vals):
    shape = vals.shape
    for i in range(shape[0]):
        for j in range(shape[1] - 1):
            if vals[i, j + 1] < vals[i, j]:
                vals[i, j + 1] = vals[i, j]
    return vals


def ci(y, n_trials):
    return 1.96 * y.std(axis=0) / (np.sqrt(n_trials))

def hv(y, ref_point):
    hvs = []
    shape = y.shape
    for i in range(shape[0]):
        hv = []
        for j in range(1, shape[1] + 1):
            # compute hv for point up to now.
            # print(y[i, :j, :2])
            bd = DominatedPartitioning(ref_point=ref_point, Y=torch.tensor(y[i, :j, :2]))
            # print(y[i, :j, :2])
            volume_true = bd.compute_hypervolume().item()
            hv.append(volume_true)
        hvs.append(hv)
    return np.array(hvs)

def cv_val(pars, scores, bounds, norm_bounds, sigma, plot_results=True):
    """"

    Code to perform model leave-one-out cross validation

    """
    X = torch.from_numpy(pars)
    # normalize input pars
    train_X = normalize(X, bounds)
    train_Y = torch.from_numpy(np.array(scores)).unsqueeze(-1)

    cv_folds = gen_loo_cv_folds(train_X=train_X, train_Y=train_Y)
    # instantiate and fit model
    cv_results = batch_cross_validation(
    model_cls=SingleTaskGP,
    mll_cls=ExactMarginalLogLikelihood,
    cv_folds=cv_folds,
    )
    posterior = cv_results.posterior
    mean = posterior.mean
    cv_error = ((cv_folds.test_Y.squeeze() - mean.squeeze()) ** 2).mean()
    print(f"Cross-validation error: {cv_error : 4.2}")
    if plot_results:
        # get lower and upper confidence bounds
        lower, upper = posterior.mvn.confidence_region()

        # scatterplot of predicted versus test
        _, axes = plt.subplots(1, 1, figsize=(6, 4))
        minval = np.min(mean.numpy().flatten() - np.min(((upper-lower)/2).numpy().flatten()))
        maxval = np.max(mean.numpy().flatten() + np.max(((upper-lower)/2).numpy().flatten()))
        plt.plot([minval, maxval], [minval, maxval], 'k', label="true objective", linewidth=2)

        axes.set_xlabel("Actual")
        axes.set_ylabel("Predicted")

        axes.errorbar(
            x=cv_folds.test_Y.numpy().flatten(),
            y=mean.numpy().flatten(),
            xerr=1.96*sigma,
            yerr=((upper-lower)/2).numpy().flatten(),
            fmt='*'
        );
    return cv_error