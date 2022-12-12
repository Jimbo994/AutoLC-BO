import numpy as np

import torch

from bo_code.MaternSortKernel import MaternSortKernel

try:
    from botorch.sampling.samplers import SobolQMCNormalSampler
except:
    from botorch.sampling.normal import SobolQMCNormalSampler

from gpytorch.kernels import ScaleKernel, MaternKernel, RBFKernel, SpectralMixtureKernel, PolynomialKernel, LinearKernel
from gpytorch.priors import GammaPrior

from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from botorch.utils.transforms import normalize, unnormalize
from botorch.models.transforms.outcome import Standardize

from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import ExpectedImprovement, UpperConfidenceBound, NoisyExpectedImprovement, qNoisyExpectedImprovement

from botorch.optim.optimize import optimize_acqf, optimize_acqf_list
from botorch.utils.multi_objective.box_decompositions.non_dominated import NondominatedPartitioning
from botorch.acquisition.multi_objective.monte_carlo import qExpectedHypervolumeImprovement

tkwargs = (
    {  # Tkwargs is a dictionary contaning data about data type and data device
        "dtype": torch.double,
        "device": torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    }
)

def BO_round(bounds, norm_bounds, scores, pars):
    """
    This function fits a Gaussian process to a single objective.
    Then uses the qNoisyExpectedImprovement acquisition function to get experiment
    to perform next.

    @param bounds: bounds of optimizable parameters
    @param norm_bounds: normalized bounds of optimizable (i.e. [0,1]^N)
    @param scores: list of previously obtained scores, i.e. train_Y
    @param pars: list of previously performed parameters, i.e train_X
    @return: pars, the parameters to evaluate next.
    """
    X = torch.from_numpy(pars).to(**tkwargs)
    # normalize input pars
    train_X = normalize(X, bounds)
    train_Y = torch.from_numpy(scores).unsqueeze(-1).to(**tkwargs)
    # standardize outcomes
    dim = len(pars[0])
    if dim == 2:
        dim_half = 1
    else:
        dim_half = int(dim / 2)

    # Initialize Gaussian Process and MarginalLogLikelihood
    gp = SingleTaskGP(train_X, train_Y, outcome_transform=Standardize(m=1))

    # custom GP equipped with MaternSortKernel
    gp.covar_module = ScaleKernel(
        MaternSortKernel(
            nu=2.5,
            # power=4,
            # num_mixtures=5,
            ard_num_dims=dim,
            batch_shape=gp._aug_batch_shape,
            lengthscale_prior=GammaPrior(3.0, 6.0),
        ),
        batch_shape=gp._aug_batch_shape,
        outputscale_prior=GammaPrior(2.0, 0.15),
    )

    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)

    # Fit parameters
    fit_gpytorch_model(mll)

    # print('outputscale ', gp.covar_module.outputscale.detach().numpy())
    # print('noise ', gp.likelihood.noise_covar.noise.detach().numpy())
    # print('lenghtscales ', gp.covar_module.base_kernel.lengthscale.detach().numpy())
    # # Initialize Acquisition functions
    sampler = SobolQMCNormalSampler(1024)
    qNEI = qNoisyExpectedImprovement(gp, train_X, sampler)

    # Optimize acq. fun.
    candidate, acq_value = optimize_acqf(
        qNEI, bounds=norm_bounds, q=1, num_restarts=20, raw_samples=512, sequential=False
    )

    candidate[0, :dim_half], _ = torch.sort(candidate[0, :dim_half])
    candidate[0, dim_half:], _ = torch.sort(candidate[0, dim_half:])
    pars = unnormalize(candidate, bounds).numpy()

    timepoints = pars[0, dim_half:]
    # In some cases timepoints are similar and chemstation and rm code does not like this,
    # this is an easy fix, better would be to have inequality constraints.
    for i in range(len(timepoints) - 1):
        if np.allclose(timepoints[i], timepoints[i + 1], rtol=0, atol=1e-2):
            timepoints[i + 1] += np.random.uniform(low=0.01, high=0.1)
    pars[0, dim_half:] = timepoints
    return pars


def MOBO_round(bounds, norm_bounds, scores, pars, ref_point):
    """
    This function fits indepented Gaussian Process to multiple objectives.
    Then uses the qExpectedHyperVolumeImprovement acquisition function to get experiment
    to perform next.
    @param bounds: bounds of optimizable parameters
    @param norm_bounds: normalized bounds of optimizable (i.e. [0,1]^N)
    @param scores: list of previously obtained scores, i.e. train_Y
    @param pars: list of previously performed parameters, i.e train_X
    @return: pars, the parameters to evaluate next.
    @param ref_point: Reference point from which to measure hypervolume
    @return:
    """
    X = torch.from_numpy(pars)
    # normalize input pars
    train_X = normalize(X, bounds)
    train_Y = torch.from_numpy(scores)
    shape = train_Y.shape
    # standardize outcomes
    # train_Y = standardize(Y)
    dim = len(pars[0])
    if dim == 2:
        dim_half = 1
    else:
        dim_half = int(dim / 2)
    # Initialize Gaussian Process and MarginalLogLikelihood
    gp = SingleTaskGP(train_X, train_Y, outcome_transform=Standardize(m=shape[1]))

    # custom GP equipped with MaternSortKernel
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

    # print('outputscale ', gp.covar_module.outputscale.detach().numpy())
    # print('noise ', gp.likelihood.noise_covar.noise.detach().numpy())
    # print('lenghtscales ', gp.covar_module.base_kernel.lengthscale.detach().numpy())
    # # Initialize Acquisition functions
    qehvi_sampler = SobolQMCNormalSampler(128)
    partitioning = NondominatedPartitioning(ref_point=ref_point, Y=train_Y)
    qEHVI = qExpectedHypervolumeImprovement(
        model=gp,
        ref_point=ref_point,  # use known reference point
        partitioning=partitioning,
        sampler=qehvi_sampler,
    )

    # Optimize acq. fun.
    candidate, acq_value = optimize_acqf(
        qEHVI, bounds=norm_bounds, q=1, num_restarts=20, raw_samples=512, sequential=True
    )

    candidate[0, :dim_half], _ = torch.sort(candidate[0, :dim_half])
    candidate[0, dim_half:], _ = torch.sort(candidate[0, dim_half:])
    pars = unnormalize(candidate, bounds).numpy()
    timepoints = pars[0, dim_half:]
    # In some cases timepoints are similar and chemstation and rm code does not like this,
    # this is an easy fix, better would be to have inequality constraints.
    for i in range(len(timepoints) - 1):
        if np.allclose(timepoints[i], timepoints[i + 1], rtol=0, atol=1e-2):
            timepoints[i + 1] += np.random.uniform(low=0.01, high=0.1)
    pars[0, dim_half:] = timepoints
    return pars