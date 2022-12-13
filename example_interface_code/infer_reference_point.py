from typing import Optional

import torch
from botorch.exceptions import BotorchError
from torch import Tensor

def infer_reference_point(
    pareto_Y: Tensor,
    max_ref_point: Optional[Tensor] = None,
    scale: float = 0.1,
    scale_max_ref_point: bool = False,
) -> Tensor:
    r"""Get reference point for hypervolume computations.

    This sets the reference point to be `ref_point = nadir - 0.1 * range`
    when there is no pareto_Y that is better than the reference point.

    [Ishibuchi2011]_ find 0.1 to be a robust multiplier for scaling the
    nadir point.

    Note: this assumes maximization of all objectives.

    Args:
        pareto_Y: A `n x m`-dim tensor of Pareto-optimal points.
        max_ref_point: A `m` dim tensor indicating the maximum reference point.
        scale: A multiplier used to scale back the reference point based on the
            range of each objective.
        scale_max_ref_point: A boolean indicating whether to apply scaling to
            the max_ref_point based on the range of each objective.

    Returns:
        A `m`-dim tensor containing the reference point.
    """

    MIN_Y_RANGE = 1e-7

    if pareto_Y.shape[0] == 0:
        # print('case a')
        if max_ref_point is None:
            raise BotorchError("Empty pareto set and no max ref point provided")
        if scale_max_ref_point:
            return max_ref_point - scale * max_ref_point.abs()
        return max_ref_point
    if max_ref_point is not None:
        # print('case b')
        better_than_ref = (pareto_Y > max_ref_point).all(dim=-1)
    else:
        # print('case c')
        better_than_ref = torch.full(
            pareto_Y.shape[:1], 1, dtype=bool, device=pareto_Y.device
        )
        # print(better_than_ref)
    if max_ref_point is not None and better_than_ref.any():
        # print('case d')
        Y_range = pareto_Y[better_than_ref].max(dim=0).values - max_ref_point
        if scale_max_ref_point:
            return max_ref_point - scale * Y_range
        return max_ref_point
    elif pareto_Y.shape[0] == 1:
        # print('case e')
        # no points better than max_ref_point and only a single observation
        # subtract MIN_Y_RANGE to handle the case that pareto_Y is a singleton
        # with objective value of 0.
        return (pareto_Y - scale * pareto_Y.abs().clamp_min(MIN_Y_RANGE)).view(-1)
    # no points better than max_ref_point and multiple observations
    # make sure that each dimension of the nadir point is no greater than
    # the max_ref_point
    nadir = pareto_Y.min(dim=0).values
    # print('nadir', nadir)
    if max_ref_point is not None:
        # print('case e')
        nadir = torch.min(nadir, max_ref_point)
    ideal = pareto_Y.max(dim=0).values
    # print('ideal', ideal)
    # handle case where all values for one objective are the same
    Y_range = (ideal - nadir).clamp_min(MIN_Y_RANGE)
    # print('Y_range', Y_range)
    return nadir - scale * Y_range
