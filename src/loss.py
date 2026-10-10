import torch
from torch import nn
import torch.nn.functional as F


def smoothl1(inputs, targets, alpha=None):
    loss_f = nn.SmoothL1Loss(reduce=False)
    weight = torch.ones(inputs.shape, dtype=torch.float).to(inputs.device)
    if alpha is not None:
        weight[targets != 0] = float(alpha)
    loss = torch.mean(loss_f(inputs, targets) * weight)
    return loss


def bce(inputs, targets):
    # inputs are raw logits — BCEWithLogitsLoss is AMP-safe.
    return F.binary_cross_entropy_with_logits(inputs, targets)


def FL(inputs, targets, alpha, gamma):
    """Focal loss for pitch detection.

    inputs  : raw logits  (B, T, N_CLASS)  — model output before sigmoid
    targets : binary labels (same shape)
    alpha   : positive-class weight
    gamma   : focusing exponent (0 = standard weighted BCE)

    Using binary_cross_entropy_with_logits instead of binary_cross_entropy
    because:
      1. AMP-safe — numerically stable in fp16.
      2. Fuses sigmoid + log into one numerically stable op (log-sum-exp trick),
         avoiding log(0) that fp16 sigmoid can produce near the extremes.
    """
    weight = torch.ones(inputs.shape, dtype=torch.float, device=inputs.device)
    weight[targets == 1] = float(alpha)

    # Per-element unweighted and weighted BCE on logits
    loss   = F.binary_cross_entropy_with_logits(inputs, targets,
                                                 reduction='none')
    loss_w = F.binary_cross_entropy_with_logits(inputs, targets,
                                                 weight=weight,
                                                 reduction='none')

    # Focal modulation: pt = exp(-BCE) ≡ the model's confidence on each label
    pt           = torch.exp(-loss)
    weight_gamma = (1 - pt) ** gamma

    return torch.mean(weight_gamma * loss_w)
