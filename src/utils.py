import os
import sys
import random
from functools import reduce
import numpy as np
import torch


from torch.nn.modules.module import _addindent


def set_seed(seed: int, deterministic: bool = True):
    """
    Fixes every RNG axis this codebase touches: Python `random`, NumPy,
    PyTorch CPU + all CUDA devices. This alone fixes model init weights
    and data order — independent of the cudnn flags below.

    deterministic=True forces cuDNN to a fixed (not necessarily fastest)
    conv algorithm, for bit-exact repro across runs. deterministic=False
    lets cuDNN autotune the fastest kernel per shape instead — this does
    NOT change achievable accuracy, only run time; results differ only in
    the last few floating-point digits since kernel summation order isn't
    fixed. There is no accuracy trade-off in either direction.

    cudnn.deterministic alone does not cover every op: some ATen kernels
    have nondeterministic CUDA implementations that only
    use_deterministic_algorithms() rejects, and cuBLAS GEMMs need the
    workspace config set to be reproducible. Both are handled below so
    deterministic=True delivers actual bit-exactness.
    """
    # Must precede the torch.cuda.* calls below: cuBLAS reads this only when
    # its handle is created, and manual_seed() is enough to initialize CUDA.
    # run_experiments.py also sets it in the subprocess env, which is the
    # airtight path; this covers direct `python train_single.py` invocation.
    if deterministic:
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    torch.backends.cudnn.deterministic = deterministic
    torch.backends.cudnn.benchmark = not deterministic

    if deterministic:
        # warn_only: an op with no deterministic kernel warns instead of
        # raising, so a long sweep never dies mid-run over one such op.
        torch.use_deterministic_algorithms(True, warn_only=True)


def seed_worker(worker_id):
    """Pass as DataLoader(..., worker_init_fn=seed_worker) so each of the
    num_workers subprocesses gets a deterministic (but distinct) seed,
    derived from torch's initial_seed rather than colliding across workers."""
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def cycle(iterable):
    while True:
        for item in iterable:
            yield item


def summary(model, file=sys.stdout):
    def repr(model):
        # We treat the extra repr like the sub-module, one item per line
        extra_lines = []
        extra_repr = model.extra_repr()
        # empty string will be split into list ['']
        if extra_repr:
            extra_lines = extra_repr.split('\n')
        child_lines = []
        total_params = 0
        for key, module in model._modules.items():
            mod_str, num_params = repr(module)
            mod_str = _addindent(mod_str, 2)
            child_lines.append('(' + key + '): ' + mod_str)
            total_params += num_params
        lines = extra_lines + child_lines

        for name, p in model._parameters.items():
            if hasattr(p, 'shape'):
                total_params += reduce(lambda x, y: x * y, p.shape)

        main_str = model._get_name() + '('
        if lines:
            # simple one-liner info, which most builtin Modules will use
            if len(extra_lines) == 1 and not child_lines:
                main_str += extra_lines[0]
            else:
                main_str += '\n  ' + '\n  '.join(lines) + '\n'

        main_str += ')'
        if file is sys.stdout:
            main_str += ', \033[92m{:,}\033[0m params'.format(total_params)
        else:
            main_str += ', {:,} params'.format(total_params)
        return main_str, total_params

    string, count = repr(model)
    if file is not None:
        if isinstance(file, str):
            file = open(file, 'w')
        print(string, file=file)
        file.flush()

    return count


def to_local_average_cents(salience, center=None, thred=0.0):
    """
    find the weighted average cents near the argmax bin
    """

    if not hasattr(to_local_average_cents, 'cents_mapping'):
        # the bin number-to-cents mapping
        to_local_average_cents.cents_mapping = (
                np.linspace(0, 7180, 360) + 1997.3794084376191)

    if salience.ndim == 1:
        if center is None:
            center = int(np.argmax(salience))
        start = max(0, center - 4)
        end = min(len(salience), center + 5)
        salience = salience[start:end]
        product_sum = np.sum(
            salience * to_local_average_cents.cents_mapping[start:end])
        weight_sum = np.sum(salience)
        return product_sum / weight_sum if np.max(salience) > thred else 0
    if salience.ndim == 2:
        return np.array([to_local_average_cents(salience[i, :], None, thred) for i in
                         range(salience.shape[0])])

    raise Exception("label should be either 1d or 2d ndarray")
