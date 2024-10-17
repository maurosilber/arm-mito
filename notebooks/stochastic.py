import xarray
from numpy.typing import ArrayLike
from simbio import Constant, Parameter, Species
from simbio_rebop.converter import create_rebop, to_rebop_loopy

from mito import ARM


def create(
    *,
    volume: float,
    L: float,
    Intrinsic: float,
    loop_values: dict[Species | Parameter | Constant, ArrayLike],
):
    reactions, y = to_rebop_loopy(
        ARM,
        ARM.mitocondria,
        values_main={
            ARM.L_concentration: L,
            ARM.IntrinsicStimuli_concentration: Intrinsic,
            ARM.volume: volume,
        },
        values_loop=loop_values,
    )
    runner = create_rebop(reactions)
    return runner, y


def run(
    seed: int,
    *,
    t_max: int,
    steps: int,
    create,
    save: list | None = [ARM.cytoplasm.C3_A, ARM.cytoplasm.C8_A, ARM.cytoplasm.Apop],
) -> xarray.Dataset:
    runner, y = create()
    if save is not None:
        save = list(map(str, save))
    return runner.run(y, tmax=t_max, nb_steps=steps, seed=seed, save=save, sparse=True)
