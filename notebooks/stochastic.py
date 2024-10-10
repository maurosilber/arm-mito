import xarray
from simbio_rebop.converter import create_rebop, to_rebop_loopy

from mito import ARM


def create(*, N: int, volume: float, L, Intrinsic):
    reactions, y = to_rebop_loopy(
        ARM,
        ARM.mitocondria,
        N=N,
        values={
            ARM.L_concentration: L,
            ARM.IntrinsicStimuli_concentration: Intrinsic,
            ARM.volume: volume,
            ARM.mitocondria_volume_fraction: 0.07 / N,
        },
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
