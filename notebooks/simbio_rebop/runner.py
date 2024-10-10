from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from converter import create_rebop, to_rebop_loopy

from mito import ARM

t_max = 15 * 3_600  # seconds
save = list(map(str, [ARM.cytoplasm.C3_A, ARM.cytoplasm.C8_A, ARM.cytoplasm.Apop]))


@dataclass
class Model:
    path: Path
    N: int
    volume: float

    def __post_init__(self):
        self.path.mkdir(parents=True, exist_ok=True)

    @cached_property
    def rebop(self):
        reactions, self.y = to_rebop_loopy(
            ARM,
            ARM.mitocondria,
            N=self.N,
            values={
                ARM.L_concentration: 1000,
                ARM.volume: self.volume,
                ARM.mitocondria_volume_fraction: 0.07 / self.N,
            },
        )
        self.save = save
        return create_rebop(reactions)

    def run(self, seed: int):
        p = self.path / f"{seed}.parquet"
        if p.exists():
            return
        df = self.rebop.run(
            self.y,
            tmax=t_max,
            nb_steps=t_max,
            seed=seed,
            save=self.save,
            sparse=True,
        ).to_dataframe()
        df[sorted(df.columns)].to_parquet(p, compression="zstd")


root = Path("results")


def batch(x: tuple[int, str, int]):
    seed, volume, N = x
    runner = Model(
        root / f"ARM{N}" / volume,
        N,
        float(volume),
    )
    runner.run(seed)


if __name__ == "__main__":
    import itertools
    import math
    import os

    import numpy as np
    import typer
    from tqdm.contrib.concurrent import process_map

    def main(*, seeds: tuple[int, int], mito: list[int], workers: int | None = None):
        if workers is None:
            workers = os.cpu_count()

        seed_range = range(*seeds)
        volumes = np.geomspace(0.1, 1, 6)[-1:]
        volumes = [f"{x:.3f}" for x in volumes]

        factors = (seed_range, volumes, mito)
        process_map(
            batch,
            itertools.product(*factors),
            total=math.prod(map(len, factors)),
            miniters=1,
            max_workers=workers,
        )

    typer.run(main)
