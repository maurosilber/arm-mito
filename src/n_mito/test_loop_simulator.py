import numpy as np
import pandas as pd
from poincare.solvers import LSODA
from pytest import mark
from simbio import Compartment, Simulator

from ..mito import ARM_Cito, Mitochondria
from .loop_simulator import LoopSimulator


class Manual(Compartment):
    cytoplasm = ARM_Cito()
    mitochondria_0 = Mitochondria(
        CytoC_C=cytoplasm.CytoC_C,
        Smac_C=cytoplasm.Smac_C,
        Bax_A=cytoplasm.Bax_A,
    )
    mitochondria_1 = Mitochondria(
        CytoC_C=cytoplasm.CytoC_C,
        Smac_C=cytoplasm.Smac_C,
        Bax_A=cytoplasm.Bax_A,
    )


@mark.parametrize(
    ["values", "main_values", "loop_values"],
    [
        (
            {Manual.cytoplasm.L: 1_000},
            {ARM_Cito.L: 1_000},
            {Mitochondria.Mito_A: [0, 0]},
        ),
        (
            {
                Manual.cytoplasm.L: 1_000,
                Manual.mitochondria_0.Bcl2: 1e4,
                Manual.mitochondria_1.Bcl2: 3e4,
            },
            {ARM_Cito.L: 1_000},
            {Mitochondria.Bcl2: [1e4, 3e4]},
        ),
        (
            {
                Manual.cytoplasm.L: 1_000,
                Manual.mitochondria_0.volume: 0.05,
                Manual.mitochondria_1.volume: 0.10,
            },
            {ARM_Cito.L: 1_000},
            {Mitochondria.volume: [0.05, 0.10]},
        ),
    ],
)
def test_with_ARM(values, main_values, loop_values):
    sim = Simulator(Manual)
    sim_loop = LoopSimulator(
        ARM_Cito,
        Mitochondria(
            CytoC_C=ARM_Cito.CytoC_C,
            Smac_C=ARM_Cito.Smac_C,
            Bax_A=ARM_Cito.Bax_A,
        ),
    )

    t = np.linspace(0, 30_000, 1_000)
    solver = LSODA(rtol=1e-6, atol=1e-6)
    df = sim.solve(
        save_at=t,
        solver=solver,
        values=values,
    )
    df_loop = sim_loop.solve(
        save_at=t,
        solver=solver,
        main_values=main_values,
        loop_values=loop_values,
        loop_output="index_as_suffix",
    )

    # Match names
    def renamer(x: str, /):
        if x.startswith("cytoplasm."):
            return x.removeprefix("cytoplasm.")

        x = x.removeprefix("mitochondria_")
        n, _, name = x.partition(".")
        return f"{name}_{n}"

    df = df.rename(columns=renamer, copy=False)
    df_loop = df_loop.rename(columns=str, copy=False)
    df_loop.index.rename("time", inplace=True)

    pd.testing.assert_frame_equal(
        df,
        df_loop,
        rtol=1e-4,
        atol=1e-4,
    )
