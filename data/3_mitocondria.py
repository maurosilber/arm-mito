import numpy as np
import xarray
from tqdm import tqdm

from mito import ARM_Cito, Mitochondria
from n_mito.loop_simulator import LoopSimulator

t = np.arange(0, 86_400, 30)

loop_sim = LoopSimulator(
    ARM_Cito,
    Mitochondria(
        volume_cell=ARM_Cito.volume,
        CytoC_C=ARM_Cito.CytoC_C,
        Smac_C=ARM_Cito.Smac_C,
        Bax_A=ARM_Cito.Bax_A,
    ),
)


def solve(
    *,
    N_mitocondria: int = 1,
    cell_volume: float = 1,
    global_mitocondria_fraction: float = 0.07,
    L_concentration: float,
    IntrinsicStimuli_concentration: float,
):
    df = loop_sim.solve(
        save_at=t,
        main_values={
            ARM_Cito.volume: cell_volume,
            ARM_Cito.L: L_concentration * cell_volume,
            ARM_Cito.IntrinsicStimuli: IntrinsicStimuli_concentration * cell_volume,
        },
        loop_values={
            Mitochondria.volume: cell_volume
            * global_mitocondria_fraction
            / N_mitocondria,
            Mitochondria.Bax4: np.zeros(N_mitocondria),  # just to create N submodules
        },
    )

    cols = {
        ARM_Cito.caspam.sCas3.monomer.variable: "sCas3",
        ARM_Cito.caspam.sCas8.monomer.variable: "sCas8",
        ARM_Cito.caspam.sCas9.monomer.variable: "sCas9",
        ARM_Cito.Apop.variable: "Apop",
        ARM_Cito.C3_A.variable: "C3",
        ARM_Cito.C8_A.variable: "C8",
    }
    df.index.rename("time", inplace=True)
    df = df[cols.keys()].rename(columns=cols)
    return df


N = np.unique(np.geomspace(1, 100, 50, dtype=int))
stimuli = {
    "extrinsic": {"L_concentration": 1000, "IntrinsicStimuli_concentration": 0},
    "intrinsic": {"L_concentration": 0, "IntrinsicStimuli_concentration": 100},
}
df = xarray.combine_by_coords(
    [
        solve(N_mitocondria=n, **stimulus_values)
        .to_xarray()
        .expand_dims(
            {
                "stimulus": [stimulus_name],
                "N_mitocondrias": [n],
            }
        )
        for n in tqdm(N)
        for stimulus_name, stimulus_values in stimuli.items()
    ]
)
df.to_zarr("3_mitochondria.zarr")
