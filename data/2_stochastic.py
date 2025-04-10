import numpy as np
import xarray
from simrebop import RebopSimulator

from mito import ARM

T = 86_400
dt = 30

stimuli = {
    "extrinsic": {ARM.L_concentration: 1000, ARM.IntrinsicStimuli_concentration: 0},
    "intrinsic": {ARM.L_concentration: 0, ARM.IntrinsicStimuli_concentration: 100},
}

volumes = np.geomspace(0.01, 1, 30)
seeds = range(100)

if __name__ == "__main__":
    from tqdm import tqdm

    sim = RebopSimulator(ARM)
    df = [
        sim.solve_many(
            values={
                **stimulus_values,
                ARM.volume: v,
            },
            upto_t=T,
            n_points=T // dt,
            var_names=[
                ARM.cytoplasm.Apop,
                ARM.cytoplasm.C3_A,
                ARM.cytoplasm.C8_A,
                ARM.cytoplasm.caspam.sCas3.monomer,
                ARM.cytoplasm.caspam.sCas8.monomer,
                ARM.cytoplasm.caspam.sCas9.monomer,
            ],
            seeds=seeds,
        ).expand_dims(
            {
                "stimulus": [stimulus_name],
                "volume": [v],
            }
        )
        for stimulus_name, stimulus_values in stimuli.items()
        for v in tqdm(volumes)
    ]

    df = xarray.combine_by_coords(df).rename_vars(
        {
            "cytoplasm.C3_A": "C3_A",
            "cytoplasm.caspam.sCas8.monomer": "sCas8",
            "cytoplasm.caspam.sCas9.monomer": "sCas9",
            "cytoplasm.caspam.sCas3.monomer": "sCas3",
            "cytoplasm.C8_A": "C8_A",
            "cytoplasm.Apop": "Apop",
        }
    )
    df.to_zarr("2_stochastic.zarr", mode="w")
