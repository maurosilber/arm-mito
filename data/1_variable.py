import numpy as np
import xarray
from simbio import Simulator
from symbolite.core import evaluate
from tqdm import tqdm

from mito import ARM

sim = Simulator(
    ARM,
    transform={
        "Apop": ARM.cytoplasm.Apop,
        "C3": ARM.cytoplasm.C3_A,
        "C8": ARM.cytoplasm.C8_A,
        "sCas3": ARM.cytoplasm.caspam.sCas3.monomer,
        "sCas8": ARM.cytoplasm.caspam.sCas8.monomer,
        "sCas9": ARM.cytoplasm.caspam.sCas9.monomer,
    },
    backend="numba",
)


def get_from_name(model, name):
    for k in name.split("."):
        model = getattr(model, k)
    return model


initials = ARM.variables.map(evaluate).pipe(lambda x: x[x > 0])
initials.index = initials.index.map(lambda x: get_from_name(ARM, x))
initials = initials.to_dict()
initials[ARM.cytoplasm.caspam_0] = evaluate(ARM.cytoplasm.caspam_0)

stimuli = {
    "extrinsic": {ARM.L_concentration: 1000, ARM.IntrinsicStimuli_concentration: 0},
    "intrinsic": {ARM.L_concentration: 0, ARM.IntrinsicStimuli_concentration: 100},
}

t = np.arange(0, 86_400, 30)
perturbation = np.logspace(-5, 5, 100)
df = xarray.combine_by_coords(
    [
        sim.solve(save_at=t, values=stimulus_values | {k: v * pert})
        .to_xarray()
        .expand_dims(
            {
                "stimulus": [stimulus_name],
                "perturbation": [pert],
                "variable": [str(k)],
            }
        )
        for pert in tqdm(perturbation)
        for k, v in tqdm(initials.items(), total=len(initials), leave=False)
        for stimulus_name, stimulus_values in stimuli.items()
    ]
)

df.to_zarr("1_variable.zarr", mode="w")
