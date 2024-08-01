import numpy as np
import pandas as pd
import pysb
import simbio

TOL = 1e-6


def run_pysb(model: pysb.Model, times: np.ndarray):
    import sys
    import types

    # Mock distutils module required by PySB simulator
    sys.modules["distutils"] = types.ModuleType("distutils")
    import pysb.simulator

    sim = pysb.simulator.ScipyOdeSimulator(
        model,
        compiler="python",
        integrator="lsoda",
        integrator_options={"atol": TOL, "rtol": TOL},
    )
    result = sim.run(times)
    return pd.DataFrame(
        data=result.species,
        columns=model.species,
        index=pd.Index(times, name="time"),
    )


def run_simbio(sim: simbio.Simulator, times: np.ndarray):
    from poincare import solvers

    return sim.solve(
        save_at=times,
        solver=solvers.LSODA(atol=TOL, rtol=TOL),
    )


def run_models(
    *,
    simbio_model: simbio.Compartment,
    pysb_model: pysb.Model,
    times: np.ndarray,
    mapping: dict[simbio.Species, str],
):
    df_pysb = run_pysb(pysb_model, times)
    df_simbio = run_simbio(simbio_model, times)

    df_pysb = df_pysb.rename(columns=str).rename(
        columns={str(v): str(k) for k, v in mapping.items()},
    )

    names = [str(k) for k in mapping]
    return df_simbio[names], df_pysb[names]
