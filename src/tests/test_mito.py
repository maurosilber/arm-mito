import numpy as np
import pandas as pd
from poincare import solvers
from pytest import mark
from simbio import Simulator

from ..corbat import IntrinsicARM as ARM
from ..mito import ARM as ARM_Mito


@mark.parametrize("stimuli", ["extrinsic", "intrinsic"])
def test_mito(stimuli):
    mapping = {
        ARM.IntrinsicStimuli: ARM_Mito.IntrinsicStimuli,
        ARM.arm.L: ARM_Mito.L,
        ARM.arm.L: ARM_Mito.L,
        ARM.arm.R: ARM_Mito.R,
        ARM.arm.DISC: ARM_Mito.DISC,
        ARM.arm.flip: ARM_Mito.flip,
        ARM.arm.C8_pro: ARM_Mito.C8_pro,
        ARM.arm.C8_A: ARM_Mito.C8_A,
        ARM.arm.BAR: ARM_Mito.BAR,
        ARM.arm.Bid_U: ARM_Mito.Bid_U,
        ARM.arm.Bid_T: ARM_Mito.Bid_T,
        ARM.arm.Bid_M: ARM_Mito.Bid_M,
        ARM.arm.Bax_C: ARM_Mito.Bax_C,
        ARM.arm.Bax_M: ARM_Mito.mitocondria.Bax,
        ARM.arm.Bax_A: ARM_Mito.Bax_A,
        ARM.arm.Bcl2: ARM_Mito.mitocondria.Bcl2,
        ARM.arm.CytoC_M: ARM_Mito.mitocondria.CytoC_M,
        ARM.arm.CytoC_C: ARM_Mito.CytoC_C,
        ARM.arm.CytoC_A: ARM_Mito.CytoC_A,
        ARM.arm.Smac_M: ARM_Mito.mitocondria.Smac_M,
        ARM.arm.Smac_C: ARM_Mito.Smac_C,
        ARM.arm.Smac_A: ARM_Mito.Smac_A,
        ARM.arm.Apaf_I: ARM_Mito.Apaf_I,
        ARM.arm.Apaf_A: ARM_Mito.Apaf_A,
        ARM.arm.Apop: ARM_Mito.Apop,
        ARM.arm.C3_pro: ARM_Mito.C3_pro,
        ARM.arm.C3_A: ARM_Mito.C3_A,
        ARM.arm.C3_ub: ARM_Mito.C3_ub,
        ARM.arm.C6_pro: ARM_Mito.C6_pro,
        ARM.arm.C6_A: ARM_Mito.C6_A,
        ARM.arm.PARP_U: ARM_Mito.PARP_U,
        ARM.arm.PARP_C: ARM_Mito.PARP_C,
        ARM.arm.XIAP: ARM_Mito.XIAP,
        ARM.arm.Bax_M2: ARM_Mito.mitocondria.Bax2,
        ARM.arm.Bax_M4: ARM_Mito.mitocondria.Bax4,
        ARM.arm.Mito_I: ARM_Mito.mitocondria.Mito_I,
        ARM.arm.Mito_A: ARM_Mito.mitocondria.Mito_A,
        ARM.arm.Bcl2c: ARM_Mito.Bcl2c,
        ARM.arm.r_Bax2_Bcl2.AB: ARM_Mito.mitocondria.r_Bax2_Bcl2.AB,
        ARM.arm.r_Bax4_Bcl2.AB: ARM_Mito.mitocondria.r_Bax4_Bcl2.AB,
        ARM.arm.r_Bax4_Mito.AB: ARM_Mito.mitocondria.r_Bax4_Mito.AB,
        ARM.arm.r_Bax_Bcl2.AB: ARM_Mito.mitocondria.r_Bax_Bcl2.AB,
        ARM.arm.r_CytoC_pore.ES: ARM_Mito.mitocondria.r_CytoC_pore.ES,
        ARM.arm.r_Smac_pore.ES: ARM_Mito.mitocondria.r_Smac_pore.ES,
        ARM.arm.r_ApafA_C3pro.ES: ARM_Mito.r_ApafA_C3pro.ES,
        ARM.arm.r_ApafA_XIAP.AB: ARM_Mito.r_ApafA_XIAP.AB,
        ARM.arm.r_Apop_C3.ES: ARM_Mito.r_Apop_C3.ES,
        ARM.arm.r_BAR_C8.AB: ARM_Mito.r_BAR_C8.AB,
        ARM.arm.r_Bcl2_Bid.AB: ARM_Mito.r_Bcl2_Bid.AB,
        ARM.arm.r_Bid_Bax.ES: ARM_Mito.r_Bid_Bax.ES,
        ARM.arm.r_C3A_ApafA.ES: ARM_Mito.r_C3A_ApafA.ES,
        ARM.arm.r_C3_C6.ES: ARM_Mito.r_C3_C6.ES,
        ARM.arm.r_C3_PARP.ES: ARM_Mito.r_C3_PARP.ES,
        ARM.arm.r_C6_C8.ES: ARM_Mito.r_C6_C8.ES,
        ARM.arm.r_C8_Bid.ES: ARM_Mito.r_C8_Bid.ES,
        ARM.arm.r_C8_C3.ES: ARM_Mito.r_C8_C3.ES,
        ARM.arm.r_CytoC_Apaf.ES: ARM_Mito.r_CytoC_Apaf.ES,
        ARM.arm.r_DISC_C8.ES: ARM_Mito.r_DISC_C8.ES,
        ARM.arm.r_DISC_flip.AB: ARM_Mito.r_DISC_flip.AB,
        ARM.arm.r_L_R.AB: ARM_Mito.r_L_R.AB,
        ARM.arm.r_Smac_XIAP.AB: ARM_Mito.r_Smac_XIAP.AB,
        ARM.arm.r_XIAP_C3.ES: ARM_Mito.r_XIAP_C3.ES,
        ARM.r_intrinsic.ES: ARM_Mito.r_intrinsic.ES,
    }

    t = np.linspace(0, 30_000, 1_000)
    solver = solvers.LSODA(atol=1e-6, rtol=1e-6)
    if stimuli == "extrinsic":
        extrinsic, intrinsic = 1000, 0
    elif stimuli == "intrinsic":
        extrinsic, intrinsic = 0, 200
    else:
        raise ValueError(stimuli)
    df_ARM = Simulator(ARM).solve(
        solver=solver,
        save_at=t,
        values={ARM.arm.L: extrinsic, ARM.IntrinsicStimuli: intrinsic},
    )
    df_Mito = Simulator(ARM_Mito).solve(
        solver=solver,
        save_at=t,
        values={ARM_Mito.L: extrinsic, ARM_Mito.IntrinsicStimuli: intrinsic},
    )

    df_ARM.rename(inplace=True, columns={str(k): str(v) for k, v in mapping.items()})
    pd.testing.assert_frame_equal(
        df_ARM,
        df_Mito,
        atol=1e-4,
        rtol=1e-4,
        check_like=True,  # ignore order
    )
