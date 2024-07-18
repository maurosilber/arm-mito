import numpy as np
import pandas as pd
import pysb
import pysb.util
from earm import albeck_modules
from simbio import Simulator

from . import common


def test_albeck_receptor():
    from ..albeck import Receptor

    model = pysb.Model()
    albeck_modules.ligand_to_c8_monomers()
    Bid = pysb.Monomer("Bid", ["bf", "state"], {"state": ["U", "T"]})
    pysb.Initial(Bid(state="U", bf=None), pysb.Parameter("Bid_0", 1e5))
    albeck_modules.rec_to_bid()

    mapping = {
        L(bf=None): Receptor.L,
        R(bf=None): Receptor.R,
        flip(bf=None): Receptor.flip,
        C8(bf=None, state="pro"): Receptor.C8_pro,
        BAR(bf=None): Receptor.BAR,
        L(bf=1) % R(bf=1): Receptor.r_L_R.AB,
        DISC(bf=None): Receptor.DISC,
        Bid(bf=1, state="U") % C8(bf=1, state="A"): Receptor.r_C8_Bid.ES,
        DISC(bf=1) % flip(bf=1): Receptor.r_DISC_flip.AB,
        C8(bf=None, state="A"): Receptor.C8_A,
        BAR(bf=1) % C8(bf=1, state="A"): Receptor.r_BAR_C8.AB,
        Bid(bf=None, state="U"): Receptor.Bid_U,
        Bid(bf=None, state="T"): Receptor.Bid_T,
        C8(bf=1, state="pro") % DISC(bf=1): Receptor.r_DISC_C8.ES,
    }

    times = np.linspace(0, 10_000, 1000)
    sim = Simulator(Receptor)
    df_simbio = common.run_simbio(sim, times)
    df_pysb = common.run_pysb(model, times)
    df_pysb = df_pysb.rename(columns=str).rename(
        columns={str(k): str(v) for k, v in mapping.items()},
    )
    pd.testing.assert_frame_equal(
        df_simbio,
        df_pysb,
        atol=1e-5,
        rtol=1e-5,
        check_like=True,  # ignore order
    )


def test_albeck11b():
    from earm.albeck_11b import model

    from ..albeck import Albeck11b as Model

    pysb.util.alias_model_components()

    mapping = {
        L(bf=None): Model.receptor.L,
        R(bf=None): Model.receptor.R,
        flip(bf=None): Model.receptor.flip,
        C8(bf=None, state="pro"): Model.receptor.C8_pro,
        BAR(bf=None): Model.receptor.BAR,
        L(bf=1) % R(bf=1): Model.receptor.r_L_R.AB,
        DISC(bf=None): Model.receptor.DISC,
        Bid(bf=1, state="U") % C8(bf=1, state="A"): Model.receptor.r_C8_Bid.ES,
        DISC(bf=1) % flip(bf=1): Model.receptor.r_DISC_flip.AB,
        C8(bf=None, state="A"): Model.receptor.C8_A,
        BAR(bf=1) % C8(bf=1, state="A"): Model.receptor.r_BAR_C8.AB,
        Bid(bf=None, state="U"): Model.receptor.Bid_U,
        Bid(bf=None, state="T"): Model.receptor.Bid_T,
        C8(bf=1, state="pro") % DISC(bf=1): Model.receptor.r_DISC_C8.ES,
        Bax(bf=None, s1=None, s2=None, state="C"): Model.mitocondria.Bax_C,
        Bax(bf=None, s1=None, s2=None, state="A"): Model.mitocondria.Bax_A,
        Bcl2(bf=None): Model.mitocondria.Bcl2,
        CytoC(bf=None, state="M"): Model.mitocondria.CytoC_M,
        CytoC(bf=None, state="A"): Model.mitocondria.CytoC_A,
        CytoC(bf=None, state="C"): Model.mitocondria.CytoC_C,
        Smac(bf=None, state="M"): Model.mitocondria.Smac_M,
        Smac(bf=None, state="A"): Model.mitocondria.Smac_A,
        Smac(bf=None, state="C"): Model.mitocondria.Smac_C,
        Apaf(bf=None, state="I"): Model.intrinsic.Apaf_I,
        Apaf(bf=None, state="A"): Model.intrinsic.Apaf_A,
        Apop(bf=None): Model.intrinsic.Apop,
        C3(bf=None, state="pro"): Model.intrinsic.C3_pro,
        C3(bf=None, state="A"): Model.intrinsic.C3_A,
        C3(bf=None, state="ub"): Model.intrinsic.C3_ub,
        C6(bf=None, state="pro"): Model.intrinsic.C6_pro,
        C6(bf=None, state="A"): Model.intrinsic.C6_A,
        C9(bf=None): Model.intrinsic.C9,
        PARP(bf=None, state="U"): Model.intrinsic.PARP_U,
        PARP(bf=None, state="C"): Model.intrinsic.PARP_C,
        XIAP(bf=None): Model.intrinsic.XIAP,
        Apaf(bf=1, state="I")
        % CytoC(bf=1, state="A"): Model.intrinsic.r_CytoC_Apaf.ES,
        Apop(bf=1) % C3(bf=1, state="pro"): Model.intrinsic.r_Apop_C3.ES,
        Apop(bf=1) % XIAP(bf=1): Model.intrinsic.r_Apop_XIAP.AB,
        Bax(bf=1, s1=None, s2=None, state="A")
        % Bcl2(bf=1): Model.r_Bax_Bcl2.AB,
        Bax(bf=1, s1=None, s2=None, state="A")
        % CytoC(bf=1, state="M"): Model.r_CytoC_pore.ES,
        Bax(bf=1, s1=None, s2=None, state="A")
        % Smac(bf=1, state="M"): Model.r_Smac_pore.ES,
        Bax(bf=1, s1=None, s2=None, state="C")
        % Bid(bf=1, state="T"): Model.r_Bid_Bax.ES,
        C3(bf=1, state="A") % C6(bf=1, state="pro"): Model.intrinsic.r_C3_C6.ES,
        C3(bf=1, state="A") % PARP(bf=1, state="U"): Model.intrinsic.r_C3_PARP.ES,
        C3(bf=1, state="A") % XIAP(bf=1): Model.intrinsic.r_XIAP_C3.ES,
        C3(bf=1, state="pro") % C8(bf=1, state="A"): Model.intrinsic.r_C8_C3.ES,
        C6(bf=1, state="A") % C8(bf=1, state="pro"): Model.intrinsic.r_C6_C8.ES,
        Smac(bf=1, state="A") % XIAP(bf=1): Model.intrinsic.r_Smac_XIAP.AB,
    }

    times = np.linspace(0, 10_000, 1000)
    sim = Simulator(Model)
    df_simbio = common.run_simbio(sim, times)
    df_pysb = common.run_pysb(model, times)
    df_pysb = df_pysb.rename(columns=str).rename(
        columns={str(k): str(v) for k, v in mapping.items()},
    )

    pd.testing.assert_frame_equal(
        df_simbio,
        df_pysb,
        atol=1e-4,
        rtol=1e-4,
        check_like=True,  # ignore order
    )
