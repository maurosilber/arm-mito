import numpy as np
import pandas as pd
import pysb
import pysb.util
from simbio import Simulator

from . import common


def test_albeck_as_matlab():
    from caspase_model.models import albeck_as_matlab

    from ..corbat import AlbeckAsMatlab as Model

    model = albeck_as_matlab()
    pysb.util.alias_model_components()

    mapping = {
        L(bf=None): Model.L,
        R(bf=None): Model.R,
        flip(bf=None): Model.flip,
        C8(bf=None, state="pro"): Model.C8_pro,
        BAR(bf=None): Model.BAR,
        L(bf=1) % R(bf=1): Model.r_L_R.AB,
        DISC(bf=None): Model.DISC,
        Bid(bf=1, state="U") % C8(bf=1, state="A"): Model.r_C8_Bid.ES,
        DISC(bf=1) % flip(bf=1): Model.r_DISC_flip.AB,
        C8(bf=None, state="A"): Model.C8_A,
        BAR(bf=1) % C8(bf=1, state="A"): Model.r_BAR_C8.AB,
        Bid(bf=None, state="U"): Model.Bid_U,
        Bid(bf=None, state="T"): Model.Bid_T,
        C8(bf=1, state="pro") % DISC(bf=1): Model.r_DISC_C8.ES,
        Bax(bf=None, s1=None, s2=None, state="C"): Model.mitocondria.Bax_C,
        Bax(bf=None, s1=None, s2=None, state="A"): Model.mitocondria.Bax_A,
        Bcl2(bf=None): Model.mitocondria.Bcl2,
        CytoC(bf=None, state="M"): Model.mitocondria.CytoC_M,
        CytoC(bf=None, state="A"): Model.mitocondria.CytoC_A,
        CytoC(bf=None, state="C"): Model.mitocondria.CytoC_C,
        Smac(bf=None, state="M"): Model.mitocondria.Smac_M,
        Smac(bf=None, state="A"): Model.mitocondria.Smac_A,
        Smac(bf=None, state="C"): Model.mitocondria.Smac_C,
        Apaf(bf=None, state="I"): Model.Apaf_I,
        Apaf(bf=None, state="A"): Model.Apaf_A,
        Apop(bf=None): Model.Apop,
        C3(bf=None, state="pro"): Model.C3_pro,
        C3(bf=None, state="A"): Model.C3_A,
        C3(bf=None, state="ub"): Model.C3_ub,
        C6(bf=None, state="pro"): Model.C6_pro,
        C6(bf=None, state="A"): Model.C6_A,
        C9(bf=None): Model.C9,
        PARP(bf=None, state="U"): Model.PARP_U,
        PARP(bf=None, state="C"): Model.PARP_C,
        XIAP(bf=None): Model.XIAP,
        Apaf(bf=1, state="I") % CytoC(bf=1, state="A"): Model.r_CytoC_Apaf.ES,
        Apop(bf=1) % C3(bf=1, state="pro"): Model.r_Apop_C3.ES,
        Apop(bf=1) % XIAP(bf=1): Model.r_Apop_XIAP.AB,
        Bax(bf=1, s1=None, s2=None, state="M") % Bcl2(bf=1): Model.r_Bax_Bcl2.AB,
        Bax(bf=1, s1=None, s2=None, state="C")
        % Bid(bf=1, state="T"): Model.r_Bid_Bax.ES,
        C3(bf=1, state="A") % C6(bf=1, state="pro"): Model.r_C3_C6.ES,
        C3(bf=1, state="A") % PARP(bf=1, state="U"): Model.r_C3_PARP.ES,
        C3(bf=1, state="A") % XIAP(bf=1): Model.r_XIAP_C3.ES,
        C3(bf=1, state="pro") % C8(bf=1, state="A"): Model.r_C8_C3.ES,
        C6(bf=1, state="A") % C8(bf=1, state="pro"): Model.r_C6_C8.ES,
        Smac(bf=1, state="A") % XIAP(bf=1): Model.r_Smac_XIAP.AB,
        # NEw
        Bax(bf=1, s1=2, s2=3, state="M")
        % Bax(bf=None, s1=4, s2=2, state="M")
        % Bax(bf=None, s1=5, s2=4, state="M")
        % Bax(bf=None, s1=3, s2=5, state="M")
        % Bcl2(bf=1): Model.r_Bax4_Bcl2.AB,
        Bax(bf=1, s1=None, s2=2, state="M")
        % Bax(bf=None, s1=2, s2=None, state="M")
        % Bcl2(bf=1): Model.r_Bax2_Bcl2.AB,
        Bax(bf=None, s1=1, s2=2, state="M")
        % Bax(bf=None, s1=3, s2=1, state="M")
        % Bax(bf=None, s1=4, s2=3, state="M")
        % Bax(bf=None, s1=2, s2=4, state="M"): Model.Bax_M4,
        Bax(bf=None, s1=1, s2=None, state="M")
        % Bax(bf=None, s1=None, s2=1, state="M"): Model.Bax_M2,
        Bax(bf=None, s1=None, s2=None, state="M"): Model.mitocondria.Bax_M,
        Bax(bf=1, s1=2, s2=3, state="M")
        % Bax(bf=None, s1=4, s2=2, state="M")
        % Bax(bf=None, s1=5, s2=4, state="M")
        % Bax(bf=None, s1=3, s2=5, state="M")
        % Mito(bf=1, state="I"): Model.r_Bax4_Mito.AB,
        CytoC(bf=1, state="M") % Mito(bf=1, state="A"): Model.r_CytoC_pore.ES,
        Mito(bf=1, state="A") % Smac(bf=1, state="M"): Model.r_Smac_pore.ES,
        Mito(bf=None, state="A"): Model.Mito_A,
        Mito(bf=None, state="I"): Model.Mito_I,
        Bcl2c(b=None): Model.Bcl2c,
        Bcl2c(b=1) % Bid(bf=1, state="T"): Model.r_Bcl2_Bid.AB,
    }

    times = np.linspace(0, 10_000, 1000)
    sim = Simulator(Model)
    df_simbio = common.run_simbio(sim, times)
    df_pysb = common.run_pysb(model, times)
    df_pysb = df_pysb.rename(columns=str).rename(
        columns={str(k): str(v) for k, v in mapping.items()},
    )
    assert len(df_simbio.columns.symmetric_difference(df_pysb.columns)) == 0
    pd.testing.assert_frame_equal(
        df_simbio,
        df_pysb,
        atol=1e-4,
        rtol=1e-4,
        check_like=True,  # ignore order
    )


def test_corbat_2018():
    from caspase_model.models import corbat_2018

    from ..corbat import Corbat2018 as Model

    model = corbat_2018(stimuli="extrinsic", add_CASPAM=False)
    pysb.util.alias_model_components()

    mapping = {
        L(bf=None): Model.albeck.L,
        R(bf=None): Model.albeck.R,
        flip(bf=None): Model.albeck.flip,
        C8(bf=None, state="pro"): Model.albeck.C8_pro,
        BAR(bf=None): Model.albeck.BAR,
        L(bf=1) % R(bf=1): Model.albeck.r_L_R.AB,
        DISC(bf=None): Model.albeck.DISC,
        Bid(bf=1, state="U") % C8(bf=1, state="A"): Model.albeck.r_C8_Bid.ES,
        DISC(bf=1) % flip(bf=1): Model.albeck.r_DISC_flip.AB,
        C8(bf=None, state="A"): Model.albeck.C8_A,
        BAR(bf=1) % C8(bf=1, state="A"): Model.albeck.r_BAR_C8.AB,
        Bid(bf=None, state="U"): Model.albeck.Bid_U,
        Bid(bf=None, state="T"): Model.albeck.Bid_T,
        C8(bf=1, state="pro") % DISC(bf=1): Model.albeck.r_DISC_C8.ES,
        Bax(bf=None, s1=None, s2=None, state="C"): Model.albeck.mitocondria.Bax_C,
        Bax(bf=None, s1=None, s2=None, state="A"): Model.albeck.mitocondria.Bax_A,
        Bcl2(bf=None): Model.albeck.mitocondria.Bcl2,
        CytoC(bf=None, state="M"): Model.albeck.mitocondria.CytoC_M,
        CytoC(bf=None, state="A"): Model.albeck.mitocondria.CytoC_A,
        CytoC(bf=None, state="C"): Model.albeck.mitocondria.CytoC_C,
        Smac(bf=None, state="M"): Model.albeck.mitocondria.Smac_M,
        Smac(bf=None, state="A"): Model.albeck.mitocondria.Smac_A,
        Smac(bf=None, state="C"): Model.albeck.mitocondria.Smac_C,
        Apaf(bf=None, state="I"): Model.albeck.Apaf_I,
        Apaf(bf=None, state="A"): Model.albeck.Apaf_A,
        Apop(bf=None): Model.albeck.Apop,
        C3(bf=None, state="pro"): Model.albeck.C3_pro,
        C3(bf=None, state="A"): Model.albeck.C3_A,
        C3(bf=None, state="ub"): Model.albeck.C3_ub,
        C6(bf=None, state="pro"): Model.albeck.C6_pro,
        C6(bf=None, state="A"): Model.albeck.C6_A,
        C9(bf=None): Model.albeck.C9,
        PARP(bf=None, state="U"): Model.albeck.PARP_U,
        PARP(bf=None, state="C"): Model.albeck.PARP_C,
        XIAP(bf=None): Model.albeck.XIAP,
        Apaf(bf=1, state="I") % CytoC(bf=1, state="A"): Model.albeck.r_CytoC_Apaf.ES,
        Apop(bf=1) % C3(bf=1, state="pro"): Model.albeck.r_Apop_C3.ES,
        Apop(bf=1) % XIAP(bf=1): Model.albeck.r_Apop_XIAP.AB,
        Bax(bf=1, s1=None, s2=None, state="M") % Bcl2(bf=1): Model.albeck.r_Bax_Bcl2.AB,
        Bax(bf=1, s1=None, s2=None, state="C")
        % Bid(bf=1, state="T"): Model.albeck.r_Bid_Bax.ES,
        C3(bf=1, state="A") % C6(bf=1, state="pro"): Model.albeck.r_C3_C6.ES,
        C3(bf=1, state="A") % PARP(bf=1, state="U"): Model.albeck.r_C3_PARP.ES,
        C3(bf=1, state="A") % XIAP(bf=1): Model.albeck.r_XIAP_C3.ES,
        C3(bf=1, state="pro") % C8(bf=1, state="A"): Model.albeck.r_C8_C3.ES,
        C6(bf=1, state="A") % C8(bf=1, state="pro"): Model.albeck.r_C6_C8.ES,
        Smac(bf=1, state="A") % XIAP(bf=1): Model.albeck.r_Smac_XIAP.AB,
        # NEw
        Bax(bf=1, s1=2, s2=3, state="M")
        % Bax(bf=None, s1=4, s2=2, state="M")
        % Bax(bf=None, s1=5, s2=4, state="M")
        % Bax(bf=None, s1=3, s2=5, state="M")
        % Bcl2(bf=1): Model.albeck.r_Bax4_Bcl2.AB,
        Bax(bf=1, s1=None, s2=2, state="M")
        % Bax(bf=None, s1=2, s2=None, state="M")
        % Bcl2(bf=1): Model.albeck.r_Bax2_Bcl2.AB,
        Bax(bf=None, s1=1, s2=2, state="M")
        % Bax(bf=None, s1=3, s2=1, state="M")
        % Bax(bf=None, s1=4, s2=3, state="M")
        % Bax(bf=None, s1=2, s2=4, state="M"): Model.albeck.Bax_M4,
        Bax(bf=None, s1=1, s2=None, state="M")
        % Bax(bf=None, s1=None, s2=1, state="M"): Model.albeck.Bax_M2,
        Bax(bf=None, s1=None, s2=None, state="M"): Model.albeck.mitocondria.Bax_M,
        Bax(bf=1, s1=2, s2=3, state="M")
        % Bax(bf=None, s1=4, s2=2, state="M")
        % Bax(bf=None, s1=5, s2=4, state="M")
        % Bax(bf=None, s1=3, s2=5, state="M")
        % Mito(bf=1, state="I"): Model.albeck.r_Bax4_Mito.AB,
        CytoC(bf=1, state="M") % Mito(bf=1, state="A"): Model.albeck.r_CytoC_pore.ES,
        Mito(bf=1, state="A") % Smac(bf=1, state="M"): Model.albeck.r_Smac_pore.ES,
        Mito(bf=None, state="A"): Model.albeck.Mito_A,
        Mito(bf=None, state="I"): Model.albeck.Mito_I,
        Bcl2c(b=None): Model.albeck.Bcl2c,
        Bcl2c(b=1) % Bid(bf=1, state="T"): Model.albeck.r_Bcl2_Bid.AB,
    }

    times = np.linspace(0, 10_000, 1000)
    sim = Simulator(Model)
    df_simbio = common.run_simbio(sim, times)
    df_pysb = common.run_pysb(model, times)
    df_pysb = df_pysb.rename(columns=str).rename(
        columns={str(k): str(v) for k, v in mapping.items()},
    )
    assert len(df_simbio.columns.symmetric_difference(df_pysb.columns)) == 0
    pd.testing.assert_frame_equal(
        df_simbio,
        df_pysb,
        atol=1e-4,
        rtol=1e-4,
        check_like=True,  # ignore order
    )
