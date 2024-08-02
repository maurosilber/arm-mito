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
        Bax(bf=None, s1=None, s2=None, state="C"): Model.Bax_C,
        Bax(bf=None, s1=None, s2=None, state="A"): Model.Bax_A,
        Bcl2(bf=None): Model.Bcl2,
        CytoC(bf=None, state="M"): Model.CytoC_M,
        CytoC(bf=None, state="A"): Model.CytoC_A,
        CytoC(bf=None, state="C"): Model.CytoC_C,
        Smac(bf=None, state="M"): Model.Smac_M,
        Smac(bf=None, state="A"): Model.Smac_A,
        Smac(bf=None, state="C"): Model.Smac_C,
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
        Bax(bf=None, s1=None, s2=None, state="M"): Model.Bax_M,
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

    times = np.linspace(0, 30_000, 1000)
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
        Bax(bf=None, s1=None, s2=None, state="C"): Model.albeck.Bax_C,
        Bax(bf=None, s1=None, s2=None, state="A"): Model.albeck.Bax_A,
        Bcl2(bf=None): Model.albeck.Bcl2,
        CytoC(bf=None, state="M"): Model.albeck.CytoC_M,
        CytoC(bf=None, state="A"): Model.albeck.CytoC_A,
        CytoC(bf=None, state="C"): Model.albeck.CytoC_C,
        Smac(bf=None, state="M"): Model.albeck.Smac_M,
        Smac(bf=None, state="A"): Model.albeck.Smac_A,
        Smac(bf=None, state="C"): Model.albeck.Smac_C,
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
        Bax(bf=None, s1=None, s2=None, state="M"): Model.albeck.Bax_M,
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

    times = np.linspace(0, 30_000, 1000)
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


def test_arm():
    from caspase_model.models import arm

    from ..corbat import ARM as Model

    model = arm(stimuli="extrinsic", add_CASPAM=False)
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
        Bax(bf=None, s1=None, s2=None, state="C"): Model.Bax_C,
        Bax(bf=None, s1=None, s2=None, state="A"): Model.Bax_A,
        Bcl2(bf=None): Model.Bcl2,
        CytoC(bf=None, state="M"): Model.CytoC_M,
        CytoC(bf=None, state="A"): Model.CytoC_A,
        CytoC(bf=None, state="C"): Model.CytoC_C,
        Smac(bf=None, state="M"): Model.Smac_M,
        Smac(bf=None, state="A"): Model.Smac_A,
        Smac(bf=None, state="C"): Model.Smac_C,
        Apaf(bf=None, state="I"): Model.Apaf_I,
        Apaf(bf=None, state="A"): Model.Apaf_A,
        Apop(bf=None): Model.Apop,
        C3(bf=None, state="pro"): Model.C3_pro,
        C3(bf=None, state="A"): Model.C3_A,
        C3(bf=None, state="ub"): Model.C3_ub,
        C6(bf=None, state="pro"): Model.C6_pro,
        C6(bf=None, state="A"): Model.C6_A,
        PARP(bf=None, state="U"): Model.PARP_U,
        PARP(bf=None, state="C"): Model.PARP_C,
        XIAP(bf=None): Model.XIAP,
        Apaf(bf=1, state="I") % CytoC(bf=1, state="A"): Model.r_CytoC_Apaf.ES,
        Apop(bf=1) % C3(bf=1, state="pro"): Model.r_Apop_C3.ES,
        Apop(bf=1) % XIAP(bf=1): Model.r_ApafA_XIAP.AB,
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
        Bax(bf=None, s1=None, s2=None, state="M"): Model.Bax_M,
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
        Apaf(bf=1, state='A') % C3(bf=1, state='A'): Model.r_C3A_ApafA.ES,
        Apaf(bf=1, state='A') % C3(bf=1, state='pro'): Model.r_ApafA_C3pro.ES,
        Apaf(bf=1, state='A') % XIAP(bf=1): Model.r_ApafA_XIAP.AB,
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


def test_intrinsic_arm():
    from caspase_model.models import arm

    from ..corbat import IntrinsicARM as Model

    model = arm(stimuli="intrinsic", add_CASPAM=False)
    pysb.util.alias_model_components()

    mapping = {
        L(bf=None): Model.arm.L,
        R(bf=None): Model.arm.R,
        flip(bf=None): Model.arm.flip,
        C8(bf=None, state="pro"): Model.arm.C8_pro,
        BAR(bf=None): Model.arm.BAR,
        L(bf=1) % R(bf=1): Model.arm.r_L_R.AB,
        DISC(bf=None): Model.arm.DISC,
        Bid(bf=1, state="U") % C8(bf=1, state="A"): Model.arm.r_C8_Bid.ES,
        DISC(bf=1) % flip(bf=1): Model.arm.r_DISC_flip.AB,
        C8(bf=None, state="A"): Model.arm.C8_A,
        BAR(bf=1) % C8(bf=1, state="A"): Model.arm.r_BAR_C8.AB,
        Bid(bf=None, state="U"): Model.arm.Bid_U,
        Bid(bf=None, state="T"): Model.arm.Bid_T,
        C8(bf=1, state="pro") % DISC(bf=1): Model.arm.r_DISC_C8.ES,
        Bax(bf=None, s1=None, s2=None, state="C"): Model.arm.Bax_C,
        Bax(bf=None, s1=None, s2=None, state="A"): Model.arm.Bax_A,
        Bcl2(bf=None): Model.arm.Bcl2,
        CytoC(bf=None, state="M"): Model.arm.CytoC_M,
        CytoC(bf=None, state="A"): Model.arm.CytoC_A,
        CytoC(bf=None, state="C"): Model.arm.CytoC_C,
        Smac(bf=None, state="M"): Model.arm.Smac_M,
        Smac(bf=None, state="A"): Model.arm.Smac_A,
        Smac(bf=None, state="C"): Model.arm.Smac_C,
        Apaf(bf=None, state="I"): Model.arm.Apaf_I,
        Apaf(bf=None, state="A"): Model.arm.Apaf_A,
        Apop(bf=None): Model.arm.Apop,
        C3(bf=None, state="pro"): Model.arm.C3_pro,
        C3(bf=None, state="A"): Model.arm.C3_A,
        C3(bf=None, state="ub"): Model.arm.C3_ub,
        C6(bf=None, state="pro"): Model.arm.C6_pro,
        C6(bf=None, state="A"): Model.arm.C6_A,
        PARP(bf=None, state="U"): Model.arm.PARP_U,
        PARP(bf=None, state="C"): Model.arm.PARP_C,
        XIAP(bf=None): Model.arm.XIAP,
        Apaf(bf=1, state="I") % CytoC(bf=1, state="A"): Model.arm.r_CytoC_Apaf.ES,
        Apop(bf=1) % C3(bf=1, state="pro"): Model.arm.r_Apop_C3.ES,
        Apop(bf=1) % XIAP(bf=1): Model.arm.r_ApafA_XIAP.AB,
        Bax(bf=1, s1=None, s2=None, state="M") % Bcl2(bf=1): Model.arm.r_Bax_Bcl2.AB,
        Bax(bf=1, s1=None, s2=None, state="C")
        % Bid(bf=1, state="T"): Model.arm.r_Bid_Bax.ES,
        C3(bf=1, state="A") % C6(bf=1, state="pro"): Model.arm.r_C3_C6.ES,
        C3(bf=1, state="A") % PARP(bf=1, state="U"): Model.arm.r_C3_PARP.ES,
        C3(bf=1, state="A") % XIAP(bf=1): Model.arm.r_XIAP_C3.ES,
        C3(bf=1, state="pro") % C8(bf=1, state="A"): Model.arm.r_C8_C3.ES,
        C6(bf=1, state="A") % C8(bf=1, state="pro"): Model.arm.r_C6_C8.ES,
        Smac(bf=1, state="A") % XIAP(bf=1): Model.arm.r_Smac_XIAP.AB,
        # NEw
        Bax(bf=1, s1=2, s2=3, state="M")
        % Bax(bf=None, s1=4, s2=2, state="M")
        % Bax(bf=None, s1=5, s2=4, state="M")
        % Bax(bf=None, s1=3, s2=5, state="M")
        % Bcl2(bf=1): Model.arm.r_Bax4_Bcl2.AB,
        Bax(bf=1, s1=None, s2=2, state="M")
        % Bax(bf=None, s1=2, s2=None, state="M")
        % Bcl2(bf=1): Model.arm.r_Bax2_Bcl2.AB,
        Bax(bf=None, s1=1, s2=2, state="M")
        % Bax(bf=None, s1=3, s2=1, state="M")
        % Bax(bf=None, s1=4, s2=3, state="M")
        % Bax(bf=None, s1=2, s2=4, state="M"): Model.arm.Bax_M4,
        Bax(bf=None, s1=1, s2=None, state="M")
        % Bax(bf=None, s1=None, s2=1, state="M"): Model.arm.Bax_M2,
        Bax(bf=None, s1=None, s2=None, state="M"): Model.arm.Bax_M,
        Bax(bf=1, s1=2, s2=3, state="M")
        % Bax(bf=None, s1=4, s2=2, state="M")
        % Bax(bf=None, s1=5, s2=4, state="M")
        % Bax(bf=None, s1=3, s2=5, state="M")
        % Mito(bf=1, state="I"): Model.arm.r_Bax4_Mito.AB,
        CytoC(bf=1, state="M") % Mito(bf=1, state="A"): Model.arm.r_CytoC_pore.ES,
        Mito(bf=1, state="A") % Smac(bf=1, state="M"): Model.arm.r_Smac_pore.ES,
        Mito(bf=None, state="A"): Model.arm.Mito_A,
        Mito(bf=None, state="I"): Model.arm.Mito_I,
        Bcl2c(b=None): Model.arm.Bcl2c,
        Bcl2c(b=1) % Bid(bf=1, state="T"): Model.arm.r_Bcl2_Bid.AB,
        Apaf(bf=1, state='A') % C3(bf=1, state='A'): Model.arm.r_C3A_ApafA.ES,
        Apaf(bf=1, state='A') % C3(bf=1, state='pro'): Model.arm.r_ApafA_C3pro.ES,
        Apaf(bf=1, state='A') % XIAP(bf=1): Model.arm.r_ApafA_XIAP.AB,
        IntrinsicStimuli(bf=None): Model.IntrinsicStimuli,
        Bid(bf=1, state='U') % IntrinsicStimuli(bf=1): Model.r_intrinsic.ES,
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
