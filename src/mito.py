from simbio import (
    Compartment,
    Constant,
    MassAction,
    Parameter,
    Species,
    assign,
    initial,
    reactions,
)


def _add_volume_factor(model: type[Compartment], volume):
    for r in model._yield(MassAction):
        power = 1 - sum(s.stoichiometry for s in r.reactants)
        match power:
            case 0:
                pass
            case 1:
                r.rate = r.rate * volume
            case -1:
                r.rate = r.rate / volume
            case _:
                raise NotImplementedError

        r.equations = tuple(r._yield_equations())


class Mitochondria(Compartment):
    volume: Constant = assign(default=0.07, constant=True)
    _Albeck_volume_fraction = 0.07

    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    pore_transport_rate: Parameter = assign(
        default=10 / _Albeck_volume_fraction * volume
    )
    transloc_rate: Parameter = assign(default=1e-2 / _Albeck_volume_fraction * volume)

    CytoC_C: Species = initial()
    Smac_C: Species = initial()
    Bax_A: Species = initial()

    # Albeck concentrations
    Bcl2_0: Constant = assign(default=2e4 / _Albeck_volume_fraction, constant=True)
    CytoC_M_0: Constant = assign(default=1e5 / _Albeck_volume_fraction, constant=True)
    Smac_M_0: Constant = assign(default=1e5 / _Albeck_volume_fraction, constant=True)
    Mito_I_0: Constant = assign(default=5e5 / _Albeck_volume_fraction, constant=True)
    Mito_A_0: Constant = assign(default=0 / _Albeck_volume_fraction, constant=True)
    Bax_0: Constant = assign(default=0 / _Albeck_volume_fraction, constant=True)
    Bax2_0: Constant = assign(default=0 / _Albeck_volume_fraction, constant=True)
    Bax4_0: Constant = assign(default=0 / _Albeck_volume_fraction, constant=True)

    Bcl2: Species = initial(default=Bcl2_0 * volume)
    CytoC_M: Species = initial(default=CytoC_M_0 * volume)
    Smac_M: Species = initial(default=Smac_M_0 * volume)
    Mito_I: Species = initial(default=Mito_I_0 * volume)
    Mito_A: Species = initial(default=Mito_A_0 * volume)
    Bax: Species = initial(default=Bax_0 * volume)
    Bax2: Species = initial(default=Bax2_0 * volume)
    Bax4: Species = initial(default=Bax4_0 * volume)

    r_Bax_transloc = reactions.Equilibration(
        A=Bax_A,
        B=Bax,
        forward_rate=transloc_rate,
        reverse_rate=transloc_rate,
    )

    r_Bax_dimerization = reactions.Equilibration(
        A=2 * Bax,
        B=Bax2,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_Bax_tetramerization = reactions.Equilibration(
        A=2 * Bax2,
        B=Bax4,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_Bax_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax,
        B=Bcl2,
        AB=0,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_Bax2_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax2,
        B=Bcl2,
        AB=0,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_Bax4_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax4,
        B=Bcl2,
        AB=0,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_Bax4_Mito = reactions.CatalyzeConvert(
        A=Bax4,
        B=Mito_I,
        AB=0,
        P=Mito_A,
        forward_rate=KF / volume,
        reverse_rate=KR,
        conversion_rate=KC,
    )
    r_Smac_pore = reactions.MichaelisMenten(
        E=Mito_A,
        S=Smac_M,
        ES=0,
        P=Smac_C,
        forward_rate=2 * KF / volume,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )
    r_CytoC_pore = reactions.MichaelisMenten(
        E=Mito_A,
        S=CytoC_M,
        ES=0,
        P=CytoC_C,
        forward_rate=2 * KF / volume,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )


class ARM_Cito(Compartment):
    volume: Constant = assign(default=1, constant=True)

    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    transloc_rate: Parameter = assign(default=1e-2)

    L: Species = initial(default=0 * volume)
    R: Species = initial(default=200 * volume)
    DISC: Species = initial(default=0 * volume)
    flip: Species = initial(default=100 * volume)
    C8_pro: Species = initial(default=20_000 * volume)
    C8_A: Species = initial(default=0 * volume)
    BAR: Species = initial(default=1_000 * volume)
    Bid_U: Species = initial(default=4e4 * volume)
    Bid_T: Species = initial(default=0 * volume)
    Bax_C: Species = initial(default=1e5 * volume)
    Bax_A: Species = initial(default=0 * volume)
    CytoC_C: Species = initial(default=0 * volume)
    CytoC_A: Species = initial(default=0 * volume)
    Smac_C: Species = initial(default=0 * volume)
    Smac_A: Species = initial(default=0 * volume)
    Apaf_I: Species = initial(default=1e3 * volume)
    Apaf_A: Species = initial(default=0 * volume)
    Apop: Species = initial(default=0 * volume)
    C3_pro: Species = initial(default=1e4 * volume)
    C3_A: Species = initial(default=0 * volume)
    C3_ub: Species = initial(default=0 * volume)
    C6_pro: Species = initial(default=1e4 * volume)
    C6_A: Species = initial(default=0 * volume)
    PARP_U: Species = initial(default=1e6 * volume)
    PARP_C: Species = initial(default=0 * volume)
    XIAP: Species = initial(default=1e4 * volume)
    Bcl2c: Species = initial(default=2e4 * volume)
    IntrinsicStimuli: Species = initial(default=0 * volume)

    r_Smac_transloc = reactions.Equilibration(
        A=Smac_C,
        B=Smac_A,
        forward_rate=transloc_rate,
        reverse_rate=transloc_rate,
    )
    r_CytoC_transloc = reactions.Equilibration(
        A=CytoC_C,
        B=CytoC_A,
        forward_rate=transloc_rate,
        reverse_rate=transloc_rate,
    )

    # Reactions
    r_L_R = reactions.CatalyzeConvert(
        A=L,
        B=R,
        AB=0,
        P=DISC,
        forward_rate=4e-7 / volume,
        reverse_rate=KR,
        conversion_rate=1e-5,
    )
    r_DISC_C8 = reactions.MichaelisMenten(
        E=DISC,
        S=C8_pro,
        ES=0,
        P=C8_A,
        forward_rate=KF / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_C8_Bid = reactions.MichaelisMenten(
        E=C8_A,
        S=Bid_U,
        ES=0,
        P=Bid_T,
        forward_rate=1e-7 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_DISC_flip = reactions.ReversibleSynthesis(
        A=DISC,
        B=flip,
        AB=0,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_BAR_C8 = reactions.ReversibleSynthesis(
        A=BAR,
        B=C8_A,
        AB=0,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_CytoC_Apaf = reactions.MichaelisMenten(
        E=CytoC_A,
        S=Apaf_I,
        ES=0,
        P=Apaf_A,
        forward_rate=5e-7 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_ApafA_C3pro = reactions.MichaelisMenten(
        E=Apaf_A,
        S=C3_pro,
        ES=0,
        P=C3_A,
        forward_rate=5e-09 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_C3A_ApafA = reactions.MichaelisMenten(
        E=C3_A,
        S=Apaf_A,
        ES=0,
        P=Apop,
        forward_rate=1.3e-06 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_Apop_C3 = reactions.MichaelisMenten(
        E=Apop,
        S=C3_pro,
        ES=0,
        P=C3_A,
        forward_rate=5e-9 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_ApafA_XIAP = reactions.ReversibleSynthesis(
        A=Apaf_A,
        B=XIAP,
        AB=0,
        forward_rate=2e-6 / volume,
        reverse_rate=KR,
    )
    r_Smac_XIAP = reactions.ReversibleSynthesis(
        A=Smac_A,
        B=XIAP,
        AB=0,
        forward_rate=7e-6 / volume,
        reverse_rate=KR,
    )
    r_C8_C3 = reactions.MichaelisMenten(
        E=C8_A,
        S=C3_pro,
        ES=0,
        P=C3_A,
        forward_rate=1e-7 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_XIAP_C3 = reactions.MichaelisMenten(
        E=XIAP,
        S=C3_A,
        ES=0,
        P=C3_ub,
        forward_rate=2e-6 / volume,
        reverse_rate=KR,
        catalytic_rate=1e-1,
    )
    r_C3_PARP = reactions.MichaelisMenten(
        E=C3_A,
        S=PARP_U,
        ES=0,
        P=PARP_C,
        forward_rate=KF / volume,
        reverse_rate=1e-2,
        catalytic_rate=KC,
    )
    r_C3_C6 = reactions.MichaelisMenten(
        E=C3_A,
        S=C6_pro,
        ES=0,
        P=C6_A,
        forward_rate=KF / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_C6_C8 = reactions.MichaelisMenten(
        E=C6_A,
        S=C8_pro,
        ES=0,
        P=C8_A,
        forward_rate=3e-8 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_Bid_Bax = reactions.MichaelisMenten(
        E=Bid_T,
        S=Bax_C,
        ES=0,
        P=Bax_A,
        forward_rate=1e-7 / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_Bcl2_Bid = reactions.ReversibleSynthesis(
        A=Bid_T,
        B=Bcl2c,
        AB=0,
        forward_rate=KF / volume,
        reverse_rate=KR,
    )
    r_intrinsic = reactions.MichaelisMenten(
        E=IntrinsicStimuli,
        S=Bid_U,
        ES=0,
        P=Bid_T,
        forward_rate=KF / volume,
        reverse_rate=KR,
        catalytic_rate=KC,
    )


class ARM(Compartment):
    volume: Constant = assign(default=1, constant=True)
    mitocondria_volume_fraction: Constant = assign(default=0.07, constant=True)
    mitocondria_volume: Constant = assign(
        default=volume * mitocondria_volume_fraction, constant=True
    )

    L_concentration: Constant = assign(default=1000, constant=True)
    IntrinsicStimuli_concentration: Constant = assign(default=0, constant=True)
    L: Constant = assign(default=L_concentration * volume, constant=True)
    IntrinsicStimuli: Constant = assign(
        default=IntrinsicStimuli_concentration * volume, constant=True
    )

    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    transloc_rate: Parameter = assign(default=1e-2)
    cytoplasm = ARM_Cito(
        volume=volume,
        KF=KF,
        KR=KR,
        KC=KC,
        transloc_rate=transloc_rate,
        L=L,
        IntrinsicStimuli=IntrinsicStimuli,
    )
    mitocondria = Mitochondria(
        volume=mitocondria_volume,
        KF=KF,
        KR=KR,
        KC=KC,
        transloc_rate=transloc_rate,
        CytoC_C=cytoplasm.CytoC_C,
        Smac_C=cytoplasm.Smac_C,
        Bax_A=cytoplasm.Bax_A,
    )
