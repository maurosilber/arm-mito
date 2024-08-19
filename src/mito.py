from simbio import Compartment, Parameter, Species, assign, initial, reactions


class Mitochondria(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    area: Parameter = assign(default=0.07)
    pore_transport_rate: Parameter = assign(default=10)
    transloc_rate: Parameter = assign(default=1e-2)

    CytoC_C: Species = initial()
    Smac_C: Species = initial()
    Bax_A: Species = initial()

    Bcl2: Species = initial(default=2e4)
    CytoC_M: Species = initial(default=1e5)
    Smac_M: Species = initial(default=1e5)
    Mito_I: Species = initial(default=5e5)
    Mito_A: Species = initial(default=0)
    Bax: Species = initial(default=0)
    Bax2: Species = initial(default=0)
    Bax4: Species = initial(default=0)

    r_Bax_transloc = reactions.Equilibration(
        A=Bax_A,
        B=Bax,
        forward_rate=transloc_rate,
        reverse_rate=transloc_rate,
    )

    r_Bax_dimerization = reactions.Equilibration(
        A=2 * Bax,
        B=Bax2,
        forward_rate=KF / area,
        reverse_rate=KR,
    )
    r_Bax_tetramerization = reactions.Equilibration(
        A=2 * Bax2,
        B=Bax4,
        forward_rate=KF / area,
        reverse_rate=KR,
    )
    r_Bax_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax,
        B=Bcl2,
        AB=0,
        forward_rate=KF / area,
        reverse_rate=KR,
    )
    r_Bax2_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax2,
        B=Bcl2,
        AB=0,
        forward_rate=KF / area,
        reverse_rate=KR,
    )
    r_Bax4_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax4,
        B=Bcl2,
        AB=0,
        forward_rate=KF / area,
        reverse_rate=KR,
    )
    r_Bax4_Mito = reactions.CatalyzeConvert(
        A=Bax4,
        B=Mito_I,
        AB=0,
        P=Mito_A,
        forward_rate=KF / area,
        reverse_rate=KR,
        conversion_rate=KC,
    )
    r_Smac_pore = reactions.MichaelisMenten(
        E=Mito_A,
        S=Smac_M,
        ES=0,
        P=Smac_C,
        forward_rate=2 * KF / area,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )
    r_CytoC_pore = reactions.MichaelisMenten(
        E=Mito_A,
        S=CytoC_M,
        ES=0,
        P=CytoC_C,
        forward_rate=2 * KF / area,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )


class ARM_Cito(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    transloc_rate: Parameter = assign(default=1e-2)

    L: Species = initial(default=0)
    R: Species = initial(default=200)
    DISC: Species = initial(default=0)
    flip: Species = initial(default=100)
    C8_pro: Species = initial(default=20_000)
    C8_A: Species = initial(default=0)
    BAR: Species = initial(default=1_000)
    Bid_U: Species = initial(default=4e4)
    Bid_T: Species = initial(default=0)
    Bid_M: Species = initial(default=0)
    Bax_C: Species = initial(default=1e5)
    Bax_A: Species = initial(default=0)
    CytoC_C: Species = initial(default=0)
    CytoC_A: Species = initial(default=0)
    Smac_C: Species = initial(default=0)
    Smac_A: Species = initial(default=0)
    Apaf_I: Species = initial(default=1e3)
    Apaf_A: Species = initial(default=0)
    Apop: Species = initial(default=0)
    C3_pro: Species = initial(default=1e4)
    C3_A: Species = initial(default=0)
    C3_ub: Species = initial(default=0)
    C6_pro: Species = initial(default=1e4)
    C6_A: Species = initial(default=0)
    PARP_U: Species = initial(default=1e6)
    PARP_C: Species = initial(default=0)
    XIAP: Species = initial(default=1e4)
    Bcl2c: Species = initial(default=2e4)
    IntrinsicStimuli: Species = initial(default=0)

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
        forward_rate=4e-7,
        reverse_rate=KR,
        conversion_rate=1e-5,
    )
    r_DISC_C8 = reactions.MichaelisMenten(
        E=DISC,
        S=C8_pro,
        ES=0,
        P=C8_A,
        forward_rate=KF,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_C8_Bid = reactions.MichaelisMenten(
        E=C8_A,
        S=Bid_U,
        ES=0,
        P=Bid_T,
        forward_rate=1e-7,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_DISC_flip = reactions.ReversibleSynthesis(
        A=DISC,
        B=flip,
        AB=0,
        forward_rate=KF,
        reverse_rate=KR,
    )
    r_BAR_C8 = reactions.ReversibleSynthesis(
        A=BAR,
        B=C8_A,
        AB=0,
        forward_rate=KF,
        reverse_rate=KR,
    )
    r_CytoC_Apaf = reactions.MichaelisMenten(
        E=CytoC_A,
        S=Apaf_I,
        ES=0,
        P=Apaf_A,
        forward_rate=5e-7,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_ApafA_C3pro = reactions.MichaelisMenten(
        E=Apaf_A,
        S=C3_pro,
        ES=0,
        P=C3_A,
        forward_rate=5e-09,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_C3A_ApafA = reactions.MichaelisMenten(
        E=C3_A,
        S=Apaf_A,
        ES=0,
        P=Apop,
        forward_rate=1.3e-06,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_Apop_C3 = reactions.MichaelisMenten(
        E=Apop,
        S=C3_pro,
        ES=0,
        P=C3_A,
        forward_rate=5e-9,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_ApafA_XIAP = reactions.ReversibleSynthesis(
        A=Apaf_A,
        B=XIAP,
        AB=0,
        forward_rate=2e-6,
        reverse_rate=KR,
    )
    r_Smac_XIAP = reactions.ReversibleSynthesis(
        A=Smac_A,
        B=XIAP,
        AB=0,
        forward_rate=7e-6,
        reverse_rate=KR,
    )
    r_C8_C3 = reactions.MichaelisMenten(
        E=C8_A,
        S=C3_pro,
        ES=0,
        P=C3_A,
        forward_rate=1e-7,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_XIAP_C3 = reactions.MichaelisMenten(
        E=XIAP,
        S=C3_A,
        ES=0,
        P=C3_ub,
        forward_rate=2e-6,
        reverse_rate=KR,
        catalytic_rate=1e-1,
    )
    r_C3_PARP = reactions.MichaelisMenten(
        E=C3_A,
        S=PARP_U,
        ES=0,
        P=PARP_C,
        forward_rate=KF,
        reverse_rate=1e-2,
        catalytic_rate=KC,
    )
    r_C3_C6 = reactions.MichaelisMenten(
        E=C3_A,
        S=C6_pro,
        ES=0,
        P=C6_A,
        forward_rate=KF,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_C6_C8 = reactions.MichaelisMenten(
        E=C6_A,
        S=C8_pro,
        ES=0,
        P=C8_A,
        forward_rate=3e-8,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_Bid_Bax = reactions.MichaelisMenten(
        E=Bid_T,
        S=Bax_C,
        ES=0,
        P=Bax_A,
        forward_rate=1e-7,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_Bcl2_Bid = reactions.ReversibleSynthesis(
        A=Bid_T,
        B=Bcl2c,
        AB=0,
        forward_rate=KF,
        reverse_rate=KR,
    )
    r_intrinsic = reactions.MichaelisMenten(
        E=IntrinsicStimuli,
        S=Bid_U,
        ES=0,
        P=Bid_T,
        forward_rate=KF,
        reverse_rate=KR,
        catalytic_rate=KC,
    )


class ARM(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    transloc_rate: Parameter = assign(default=1e-2)
    cytoplasm = ARM_Cito()
    mitocondria = Mitochondria(
        KF=KF,
        KR=KR,
        KC=KC,
        transloc_rate=transloc_rate,
        CytoC_C=cytoplasm.CytoC_C,
        Smac_C=cytoplasm.Smac_C,
        Bax_A=cytoplasm.Bax_A,
    )
