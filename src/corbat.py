from simbio import Compartment, Parameter, Species, assign, initial, reactions

from .albeck import MOMP, Intrinsic


class AlbeckAsMatlab(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    mitocondria_volume_fraction: Parameter = assign(default=0.07)

    L: Species = initial(default=3_000)
    R: Species = initial(default=200)
    DISC: Species = initial(default=0)
    flip: Species = initial(default=100)
    C8_pro: Species = initial(default=20_000)
    C8_A: Species = initial(default=0)
    BAR: Species = initial(default=1_000)

    Bid_U: Species = initial(default=4e4)
    Bid_T: Species = initial(default=0)

    # =====================
    # tBID Activation Rules
    # ---------------------
    #        L + R <--> L:R --> DISC
    #        pC8 + DISC <--> DISC:pC8 --> C8 + DISC
    #        Bid + C8 <--> Bid:C8 --> tBid + C8
    # ---------------------

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

    # ---------------------
    # Inhibition Rules
    # ---------------------
    #        flip + DISC <-->  flip:DISC
    #        C8 + BAR <--> BAR:C8
    # ---------------------

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

    mitocondria = MOMP(
        KF=KF,
        KR=KR,
        KC=KC,
        Bid_U=Bid_U,
        Bid_T=Bid_T,
    )
    intrinsic = Intrinsic(
        KF=KF,
        KR=KR,
        KC=KC,
        CytoC_A=mitocondria.CytoC_A,
        Smac_A=mitocondria.Smac_A,
        C8_pro=C8_pro,
        C8_A=C8_A,
    )
    # MOMP Mechanism
    r_Bid_Bax = reactions.MichaelisMenten(
        E=mitocondria.Bid_T,
        S=mitocondria.Bax_C,
        ES=0,
        P=mitocondria.Bax_A,
        forward_rate=1e-7,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    # Bax dimerizes/tetramerizes
    Bax_M2: Species = initial(default=0)
    Bax_M4: Species = initial(default=0)
    r_Bax_dimerization = reactions.Equilibration(
        A=2 * mitocondria.Bax_M,
        B=Bax_M2,
        forward_rate=KF / mitocondria_volume_fraction,
        reverse_rate=KR,
    )
    r_Bax_tetramerization = reactions.Equilibration(
        A=2 * Bax_M2,
        B=Bax_M4,
        forward_rate=KF / mitocondria_volume_fraction,
        reverse_rate=KR,
    )

    # Bcl2 inhibits Bax, Bax2, and Bax4
    r_Bax_Bcl2 = reactions.ReversibleSynthesis(
        A=mitocondria.Bax_M,
        B=mitocondria.Bcl2,
        AB=0,
        forward_rate=KF / mitocondria_volume_fraction,
        reverse_rate=KR,
    )
    r_Bax2_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax_M2,
        B=mitocondria.Bcl2,
        AB=0,
        forward_rate=KF / mitocondria_volume_fraction,
        reverse_rate=KR,
    )
    r_Bax4_Bcl2 = reactions.ReversibleSynthesis(
        A=Bax_M4,
        B=mitocondria.Bcl2,
        AB=0,
        forward_rate=KF / mitocondria_volume_fraction,
        reverse_rate=KR,
    )

    r_Bax_transloc = reactions.Equilibration(
        A=mitocondria.Bax_A,
        B=mitocondria.Bax_M,
        forward_rate=mitocondria.transloc_rate,
        reverse_rate=mitocondria.transloc_rate,
    )

    Mito_I: Species = initial(default=5e5)
    Mito_A: Species = initial(default=0)

    r_Bax4_Mito = reactions.CatalyzeConvert(
        A=Bax_M4,
        B=Mito_I,
        AB=0,
        P=Mito_A,
        forward_rate=KF / mitocondria_volume_fraction,
        reverse_rate=KR,
        conversion_rate=KC,
    )

    # Pore transport
    pore_transport_rate: Parameter = assign(default=10)
    r_Smac_pore = reactions.MichaelisMenten(
        E=Mito_A,
        S=mitocondria.Smac_M,
        ES=0,
        P=mitocondria.Smac_C,
        forward_rate=2 * KF / mitocondria_volume_fraction,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )
    r_CytoC_pore = reactions.MichaelisMenten(
        E=Mito_A,
        S=mitocondria.CytoC_M,
        ES=0,
        P=mitocondria.CytoC_C,
        forward_rate=2 * KF / mitocondria_volume_fraction,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )

    # Add citoplasmic Bcl2 as it was in Albeck's model because it's absent in
    # EARM implementation.
    Bcl2c: Species = initial(default=2e4)
    r_Bcl2_Bid = reactions.ReversibleSynthesis(
        A=Bid_T,
        B=Bcl2c,
        AB=0,
        forward_rate=KF,
        reverse_rate=KR,
    )
