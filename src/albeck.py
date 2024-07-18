from simbio import Compartment, Parameter, Species, assign, initial, reactions


class Receptor(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)

    L: Species = initial(default=3_000)
    R: Species = initial(default=200)
    DISC: Species = initial(default=0)
    flip: Species = initial(default=100)
    C8_pro: Species = initial(default=20_000)
    C8_A: Species = initial(default=0)
    BAR: Species = initial(default=1_000)

    Bid_U: Species = initial(default=1e5)
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
        forward_rate=KF,
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


class Intrinsic(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)

    CytoC_A: Species = initial()
    Smac_A: Species = initial()
    C8_pro: Species = initial()
    C8_A: Species = initial()

    Apaf_I: Species = initial(default=1e5)
    Apaf_A: Species = initial(default=0)
    Apop: Species = initial(default=0)  # Apoptosome (activated Apaf-1 + caspase 9)
    C3_pro: Species = initial(default=1e4)
    C3_A: Species = initial(default=0)
    C3_ub: Species = initial(default=0)
    C6_pro: Species = initial(default=1e4)
    C6_A: Species = initial(default=0)
    C9: Species = initial(default=1e5)
    PARP_U: Species = initial(default=1e6)  # Uncleaved
    PARP_C: Species = initial(default=0)  # Cleaved
    XIAP: Species = initial(default=1e5)  # # X-linked Inhibitor of Apoptosis Protein

    # Apoptosome formation
    # --------------------
    #   Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
    #   aApaf + pC9 <-->  Apop
    #   Apop + pC3 <-->  Apop:pC3 --> Apop + C3

    r_CytoC_Apaf = reactions.MichaelisMenten(
        E=CytoC_A,
        S=Apaf_I,
        ES=0,
        P=Apaf_A,
        forward_rate=5e-7,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    r_Apaf_C9 = reactions.ReversibleSynthesis(
        A=Apaf_A,
        B=C9,
        AB=Apop,
        forward_rate=5e-8,
        reverse_rate=KR,
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

    # Apoptosome-related inhibitors
    # -----------------------------
    #   Apop + XIAP <-->  Apop:XIAP
    #   cSmac + XIAP <-->  cSmac:XIAP

    r_Apop_XIAP = reactions.ReversibleSynthesis(
        A=Apop,
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

    # Caspase reactions
    # -----------------
    # Includes effectors, inhibitors, and feedback initiators:
    #
    #   pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
    #   pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
    #   XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
    #   PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
    #   pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS

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


class MOMP(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)

    # Activators
    Bid_U: Species = initial(default=1e5)  # Untruncated
    Bid_T: Species = initial(default=0)  # Truncated
    Bid_M: Species = initial(default=0)  # truncated andMitocondrial

    # Effectors
    Bax_C: Species = initial(default=1e5)  # Cytoplasmic
    Bax_M: Species = initial(default=0)  # Mitochondrial
    Bax_A: Species = initial(default=0)  # Active

    # Anti-Apoptotics
    Bcl2: Species = initial(default=2e4)

    # Cytochrome C and Smac
    CytoC_M: Species = initial(default=1e6)
    CytoC_C: Species = initial(default=0)
    CytoC_A: Species = initial(default=0)
    Smac_M: Species = initial(default=1e6)
    Smac_C: Species = initial(default=0)
    Smac_A: Species = initial(default=0)

    # CytoC and Smac activation after release
    transloc_rate: Parameter = assign(default=1e-2)
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


class Albeck11b(Compartment):
    KF: Parameter = assign(default=1e-6)
    KR: Parameter = assign(default=1e-3)
    KC: Parameter = assign(default=1)
    receptor = Receptor(KF=KF, KR=KR, KC=KC)
    mitocondria = MOMP(
        KF=KF,
        KR=KR,
        KC=KC,
        Bid_U=receptor.Bid_U,
        Bid_T=receptor.Bid_T,
    )
    intrinsic = Intrinsic(
        KF=KF,
        KR=KR,
        KC=KC,
        CytoC_A=mitocondria.CytoC_A,
        Smac_A=mitocondria.Smac_A,
        C8_pro=receptor.C8_pro,
        C8_A=receptor.C8_A,
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
    r_Bax_Bcl2 = reactions.ReversibleSynthesis(
        A=mitocondria.Bax_A,
        B=mitocondria.Bcl2,
        AB=0,
        forward_rate=KF,
        reverse_rate=KR,
    )
    # Pore transport
    pore_transport_rate: Parameter = assign(default=10)
    r_Smac_pore = reactions.MichaelisMenten(
        E=mitocondria.Bax_A,
        S=mitocondria.Smac_M,
        ES=0,
        P=mitocondria.Smac_C,
        forward_rate=KF,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )
    r_CytoC_pore = reactions.MichaelisMenten(
        E=mitocondria.Bax_A,
        S=mitocondria.CytoC_M,
        ES=0,
        P=mitocondria.CytoC_C,
        forward_rate=KF,
        reverse_rate=KR,
        catalytic_rate=pore_transport_rate,
    )
