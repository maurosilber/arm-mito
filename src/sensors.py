from poincare import assign
from simbio import Compartment, Constant, Parameter, Species, initial, reactions


class Sensor(Compartment):
    dimer: Species = initial()
    monomer: Species = initial(default=0)


class CASPAM(Compartment):
    concentration: Constant = assign(constant=True)
    sCas3: Sensor = Sensor(dimer=concentration)
    sCas8: Sensor = Sensor(dimer=concentration)
    sCas9: Sensor = Sensor(dimer=concentration)


class SensorReaction(Compartment):
    sensor: Sensor = Sensor(dimer=0)
    enzyme: Species = initial()
    forward_rate: Parameter = assign()
    reverse_rate: Parameter = assign()
    catalytic_rate: Parameter = assign()

    r = reactions.MichaelisMenten(
        E=enzyme,
        S=sensor.dimer,
        ES=0,
        P=2 * sensor.monomer,
        forward_rate=forward_rate,
        reverse_rate=reverse_rate,
        catalytic_rate=catalytic_rate,
    )
