from typing import Protocol, cast

import rebop
import xarray
from simbio import Compartment, MassAction, Simulator
from symbolite.core import evaluate
from symbolite.impl import libstd


class Rebop(Protocol):
    def add_reaction(self, rate: float, reactants: list[str], products: list[str]): ...
    def run(
        self,
        values: dict[str, int],
        /,
        tmax: float,
        nb_steps: int,
        seed: int | None = None,
    ) -> xarray.Dataset: ...


def to_rebop(model: type[Compartment], /, values: dict = {}):
    sim = Simulator(model)
    problem = sim.create_problem(values)
    p = dict(zip(sim.compiled.parameters, problem.p))
    y = {str(k): int(v) for k, v in zip(sim.compiled.variables, problem.y)}
    reactions = []
    for r in model._yield(MassAction):
        rate = evaluate(r.rate.subs(p), libsl=libstd)
        reactants = [
            str(s.variable) for s in r.reactants for _ in range(s.stoichiometry)
        ]
        products = [str(s.variable) for s in r.products for _ in range(s.stoichiometry)]
        reactions.append((rate, reactants, products))

    return reactions, y


def to_rebop_loopy(
    model: type[Compartment],
    loop: Compartment,
    /,
    N: int,
    values: dict = {},
):
    sim = Simulator(model)
    problem = sim.create_problem(values)
    p = dict(zip(sim.compiled.parameters, problem.p))
    y = {}
    for k, v in zip(sim.compiled.variables, problem.y):
        k = str(k)
        if k.startswith(loop.name):
            v = int(v / N)
            for i in range(N):
                y[k.replace(loop.name, f"{loop.name}_{i}")] = v
        else:
            y[k] = int(v)

    reactions = []
    for r in model._yield(MassAction):
        rate = evaluate(r.rate.subs(p), libsl=libstd)
        reactants = [
            str(s.variable) for s in r.reactants for _ in range(s.stoichiometry)
        ]
        products = [str(s.variable) for s in r.products for _ in range(s.stoichiometry)]
        if any(r.startswith(loop.name) for r in [*reactants, *products]):
            for i in range(N):
                reactions.append(
                    (
                        rate,
                        [r.replace(loop.name, f"{loop.name}_{i}") for r in reactants],
                        [r.replace(loop.name, f"{loop.name}_{i}") for r in products],
                    )
                )
        else:
            reactions.append((rate, reactants, products))

    return reactions, y


def create_rebop(reactions):
    runner = cast(Rebop, rebop.Gillespie())
    for r in sorted(reactions):
        runner.add_reaction(*r)
    return runner
