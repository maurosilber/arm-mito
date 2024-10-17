import functools
import itertools
from typing import Protocol, cast

import numpy as np
import rebop
import xarray
from numpy.typing import ArrayLike
from poincare import Variable
from simbio import Compartment, Constant, MassAction, Parameter, Simulator, Species
from symbolite.core import evaluate
from symbolite.impl import libstd

type VALUES = Species | Parameter | Constant


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


def to_rebop(model: type[Compartment], /, values: dict[VALUES, float] = {}):
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
    values_main: dict[VALUES, float] = {},
    values_loop: dict[VALUES, ArrayLike] = {},
):
    @functools.cache
    def is_loop_var(x: VALUES, /) -> bool:
        parent = x.parent
        while parent is not None and parent is not model:
            if parent is loop:
                return True
            parent = parent.parent
        return False

    def is_loop_reaction(reaction: MassAction, /) -> bool:
        components = itertools.chain(reaction.reactants, reaction.products)
        return any(map(is_loop_var, components))

    def variable_to_str(
        variable: Variable | Parameter | Constant, loop_index: int, /
    ) -> str:
        name = str(variable)
        if is_loop_var(variable):
            name = name.replace(loop.name, f"{loop.name}_{loop_index}")
        return name

    sim = Simulator(model)
    y: dict[str, int] = {}
    reactions = []

    values = values_main.copy()
    for loop_ix, loop_vals in enumerate(
        zip(*np.broadcast_arrays(*values_loop.values()))
    ):
        values.update(zip(values_loop.keys(), loop_vals))
        problem = sim.create_problem(values)

        y.update(
            zip(
                (variable_to_str(v, loop_ix) for v in sim.compiled.variables),
                map(int, problem.y),
            )
        )
        p = dict(zip(sim.compiled.parameters, problem.p))

        for r in model._yield(MassAction):
            if loop_ix > 0 and not is_loop_reaction(r):
                continue

            rate = evaluate(r.rate.subs(p), libsl=libstd)
            reactants = [
                variable_to_str(s.variable, loop_ix)
                for s in r.reactants
                for _ in range(s.stoichiometry)
            ]
            products = [
                variable_to_str(s.variable, loop_ix)
                for s in r.products
                for _ in range(s.stoichiometry)
            ]
            reactions.append((rate, reactants, products))

    return reactions, y


def create_rebop(reactions):
    runner = cast(Rebop, rebop.Gillespie())
    for r in sorted(reactions):
        runner.add_reaction(*r)
    return runner
