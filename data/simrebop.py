import dataclasses
from collections.abc import Iterable, Mapping, Sequence

import pandas as pd
import xarray
from rebop import Gillespie
from rebop.gillespie import RNGLike, SeedLike
from simbio import Compartment, MassAction, RateLaw, Simulator, Species


class RebopSimulator:
    def __init__(
        self,
        model: type[Compartment],
        /,
    ):
        if not all(isinstance(r, MassAction) for r in model._yield(RateLaw)):
            # TODO: check non-reaction equations?
            raise NotImplementedError("only MassAction reactions are implemented")

        self.model = model
        self._sim = Simulator(model)
        self._sim.compiled = dataclasses.replace(self._sim.compiled, func=None)

        self._map = pd.DataFrame({"variable": self._sim.compiled.variables})
        self._map["name"] = self._map["variable"].map(str)
        self._map["rebop_name"] = self._map["name"].str.replace(".", "__")

        self._name_map = self._map.set_index("rebop_name")["name"].to_dict()
        self._variable_map = self._map.set_index("variable")["rebop_name"].to_dict()

    def _yield_species(self, species: Iterable[Species], /):
        for s in species:
            if not s.stoichiometry.is_integer():
                raise NotImplementedError

            name = self._variable_map[s.variable]
            for _ in range(int(s.stoichiometry)):
                yield name

    def _rename(self, df: xarray.Dataset, /):
        return df.rename({k: self._name_map[k] for k in df.keys()})

    def create_problem(
        self,
        values: Mapping = {},
        *,
        upto_t: float,
        n_points: int | None = None,
        sparse: bool = True,
        var_names: Sequence[Species] | None = None,
    ):
        problem = self._sim.create_problem(values)
        init = {k: int(v) for k, v in zip(self._variable_map.values(), problem.y)}
        parameters = dict(zip(self._sim.compiled.parameters, problem.p))
        reactions = [
            Reaction(
                rate=parameters[r.rate],
                reactants=list(self._yield_species(r.reactants)),
                products=list(self._yield_species(r.products)),
            )
            for r in self.model._yield(MassAction)
        ]
        return Problem(
            reactions,
            init,
            upto_t=upto_t,
            n_points=n_points if n_points is not None else 0,
            sparse=sparse,
            var_names=[self._variable_map[v.variable] for v in var_names]
            if var_names is not None
            else None,
        )

    def solve(
        self,
        values: Mapping = {},
        *,
        upto_t: float,
        n_points: int | None = None,
        sparse: bool = True,
        var_names: Sequence[Species] | None = None,
        rng: RNGLike | SeedLike | None = None,
    ):
        problem = self.create_problem(
            values=values,
            upto_t=upto_t,
            n_points=n_points,
            sparse=sparse,
            var_names=var_names,
        )
        df = problem.run(rng)
        return self._rename(df)

    def solve_many(
        self,
        values: Mapping = {},
        *,
        upto_t: float,
        n_points: int | None = None,
        sparse: bool = True,
        var_names: Sequence[Species] | None = None,
        seeds: Sequence[int],
        progress: bool = False,
    ):
        problem = self.create_problem(
            values=values,
            upto_t=upto_t,
            n_points=n_points,
            sparse=sparse,
            var_names=var_names,
        )

        if progress:
            from tqdm.contrib.concurrent import process_map

            solutions: list[xarray.Dataset] = process_map(problem.run, seeds)
        else:
            from concurrent.futures.process import ProcessPoolExecutor

            with ProcessPoolExecutor() as p:
                solutions = list(p.map(problem.run, seeds))
        df = xarray.concat(solutions, dim="seed").assign_coords({"seed": seeds})
        return self._rename(df)


@dataclasses.dataclass(frozen=True)
class Reaction:
    rate: float
    reactants: Sequence[str]
    products: Sequence[str]


@dataclasses.dataclass(frozen=True)
class Problem:
    reactions: Sequence[Reaction]
    init: Mapping[str, int]
    upto_t: float
    n_points: int
    sparse: bool = True
    var_names: Sequence[str] | None = None

    def run(self, rng: RNGLike | SeedLike | None = None):
        runner = Gillespie()
        for r in self.reactions:
            runner.add_reaction(r.rate, r.reactants, r.products)
        return runner.run(
            self.init,
            tmax=self.upto_t,
            nb_steps=self.n_points,
            rng=rng,
            sparse=self.sparse,
            var_names=self.var_names,
        )
