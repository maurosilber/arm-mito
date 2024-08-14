import numpy as np
from pandas.testing import assert_frame_equal
from simbio import Compartment, Simulator, Species, initial, reactions

from ..loop_simulator import LoopSimulator


def test_main_disconnected_from_loop():
    class Main(Compartment):
        x: Species = initial(default=0)
        r = reactions.Creation(A=x, rate=1)

    class Loop(Compartment):
        y: Species = initial(default=0)
        r = reactions.Creation(A=y, rate=1)

    sim = LoopSimulator(main=Main, loop=Loop, main_variables_in_loop=[])

    t = np.linspace(0, 1, 10)
    main_values = {Main.x.variable: 0.0}
    df = sim.solve(
        main_values=main_values,
        loop_values={Loop.y.variable: [1, 2]},
        save_at=t,
    )

    assert Main.x.variable in df
    assert Loop.y.variable not in df

    df_main = Simulator(Main).solve(values=main_values, save_at=t)
    assert_frame_equal(
        df.rename(columns=str).rename_axis(index="time"),
        df_main,
    )
