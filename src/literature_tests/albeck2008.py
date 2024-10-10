"""
  Title: Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic Cell Death
Authors: John G Albeck,  John M Burke,  Sabrina L Spencer,  Douglas A Lauffenburger,  Peter K Sorger
    DOI: https://doi.org/10.1371/journal.pbio.0060299
"""

import pandas as pd
import pint_pandas


def test_onset():
    """Delay until apoptosis onset.

    Defined as the half-maximum of cleavage.

    Table 1.
    """
    data = pd.DataFrame(
        {
            "drug": ["TRAIL", "TRAIL", "TRAIL", "TRAIL", "TRAIL", "TNF"],
            "treatment": pd.Series([1000, 250, 50, 10, 2, 100], dtype="pint[ng / ml]"),
            "time_mean": pd.Series([140, 180, 240, 360, 660, 460], dtype="pint[min]"),
            "time_std": pd.Series([32, 32, 36, 79, 170, 190], dtype="pint[min]"),
        }
    )


def test_switch_width():
    """Width of the cleavage curve.

    Defined as:

        c(t) \propto 1 - 1 / (1 + e^[(t - t_d) / (4 t_s)])

    Table 1.
    """
    data = pd.DataFrame(
        {
            "drug": ["TRAIL", "TRAIL", "TRAIL", "TRAIL", "TRAIL", "TNF"],
            "treatment": pd.Series([1000, 250, 50, 10, 2, 100], dtype="pint[ng / ml]"),
            "time_mean": pd.Series([22, 24, 27, 22, 19, 22], dtype="pint[min]"),
            "time_std": pd.Series([9.5, 9.5, 13, 7.7, 10, 27], dtype="pint[min]"),
        }
    )
