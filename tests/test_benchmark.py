"""The Monte Carlo benchmark must run with any multiprocessing start method."""

import importlib
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.mark.parametrize("method", ["spawn", "forkserver"])
def test_benchmark_workers_get_their_state(method):
    if method not in mp.get_all_start_methods():
        pytest.skip(f"{method} is not available on this platform")
    sys.path.insert(0, str(EXAMPLES))
    try:
        bench = importlib.import_module("07_benchmark_interferometry")
        # a record without a wake: both methods must report "no fix", not crash
        state = {"ferry": {"flat": (None, np.zeros((9, bench.F_TM.size)))}}
        with mp.get_context(method).Pool(1, initializer=bench._init_worker, initargs=(state,)) as pool:
            rows = pool.map(bench._dispatch, [(bench.run_ferry, ("flat", 0.0, 0))])[0]
    finally:
        sys.path.remove(str(EXAMPLES))
    assert [r[2] for r in rows] == ["kelvin", "interferometry"]
    for r in rows:
        assert r[5] == 0 and r[9].startswith("no fix"), r
