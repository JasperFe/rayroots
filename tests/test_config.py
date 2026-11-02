"""The yaml defaults must equal the literal values in the frozen modules, and apply() must leave them unchanged."""
import numpy as np

from rayroots.plot.roots import config, connect as cn, cuts, flip, ground, v2


def snapshot():
    out = {}
    for m in (cn, cuts, flip, ground, v2):
        for k, v in vars(m).items():
            if k.isupper() and isinstance(v, (int, float, tuple, dict, np.ndarray)):
                out[(m.__name__, k)] = v.tolist() if isinstance(v, np.ndarray) else v
    return out


def test_defaults_equal_frozen_literals():
    before = snapshot()
    config.apply(config.load())
    after = snapshot()
    assert before == after


def test_frozen_flag():
    cfg = config.load()
    assert config.is_frozen(cfg)
    cfg["h_cut"] = 3.0
    assert not config.is_frozen(cfg)


def test_colour_roundtrip():
    from rayroots.plot.roots import plotdata as cm
    ids = np.array([-1, 0, 1, 262, 70000])
    r, g, b = cm.encode_ids(ids)
    assert (ground.decode_ids(r, g, b) == ids).all()
