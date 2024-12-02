# ruff: noqa: F401, F821, F841
import equinox as eqx
import jax.nn
import jax.random as jr
import jax
from amici.jax.nn import Flatten


class TPL_MODEL_ID(eqx.Module):
    layers: dict
    inputs: list[str]

    def __init__(self, key):
        super().__init__()
        keys = jr.split(key, TPL_N_LAYERS)
        self.layers = {TPL_LAYERS}
        self.inputs = [TPL_INPUT]

    def forward(self, input, inference=False, key=None):
        TPL_FORWARD
        return output


net = TPL_MODEL_ID
