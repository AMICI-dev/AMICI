# ruff: noqa: F401, F821, F841
import amici
import equinox as eqx
import jax
import jax.random as jr


class net1(eqx.Module):
    layers: dict
    inputs: list[str]
    outputs: list[str]

    def __init__(self, key):
        super().__init__()
        keys = jr.split(key, 3)
        self.layers = {
            "layer1": eqx.nn.Linear(
                in_features=2, out_features=5, use_bias=True, key=keys[0]
            ),
            "layer2": eqx.nn.Linear(
                in_features=5, out_features=5, use_bias=True, key=keys[1]
            ),
            "layer3": eqx.nn.Linear(
                in_features=5, out_features=1, use_bias=True, key=keys[2]
            ),
        }
        self.inputs = ["input0"]
        self.outputs = ["layer3"]

    def forward(self, input, key=None):
        net_input = input
        layer1 = (
            jax.vmap(self.layers["layer1"])
            if len(net_input.shape) == 2
            else self.layers["layer1"]
        )(net_input)
        tanh = jax.nn.tanh(layer1)
        layer2 = (
            jax.vmap(self.layers["layer2"])
            if len(tanh.shape) == 2
            else self.layers["layer2"]
        )(tanh)
        tanh_1 = jax.nn.tanh(layer2)
        layer3 = (
            jax.vmap(self.layers["layer3"])
            if len(tanh_1.shape) == 2
            else self.layers["layer3"]
        )(tanh_1)
        output = layer3
        return output


net = net1
