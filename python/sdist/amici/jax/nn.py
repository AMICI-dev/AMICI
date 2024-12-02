from pathlib import Path

from petab_sciml import MLModel, Layer, Node
import equinox as eqx
import jax.numpy as jnp

from amici._codegen.template import apply_template
from amici import amiciModulePath


class Flatten(eqx.Module):
    start_dim: int
    end_dim: int

    def __init__(self, start_dim: int, end_dim: int):
        super().__init__()
        self.start_dim = start_dim
        self.end_dim = end_dim

    def __call__(self, x):
        if self.end_dim == -1:
            return jnp.reshape(x, x.shape[: self.start_dim] + (-1,))
        else:
            return jnp.reshape(
                x, x.shape[: self.start_dim] + (-1,) + x.shape[self.end_dim :]
            )


def generate_equinox(ml_model: MLModel, filename: Path | str):
    filename = Path(filename)
    layer_indent = 12
    node_indent = 8

    layers = {layer.layer_id: layer for layer in ml_model.layers}

    tpl_data = {
        "MODEL_ID": ml_model.mlmodel_id,
        "LAYERS": ",\n".join(
            [
                _generate_layer(layer, layer_indent, ilayer)
                for ilayer, layer in enumerate(ml_model.layers)
            ]
        )[layer_indent:],
        "FORWARD": "\n".join(
            [
                _generate_forward(
                    node,
                    node_indent,
                    layers.get(
                        node.target,
                        Layer(layer_id="dummy", layer_type="Linear"),
                    ).layer_type,
                )
                for node in ml_model.forward
            ]
        )[node_indent:],
        "INPUT": ", ".join([f"'{inp.input_id}'" for inp in ml_model.inputs]),
        "N_LAYERS": len(ml_model.layers),
    }

    filename.parent.mkdir(parents=True, exist_ok=True)

    apply_template(
        Path(amiciModulePath) / "jax" / "nn.template.py",
        filename,
        tpl_data,
    )


def _process_argval(v):
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, bool):
        return str(v)
    return str(v)


def _generate_layer(layer: Layer, indent: int, ilayer: int) -> str:
    layer_map = {
        "InstanceNorm1d": "eqx.nn.LayerNorm",
        "InstanceNorm2d": "eqx.nn.LayerNorm",
        "InstanceNorm3d": "eqx.nn.LayerNorm",
        "Dropout1d": "eqx.nn.Dropout",
        "Dropout2d": "eqx.nn.Dropout",
        "Flatten": "Flatten",
    }
    kwarg_map = {
        "Linear": {
            "bias": "use_bias",
        },
        "Conv1d": {
            "bias": "use_bias",
        },
        "Conv2d": {
            "bias": "use_bias",
        },
        "InstanceNorm1d": {
            "affine": "elementwise_affine",
            "num_features": "shape",
        },
        "InstanceNorm2d": {
            "affine": "elementwise_affine",
            "num_features": "shape",
        },
        "InstanceNorm3d": {
            "affine": "elementwise_affine",
            "num_features": "shape",
        },
    }
    kwarg_ignore = {
        "InstanceNorm1d": ("track_running_stats", "momentum"),
        "InstanceNorm2d": ("track_running_stats", "momentum"),
        "InstanceNorm3d": ("track_running_stats", "momentum"),
        "Dropout1d": ("inplace",),
        "Dropout2d": ("inplace",),
    }
    kwargs = [
        f"{kwarg_map.get(layer.layer_type, {}).get(k, k)}={_process_argval(v)}"
        for k, v in layer.args.items()
        if k not in kwarg_ignore.get(layer.layer_type, ())
    ]
    # add key for initialization
    if layer.layer_type in ("Linear", "Conv1d", "Conv2d", "Conv3d"):
        kwargs += [f"key=keys[{ilayer}]"]
    type_str = layer_map.get(layer.layer_type, f"eqx.nn.{layer.layer_type}")
    layer_str = f"{type_str}({', '.join(kwargs)})"
    if layer.layer_type.startswith(("InstanceNorm",)):
        if layer.layer_type.endswith(("1d", "2d", "3d")):
            layer_str = f"jax.vmap({layer_str}, in_axes=1, out_axes=1)"
        if layer.layer_type.endswith(("2d", "3d")):
            layer_str = f"jax.vmap({layer_str}, in_axes=2, out_axes=2)"
        if layer.layer_type.endswith("3d"):
            layer_str = f"jax.vmap({layer_str}, in_axes=3, out_axes=3)"
    return f"{' ' * indent}'{layer.layer_id}': {layer_str}"


def _generate_forward(node: Node, indent, layer_type=str) -> str:
    if node.op == "placeholder":
        # TODO: inconsistent target vs name
        return f"{' ' * indent}{node.name} = input"

    if node.op == "call_module":
        fun_str = f"self.layers['{node.target}']"
        if layer_type.startswith(("InstanceNorm", "Conv", "Linear")):
            if layer_type == "Linear":
                dims = 1
            if layer_type.endswith(("1d",)):
                dims = 2
            elif layer_type.endswith(("2d",)):
                dims = 3
            elif layer_type.endswith("3d"):
                dims = 4
            fun_str = f"(jax.vmap({fun_str}) if len({node.args[0]}.shape) == {dims + 1} else {fun_str})"

    if node.op in ("call_function", "call_method"):
        map_fun = {
            "hardtanh": "jax.nn.hard_tanh",
        }
        if node.target == "hardtanh":
            if node.kwargs.pop("min_val", -1.0) != -1.0:
                raise NotImplementedError(
                    "min_val != -1.0 not supported for hardtanh"
                )
            if node.kwargs.pop("max_val", 1.0) != 1.0:
                raise NotImplementedError(
                    "max_val != 1.0 not supported for hardtanh"
                )
        fun_str = map_fun.get(node.target, f"jax.nn.{node.target}")

    args = ", ".join([f"{arg}" for arg in node.args])
    kwargs = [
        f"{k}={v}" for k, v in node.kwargs.items() if k not in ("inplace",)
    ]
    if layer_type.startswith(("Dropout",)):
        kwargs += ["inference=inference", "key=key"]
    kwargs_str = ", ".join(kwargs)
    if node.op in ("call_module", "call_function", "call_method"):
        return f"{' ' * indent}{node.name} = {fun_str}({args + ', ' + kwargs_str})"
    if node.op == "output":
        return f"{' ' * indent}{node.target} = {args}"
