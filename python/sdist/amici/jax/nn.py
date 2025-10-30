from pathlib import Path

import equinox as eqx
import jax.numpy as jnp

from amici import amiciModulePath
from amici._codegen.template import apply_template


class Flatten(eqx.Module):
    """Custom implementation of a torch.flatten layer for Equinox."""

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


def tanhshrink(x: jnp.ndarray) -> jnp.ndarray:
    """Custom implementation of the torch.nn.Tanhshrink activation function for JAX."""
    return x - jnp.tanh(x)


def generate_equinox(
    nn_model: "NNModel",  # noqa: F821
    filename: Path | str,
    frozen_layers: dict[str, bool] | None = None,
) -> None:
    """
    Generate Equinox model file from petab_sciml neural network object.

    :param nn_model:
        Neural network model in petab_sciml format
    :param filename:
        output filename for generated Equinox model
    :param frozen_layers:
        list of layer names to freeze during training
    :return:

    """
    # TODO: move to top level import and replace forward type definitions
    from petab_sciml import Layer

    if frozen_layers is None:
        frozen_layers = {}

    filename = Path(filename)
    layer_indent = 12
    node_indent = 8

    layers = {layer.layer_id: layer for layer in nn_model.layers}

    tpl_data = {
        "MODEL_ID": nn_model.nn_model_id,
        "LAYERS": ",\n".join(
            [
                _generate_layer(layer, layer_indent, ilayer)
                for ilayer, layer in enumerate(nn_model.layers)
            ]
        )[layer_indent:],
        "FORWARD": "\n".join(
            [
                _generate_forward(
                    node,
                    node_indent,
                    frozen_layers,
                    layers.get(
                        node.target,
                        Layer(layer_id="dummy", layer_type="Linear"),
                    ).layer_type,
                )
                for node in nn_model.forward
            ]
        )[node_indent:],
        "INPUT": ", ".join([f"'{inp.input_id}'" for inp in nn_model.inputs]),
        "OUTPUT": ", ".join(
            [
                f"'{arg}'"
                for arg in next(
                    node for node in nn_model.forward if node.op == "output"
                ).args
            ]
        ),
        "N_LAYERS": len(nn_model.layers),
    }

    filename.parent.mkdir(parents=True, exist_ok=True)

    apply_template(
        Path(amiciModulePath) / "jax" / "nn.template.py",
        filename,
        tpl_data,
    )


def _process_argval(v):
    """
    Process argument value for layer instantiation string
    """
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, bool):
        return str(v)
    return str(v)


def _generate_layer(layer: "Layer", indent: int, ilayer: int) -> str:  # noqa: F821
    """
    Generate layer definition string for a given layer

    :param layer:
        petab_sciml Layer object
    :param indent:
        indentation level for generated string
    :param ilayer:
        layer index for key generation

    :return:
        string defining the layer in equinox syntax
    """
    if layer.layer_type.startswith(
        ("BatchNorm", "AlphaDropout", "InstanceNorm")
    ):
        raise NotImplementedError(
            f"{layer.layer_type} layers currently not supported"
        )
    if layer.layer_type.startswith("MaxPool") and "dilation" in layer.args:
        raise NotImplementedError("MaxPool layers with dilation not supported")
    if layer.layer_type.startswith("Dropout") and "inplace" in layer.args:
        raise NotImplementedError("Dropout layers with inplace not supported")
    if layer.layer_type == "Bilinear":
        raise NotImplementedError("Bilinear layers not supported")

    # mapping of layer names in sciml yaml format to equinox/custom amici implementations
    layer_map = {
        "Dropout1d": "eqx.nn.Dropout",
        "Dropout2d": "eqx.nn.Dropout",
        "Flatten": "amici.jax.Flatten",
    }

    # mapping of keyword argument names in sciml yaml format to equinox/custom amici implementations
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
        "LayerNorm": {
            "elementwise_affine": "use_bias",  # Deprecation warning - replace LayerNorm(elementwise_affine) with LayerNorm(use_bias)
            "normalized_shape": "shape",
        },
    }
    # list of keyword arguments to ignore when generating layer, as they are not supported in equinox (see above)
    kwarg_ignore = {
        "Dropout1d": ("inplace",),
        "Dropout2d": ("inplace",),
    }
    # construct argument string for layer instantiation
    kwargs = [
        f"{kwarg_map.get(layer.layer_type, {}).get(k, k)}={_process_argval(v)}"
        for k, v in layer.args.items()
        if k not in kwarg_ignore.get(layer.layer_type, ())
    ]
    # add key for initialization
    if layer.layer_type in (
        "Linear",
        "Conv1d",
        "Conv2d",
        "Conv3d",
        "ConvTranspose1d",
        "ConvTranspose2d",
        "ConvTranspose3d",
    ):
        kwargs += [f"key=keys[{ilayer}]"]
    type_str = layer_map.get(layer.layer_type, f"eqx.nn.{layer.layer_type}")
    layer_str = f"{type_str}({', '.join(kwargs)})"
    return f"{' ' * indent}'{layer.layer_id}': {layer_str}"


def _format_function_call(
    var_name: str, fun_str: str, args: list, kwargs: list[str], indent: int
) -> str:
    """
    Utility function to format a function call assignment string.

    :param var_name:
        name of the variable to assign the result to
    :param fun_str:
        string representation of the function to call
    :param args:
        list of positional arguments
    :param kwargs:
        list of keyword arguments as strings
    :param indent:
        indentation level for generated string

    :return:
        formatted string representing the function call assignment
    """
    args_str = ", ".join([f"{arg}" for arg in args])
    kwargs_str = ", ".join(kwargs)
    all_args = ", ".join(filter(None, [args_str, kwargs_str]))
    return f"{' ' * indent}{var_name} = {fun_str}({all_args})"


def _process_layer_call(
    node: "Node",  # noqa: F821
    layer_type: str,
    frozen_layers: dict[str, bool],
) -> tuple[str, str]:
    """
    Process a layer (call_module) node and return function string and optional tree string.

    :param node:
        petab sciml Node object representing a layer call
    :param layer_type:
        petab sciml layer type of the node
    :param frozen_layers:
        dict of layer names to boolean indicating whether layer is frozen

    :return:
        tuple of (function_string, tree_string) where tree_string is empty if no tree is needed
    """
    fun_str = f"self.layers['{node.target}']"
    tree_string = ""

    # Handle frozen layers
    if node.name in frozen_layers:
        if frozen_layers[node.name]:
            arr_attr = frozen_layers[node.name]
            get_lambda = f"lambda layer: getattr(layer, '{arr_attr}')"
            replacer = "replace_fn = lambda arr: jax.lax.stop_gradient(arr)"
            tree_string = f"tree_{node.name} = eqx.tree_at({get_lambda}, {fun_str}, {replacer})"
            fun_str = f"tree_{node.name}"
        else:
            fun_str = f"jax.lax.stop_gradient({fun_str})"

    # Handle vmap for certain layer types
    if layer_type.startswith(("Conv", "Linear", "LayerNorm")):
        if layer_type in ("LayerNorm",):
            dims = f"len({fun_str}.shape)+1"
        elif layer_type == "Linear":
            dims = 2
        elif layer_type.endswith("1d"):
            dims = 3
        elif layer_type.endswith("2d"):
            dims = 4
        elif layer_type.endswith("3d"):
            dims = 5
        fun_str = f"(jax.vmap({fun_str}) if len({node.args[0]}.shape) == {dims} else {fun_str})"

    return fun_str, tree_string


def _process_activation_call(node: "Node") -> str:  # noqa: F821
    """
    Process an activation function (call_function/call_method) node and return function string.

    :param node:
        petab sciml Node object representing an activation function call

    :return:
        string representation of the activation function
    """
    # Mapping of function names in sciml yaml format to equinox/custom amici implementations
    activation_map = {
        "hardtanh": "jax.nn.hard_tanh",
        "hardsigmoid": "jax.nn.hard_sigmoid",
        "hardswish": "jax.nn.hard_swish",
        "tanhshrink": "amici.jax.tanhshrink",
        "softsign": "jax.nn.soft_sign",
    }

    # Validate hardtanh parameters
    if node.target == "hardtanh":
        if node.kwargs.pop("min_val", -1.0) != -1.0:
            raise NotImplementedError(
                "min_val != -1.0 not supported for hardtanh"
            )
        if node.kwargs.pop("max_val", 1.0) != 1.0:
            raise NotImplementedError(
                "max_val != 1.0 not supported for hardtanh"
            )

    return activation_map.get(node.target, f"jax.nn.{node.target}")


def _generate_forward(
    node: "Node",  # noqa: F821
    indent,
    frozen_layers: dict[str, bool] | None = None,
    layer_type: str = "",
) -> str:
    """
    Generate forward pass line for a given node

    :param node:
        petab sciml Node object representing a step in the forward pass
    :param indent:
        indentation level for generated string
    :param frozen_layers:
        dict of layer names to boolean indicating whether layer is frozen
    :param layer_type:
        petab sciml layer type of the node (only relevant for call_module nodes)

    :return:
        string defining the forward pass implementation for the given node in equinox syntax
    """
    if frozen_layers is None:
        frozen_layers = {}

    # Handle placeholder nodes
    if node.op == "placeholder":
        # TODO: inconsistent target vs name
        return f"{' ' * indent}{node.name} = input"

    # Handle output nodes
    if node.op == "output":
        args_str = ", ".join([f"{arg}" for arg in node.args])
        return f"{' ' * indent}{node.target} = {args_str}"

    # Process layer calls
    tree_string = ""
    if node.op == "call_module":
        fun_str, tree_string = _process_layer_call(
            node, layer_type, frozen_layers
        )

    # Process activation function calls
    if node.op in ("call_function", "call_method"):
        fun_str = _process_activation_call(node)

    # Build kwargs list, filtering out unsupported arguments
    kwargs = [
        f"{k}={item}"
        for k, item in node.kwargs.items()
        if k not in ("inplace",)
    ]

    # Add key parameter for Dropout layers
    if layer_type.startswith("Dropout"):
        kwargs += ["key=key"]

    # Format the function call
    if node.op in ("call_module", "call_function", "call_method"):
        result = _format_function_call(
            node.name, fun_str, node.args, kwargs, indent
        )
        # Prepend tree_string if needed for frozen layers
        if tree_string:
            return f"{' ' * indent}{tree_string}\n{result}"
        return result

    raise NotImplementedError(f"Operation {node.op} not supported")
