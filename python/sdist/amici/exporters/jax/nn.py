from pathlib import Path
from shutil import copyfile

import equinox as eqx
import jax
import jax.numpy as jnp

from amici import amiciModulePath
from amici.exporters.template import apply_template


class BatchNorm(eqx.Module):
    """Custom implementation of PyTorch BatchNorm1d/2d/3d for Equinox.

    Computes batch normalisation using statistics computed from the current
    input batch.

    **Note**: Unlike PyTorch's eval-mode BatchNorm, running statistics are
    **not** supported.

    Parameters
    ----------
    num_features : int
        Number of features/channels ``C``.
    eps : float, optional
        Value added to the denominator for numerical stability.
        Default: ``1e-5``.
    affine : bool, optional
        If ``True``, learnable ``weight`` and ``bias`` parameters are
        included.  Default: ``True``.
    bias : bool, optional
        If ``True`` and ``affine=True``, a learnable bias term is included.
        Default: ``True``.
    """

    weight: jnp.ndarray | None
    bias: jnp.ndarray | None
    eps: float
    num_features: int

    def __init__(
        self,
        num_features: int,
        eps: float = 1e-5,
        affine: bool = True,
        bias: bool = True,
    ):
        self.num_features = num_features
        self.eps = eps
        self.weight = jnp.ones(num_features) if affine else None
        self.bias = jnp.zeros(num_features) if (affine and bias) else None

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        """Apply batch normalisation to ``x``.

        Parameters
        ----------
        x : jnp.ndarray
            Input array of shape ``(N, C)``, ``(N, C, L)``,
            ``(N, C, H, W)``, or ``(N, C, D, H, W)``.

        Returns
        -------
        jnp.ndarray
            Normalised array of the same shape as ``x``.
        """
        n_spatial = x.ndim - 2  # 0 for (N, C), >=1 for spatial inputs
        reduce_axes = (0,) + tuple(range(2, x.ndim)) if n_spatial > 0 else (0,)
        mean = jnp.mean(x, axis=reduce_axes, keepdims=True)
        var = jnp.var(x, axis=reduce_axes, keepdims=True)

        x_norm = (x - mean) / jnp.sqrt(var + self.eps)

        if self.weight is not None:
            w = (
                self.weight.reshape((self.num_features,) + (1,) * n_spatial)
                if n_spatial > 0
                else self.weight
            )
            x_norm = x_norm * w

        if self.bias is not None:
            b = (
                self.bias.reshape((self.num_features,) + (1,) * n_spatial)
                if n_spatial > 0
                else self.bias
            )
            x_norm = x_norm + b

        return x_norm


class InstanceNorm(eqx.Module):
    """Custom implementation of PyTorch InstanceNorm1d/2d/3d for Equinox.

    Applies Instance Normalisation over a batched input with optional learnable
    affine parameters.  Statistics are computed per instance per channel.

    **Note**: Unlike PyTorch's eval-mode InstanceNorm, running statistics are
    **not** supported.

    Parameters
    ----------
    num_features : int
        Number of features/channels ``C``.
    eps : float, optional
        Value added to the denominator for numerical stability.
        Default: ``1e-5``.
    affine : bool, optional
        If ``True``, learnable ``weight`` and ``bias`` parameters are included,
        initialised to ones and zeros respectively.
        Default: ``False``.
    bias : bool, optional
        If ``True`` and ``affine=True``, a learnable bias term is included.
        Default: ``True``.
    """

    weight: jnp.ndarray | None
    bias: jnp.ndarray | None
    eps: float
    num_features: int

    def __init__(
        self,
        num_features: int,
        eps: float = 1e-5,
        affine: bool = False,
        bias: bool = True,
    ):
        self.num_features = num_features
        self.eps = eps
        self.weight = jnp.ones(num_features) if affine else None
        self.bias = jnp.zeros(num_features) if (affine and bias) else None

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        """Apply instance normalisation to ``x``.

        Parameters
        ----------
        x : jnp.ndarray
            Batched input of shape ``(N, C, *)`` or unbatched ``(C, *)``.

        Returns
        -------
        jnp.ndarray
            Normalised array of the same shape as ``x``.
        """
        is_batched = x.ndim >= 3

        if is_batched:
            n_spatial = x.ndim - 2
            reduce_axes = tuple(range(2, x.ndim))
            mean = jnp.mean(x, axis=reduce_axes, keepdims=True)
            var = jnp.var(x, axis=reduce_axes, keepdims=True)
            x_norm = (x - mean) / jnp.sqrt(var + self.eps)

            if self.weight is not None:
                w = self.weight.reshape(
                    (1, self.num_features) + (1,) * n_spatial
                )
                x_norm = x_norm * w
            if self.bias is not None:
                b = self.bias.reshape(
                    (1, self.num_features) + (1,) * n_spatial
                )
                x_norm = x_norm + b
        else:
            n_spatial = x.ndim - 1
            reduce_axes = tuple(range(1, x.ndim))
            mean = jnp.mean(x, axis=reduce_axes, keepdims=True)
            var = jnp.var(x, axis=reduce_axes, keepdims=True)
            x_norm = (x - mean) / jnp.sqrt(var + self.eps)

            if self.weight is not None:
                w = self.weight.reshape(
                    (self.num_features,) + (1,) * n_spatial
                )
                x_norm = x_norm * w
            if self.bias is not None:
                b = self.bias.reshape((self.num_features,) + (1,) * n_spatial)
                x_norm = x_norm + b

        return x_norm


class AlphaDropout(eqx.Module):
    """Custom implementation of PyTorch AlphaDropout for Equinox.

    Applies Alpha Dropout over the input, maintaining the self-normalizing
    property for inputs with zero mean and unit standard deviation.

    During training, randomly masks elements of the input with probability
    ``p``, replacing dropped elements with the saturated negative SELU value
    ``alpha' = -scale * alpha ≈ -1.7581``, then applies an affine
    transformation to restore zero mean and unit standard deviation.

    **Note**: ``inplace`` mode is not supported.

    Parameters
    ----------
    p : float, optional
        Probability of an element being dropped. Default: ``0.5``.
    inference : bool, optional
        If ``True``, acts as an identity function (eval mode).
    """

    p: float
    inference: bool

    def __init__(self, p: float = 0.5, inference: bool = False):
        self.p = p
        self.inference = inference

    def __call__(self, x: jnp.ndarray, *, key=None) -> jnp.ndarray:
        """Apply Alpha Dropout to ``x``.

        Parameters
        ----------
        x : jnp.ndarray
            Input array of any shape.
        key : jax.random.PRNGKey, optional
            Random key required during training. Ignored in inference mode.

        Returns
        -------
        jnp.ndarray
            Output array of the same shape as ``x``.
        """
        if self.inference or self.p == 0.0:
            return x

        if key is None:
            raise RuntimeError(
                "Dropout requires a key when running in non-deterministic mode."
            )

        # alpha' = -scale * alpha (saturated negative SELU value)
        alpha_prime = -1.0507009873554805 * 1.6732632423543772
        keep_prob = 1.0 - self.p

        # Bernoulli mask: True = keep original, False = replace with alpha'
        mask = jax.random.bernoulli(key, keep_prob, x.shape)
        x_dropped = jnp.where(mask, x, alpha_prime)

        # Affine transform to restore zero mean and unit variance:
        #   E[x_dropped]   = p * alpha'
        #   Var[x_dropped] = keep_prob * (1 + p * alpha'^2)
        a = (keep_prob * (1.0 + self.p * alpha_prime**2)) ** (-0.5)
        b = -a * self.p * alpha_prime

        return a * x_dropped + b


class Bilinear(eqx.Module):
    """Custom implementation of PyTorch Bilinear for Equinox.

    Applies a bilinear transformation to the incoming data:
    ``y = x1^T A x2 + b``.

    Parameters
    ----------
    in1_features : int
        Size of each first input sample.
    in2_features : int
        Size of each second input sample.
    out_features : int
        Size of each output sample.
    bias : bool, optional
        If ``False``, the layer will not learn an additive bias.
        Default: ``True``.
    key : jax.random.PRNGKey
        Random key used to initialise ``weight`` and ``bias``.
    """

    weight: jnp.ndarray
    bias: jnp.ndarray | None
    in1_features: int
    in2_features: int
    out_features: int

    def __init__(
        self,
        in1_features: int,
        in2_features: int,
        out_features: int,
        bias: bool = True,
        *,
        key: "jax.random.PRNGKey",
    ):
        import math

        self.in1_features = in1_features
        self.in2_features = in2_features
        self.out_features = out_features
        k = 1.0 / math.sqrt(in1_features)
        w_key, b_key = jax.random.split(key)
        self.weight = jax.random.uniform(
            w_key,
            (out_features, in1_features, in2_features),
            minval=-k,
            maxval=k,
        )
        self.bias = (
            jax.random.uniform(b_key, (out_features,), minval=-k, maxval=k)
            if bias
            else None
        )

    def __call__(self, x1: jnp.ndarray, x2: jnp.ndarray) -> jnp.ndarray:
        """Apply bilinear transformation: ``y = x1^T A x2 + b``.

        Parameters
        ----------
        x1 : jnp.ndarray
            First input of shape ``(*, in1_features)``.
        x2 : jnp.ndarray
            Second input of shape ``(*, in2_features)``.

        Returns
        -------
        jnp.ndarray
            Output of shape ``(*, out_features)``.
        """
        # y_o = sum_i sum_j x1_i * W_oij * x2_j  (per batch element)
        out = jnp.einsum("...i,oij,...j->...o", x1, self.weight, x2)
        if self.bias is not None:
            out = out + self.bias
        return out


class Flatten(eqx.Module):
    """Custom implementation of a `torch.flatten` layer for Equinox."""

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


def cat(tensors, axis: int = 0):
    """Alias for torch.cat using JAX's concatenate/stack function.

    Handles both regular arrays and zero-dimensional (scalar) arrays by
    using stack instead of concatenate for 0D arrays.

    :param tensors:
        List of arrays to concatenate
    :param axis:
        Dimension along which to concatenate (default: 0)

    :return:
        Concatenated array
    """
    # Check if all tensors are 0-dimensional (scalars)
    if all(jnp.ndim(t) == 0 for t in tensors):
        # For 0D arrays, use stack instead of concatenate
        return jnp.stack(tensors, axis=axis)
    return jnp.concatenate(tensors, axis=axis)


def generate_equinox(
    nn_model: "NNModel | str",  # noqa: F821
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
    """
    # TODO: move to top level import and replace forward type definitions
    from petab_sciml import Layer

    if isinstance(nn_model, Path):
        copyfile(nn_model, filename)
        return

    if frozen_layers is None:
        frozen_layers = {}

    filename = Path(filename)
    layer_indent = 12
    node_indent = 8

    layers = {layer.layer_id: layer for layer in nn_model.layers}

    # Collect placeholder nodes to determine input handling
    placeholder_nodes = [
        node for node in nn_model.forward if node.op == "placeholder"
    ]
    input_names = [node.name for node in placeholder_nodes]

    # Generate input unpacking line
    if len(input_names) == 1:
        input_unpack = f"{input_names[0]} = input"
    else:
        input_unpack = f"{', '.join(input_names)} = input"

    # Generate forward pass lines (excluding placeholder nodes)
    forward_lines = [
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
    # Filter out empty lines from placeholder processing
    forward_lines = [line for line in forward_lines if line]
    # Prepend input unpacking
    forward_code = f"{' ' * node_indent}{input_unpack}\n" + "\n".join(
        forward_lines
    )

    tpl_data = {
        "MODEL_ID": nn_model.nn_model_id,
        "LAYERS": ",\n".join(
            [
                _generate_layer(layer, layer_indent, ilayer)
                for ilayer, layer in enumerate(nn_model.layers)
            ]
        )[layer_indent:],
        "FORWARD": forward_code[node_indent:],
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
        Path(amiciModulePath) / "exporters" / "jax" / "nn.template.py",
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
    if (
        layer.layer_type.startswith("MaxPool")
        and layer.args["dilation"]
        and layer.args["dilation"]
    ):
        raise NotImplementedError("MaxPool layers with dilation not supported")
    if (
        layer.layer_type.startswith("Dropout")
        and "inplace" in layer.args
        and layer.args["inplace"]
    ):
        raise NotImplementedError("Dropout layers with inplace not supported")
    # mapping of layer names in sciml yaml format to equinox/custom amici implementations
    layer_map = {
        "BatchNorm1d": "amici.exporters.jax.BatchNorm",
        "BatchNorm2d": "amici.exporters.jax.BatchNorm",
        "BatchNorm3d": "amici.exporters.jax.BatchNorm",
        "InstanceNorm1d": "amici.exporters.jax.InstanceNorm",
        "InstanceNorm2d": "amici.exporters.jax.InstanceNorm",
        "InstanceNorm3d": "amici.exporters.jax.InstanceNorm",
        "AlphaDropout": "amici.exporters.jax.AlphaDropout",
        "Bilinear": "amici.exporters.jax.Bilinear",
        "Dropout1d": "eqx.nn.Dropout",
        "Dropout2d": "eqx.nn.Dropout",
        "Flatten": "amici.exporters.jax.Flatten",
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
        "AlphaDropout": ("inplace",),
        "Dropout": ("inplace",),
        "Dropout1d": ("inplace",),
        "Dropout2d": ("inplace",),
        "LayerNorm": ("bias",),
        "BatchNorm1d": ("track_running_stats", "momentum"),
        "BatchNorm2d": ("track_running_stats", "momentum"),
        "BatchNorm3d": ("track_running_stats", "momentum"),
        "InstanceNorm1d": ("track_running_stats", "momentum"),
        "InstanceNorm2d": ("track_running_stats", "momentum"),
        "InstanceNorm3d": ("track_running_stats", "momentum"),
    }
    # construct argument string for layer instantiation
    kwargs = [
        f"{kwarg_map.get(layer.layer_type, {}).get(k, k)}={_process_argval(v)}"
        for k, v in layer.args.items()
        if k not in kwarg_ignore.get(layer.layer_type, ())
    ]
    # add key for initialization
    if layer.layer_type in (
        "Bilinear",
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
        "tanhshrink": "amici.exporters.jax.tanhshrink",
        "softsign": "jax.nn.soft_sign",
        "cat": "amici.exporters.jax.cat",
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

    # Handle kwarg aliasing for cat (dim -> axis)
    if node.target == "cat":
        if "dim" in node.kwargs:
            node.kwargs["axis"] = node.kwargs.pop("dim")
        # Convert list of variable names to proper bracket-enclosed list
        if isinstance(node.args[0], list):
            # node.args[0] is a list like ['net_input1', 'net_input2']
            # We need to convert it to a single string representing the list: [net_input1, net_input2]
            node.args = tuple(
                ["[" + ", ".join(node.args[0]) + "]"] + list(node.args[1:])
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

    # Handle placeholder nodes - skip individual processing, handled collectively in generate_equinox
    if node.op == "placeholder":
        return ""

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
    if layer_type.startswith("Dropout") or layer_type == "AlphaDropout":
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
