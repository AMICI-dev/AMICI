"""Tests for SBML/SciML functionality, including JAX neural network code generation."""

from unittest.mock import Mock

import pytest

pytest.importorskip("jax")
pytest.importorskip("equinox")

import pytest
from amici.exporters.jax.nn import (
    _format_function_call,
    _generate_forward,
    _generate_layer,
    _process_activation_call,
    _process_layer_call,
)
from amici.importers.utils import symbol_with_assumptions


class TestFormatFunctionCall:
    """Test the utility function for formatting function calls."""

    def test_format_with_args_only(self):
        """Test formatting with only positional arguments."""
        result = _format_function_call(
            var_name="output",
            fun_str="my_function",
            args=["x", "y"],
            kwargs=[],
            indent=4,
        )
        assert result == "    output = my_function(x, y)"

    def test_format_with_kwargs_only(self):
        """Test formatting with only keyword arguments."""
        result = _format_function_call(
            var_name="output",
            fun_str="my_function",
            args=[],
            kwargs=["a=1", "b=2"],
            indent=4,
        )
        assert result == "    output = my_function(a=1, b=2)"

    def test_format_with_args_and_kwargs(self):
        """Test formatting with both positional and keyword arguments."""
        result = _format_function_call(
            var_name="result",
            fun_str="jax.nn.relu",
            args=["input_tensor"],
            kwargs=["axis=1"],
            indent=8,
        )
        assert result == "        result = jax.nn.relu(input_tensor, axis=1)"

    def test_format_with_no_args(self):
        """Test formatting with no arguments."""
        result = _format_function_call(
            var_name="output",
            fun_str="get_value",
            args=[],
            kwargs=[],
            indent=0,
        )
        assert result == "output = get_value()"

    def test_format_with_zero_indent(self):
        """Test formatting with zero indentation."""
        result = _format_function_call(
            var_name="x",
            fun_str="func",
            args=["a"],
            kwargs=["b=2"],
            indent=0,
        )
        assert result == "x = func(a, b=2)"


class TestGenerateLayerCall:
    def test_flatten(self):
        """Test generation for Flatten layer."""
        flatten = Mock()
        flatten.layer_id = "flat"
        flatten.layer_type = "Flatten"
        flatten.args = {}

        layer_str = _generate_layer(flatten, indent=4, ilayer=0)

        assert "amici.exporters.jax.Flatten" in layer_str

    def test_batchnorm(self):
        """Test generation for BatchNorm layer."""
        batchnorm = Mock()
        batchnorm.layer_id = "bn1"
        batchnorm.layer_type = "BatchNorm1d"
        batchnorm.args = {
            "num_features": 64,
            "eps": 1e-5,
            "momentum": 0.1,
            "affine": True,
            "track_running_stats": True,
        }

        layer_str = _generate_layer(batchnorm, indent=4, ilayer=0)

        assert "amici.exporters.jax.BatchNorm" in layer_str
        assert "track_running_stats" not in layer_str
        assert "momentum" not in layer_str
        assert "num_features" in layer_str
        assert "eps" in layer_str
        assert "affine" in layer_str

    def test_instancenorm(self):
        """Test generation for InstanceNorm layer."""
        instancenorm = Mock()
        instancenorm.layer_id = "in1"
        instancenorm.layer_type = "InstanceNorm1d"
        instancenorm.args = {
            "num_features": 64,
            "eps": 1e-5,
            "momentum": 0.1,
            "affine": True,
            "track_running_stats": True,
        }

        layer_str = _generate_layer(instancenorm, indent=4, ilayer=0)

        assert "amici.exporters.jax.InstanceNorm" in layer_str
        assert "track_running_stats" not in layer_str
        assert "momentum" not in layer_str
        assert "num_features" in layer_str
        assert "eps" in layer_str
        assert "affine" in layer_str

    def test_dropout(self):
        """Test generation for Dropout layer."""
        dropout = Mock()
        dropout.layer_id = "drop"
        dropout.layer_type = "Dropout1d"
        dropout.args = {}

        layer_str = _generate_layer(dropout, indent=4, ilayer=0)

        assert "eqx.nn.Dropout" in layer_str
        assert "inplace" not in layer_str

    def test_dropout2d(self):
        """Test generation for Dropout2d layer."""
        dropout2d = Mock()
        dropout2d.layer_id = "drop2d"
        dropout2d.layer_type = "Dropout2d"
        dropout2d.args = {}

        layer_str = _generate_layer(dropout2d, indent=4, ilayer=0)

        assert "eqx.nn.Dropout" in layer_str
        assert "inplace" not in layer_str

    def test_alphadropout(self):
        """Test generation for AlphaDropout layer."""
        alphadropout = Mock()
        alphadropout.layer_id = "adrop"
        alphadropout.layer_type = "AlphaDropout"
        alphadropout.args = {"p": 0.5, "inplace": False}

        layer_str = _generate_layer(alphadropout, indent=4, ilayer=0)

        assert "amici.exporters.jax.AlphaDropout" in layer_str
        assert "p=0.5" in layer_str
        assert "inplace" not in layer_str

    def test_bilinear(self):
        """Test generation for Bilinear layer."""
        bilinear = Mock()
        bilinear.layer_id = "bil"
        bilinear.layer_type = "Bilinear"
        bilinear.args = {
            "in1_features": 128,
            "in2_features": 64,
            "out_features": 32,
            "bias": True,
        }

        layer_str = _generate_layer(bilinear, indent=4, ilayer=0)

        assert "amici.exporters.jax.Bilinear" in layer_str
        assert "in1_features" in layer_str
        assert "in2_features" in layer_str
        assert "out_features" in layer_str
        assert "bias" in layer_str

    def test_maxpool_not_implemented(self):
        """Test that MaxPool layer with dilation raises NotImplementedError."""
        maxpool = Mock()
        maxpool.layer_id = "mp"
        maxpool.layer_type = "MaxPool2d"
        maxpool.args = {"dilation": 1}

        with pytest.raises(
            NotImplementedError,
            match="MaxPool layers with dilation not supported",
        ):
            _generate_layer(maxpool, indent=4, ilayer=0)

    def test_dropout_not_implemented(self):
        """Test that Dropout layer with inplace raises NotImplementedError."""
        dropout = Mock()
        dropout.layer_id = "drop"
        dropout.layer_type = "Dropout1d"
        dropout.args = {"inplace": True}

        with pytest.raises(
            NotImplementedError,
            match="Dropout layers with inplace not supported",
        ):
            _generate_layer(dropout, indent=4, ilayer=0)


class TestProcessLayerCall:
    """Test layer-specific processing logic."""

    def test_simple_layer_no_freezing(self):
        """Test processing a simple layer without freezing."""
        node = Mock()
        node.target = "layer1"
        node.name = "conv1"
        node.args = ["input"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Conv2d", frozen_layers={}
        )

        assert fun_str.startswith("(jax.vmap(self.layers['layer1'])")
        assert tree_string == ""

    def test_frozen_layer_with_attribute(self):
        """Test processing a frozen layer with specific attribute."""
        node = Mock()
        node.target = "layer1"
        node.name = "conv1"
        node.args = ["input"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Conv2d", frozen_layers={"conv1": "weight"}
        )

        assert "tree_conv1" in fun_str
        assert "tree_conv1 = eqx.tree_at(" in tree_string
        assert "'weight'" in tree_string

    def test_frozen_layer_full_stop_gradient(self):
        """Test processing a fully frozen layer."""
        node = Mock()
        node.target = "layer1"
        node.name = "linear1"
        node.args = ["input"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Linear", frozen_layers={"linear1": False}
        )

        assert "jax.lax.stop_gradient(self.layers['layer1'])" in fun_str
        assert tree_string == ""

    def test_linear_layer_vmap(self):
        """Test that Linear layer gets vmap wrapper."""
        node = Mock()
        node.target = "fc1"
        node.name = "fc1"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Linear", frozen_layers={}
        )

        assert "jax.vmap" in fun_str
        assert "if len(x.shape) == 2" in fun_str

    def test_conv1d_layer_vmap(self):
        """Test that Conv1d layer gets vmap wrapper with correct dimensions."""
        node = Mock()
        node.target = "conv"
        node.name = "conv"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Conv1d", frozen_layers={}
        )

        assert "jax.vmap" in fun_str
        assert "if len(x.shape) == 3" in fun_str

    def test_conv2d_layer_vmap(self):
        """Test that Conv2d layer gets vmap wrapper with correct dimensions."""
        node = Mock()
        node.target = "conv"
        node.name = "conv"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Conv2d", frozen_layers={}
        )

        assert "jax.vmap" in fun_str
        assert "if len(x.shape) == 4" in fun_str

    def test_layernorm_vmap(self):
        """Test that LayerNorm layer gets vmap wrapper."""
        node = Mock()
        node.target = "norm"
        node.name = "norm"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="LayerNorm", frozen_layers={}
        )

        assert "jax.vmap" in fun_str
        assert "len(self.layers['norm'].shape)+1" in fun_str

    def test_non_vmap_layer(self):
        """Test layer that doesn't require vmap."""
        node = Mock()
        node.target = "dropout"
        node.name = "dropout"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Dropout", frozen_layers={}
        )

        assert "jax.vmap" not in fun_str
        assert fun_str == "self.layers['dropout']"

    def test_batchnorm1d(self):
        """Test process BatchNorm1d layer."""
        node = Mock()
        node.target = "bn1"
        node.name = "bn1"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="BatchNorm1d", frozen_layers={}
        )

        assert "self.layers['bn1']" in fun_str

    def test_batchnorm2d(self):
        """Test process BatchNorm2d layer."""
        node = Mock()
        node.target = "bn2"
        node.name = "bn2"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="BatchNorm2d", frozen_layers={}
        )

        assert "self.layers['bn2']" in fun_str

    def test_instancenorm1d(self):
        """Test process InstanceNorm1d layer."""
        node = Mock()
        node.target = "in1"
        node.name = "in1"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="InstanceNorm1d", frozen_layers={}
        )

        assert "self.layers['in1']" in fun_str

    def test_instancenorm2d(self):
        """Test process InstanceNorm2d layer."""
        node = Mock()
        node.target = "in2"
        node.name = "in2"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="InstanceNorm2d", frozen_layers={}
        )

        assert "self.layers['in2']" in fun_str

    def test_alphadropout(self):
        """Test process AlphaDropout layer."""
        node = Mock()
        node.target = "adrop"
        node.name = "adrop"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="AlphaDropout", frozen_layers={}
        )

        assert "self.layers['adrop']" in fun_str

    def test_bilinear(self):
        """Test process Bilinear layer."""
        node = Mock()
        node.target = "bil"
        node.name = "bil"
        node.args = ["x1", "x2"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Bilinear", frozen_layers={}
        )

        assert "self.layers['bil']" in fun_str

    def test_flatten(self):
        """Test process Flatten layer."""
        node = Mock()
        node.target = "flat"
        node.name = "flat"
        node.args = ["x"]

        fun_str, tree_string = _process_layer_call(
            node, layer_type="Flatten", frozen_layers={}
        )

        assert "self.layers['flat']" in fun_str


class TestProcessActivationCall:
    """Test activation function processing logic."""

    def test_standard_activation(self):
        """Test standard JAX activation function."""
        node = Mock()
        node.target = "relu"
        node.kwargs = {}

        fun_str = _process_activation_call(node)
        assert fun_str == "jax.nn.relu"

    def test_mapped_activation_hardtanh(self):
        """Test hardtanh activation with custom mapping."""
        node = Mock()
        node.target = "hardtanh"
        node.kwargs = {}

        fun_str = _process_activation_call(node)
        assert fun_str == "jax.nn.hard_tanh"

    def test_mapped_activation_hardsigmoid(self):
        """Test hardsigmoid activation with custom mapping."""
        node = Mock()
        node.target = "hardsigmoid"
        node.kwargs = {}

        fun_str = _process_activation_call(node)
        assert fun_str == "jax.nn.hard_sigmoid"

    def test_mapped_activation_tanhshrink(self):
        """Test tanhshrink activation with custom mapping."""
        node = Mock()
        node.target = "tanhshrink"
        node.kwargs = {}

        fun_str = _process_activation_call(node)
        assert fun_str == "amici.exporters.jax.tanhshrink"

    def test_hardtanh_valid_params(self):
        """Test hardtanh with valid default parameters."""
        node = Mock()
        node.target = "hardtanh"
        node.kwargs = {"min_val": -1.0, "max_val": 1.0}

        fun_str = _process_activation_call(node)
        assert fun_str == "jax.nn.hard_tanh"

    def test_hardtanh_invalid_min_val(self):
        """Test hardtanh raises error for non-default min_val."""
        node = Mock()
        node.target = "hardtanh"
        node.kwargs = {"min_val": -2.0}

        with pytest.raises(NotImplementedError, match="min_val != -1.0"):
            _process_activation_call(node)

    def test_hardtanh_invalid_max_val(self):
        """Test hardtanh raises error for non-default max_val."""
        node = Mock()
        node.target = "hardtanh"
        node.kwargs = {"max_val": 2.0}

        with pytest.raises(NotImplementedError, match="max_val != 1.0"):
            _process_activation_call(node)


class TestGenerateForward:
    """Test the main forward pass generation function."""

    def test_placeholder_node(self):
        """Test generation for placeholder nodes."""
        node = Mock()
        node.op = "placeholder"
        node.name = "input_x"

        result = _generate_forward(node, indent=4)
        assert result == ""

    def test_output_node(self):
        """Test generation for output nodes."""
        node = Mock()
        node.op = "output"
        node.target = "output"
        node.args = ["y1", "y2"]

        result = _generate_forward(node, indent=8)
        assert result == "        output = y1, y2"

    def test_call_module_simple(self):
        """Test generation for simple module call."""
        node = Mock()
        node.op = "call_module"
        node.name = "x1"
        node.target = "layer1"
        node.args = ["input"]
        node.kwargs = {}

        result = _generate_forward(
            node, indent=4, frozen_layers={}, layer_type="Dropout"
        )
        assert "x1 = self.layers['layer1'](input, key=key)" in result

    def test_call_function_activation(self):
        """Test generation for activation function call."""
        node = Mock()
        node.op = "call_function"
        node.name = "act1"
        node.target = "relu"
        node.args = ["x"]
        node.kwargs = {}

        result = _generate_forward(
            node, indent=4, frozen_layers={}, layer_type=""
        )
        assert result == "    act1 = jax.nn.relu(x)"

    def test_call_module_with_frozen_layer(self):
        """Test generation for frozen layer with tree_string."""
        node = Mock()
        node.op = "call_module"
        node.name = "conv1"
        node.target = "layer1"
        node.args = ["input"]
        node.kwargs = {}

        result = _generate_forward(
            node,
            indent=4,
            frozen_layers={"conv1": "weight"},
            layer_type="Conv2d",
        )

        assert "tree_conv1 = eqx.tree_at(" in result
        assert "conv1 = " in result
        assert "\n" in result  # Should have tree_string on separate line

    def test_unsupported_operation(self):
        """Test that unsupported operations raise NotImplementedError."""
        node = Mock()
        node.op = "unknown_op"
        node.kwargs = {}

        with pytest.raises(
            NotImplementedError, match="Operation unknown_op not supported"
        ):
            _generate_forward(node, indent=4)

    def test_kwargs_filtering(self):
        """Test that 'inplace' kwarg is filtered out."""
        node = Mock()
        node.op = "call_function"
        node.name = "act1"
        node.target = "relu"
        node.args = ["x"]
        node.kwargs = {"inplace": True, "other": "value"}

        result = _generate_forward(node, indent=4, layer_type="")
        assert "inplace" not in result
        assert "other=value" in result

    def test_dropout_layer_adds_key(self):
        """Test that Dropout layers get key parameter added."""
        node = Mock()
        node.op = "call_module"
        node.name = "drop1"
        node.target = "dropout1"
        node.args = ["x"]
        node.kwargs = {}

        result = _generate_forward(
            node, indent=4, frozen_layers={}, layer_type="Dropout1d"
        )
        assert "key=key" in result

    def test_alphadropout_layer_adds_key(self):
        """Test that AlphaDropout layers get key parameter added."""
        node = Mock()
        node.op = "call_module"
        node.name = "aldrop1"
        node.target = "aldropout1"
        node.args = ["x"]
        node.kwargs = {}

        result = _generate_forward(
            node, indent=4, frozen_layers={}, layer_type="AlphaDropout"
        )
        assert "key=key" in result


class TestProcessHybridizationErrors:
    """Test the improved error messages in _process_hybridization."""

    @pytest.fixture
    def mock_de_model(self):
        """Create a mock DEModel instance for testing."""
        import sympy as sp
        from amici._symbolic.de_model import DEModel
        from amici._symbolic.de_model_components import (
            DifferentialState,
            Expression,
            FreeParameter,
            Observable,
        )

        model = DEModel()

        # Add some parameters
        model._free_parameters = [
            FreeParameter(
                symbol_with_assumptions("p1"), "param1", sp.Float(1.0)
            ),
            FreeParameter(
                symbol_with_assumptions("p2"), "param2", sp.Float(2.0)
            ),
        ]

        # Add some expressions
        model._expressions = [
            Expression(
                symbol_with_assumptions("expr1"), "expression1", sp.Float(0.5)
            ),
            Expression(
                symbol_with_assumptions("expr2"), "expression2", sp.Float(0.7)
            ),
        ]

        # Add some differential states
        model._differential_states = [
            DifferentialState(
                symbol_with_assumptions("x1"),
                "state1",
                sp.Float(0.0),
                sp.Float(0.1),
            ),
            DifferentialState(
                symbol_with_assumptions("x2"),
                "state2",
                sp.Float(0.0),
                sp.Float(0.2),
            ),
        ]

        # Add some observables
        model._observables = [
            Observable(
                symbol_with_assumptions("obs1"),
                "observable1",
                symbol_with_assumptions("x1"),
            ),
            Observable(
                symbol_with_assumptions("obs2"),
                "observable2",
                symbol_with_assumptions("x2"),
            ),
        ]

        return model

    def test_missing_input_variables(self, mock_de_model):
        """Test error message for missing input variables."""
        hybridization = {
            "neural_net1": {
                "pre_initialization": False,
                "input_vars": [
                    "p1",
                    "p2",
                    "p_missing",
                ],  # p_missing doesn't exist
                "output_vars": {"expr1": 0},
                "observable_vars": {},
            }
        }

        with pytest.raises(ValueError) as exc_info:
            mock_de_model._process_hybridization(hybridization)

        error_msg = str(exc_info.value)
        assert (
            "Could not find all input variables for neural network neural_net1"
            in error_msg
        )
        assert "Missing variables:" in error_msg
        assert "p_missing" in error_msg

    def test_missing_output_variables(self, mock_de_model):
        """Test error message for missing output variables."""
        hybridization = {
            "neural_net2": {
                "pre_initialization": False,
                "input_vars": ["p1", "p2"],
                "output_vars": {
                    "expr1": 0,
                    "expr_missing": 1,
                    "expr_also_missing": 2,
                },
                "observable_vars": {},
            }
        }

        with pytest.raises(ValueError) as exc_info:
            mock_de_model._process_hybridization(hybridization)

        error_msg = str(exc_info.value)
        assert (
            "Could not find all output variables for neural network neural_net2"
            in error_msg
        )
        assert "Missing variables:" in error_msg
        # Check that missing variables are in the message
        assert "expr_missing" in error_msg or "expr_also_missing" in error_msg

    def test_missing_observable_variables(self, mock_de_model):
        """Test error message for missing observable variables."""
        hybridization = {
            "neural_net3": {
                "pre_initialization": False,
                "input_vars": ["p1"],
                "output_vars": {"expr1": 0},
                "observable_vars": {
                    "obs1": {
                        "index": 0,
                        "formula": "net1_output1",
                        "petab_id": "net1_output1",
                    },
                    "obs_missing": {
                        "index": 1,
                        "formula": "net1_output2",
                        "petab_id": "net1_output2",
                    },
                },
            }
        }

        with pytest.raises(ValueError) as exc_info:
            mock_de_model._process_hybridization(hybridization)

        error_msg = str(exc_info.value)
        assert (
            "Could not find all observable variables for neural network neural_net3"
            in error_msg
        )
        assert "Missing variables:" in error_msg
        assert "obs_missing" in error_msg

    def test_valid_hybridization_no_error(self, mock_de_model):
        """Test that valid hybridization doesn't raise errors."""
        hybridization = {
            "valid_net": {
                "pre_initialization": False,
                "input_vars": ["p1", "p2"],
                "output_vars": {"expr1": 0},
                "observable_vars": {
                    "obs1": {
                        "index": 0,
                        "formula": "net1_output1",
                        "petab_id": "net1_output1",
                    }
                },
            }
        }

        # Should not raise any errors
        mock_de_model._process_hybridization(hybridization)
