"""Tests for SBML/SciML functionality, including JAX neural network code generation."""

from unittest.mock import Mock

import pytest

pytest.importorskip("jax")
pytest.importorskip("equinox")

import pytest
from amici.jax.nn import (
    _format_function_call,
    _generate_forward,
    _process_activation_call,
    _process_layer_call,
)


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
        assert fun_str == "amici.jax.tanhshrink"

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


class TestProcessHybridizationErrors:
    """Test the improved error messages in _process_hybridization."""

    @pytest.fixture
    def mock_de_model(self):
        """Create a mock DEModel instance for testing."""
        import sympy as sp
        from amici.de_model import DEModel
        from amici.de_model_components import (
            DifferentialState,
            Expression,
            FreeParameter,
            Observable,
        )

        model = DEModel()

        # Add some parameters
        model._parameters = [
            FreeParameter(sp.Symbol("p1"), "param1", sp.Float(1.0)),
            FreeParameter(sp.Symbol("p2"), "param2", sp.Float(2.0)),
        ]

        # Add some expressions
        model._expressions = [
            Expression(sp.Symbol("expr1"), "expression1", sp.Float(0.5)),
            Expression(sp.Symbol("expr2"), "expression2", sp.Float(0.7)),
        ]

        # Add some differential states
        model._differential_states = [
            DifferentialState(
                sp.Symbol("x1"), "state1", sp.Float(0.0), sp.Float(0.1)
            ),
            DifferentialState(
                sp.Symbol("x2"), "state2", sp.Float(0.0), sp.Float(0.2)
            ),
        ]

        # Add some observables
        model._observables = [
            Observable(sp.Symbol("obs1"), "observable1", sp.Symbol("x1")),
            Observable(sp.Symbol("obs2"), "observable2", sp.Symbol("x2")),
        ]

        return model

    def test_missing_input_variables(self, mock_de_model):
        """Test error message for missing input variables."""
        hybridization = {
            "neural_net1": {
                "static": False,
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
                "static": False,
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
                "static": False,
                "input_vars": ["p1"],
                "output_vars": {"expr1": 0},
                "observable_vars": {"obs1": 0, "obs_missing": 1},
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
                "static": False,
                "input_vars": ["p1", "p2"],
                "output_vars": {"expr1": 0},
                "observable_vars": {"obs1": 0},
            }
        }

        # Should not raise any errors
        mock_de_model._process_hybridization(hybridization)
