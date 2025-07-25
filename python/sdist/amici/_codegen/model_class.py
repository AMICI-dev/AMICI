"""Function for generating the ``amici::Model`` subclass for an amici model."""

from __future__ import annotations

from .cxx_functions import functions, multiobs_functions


def get_function_extern_declaration(fun: str, name: str, ode: bool) -> str:
    """
    Constructs the extern function declaration for a given function

    :param fun:
        function name
    :param name:
        model name
    :param ode:
        whether to generate declaration for DAE or ODE

    :return:
        C++ function definition string
    """
    f = functions[fun]
    return f"extern {f.return_type} {fun}_{name}({f.arguments(ode)});"


def get_sunindex_extern_declaration(
    fun: str, name: str, indextype: str
) -> str:
    """
    Constructs the function declaration for an index function of a given
    function

    :param fun:
        function name

    :param name:
        model name

    :param indextype:
        index function {'colptrs', 'rowvals'}

    :return:
        C++ function declaration string
    """
    index_arg = ", int index" if fun in multiobs_functions else ""
    return (
        f"extern void {fun}_{indextype}_{name}"
        f"(SUNMatrixWrapper &{indextype}{index_arg});"
    )


def get_model_override_implementation(
    fun: str, name: str, ode: bool, nobody: bool = False
) -> str:
    """
    Constructs ``amici::Model::*`` override implementation for a given function

    :param fun:
        function name

    :param name:
        model name

    :param nobody:
        whether the function has a nontrivial implementation

    :return:
        C++ function implementation string
    """
    func_info = functions[fun]
    body = (
        f" return {func_info.default_return_value}; "
        if nobody and func_info.default_return_value
        else ""
        if nobody
        else "\n{ind8}{maybe_return}{fun}_{name}({eval_signature});\n{ind4}".format(
            ind4=" " * 4,
            ind8=" " * 8,
            maybe_return="" if func_info.return_type == "void" else "return ",
            fun=fun,
            name=name,
            eval_signature=remove_argument_types(func_info.arguments(ode)),
        )
    )
    return f"{func_info.return_type} f{fun}({func_info.arguments(ode)}) override {{{body}}}\n"


def get_sunindex_override_implementation(
    fun: str, name: str, indextype: str, nobody: bool = False
) -> str:
    """
    Constructs the ``amici::Model`` function implementation for an index
    function of a given function

    :param fun:
        function name

    :param name:
        model name

    :param indextype:
        index function {'colptrs', 'rowvals'}

    :param nobody:
        whether the corresponding function has a nontrivial implementation

    :return:
        C++ function implementation string
    """
    index_arg = ", int index" if fun in multiobs_functions else ""
    index_arg_eval = ", index" if fun in multiobs_functions else ""

    impl = "void f{fun}_{indextype}({signature}) override {{"

    if nobody:
        impl += "}}\n"
    else:
        impl += (
            "\n{ind8}{fun}_{indextype}_{name}({eval_signature});\n{ind4}}}\n"
        )

    return impl.format(
        ind4=" " * 4,
        ind8=" " * 8,
        fun=fun,
        indextype=indextype,
        name=name,
        signature=f"SUNMatrixWrapper &{indextype}{index_arg}",
        eval_signature=f"{indextype}{index_arg_eval}",
    )


def remove_argument_types(signature: str) -> str:
    """
    Strips argument types from a function signature

    :param signature:
        function signature

    :return:
        string that can be used to construct function calls with the same
        variable names and ordering as in the function signature
    """
    # remove * prefix for pointers (pointer must always be removed before
    # values otherwise we will inadvertently dereference values,
    # same applies for const specifications)
    #
    # always add whitespace after type definition for cosmetic reasons
    known_types = [
        "const realtype *",
        "const double *",
        "const realtype ",
        "double *",
        "realtype *",
        "const int ",
        "int ",
        "bool ",
        "SUNMatrixContent_Sparse ",
        "gsl::span<const int>",
    ]

    for type_str in known_types:
        signature = signature.replace(type_str, "")

    return signature
