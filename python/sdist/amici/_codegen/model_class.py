"""Function for generating the ``amici::Model`` subclass for an amici model."""

from __future__ import annotations

from .cxx_functions import functions, multiobs_functions

from ..de_model_components import Event


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
        ""
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
    return "{return_type} f{fun}({signature}) override {{{body}}}\n".format(
        return_type=func_info.return_type,
        fun=fun,
        signature=func_info.arguments(ode),
        body=body,
    )


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


def get_state_independent_event_intializer(events: list[Event]) -> str:
    """Get initializer list for state independent events in amici::Model."""
    map_time_to_event_idx = {}
    for event_idx, event in enumerate(events):
        if not event.triggers_at_fixed_timepoint():
            continue
        trigger_time = float(event.get_trigger_time())
        try:
            map_time_to_event_idx[trigger_time].append(event_idx)
        except KeyError:
            map_time_to_event_idx[trigger_time] = [event_idx]

    def vector_initializer(v):
        """std::vector initializer list with elements from `v`"""
        return f"{{{', '.join(map(str, v))}}}"

    return ", ".join(
        f"{{{trigger_time}, {vector_initializer(event_idxs)}}}"
        for trigger_time, event_idxs in map_time_to_event_idx.items()
    )
