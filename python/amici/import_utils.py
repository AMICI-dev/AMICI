"""Miscellaneous functions related to model import, independent of any specific
 model format"""

from typing import Dict, Union, Optional, Callable

import sympy as sp
from toposort import toposort

SymbolDef = Dict[sp.Symbol, Union[Dict[str, sp.Expr], sp.Expr]]


def noise_distribution_to_cost_function(
        noise_distribution: str
) -> Callable[[str], str]:
    """
    Parse noise distribution string to a cost function definition amici can
    work with.

    The noise distributions listed in the following are supported. :math:`m`
    denotes the measurement, :math:`y` the simulation, and :math:`\\sigma` a
    distribution scale parameter
    (currently, AMICI only supports a single distribution parameter).

    - `'normal'`, `'lin-normal'`: A normal distribution:

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma}\\
         exp\\left(-\\frac{(m-y)^2}{2\\sigma^2}\\right)

    - `'log-normal'`: A log-normal distribution (i.e. log(m) is
      normally distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma m}\\
         exp\\left(-\\frac{(\\log m - \\log y)^2}{2\\sigma^2}\\right)

    - `'log10-normal'`: A log10-normal distribution (i.e. log10(m) is
      normally distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma m \\log(10)}\\
         exp\\left(-\\frac{(\\log_{10} m - \\log_{10} y)^2}{2\\sigma^2}\\right)

    - `'laplace'`, `'lin-laplace'`: A laplace distribution:

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma}
         \\exp\\left(-\\frac{|m-y|}{\\sigma}\\right)

    - `'log-laplace'`: A log-Laplace distribution (i.e. log(m) is Laplace
      distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma m}
         \\exp\\left(-\\frac{|\\log m - \\log y|}{\\sigma}\\right)

    - `'log10-laplace'`: A log10-Laplace distribution (i.e. log10(m) is
      Laplace distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma m \\log(10)}
         \\exp\\left(-\\frac{|\\log_{10} m - \\log_{10} y|}{\\sigma}\\right)

    - `'binomial'`, `'lin-binomial'`: A (continuation of a discrete) binomial
      distribution, parameterized via the success probability
      :math:`p=\\sigma`:

      .. math::
         \\pi(m|y,\\sigma) = \\operatorname{Heaviside}(y-m) \\cdot
                \\frac{\\Gamma(y+1)}{\\Gamma(m+1) \\Gamma(y-m+1)}
                \\sigma^m (1-\\sigma)^{(y-m)}

    - `'negative-binomial'`, `'lin-negative-binomial'`: A (continuation of a
      discrete) negative binomial distribution, with with `mean = y`,
      parameterized via success probability `p`:

      .. math::

         \\pi(m|y,\\sigma) = \\frac{\\Gamma(m+r)}{\\Gamma(m+1) \\Gamma(r)}
            (1-\\sigma)^m \\sigma^r

      where

      .. math::
         r = \\frac{1-\\sigma}{\\sigma} y

    The distributions above are for a single data point.
    For a collection :math:`D=\\{m_i\\}_i` of data points and corresponding
    simulations :math:`Y=\\{y_i\\}_i` and noise parameters
    :math:`\\Sigma=\\{\\sigma_i\\}_i`, AMICI assumes independence,
    i.e. the full distributions is

    .. math::
       \\pi(D|Y,\\Sigma) = \\prod_i\\pi(m_i|y_i,\\sigma_i)

    AMICI uses the logarithm :math:`\\log(\\pi(m|y,\\sigma)`.

    In addition to the above mentioned distributions, it is also possible to
    pass a function taking a symbol string and returning a log-distribution
    string with variables '{str_symbol}', 'm{str_symbol}', 'sigma{str_symbol}'
    for y, m, sigma, respectively.

    :param noise_distribution: An identifier specifying a noise model.
        Possible values are

        {`'normal'`, `'lin-normal'`, `'log-normal'`, `'log10-normal'`,
        `'laplace'`, `'lin-laplace'`, `'log-laplace'`, `'log10-laplace'`,
        `'binomial'`, `'lin-binomial'`, `'negative-binomial'`,
        `'lin-negative-binomial'`, `<Callable>`}

        For the meaning of the values see above.

    :return: A function that takes a strSymbol and then creates a cost
        function string (negative log-likelihood) from it, which can be
        sympified.
    """

    if isinstance(noise_distribution, Callable):
        return noise_distribution

    if noise_distribution in ['normal', 'lin-normal']:
        y_string = '0.5*log(2*pi*{sigma}**2) + 0.5*(({y} - {m}) / {sigma})**2'
    elif noise_distribution == 'log-normal':
        y_string = '0.5*log(2*pi*{sigma}**2*{m}**2) ' \
                   '+ 0.5*((log({y}) - log({m})) / {sigma})**2'
    elif noise_distribution == 'log10-normal':
        y_string = '0.5*log(2*pi*{sigma}**2*{m}**2*log(10)**2) ' \
                   '+ 0.5*((log({y}, 10) - log({m}, 10)) / {sigma})**2'
    elif noise_distribution in ['laplace', 'lin-laplace']:
        y_string = 'log(2*{sigma}) + Abs({y} - {m}) / {sigma}'
    elif noise_distribution == 'log-laplace':
        y_string = 'log(2*{sigma}*{m}) + Abs(log({y}) - log({m})) / {sigma}'
    elif noise_distribution == 'log10-laplace':
        y_string = 'log(2*{sigma}*{m}*log(10)) ' \
                   '+ Abs(log({y}, 10) - log({m}, 10)) / {sigma}'
    elif noise_distribution in ['binomial', 'lin-binomial']:
        # Binomial noise model parameterized via success probability p
        y_string = '- log(Heaviside({y} - {m})) - loggamma({y}+1) ' \
                   '+ loggamma({m}+1) + loggamma({y}-{m}+1) ' \
                   '- {m} * log({sigma}) - ({y} - {m}) * log(1-{sigma})'
    elif noise_distribution in ['negative-binomial', 'lin-negative-binomial']:
        # Negative binomial noise model of the number of successes m
        # (data) before r=(1-sigma)/sigma * y failures occur,
        # with mean number of successes y (simulation),
        # parameterized via success probability p = sigma.
        r = '{y} * (1-{sigma}) / {sigma}'
        y_string = f'- loggamma({{m}}+{r}) + loggamma({{m}}+1) ' \
                   f'+ loggamma({r}) - {r} * log(1-{{sigma}}) ' \
                   f'- {{m}} * log({{sigma}})'
    else:
        raise ValueError(
            f"Cost identifier {noise_distribution} not recognized.")

    def nllh_y_string(str_symbol):
        y, m, sigma = _get_str_symbol_identifiers(str_symbol)
        return y_string.format(y=y, m=m, sigma=sigma)

    return nllh_y_string


def _get_str_symbol_identifiers(str_symbol: str) -> tuple:
    """Get identifiers for simulation, measurement, and sigma."""
    y, m, sigma = f"{str_symbol}", f"m{str_symbol}", f"sigma{str_symbol}"
    return y, m, sigma


def smart_subs_dict(sym: sp.Expr,
                    subs: SymbolDef,
                    field: Optional[str] = None,
                    reverse: bool = True) -> sp.Expr:
    """
    Subsitutes expressions completely flattening them out. Requires
    sorting of expressions with toposort.

    :param sym:
        Symbolic expression in which expressions will be substituted

    :param subs:
        Substitutions

    :param field:
        Field of substitution expressions in subs.values(), if applicable

    :param reverse:
        Whether ordering in subs should be reversed. Note that substitution
        requires the reverse order of what is required for evaluation.

    :return:
        Substituted symbolic expression
    """
    s = [
        (eid, expr[field] if field is not None else expr)
        for eid, expr in subs.items()
    ]
    if reverse:
        s.reverse()
    for substitution in s:
        # note that substitution may change free symbols, so we have to do
        # this recursively
        if substitution[0] in sym.free_symbols:
            sym = sym.subs(*substitution)
    return sym


def smart_subs(element: sp.Expr, old: sp.Symbol, new: sp.Expr) -> sp.Expr:
    """
    Optimized substitution that checks whether anything needs to be done first

    :param element:
        substitution target

    :param old:
        to be substituted

    :param new:
        subsitution value

    :return:
        substituted expression
    """

    if old in element.free_symbols:
        return element.subs(old, new)

    return element


def toposort_symbols(symbols: SymbolDef,
                     field: Optional[str] = None) -> SymbolDef:
    """
    Topologically sort symbol definitions according to their interdependency

    :param symbols:
        symbol definitions

    :param field:
        field of definition.values() that is used to compute interdependency

    :return:
        ordered symbol definitions
    """
    sorted_symbols = toposort({
        identifier: {
            s for s in (
                definition[field] if field is not None else definition
            ).free_symbols
            if s in symbols
        }
        for identifier, definition
        in symbols.items()
    })
    return {
        s: symbols[s]
        for symbol_group in sorted_symbols
        for s in symbol_group
    }
