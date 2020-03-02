{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

{% if classes %}
.. rubric:: Classes

.. autosummary::
    :toctree: .
    {% for class in classes %}
    {{ class }}
    {% endfor %}

{% endif %}

{% if functions %}
.. rubric:: Functions Summary

.. autosummary::
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}

{% if functions %}
.. rubric:: Functions

{% for function in functions %}
.. autofunction:: {{ function }}
{% endfor %}

{% endif %}