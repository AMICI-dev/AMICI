{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__
   {% if methods %}
   .. rubric:: Methods Summary
   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes
   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods_full %}
   {% if methods %}
   .. rubric:: Methods
   {% for item in methods %}
   .. automethod:: {{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}


