{{ name }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ name }}
   :members:
   :show-inheritance:
   :special-members: __init__

   {% block methods %}

   {% if methods %}
   .. rubric:: Methods Overview

   .. autosummary::
      :nosignatures:
   
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes Overview

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block transition_to_details %}
   {% if methods or attributes %}
   .. raw:: html
      
      <hr>
      <h2>{{name}} Method & Attribute Details</h2>

   {% endif %}
   {% endblock %}
