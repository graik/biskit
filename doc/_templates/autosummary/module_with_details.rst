{{ name }}
{{ underline }}

.. automodule:: {{ fullname }}
   :members:
   :exclude-members: Test

   {% block functions %}
   {% if functions %}
   .. rubric:: Functions Overview

   .. autosummary::
      :nosignatures:
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: Classes Overview

   .. autosummary::
      :nosignatures:
   {% for item in classes %}
      {% if item != 'Test' %}
          {{ item }}
      {%- endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: Exceptions Defined in this Module

   .. autosummary::
      :nosignatures:
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   .. raw:: html
               
      <hr>
      <h2>{{name}} Module Details</h2>
