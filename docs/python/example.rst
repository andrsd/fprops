Example
=======

This code snippet evaluates the fluid properties of nitrogen at pressure = 1 MPa and
temperature = 280 K.

.. code::

   import fprops as fp
   n2 = fp.Nitrogen()
   # p = 1MPa, T = 280 K
   s = n2.p_T(1e6, 280)
   # density
   s.rho
