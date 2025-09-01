Installation
============

Requirements

* C++17 compiler
* cmake 3.16+
* `fmt`_ 11.x
* python
* pybind11 (>=2.9.2 - if building with python bindings)

.. _fmt: https://github.com/fmtlib/fmt

CMake
-----

Standard installation steps should build the library:

.. code-block:: shell

   $ cd /path/to/fprops
   $ mkdir build
   $ cd build
   $ cmake ..
   $ make
   $ make install
