Nimfem
======

.. image:: https://github.com/Psirus/Nimfem/workflows/Tests/badge.svg

A finite element library in Nim.

To run the poisson example, go to the ``example`` directory, compile it with ``nim compile -d:release poisson.nim`` and run with ``./poisson``.

For comparisons, an equivalent FEniCS file is included as ``poisson.py``.

Speed
-----

For now, Nimfem seems quite quick. This shows the runtime for different mesh sizes for the Poisson example:

.. image:: https://github.com/Psirus/Nimfem/raw/speed_comparison.py/benchmark/poisson_comparison.png

To see what was done, or run this on your machine, see ``benchmark/comparison.py``.
