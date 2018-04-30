.. QCF documentation master file, created by
   sphinx-quickstart on Fri Feb 23 15:23:59 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``qcf``: Estimating Quartet Concordance Factors
===========================================

``qcf`` is a C++ program for estimating quartet concordance factors
from multilocus sequence data. It uses the theory of phylogenetic
invariants to calculate the frequency of gene tree quartets matching
the three possible unrooted species tree topologies for a set of four
taxa. ``qcf`` can handle multiple haplotypes per taxon/population, and outputs a
CSV file that is compatible with the SNaQ method for phylogenetic network
inference implemented in the Julia package
`PhyloNetworks <http://crsl4.github.io/PhyloNetworks.jl/latest/>`__.

Documentation
-------------

.. toctree::
   :maxdepth: 1

   getting_started.rst
   tutorial.rst
   api.rst
