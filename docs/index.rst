.. QCF documentation master file, created by
   sphinx-quickstart on Fri Feb 23 15:23:59 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``qcf``: estimating quartet concordance factors
===============================================

|Build Status| |Documentation| |QCF Version|

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

.. |Build Status| image:: https://travis-ci.org/pblischak/QCF.svg?branch=master
   :target: https://travis-ci.org/pblischak/QCF

.. |Documentation| image:: http://readthedocs.org/projects/qcf/badge/?version=latest
   :target: http://qcf.readthedocs.io

.. |QCF Version| image:: https://img.shields.io/badge/qcf-v0.1.3a-orange.svg
   :target: https://github.com/pblischak/QCF/releases/tag/0.1.3a
