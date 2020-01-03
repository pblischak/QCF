.. _Tutorial:

Tutorial
========

Below we will go through an example analysis using some simulated data that is
available in the ``example/`` folder in the
`QCF GitHub repository <https://github.com/pblischak/QCF>`__.

Step 0: Get ``qcf``
-------------------

If you don't have the repo, clone it and install qcf. If you already have the QCF
repository then go ahead and skip to the next step.

.. code:: bash

  git clone https://github.com/pblischak/QCF.git
  cd QCF
  make
  make test
  sudo make install

Step 1: Look at Example Files
-----------------------------

Now we'll run the example analyses. First, we'll change into the `example/` folder in the QCF repo.

.. code:: bash

  # cd into the example folder (wherever it is on your computer)
  cd /path/to/QCF/example

  # check out what files are there
  ls

You should see the Phylip files containing the sequence data, the ``genes.txt`` file containing
a list of all of the genes, and the ``map.txt`` file.

.. code:: bash

  # Aside: Making the gene list file can be done in the folder with all of the
  # Phylip files using the following code within a terminal
  ls -1 *.phy > genes.txt

Step 2: Run ``qcf``
-------------------

The most basic analysis that we can do with ``qcf`` is to just
calculate QCF scores with the ``genes.txt`` and ``map.txt`` files.

.. code:: bash

  qcf -i genes.txt -m map.txt --prefix example1

If you want to incorporate uncertainty into the estimation of the QCF scores,
there is an option to perform bootstrap resampling of sites within each gene.
This will make the analysis slower depending on how many bootstrap replicates
you perform.

.. code:: bash

  # add -b <#> for bootstrapping
  qcf -i genes.txt -m map.txt -b 500 --prefix example2

To calculate confidence intervals for the QCF scores, we'll add the ``--printRaw``
flag when calling ``qcf``. This will generate an extra output file that can be used
with the ``qcf_boot.py`` script to

.. code:: bash

  # add the --printRaw flag
  qcf -i genes.txt -m map.txt -b 500 --prefix example3 --printRaw

Using the raw output from the previous step, we'll use the Python script
``qcf_boot.py`` to resample gene-level quartets to calculate QCF values
and their 95% confidence intervals.

.. code:: bash

  qcf_boot.py -i example3-raw.csv -b 500 --prefix resampled

This will generate the file ``resampled-boot.CFs.csv``. This file can be used to infer
a species tree with scripts from the TICR pipeline (Stenz et al. 2015), which are packaged
with QCF in the ``scripts/`` folder (should be available if you ran ``sudo make install``).

.. note:: **Citing TICR**

  If you use these scripts please be sure to cite the TICR pipeline:

  Stenz, N. W. M., B. Larget, D. A. Baum, and C. Ane. 2015.
  Exploring Tree-Like and Non-Tree-Like Patterns Using Genome Sequences:
  An Example Using the Inbreeding Plant Species *Arabidopsis thaliana* (L.) Heynh.
  *Systematic Biology* 64:809--823.

To use these scripts, you will also need to install QuartetMaxCut, which is available
`here <http://research.haifa.ac.il/~ssagi/software/QMCN.tar.gz>`__. The
`TICR README <https://github.com/nstenz/TICR>`__
has a lot of helpful information for using these scripts as well.

.. code:: bash

  # Get a tree topology using QuartetMaxCut
  # Usage:
  #   perl get-pop-tree.pl <bootstrap QCF file>
  perl get-pop-tree.pl resampled-boot.CFs.csv

  # Estimate branch lengths in coalescent units
  # Usage:
  #   Rscript --vanilla getTreeBranchLengths.R <bootstrap file prefix> <outgroup>
  Rscript --vanilla getTreeBranchLengths.R resampled-boot 6

The ``resampled-boot.CFs.csv`` file is also formatted to be analyzed using the SNaQ
species network inference method in the PhyloNetworks package. Documentation for running
SNaQ is available on the `PhyloNetworks website <https://crsl4.github.io/PhyloNetworks.jl/stable/>`__.

Analyzing Genes in Parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a large number of genes, it is possible to analyze smaller numbers of
genes separately and in parallel to make analyses more computationally efficient.
To do this, instead of listing all genes in one file, create several files listing
different groups of genes, analyze each one on its own (you can use the same map
file for each), and then combine them using the ``qcf_boot.py`` script. Because
the results are combined using ``qcf_boot.py``, each analysis will have to be run
with the ``--printRaw`` flag.

.. code:: bash

  # First we'll analyze gene set 1
  qcf -i genes1.txt -m map.txt -b 500 --prefix out1 --printRaw

  # Now gene set 2
  qcf -i genes2.txt -m map.txt -b 500 --prefix out2 --printRaw

Now we'll calculate QCFs and their confidence intervals across
the independent runs we just completed. The qcf_boot.py script
is written such that it can combine the raw data across any number
of independent runs.

.. code:: bash

  #
  qcf_boot.py -i out1-raw.csv out2-raw.csv -b 500 --prefix resampled2

If you have more than 2 input files, you can list them all after the ``-i``
flag:

.. code:: bash

  qcf_boot.py -i out1-raw.csv out2-raw.csv out3-raw.csv <...more files...> \
              -b 500 --prefix resampled3

An easy way to list them all would be to do something like this:

.. code:: bash

  qcf_boot.py -i $(ls *-raw.csv) -b 500 --prefix resampled4
