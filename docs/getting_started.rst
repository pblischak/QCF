.. _Getting_Started:

Getting Started
===============

Installation
------------

``qcf`` can be installed by cloning the code from GitHub using the following steps:

.. code:: bash

  git clone https://github.com/pblischak/QCF.git        # 1. Clone the repo from GitHub
  cd QCF/                                               # 2. cd into the QCF/ folder
  make                                                  # 3. compile the qcf executable
  make test                                             # 4. test that the executable works
  sudo make install                                     # 5. copy executable to /usr/local/bin

The ``sudo make install`` step will also cp all the files in the ``scripts/`` folder
to ``/usr/local/bin``.
Stable versions of QCF are also available on the `Releases <https://github.com/pblischak/QCF/releases>`_ page.

Input Files
-----------

Phylip Files
~~~~~~~~~~~~

Sequence data for each gene should be in its own file in Phylip format.
The setup should be the same as if you were planning to run RAxML
on each gene individually.

**Example:**

.. code::

    16  500
  sp1_1 AGTACAAGGTAGACAGTAGACG...
  sp1_2 AGTACAAGGTAGACAGTAGACG...
  sp2_1 AGTACAAGGTAGACAGTAGACG...
  .
  .
  .
  spN_3 AGTACAAGGTAGACAGTAGACG...


Gene List File
~~~~~~~~~~~~~~

The gene list file is a simple text file that has the name of each Phylip
file that is to be included in an analysis on its own line.

**Example:**

.. code::

  gene1.phy
  gene2.phy
  gene3.phy
  .
  .
  .
  geneL.phy

If this file and the gene sequence files are not in the same directory, then
you can add the relevant path information to the Phylip files here so that
the program can still find them (e.g., ``path/to/geneL.phy``).

Map File
~~~~~~~~

The mapping file maps haplotypes to sampled taxa.
The easiest way to do this is to sequentially number the haplotypes
for each gene (e.g., SpeciesName_1, SpeciesName_2, etc.).
Genes are treated as independent, so they can reuse the same
haplotype names. Also, not all genes need to have all haplotypes.
For each taxon, start with its name, followed by a colon (``:``), then the
names of the haplotypes that are present in the Phylip files containing the
sequence data, each separated by a comma (``,``). **There should be no spaces**.
This format is the same as the one used by
`ASTRAL <https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#running-on-a-multi-individual-datasets>`__.

**Example:**

.. code::

  sp1:sp1_1,sp1_2,sp1_3
  sp2:sp2_1,sp2_2
  .
  .
  .
  spN:spN_1,spN_2,spN_3,spN_4

Output Files
------------

``qcf`` by default will produce an output file that contains the estimated quartet
concordance factors in a file called ``out-qcf.CFs.csv``. If you also print the raw
quartet scores, then the program will write another file called ``out-raw.csv``. This
file contains all of the raw scores for all haplotypes for each species quartet
(it is not intended to be human readable). The ``out-raw.csv`` file is what can be
passed to the ``qcf_boot.py`` Python script to conduct bootstrap resampling for confidence
interval estimation.
