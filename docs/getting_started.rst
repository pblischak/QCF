.. _Getting_Started:

Getting Started
===============

Installation
------------

.. code:: bash

  git clone https://github.com/pblischak/QCF.git
  cd QCF
  make
  sudo make install

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
  sp1 AGTACAAGGTAGACAGTAGACG...
  sp2 AGTACAAGGTAGACAGTAGACG...
  sp3 AGTACAAGGTAGACAGTAGACG...
  .
  .
  .
  spN AGTACAAGGTAGACAGTAGACG...


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
