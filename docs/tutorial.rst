.. _Tutorial:

Tutorial
========

Below we will go through an example analysis using some simulated data that is
available in the ``example/`` folder in the
`QCF GitHub repository <https://github.com/pblischak/QCF>`__.

.. code:: bash

  # Step 0:
  # If you don't have the repo, clone it and install qcf.
  # Otherwise, skip ahead to Step 1.
  git clone https://github.com/pblischak/QCF.git
  cd QCF
  make
  sudo make install

.. code:: bash

  # Step 1:
  # Now we'll run the example analyses.
  # Change into the `example/` folder in
  # the QCF repo.
  cd /path/to/QCF/example

  # check out what files are there
  ls

.. code:: bash

  # Aside:
  # Making the gene list file can be done using
  # the following code within a terminal.
  for gene in genes/*.phy; do printf $gene"\n"; done > geneFiles.txt

.. code:: bash

  # Step 2:
  # Now we'll run qcf
  qcf -i geneFiles.txt -m map.txt --prefix example

.. code:: bash

  # Bootstrapping
  qcf -i geneFiles.txt -m map.txt -b 100 --prefix booted
