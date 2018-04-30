.. _Tutorial:

Tutorial
========

Below we will go through an example analysis using some simulated data that is
available in the ``example/`` folder in the
`QCF GitHub repository <https://github.com/pblischak/QCF>`__.

Step 0: Get ``qcf``
-------------------

.. code:: bash

  # Step 0:
  # If you don't have the repo, clone it and install qcf.
  # Otherwise, skip ahead to Step 1.
  git clone https://github.com/pblischak/QCF.git
  cd QCF
  make
  make test
  sudo make install

Step 1: Look at example files
-----------------------------

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
  ls -1 *.phy > genes.txt.txt

Step 2a: Run ``qcf``
--------------------

.. code:: bash

  # Step 2a:
  # Now we'll run qcf
  qcf -i genes.txt -m map.txt --prefix example

Step 2b: Run ``qcf`` with bootstrapping
---------------------------------------

.. code:: bash

  # Step 2b:
  # Now we'll run qcf with bootstrapping (add -b #)
  qcf -i genes.txt -m map.txt -b 100 --prefix booted
