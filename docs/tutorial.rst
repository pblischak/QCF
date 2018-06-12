.. _Tutorial:

Tutorial
========

Below we will go through an example analysis using some simulated data that is
available in the ``example/`` folder in the
`QCF GitHub repository <https://github.com/pblischak/QCF>`__.

Step 0: Get ``qcf``
-------------------

.. code:: bash

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
  ls -1 *.phy > genes.txt

Step 2: Run ``qcf``
-------------------

Option 1: Basic analysis with ``qcf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

  # Basic qcf run
  qcf -i genes.txt -m map.txt --prefix example

Option 2: Run ``qcf`` with bootstrapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

  # Now we'll run qcf with bootstrapping (add -b #)
  qcf -i genes.txt -m map.txt -b 500 --prefix booted

Option 3: Calculating QCFs and confidence intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Print raw quartet information**

.. code:: bash

  # We'll repeat Step 2b, but adding the --printRaw flag
  qcf -i genes.txt -m map.txt -b 500 --prefix booted2 --printRaw

**Perform bootstrap resampling on raw output**

.. code:: bash

  # Using the raw output from Step 3a, we'll use the Python script
  # qcf_boot.py to resample gene-level quartets to calculate QCF values
  # and their 95% confidence intervals
  qcf_boot.py -i booted2-raw.csv -b 500 --prefix resampled

Option 4: Analyzing genes in parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Perform separate runs on gene sets**

.. code:: bash

  # First we'll analyze gene set 1
  qcf -i genes1.txt -m map.txt -b 500 --prefix out1 --printRaw

  # Now gene set 2
  qcf -i genes2.txt -m map.txt -b 500 --prefix out2 --printRaw

**Processing raw output across independent runs**

.. code:: bash

  # Now we'll calculate QCFs and their confidence intervals across
  # the independent runs we just completed. The qcf_boot.py script
  # is written such that it can combine the raw data across any number
  # of independent runs.
  qcf_boot.py -i out1-raw.csv out2-raw.csv -b 500 --prefix resampled2

If you have more than 2 input files, you can list them all after the ``-i``
flag:

.. code:: bash

  qcf_boot.py -i out1-raw.csv out2-raw.csv out3-raw.csv <...more files...> \
              -b 500 --prefix resampled3

An easy way to list them all would be to do something like this:

.. code:: bash

  qcf_boot.py -i $(ls *-raw.csv) -b 500 --prefix resampled4
