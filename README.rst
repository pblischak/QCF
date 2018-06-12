|Build Status| |Documentation| |QCF Version|

``qcf``: estimating quartet concordance factors
===============================================

  **NB**: This is an alpha-test version of ``qcf``. We are still adding and
  testing new features, but if you would like give the software a try
  for yourself, please feel free to do so. We would appreciate it if you
  report any strange behavior or bad results as an issue. Thanks!

`Read the Docs <http://qcf.readthedocs.io/>`__
==============================================

Installation
~~~~~~~~~~~~

``qcf`` can be installed by cloning this repo and then compiling the main executable
using the provided ``Makefile`` (see code below). ``qcf`` is written in C++ for Unix-like
operating systems and makes use of features from the C++11 standard, which
requires a compatible compiler (GNU ``g++`` >= 4.8, Clang ``clang++`` >= 3.3).
Mac users will also need to have the Xcode
`Command Line Tools <http://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/>`__ installed.

The software also comes with several Python scripts for processing and analyzing output.
These scripts have been tested on Python versions 2.7 and 3.6, and require the ``numpy``
package, which can be installed as follows: ``pip install numpy``.

The steps to compile, test, and install ``qcf`` and the associated scripts
is given here:

.. code:: bash

  # Code to compile and install qcf
  git clone https://github.com/pblischak/QCF.git        # 1. Clone the repo from GitHub
  cd QCF                                                # 2. cd into the QCF/ folder
  make                                                  # 3. compile the qcf executable
  make test                                             # 4. test that the executable
                                                             and Python scripts are working
  sudo make install                                     # 5. copy executable and Python scripts
                                                             to /usr/local/bin

Usage
~~~~~

After compiling and installing the software, you can use the program in a terminal
window by typing ``qcf``. Below are the main commands that can be used to estimate
quartet concordance factors using ``qcf``. More details on input data formats, as well
as background on interpreting the output, can be found in the
`tutorial <http://qcf.readthedocs.io/en/latest/tutorial.html>`__ on the ReadTheDocs site.

**Basic usage**:

.. code:: bash

  # Run the data in the example folder
  cd example
  qcf -i genes.txt -m map.txt --prefix example

  # Run with 100 bootstrap reps
  qcf -i genes.txt -m mapt.txt -b 100 --prefix bootstrap

Raw quartet scores and confidence intervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In additional to calculating the QCF values for each species quartet, we have also
provided an option to print the gene-level quartet information to file to perform
per-gene bootstrap resampling for the calculation of confidence intervals on the
the species-level QCF values. This can activated by running ``qcf`` with the
``--printRaw`` flag added:

.. code::

  # Print out the raw gene-tree quartet information
  cd example/
  qcf -i genes.txt -m map.txt --printRaw --prefix example

This raw output can then be processed using the ``qcf_boot.py`` Python script,
which will conduct bootstrap resampling across genes. Options for combining
runs across different sets of genes are available as well
(see `tutorial <http://qcf.readthedocs.io/en/latest/tutorial.html>`__).

Getting Help
~~~~~~~~~~~~

If you have questions about running ``qcf``, please feel free to use the gitter chatroom to get help:

|Gitter|

If you have a problem while running ``qcf``, and you think it may be a bug, please consider filing an issue:

|QCF Issues|

.. |Build Status| image:: https://travis-ci.org/pblischak/QCF.svg?branch=master
   :target: https://travis-ci.org/pblischak/QCF

.. |Documentation| image:: http://readthedocs.org/projects/qcf/badge/?version=latest
   :target: http://qcf.readthedocs.io

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/quartet-cf/Lobby

.. |QCF Issues| image:: https://img.shields.io/badge/qcf-issues-blue.svg
   :target: https://github.com/pblischak/QCF/issues

.. |QCF Version| image:: https://img.shields.io/badge/qcf-v0.2.0a-orange.svg
   :target: https://github.com/pblischak/QCF/releases/tag/0.2.0a
