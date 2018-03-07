|Build Status| |Documentation|

QCF: Estimating Quartet Concordance Factors
===========================================

`Read the Docs <http://quartet-cf.rtfd.io/>`__
----------------------------------------------

Installation
~~~~~~~~~~~~

QCF can be installed by cloning this repo and then compiling the main executable
using the provided ``Makefile`` (see code below). QCF is written in C++ for Unix-like
operating systems and makes use of features from the C++11 standard, which
requires a compatible compiler (GNU ``g++`` >= 4.8, Clang ``clang++`` >= 3.3).
Mac users will also need to have the Xcode
`Command Line Tools <http://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/>`__ installed.

.. code:: bash

  # Cloning, compiling, testing, and installing QCF
  git clone https://github.com/pblischak/QCF.git      # 1. Clone the repo from GitHub
  cd QCF                                              # 2. cd into the QCF/ folder
  make                                                # 3. compile the qcf executable
  make test                                           # 4. test that executable works
  sudo make install                                   # 5. copy executable to /usr/local/bin

Usage
~~~~~

After compiling and installing the software, you can use the program in a terminal
window by typing ``qcf``. Below are the main commands that can be used to estimate
quartet concordance factors using QCF. More details on input data formats, as well
as background on interpreting the output, can be found in the `tutorial <>`__ on
the ReadTheDocs site.

.. code:: bash

  # Run the data in the example folder
  cd example
  qcf -i genes.txt -m map.txt --prefix example

  # Run with 100 bootstrap reps
  qcf -i genes.txt -m mapt.txt -b 100 --prefix bootstrap

Getting Help
~~~~~~~~~~~~

If you have questions about running QCF, please feel free to use the gitter chatroom to get help:

|Gitter|

If you have a problem while running QCF, and you think it may be a bug, please consider filing an issue:

|QCF Issues|

.. |Build Status| image:: https://travis-ci.com/pblischak/QCF.svg?token=3UtCuy4QMGzzqmrdSwV2&branch=master
   :target: https://travis-ci.com/pblischak/QCF

.. |Documentation| image:: http://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://quartet-cf.readthedocs.io

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/pblischak-QCF/Lobby

.. |QCF Issues| image:: https://img.shields.io/badge/QCF-Issues-blue.svg
   :target: https://github.com/pblischak/QCF/issues
