============
Installation
============

This page describes the steps necessary to install WaveDec.


Requirements
############

WaveDec needs the Python interpreter, version 3, and few additional Python packages.

If you wish :doc:`wdPicker <userguide_wdPicker>` to produce plots with Latex fonts, you will need to configure Latex environment appropriately.

Download
########

The code can be cloned from the `repository <https://github.com/StefanoMarano/WaveDec>`__ ::

  git clone https://github.com/StefanoMarano/WaveDec.git

or downloaded as `ZIP archive <https://github.com/StefanoMarano/WaveDec/archive/master.zip>`__.


GNU/Linux (Ubuntu >14.04)
*************************

* Install Python 3 and the Python package installer by running::

    sudo apt-get install python3 python3-pip

* The following Python libraries are needed: ``glob, os, errno, sys, yaml, time, csv, logging, argparse, numpy, scipy, matplotlib``

  On Ubuntu many packages are installed with the default Python 3 installation. You may need to run::

    sudo apt-get install python3-yaml python3-numpy python3-scipy python3-matplotlib

  Alternatively, you may need to try::

    sudo pip3 install <package_name>

GNU/Linux (Ubuntu <14.04)
*************************

See the instructions for Ubuntu 14.04 above and keep in mind the following:

* You will need to upgrade both your ``numpy`` and ``scipy`` installations to a recent version. Download them from `here <http://www.scipy.org/scipylib/download.html>`__ and read installation instructions `here <http://www.scipy.org/scipylib/building/linux.html>`__.
  
* On some older Ubuntu versions the ``python3-pip`` may be unavailable. See this `thread <http://askubuntu.com/questions/412178/how-to-install-pip-for-python-3-in-ubuntu-12-04-lts>`_.



Windows
*******

There are several Python distributions including many key Python packages. See `here <http://www.scipy.org/install.html#scientific-python-distributions>`__ for a list.

Consider installing WinPython with the following steps:

  * Download and install WinPython, version 3, from `here <http://winpython.sourceforge.net/>`__. This distribution already includes the ScipPy stack.

  * Download the ``zip`` archive with a recent version of PyYAML from `here <http://pyyaml.org/wiki/PyYAML>`__. Install it using the WinPython control panel (add package, install).

  * Open the WinPython command prompt. You will need to specify the full path to a script in order to execute it. For example

  .. code-block:: bash

    python C:\WaveDec\bin\WaveDec.py


Other platforms
***************

  * Download Python 3 `here <https://www.python.org/downloads/>`__.

  * Install recent versions of Numpy, Scipy, and Matplotlib. See installation of the full ScipPy stack `here <http://www.scipy.org/install.html/>`__.

  * Install `YAML <http://pyyaml.org>`_ and any other missing package.


Download
########

* Download the WaveDec files in a compressed ZIP archive from one of the following links.

* `Download WaveDec <http://mercalli.ethz.ch/~marra/WaveDec/WaveDec.zip>`_
* `Download WaveDec and WaveDecActive <http://mercalli.ethz.ch/~marra/WaveDec/WaveDecActive.zip>`_
* `Download WaveDec and WaveDecActive without documentation <http://mercalli.ethz.ch/~marra/WaveDec/WaveDecActive_nodoc.zip>`_


Installation
############

GNU/Linux (Ubuntu)
************************

* Download the WaveDec files and extract the files from the compressed archive.

* Add the folder ``WaveDec/bin`` to your search path

  * Specify where the binaries file are located by adding the following line

    .. code-block:: bash

      PATH=$PATH:/full/path/to/WaveDec/bin

    at the end of the ``.bash_profile`` file (if you are using Bash shell) or at the end of the ``.profile`` file. Both files are located in your home folder. This change will take effect at the next login.

    To make the change immediate (and not permanent) type in your shell

    .. code-block:: bash

      export PATH="$PATH:/full/path/to/WaveDec/bin"


* Make sure that ``Wavedec.py`` has execution privileges. On GNU/Linux type ``chmod u+x WaveDec.py``.

* Check whether your system is configured correctly by typing ``WaveDec.py`` at the shell prompt

  .. code-block:: bash

    $ WaveDec.py
    No configuration file (None) found. Proceeding with default values.
    No suitable input found in '/home/yourdir'

* To run WaveDec see the :doc:`user's guide <userguide_WaveDec>`.


Windows
*******

* Download the WaveDec files and extract the files from the compressed archive to the folder ``c:\WaveDec\``.

* Check whether your system is configured correctly by typing ``python c:\WaveDec\bin\WaveDec.py`` at the shell prompt

  .. code-block:: bash

    > python c:\WaveDec\bin\WaveDec.py
    No configuration file (None) found. Proceeding with default values.
    No suitable input found in 'c:\WaveDec\'

* To run WaveDec see the :doc:`user's guide <userguide_WaveDec>`.

Troubleshooting
###############

* We can check whether Python 3 is properly installed by running ``python3`` from the command line. In a working Python 3 installation, Python will start in interactive mode. It should produce an output similar to:

  .. code-block:: bash

    $ python3
    Python 3.4.0 (default, Apr 11 2014, 13:05:11) 
    [GCC 4.8.2] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

* To check whether the Python 3 packages are properly installed, let us import them as follows

  .. code-block:: python

    >>> import glob, os, errno, sys, yaml, time, csv, logging, argparse, numpy, scipy, matplotlib
    >>>

  No errors, the packages are properly installed on this system. When the packages are not correctly installed, an error is returned

  .. code-block:: python

    >>> import nonexistingpackage
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: No module named 'nonexistingpackage'

  Leave the interactive mode with ``exit()``

* Errors like the following

  .. code-block:: bash

    marra@bigstar01:~/WaveDec/bin$ ./WaveDec.py 
    Traceback (most recent call last):
      File "./WaveDec.py", line 12, in <module>
        from EstimationRoutines import *
      File "/home/marra/WaveDec/bin/EstimationRoutines.py", line 8, in <module>
        from scipy.optimize import minimize
    ImportError: cannot import name minimize

  will arise if the version of the installed packages ``numpy`` or ``scipy`` are too old. To check the currently installed versions type

  .. code-block:: python

    >>> import numpy,scipy
    >>> print(scipy.__version__)
    0.9.0
    >>> print(numpy.__version__)
    1.6.1

  The most recent versions can be downloaded from `here <http://www.scipy.org/scipylib/download.html>`__ and read installation instructions `here <http://www.scipy.org/scipylib/building/linux.html>`__.
