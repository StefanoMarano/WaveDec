A very quick start
==================

This page is about running WaveDec with the smallest amount of configuration. For production, you may wish to look into the :doc:`user's guide <userguide_WaveDec>` or :doc:`user's guide <userguide_WaveDecActive>`. For installation, you should look :doc:`here <installation>`.

With no options
###############

From the command line, simply type ``WaveDec.py`` and the program will start

.. code-block:: bash

  $ WaveDec.py

The program will run using default options.

* Input files (e.g., SAC files) are sought in the current directory.

* Output and log files will be save in the current directory.

* WaveDec will find suitable parameters for the processing on its own. But if a configuration file named ``config.yaml`` is present in the current directory, this file will be read.


With few command line options
#############################

Certain options can be set from the command line. Description of these options can be found typing ``WaveDec.py --help``. The default values are shown in brackets ``[]``.

.. code-block:: bash

  $ WaveDec.py --help
  usage: WaveDec.py [-h] [--config_file CONFIG_FILE] [--input INPUT] [--output OUTPUT] [--verbosity VERBOSITY]

  WaveDec - a tool for seismic wavefield analysis

  optional arguments:
    -h, --help            show this help message and exit
    --config_file CONFIG_FILE
                          path to YAML configuration file [config.yaml]
    --input INPUT         path to input folder (or file) [./]
    --output OUTPUT       path to output folder [./]
    --verbosity VERBOSITY
                          increase output verbosity: 0 = only warnings, 1 = info, 2 = debug. [0]

For example we may wish to specify the folder containing the SAC files to be analyzed.

.. code-block:: bash

  $ WaveDec.py --input /path/to/sac/files


More options
############

If you need to set more options please look into the :doc:`user's guide <userguide_WaveDec>` or :doc:`user's guide <userguide_WaveDecActive>`.


