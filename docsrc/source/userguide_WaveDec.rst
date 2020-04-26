=====================
User's guide: WaveDec
=====================

A typical command to run WaveDec looks as follows:

.. code-block:: bash

  $ WaveDec.py --config_file /path/to/config/file/MyConfig.yaml

As a command line option we specify the path and name of a configuration file. 
This configuration file contains all the settings that you may wish to specify, 
including input/output and processing parameters. The syntax of this configuration file is described below.

.. contents::
   :local:
   :backlinks: top

Input files
###########

.. _WaveDecConfigurationFile:

Configuration file
******************

As described in the :doc:`quick start guide <quickstart>`, few parameters can be passed via command line. However you will need to provide a configuration file to specify important processing options.

Remember that **all the parameters are optional**. If you are unsure about a specific processing option, do not set it and WaveDec will attempt to do its best.

The configuration file is written in `YAML <https://en.wikipedia.org/wiki/YAML>`__. The YAML syntax is simple and it is sufficient to look at the commented example below.



.. literalinclude:: AdditionalFiles/ExampleConfig.yaml
   :language: yaml

The example configuration file shown above is available :download:`here <AdditionalFiles/ExampleConfig.yaml>`.

The frequencies that are actually processed by WaveDec may be different from what you expected. WaveDec will always process discrete Fourier transform (DFT) frequencies. The frequency actually processed are chosen among the possible DFT frequencies in order to loosely match what is requested by the user. The possible DFT frequencies are given, for example, by the Python function `numpy.fft.fftfreq(Twindow/Ts, d=Ts) <https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html>`__.

.. warning:: The configuration file is case sensitive.

.. warning:: The syntax of the configration file is not checked properly. Many mistakes will go unoticed eg: ``F_min`` or ``fmin`` instead of ``Fmin``. Please check info printed when WaveDec starts or the log file to make sure the configuration file is read correctly.

Data
****

The program reads as input binary SAC files, both big- and little-endian. Certain fields of the SAC header must be set with the right values.

In the following table a list of the SAC header fields used by WaveDec.

==========   ====   =====================   ============================================================
Field name   Word   Value                   Comments
==========   ====   =====================   ============================================================
DELTA        0      Sampling interval [s]   | Must be the same in all SAC files of the provided input folder
                                            | The value found in the header can be overridden within the `Configuration file`_
NPTS         79     Number of samples       | Must be the same in all SAC files of the provided input folder
USER7        47     Easting [m]             | :math:`x` coordinate in meters in a Cartesian coordinate system
                                            | E.g., Swiss coordinate system, East--West
USER8        48     Northing [m]            | :math:`y` coordinate in meters in a Cartesian coordinate system
                                            | E.g., Swiss coordinate system, North-South
USER9        49     Elevation [m]           | :math:`z` coordinate, station elevation. Not in use.
KCMPNM       129    Component code          | Valid codes for the components are:
                                            | Vertical (:math:`u_z`): ``[Z, UDZ, BHZ, EHZ, HGZ, HHZ, EH1, EH4]``
                                            | East-West (:math:`u_x`): ``[E, EWX, BHE, EHE, HGE, HHE, EH3, EH6]``
                                            | North-South (:math:`u_y`): ``[N, NSY, BHN, EHN, HGN, HHN, EH2, EH5]``
                                            | Rotation around East-West (:math:`\omega_x`): ``[ROTX, RX]``
                                            | Rotation around North-South (:math:`\omega_y`): ``[ROTY, RY]``
                                            | Rotation around Vertical (:math:`\omega_z`): ``[ROTZ, RZ]``
==========   ====   =====================   ============================================================

About **polarity**. On the vertical component a positive sign of the signal implies an upward motion. 
On the east-west component an eastward motion. On the north-south component a northward motion.

More documentation on the SAC header is found `here <http://www.iris.edu/files/sac-manual/manual/file_format.html>`__.

Missing files (e.g., a sensor having recorded only one or two components) are not an issue. If the recording from a specific sensor/component is unusable, simply remove the corresponding file and leave the good ones in the ``INPUT`` folder.


Processing
##########

After reading the configuration file and the input files WaveDec will start processing the data. The configuration parameters in use are printed at statup. Check wheather your configuration file was read properly!

An estimate of the processing duration is printed and updated on the console. As soon as the processing of a signal window is completed, estimated parameters are written to file immediately in the ``OUTPUT`` folder. You may check and plot intermediate results at any time.

To terminate the processing prematurely press ``CTRL+C``.  Starting again the processing will overwrite the output files and start the processing from the beginning of the recording.

.. _LoveRayleighNoise:

Model selection: Love wave, Rayleigh wave, or nothing?
******************************************************

One strength of WaveDec is the ability to model multiple waves at the same time (See :ref:`Benefit of joint modeling Rayleigh and Love waves<secJointModeling>` for a motivation). This approach is called wavefield decomposition, hence WaveDec. The algorithm gradually increases the number of waves modeled up to a maximum number specified by ``MaxWaves``. The number of waves effectively modeled is in general different at different frequencies and in different time windows.

When decomposing the wavefield WaveDec has to choose whether to fit a Love wave, a Rayleigh wave or to stop adding waves. This is a model selection problem. The `Bayesian Information Criterion (BIC) <https://en.wikipedia.org/wiki/Bayesian_information_criterion>`_ is used to make this choice. In a nutshell, the BIC criterion selects the model with smallest BIC value. The BIC value is defined as

.. math::

   BIC = -2p(y|\theta) + N_p \ln(LK)\,,

where :math:`p(y|\theta)` is the model likelihood, :math:`N_p` is the number of model parameters, and :math:`LK` is the number of measurements. The quantity :math:`-2p(y|\theta)` describes the goodness of fit of our model. The quantity :math:`N_p \ln(LK)` is a penalty term related to the complexity of the model.

The algorithm iteratively consider the three following possibilities:

 * Should a Love wave be added to the current model? (Assuming that ``ModelLove: true``)
 * Should a Rayleigh wave be added to the current model? (Assuming that ``ModelRayleigh: true``)
 * Should we stop adding additional waves?  (Assuming that ``ModelNoise: true``) And/or was the maximum number of waves ``MaxWaves`` reached?

For each of these three possibilities, corresponding to three possible models of the wavefield, the BIC value is computed. The model corresponding to the smallest BIC value is chosen.

When doing model selection, we want to avoid underfitting or `overfitting <https://en.wikipedia.org/wiki/Overfitting>`__ the data.

Setting ``ModelNoise: false`` is not recommended as it may lead to overfitting the data. In fact, WaveDec is forced to model exactly ``MaxWaves`` waves. It will not be able to create a model with less waves.

When WaveDec is allowed to model fewer than ``MaxWaves`` (that is, with the flag ``ModelNoise: true`` --recommended choice) then the result will exhibit less outliers because most of the outliers have been rejected.

However, we observed that in certain circumstances this may lead underfitting and loosing information from weaker waves (e.g., higher modes, high frequencies, especially in active surveys). We therefore introduce another optional parameter enabling a more aggressive model selection strategy.

We introduce the parameter :math:`\gamma>0` to tune the result to achieve more or less complex models. It is controlled with the parameter ``Gamma`` in the configuration file. The new model selection strategy chooses the model with smallest

.. math::

   BIC_\gamma = -2p(y|\theta) + \gamma N_p \ln(LK)\,.

A small :math:`\gamma=0` reduces the penalty for a complex model. Model selection will be more aggressive and favour more complex models (i.e., more waves).

For :math:`\gamma=0` the ML criterion will be used for model selection. This will lead to greatly overfit the data. For :math:`\gamma=1` the standard BIC approach is followed. Setting :math:`\gamma<1` allows a more aggressive strategy compared to pure BIC, thus including more models from weaker waves but also more outliers. 

.. TIP::
   Our suggestion is to model jointly Love and Rayleigh wave (and the noise model) with the default BIC model selection criterion. Setting ``MaxWaves`` to 3 or 5 appears to work well.

.. TIP::
   If your output files are empty or present very few waves, try with a smaller ``Gamma`` (e.g. ``Gamma: 0.5``, ``Gamma: 0.1``...). In our experience, this was necessary at very few passive survey sites and at most active survey sites.

Speeding it up
**************

While it is recommended that you process multiple waves jointly, you may find it too time consuming. This section is about choosing processing options to obtain results faster. It may be a good idea for a quick preliminary analysis of your data.

* **Reduce the number of waves.** In the configuration file, set

  .. code-block:: yaml

    MaxWaves: 1

  so that a single wave is fitted for each time window / frequency. This prevents the software from iteratively estimate wave parameters. Still, for each time window / frequency both Love wave and Rayleigh wave fit are attempted and only the best one is kept.


* **Process different waves separately.** Instead of running one processing round modeling both Love waves and Rayleigh waves it may be faster to run two separate runs. You will need to create two separate configuration files.

  In the first run, only Love waves are analyzed and the configuration file looks like:

  .. code-block:: yaml

    ModelRayleighWaves: false
    ModelLoveWaves: true
    ModelNoise: true

  In the second run, only Rayleigh waves are analyzed:

  .. code-block:: yaml

    ModelRayleighWaves: true
    ModelLoveWaves: false
    ModelNoise: true


Output files
############

The output is saved to CSV files in the directory specified in the configuration file by ``OUTPUT``. The output files can be opened with a spreadsheet or with a simple text editor.

Estimated wavefield parameters
******************************

Output with the estimated parameters is written to disk after each time window is processed. Therefore you do not need to wait for the end of the elaboration to see partial results.

The output file for Rayleigh waves, ``RayleighWaves.csv``, looks like this:

.. literalinclude:: AdditionalFiles/RayleighWaves.csv
   :language: text
   :lines: 1-10

Comment lines begin with ``#``. On each row columns are separated by a tabulation character ``\t``.

Each column of the CSV file is described in the following table.

==================    =====================   =========================================================================
Column name           Unit                    Description
==================    =====================   =========================================================================
Frequency             Hertz                   | The frequency :math:`f` at which the processing was performed
                                              |
                                              | No matter what the frequency options provided in the configuration file are, the frequencies actually processed always correspond to Discrete Fourier transform (DFT) frequencies, i.e. :math:`f_n=n/T_s` where :math:`n` is a positive integer and :math:`T_s` is the sampling interval. 
Amplitude             Same unit as input      | The ML estimate of the wave amplitude :math:`\alpha`
                                              |
                                              | The estimated amplitude is always greater than zero :math:`\alpha>0`.
Wavenumber            1/meter                 | The ML estimate of the wavenumber :math:`\kappa`
                                              |
                                              | The estimated wavenumber is always greater than zero :math:`\kappa>0`.
Velocity              meter/second            | The ML estimate of the velocity of propagation, :math:`v=f/\kappa`
Azimuth               radian                  | The ML estimate of the azimuth :math:`\psi`
                                              |
                                              | It is a value in the interval :math:`\psi\in[0,2\pi]`. The azimuth is measured counterclockwise from the :math:`x` (east) axes.
                                              | E.g. :math:`\psi=0` for east, :math:`\psi=\pi/4` for north-east, :math:`\psi=\pi/2` for north.
                                              |
                                              | The wavevector is :math:`\kappa\,(\cos(\psi),\sin(\psi))` and points in the direction of propagation (not in the direction of arrival).
EllipticityAngle      radian                  | The ML estimate of the ellipticity angle :math:`\xi` of a Rayleigh wave
                                              |
                                              | It is a value in the interval :math:`\xi\in[-\pi/2,\pi/2]`. Positive for prograde particle motion, negative for retrograde motion.
                                              | The ellipticity angle is related to the H/V ratio as :math:`H/V=|\tan(\xi)|`.
==================    =====================   =========================================================================

A very similar output file is generated for Love waves (``LoveWaves.csv``).

.. caution:: The name of the output files should not be changed as it is needed by ``wdPicker.py`` when analysing the files.

Resolution limits
*****************

The file ``ArrayResolutionLimits.csv`` contains the resolutions limits of the array.

The study of resolution limits is an intricate matter. The resolution limits provided here should be considered as approximate and indicative only. Good results may show up outside the resolution limits, and bad results inside the limits. Be wise and use your expertise.

Let :math:`d_{\textrm{min}}` and :math:`d_{\textrm{max}}` be the smallest and the largest interstation distance within the array stations (:math:`d_{\textrm{max}}` is also known as array aperture).

The smallest and largest resolvable wavenumbers are defined (in 1/meter) as :math:`\kappa_{\textrm{min}}=1/d_{\textrm{max}}` and :math:`\kappa_{\textrm{max}}=0.5/d_{\textrm{min}}` (As suggested in Asten and Henstridge [1984]).

.. TIP::
   The resolution limits are computed as :math:`\kappa_{\textrm{max}}=c_1/d_{\textrm{min}}` and :math:`\kappa_{\textrm{min}}=c_2/d_{\textrm{max}}`. The value of :math:`c_1` and :math:`c_2` can be changed in the file ``wdSettings.py``.

Equivalently, in velocity at the frequency :math:`f` the limits are :math:`v_{\textrm{min}}=f/k_{\textrm{max}}` and :math:`v_{\textrm{max}}=f/k_{\textrm{min}}`.

The columns of ``ArrayResolutionLimits.csv`` have the following meaning

==================    =====================   =========================================================================
Column name           Unit                    Description
==================    =====================   =========================================================================
Frequency             Hertz                   | Frequency :math:`f`
Kmin                  1/meter                 | Smallest resolvable wavenumber :math:`\kappa_{\textrm{min}}`
Kmax                  1/meter                 | Largest resolvable wavenumber :math:`\kappa_{\textrm{max}}`
Vmin                  meter/second            | Smallest resolvable velocity :math:`v_{\textrm{min}}`
Vmax                  meter/second            | Largest resolvable velocity :math:`v_{\textrm{max}}`
==================    =====================   =========================================================================


Array Layout
************

The file ``ArrayLayout.csv`` contains the coordinates of the sensors.



The columns of ``ArrayLayout.csv`` have the following meaning

==================    =====================   =========================================================================
Column name           Unit                    Description
==================    =====================   =========================================================================
Easting               meter                   | 
Northing              meter                   | 
Elevation             meter                   | 
==================    =====================   =========================================================================

Log file
********

The file ``WaveDec.log`` is saved in the same directory specified by ``OUTPUT``. It contains all the output printed to screen from WaveDec.

Analysis of the output files
############################

The script ``wdPicker.py`` allows to visualize and perform some analysis of WaveDec output files. See :doc:`userguide_wdPicker` for details.
