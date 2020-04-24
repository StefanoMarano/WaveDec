# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
A configuration file for WaveDec, WaveDecActive and wdPicker
"""

###
### Constants used for calulation of resolution limits (Kmin and Kmax)
###
# for WaveDec
WD_C_1 = 0.5 # resolution_Kmax = WD_C_1 / dmin # As suggested in Asten and Henstridge [1984]
WD_C_2 = 1.0 # resolution_Kmin = WD_C_2 / dmax
# where dmin and dmax are the smallest and largest interstation distance
# for WaveDecActive
WDA_C_1 = 0.5 # resolution_Kmax = WDA_C_1 / d_min
WDA_C_2 = 0.5 # resolution_Kmin = WDA_C_2 / d_max
# where dmin and dmax are the smallest and largest interstation in the offset domain
# Please not that the values of resolution_Kmin and resolution_Kmax affect how the wavenumber domain grid is defined


###
### Some parameters used for producing plots with wdPicker.py
###
DEFAULT_SAVE_PLOT = True        # True: save plots to file
DEFAULT_SAVE_CSV = True         # True: save filtered and picked curves to comma separated value (CSV) file

DEFAULT_DPI=150                 # image resolution in dots per inch (DPI)
DEFAULT_IMAGE_FORMAT = 'png'    # image format.  Possible formats include: 'png', 'pdf', 'eps'
                                # See matplotlib documentation for more supported formats: 
                                # http://matplotlib.org/api/backend_bases_api.html


DEFAULT_DIR_PLOTS = 'plots'
DEFAULT_DIR_PICKED = 'filtered/picked'
DEFAULT_DIR_FILTERED = 'filtered'
DEFAULT_DIR_FILTERED_PLOTS = 'filtered/plots'


DEFAULT_USETEX = False      # True: uses LaTex font for labels. It may require the installation of additional packages
DEFAULT_FONT_SIZE_1 = 13    # Font size for xtick, ytick
DEFAULT_FONT_SIZE_2 = 14    # Font size for title, legend

DEFAULT_SWITCH_ELLIPTICITY = False # True when the vertical Z axes points downwards (e.g., SESAME synthetics)

DEFAULT_SCALE_WAVENUMBER = 0.25
DEFAULT_SCALE_ELLIPTICITY = 0.2
DEFAULT_SCALE_VELOCITY = 0.05
DEFAULT_SCALE_AZIMUTH = 0.2

DEFAULT_PLOT_ESTIMATE_DOT = True # True, plot a red dot for each estimate. False: do not plot.

DEFAULT_PLOT_FREQUENCY_LOG = False  # True: plot the frequency axis with log scale. False: linear scale
DEFAULT_PLOT_VELOCITY_LOG = False  # True: plot the velocity axis with log scale. False: linear scale


# number of points used in the search grid
#DEFAULT_NumE = 30       # Even, so that we have 0 in the final search grid
#DEFAULT_NumK = None     # If not None overrides default behaviour 2*ceil(Kmax/Kstep)
#                        # Even, so that we do not have 0 in the final search grid






###
### Do not ever change the constants below 
###

# Sensor component codes
EWX=1;      # Horizontal East-West
NSY=2;      # Horizontal North-South
UDZ=3;      # Vertical Up-Down
ROTX=4;     # Rotation around X
ROTY=5;     # Rotation around Y
ROTZ=6;     # Rotation around Z

# Wave model codes
MODEL_NOISE=1;
MODEL_VERTICAL=2;
MODEL_RAYLEIGH=3;
MODEL_LOVE=4;
MODEL_CIRCULAR_VERTICAL=5;
MODEL_CIRCULAR_RAYLEIGH=6;
MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH=7;
# MODEL_CIRCULAR_LOVE=7;

# the dictionary maps channel descriptor (as found, for example, in SAC file header)
# to the values used internally by this software
Components = {    'E':EWX,    'N':NSY,    'Z':UDZ,
                'EWX':EWX,  'NSY':NSY,  'UDZ':UDZ,
                'ROTX':ROTX,'ROTY':ROTY,'ROTZ':ROTZ,
                'RX':ROTX,'RY':ROTY,'RZ':ROTZ,
                 b'E':EWX,   b'N':NSY,   b'Z':UDZ,
               b'EWX':EWX, b'NSY':NSY, b'UDZ':UDZ,
               b'BHE':EWX, b'BHN':NSY, b'BHZ':UDZ,
               b'EHE':EWX, b'EHN':NSY, b'EHZ':UDZ,   # EH short period
               b'EH1':UDZ, b'EH2':NSY, b'EH3':EWX,   # Numbers should be used whenever the aligment is not known
               b'EH4':UDZ, b'EH5':NSY, b'EH6':EWX,   
               # # b'EH1':EWX, b'EH2':NSY, b'EHZ':UDZ, # Clashes with one line above! Used for borehole sensors. But X and Y are not actually known.
               b'HGE':EWX, b'HGN':NSY, b'HGZ':UDZ,
               b'HHE':EWX, b'HHN':NSY, b'HHZ':UDZ # Broadband
                }
