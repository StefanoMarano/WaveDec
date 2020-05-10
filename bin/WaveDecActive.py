#!/usr/bin/python3
# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################


import glob, os, errno, sys
import yaml
import time, datetime
import logging
import argparse
from numpy import shape, zeros, min, max, any, linspace, diff, logspace, pi, array, fft, ceil, unique, log10, abs, sum
from scipy.spatial.distance import cdist



from CircularEstimationRoutines import circularDecomposeWavefield
from SyntheticWavefield import readSyntheticWavefield
from ReadSAC import readSacDir
import DataUtils as db

from PlotUtils import plotArray

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from wdSettings import MODEL_NOISE, MODEL_VERTICAL, MODEL_RAYLEIGH, MODEL_LOVE, MODEL_CIRCULAR_VERTICAL, MODEL_CIRCULAR_RAYLEIGH, MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH
from wdSettings import EWX, NSY, UDZ, WDA_C_1, WDA_C_2


def main():
    
    
    ### Initialization: commandline args, config file, and other parameters
    parser = argparse.ArgumentParser(description='WaveDec - a tool for seismic wavefield analysis')
    parser.add_argument("--config_file", default='config.yaml', help="path to YAML configuration file [config.yaml]")
    parser.add_argument("--input",  help="path to input folder (or file) [./]")
    parser.add_argument("--output", help="path to output folder [./]")
    parser.add_argument('--plot', action='store_true', default=False, help="Display plots about array, traces, and shots. No processing is performed.")
    parser.add_argument("--verbosity", type = int, default=1, help="increase output verbosity: 0 = only warnings, 1 = info, 2 = debug.  [0]")
    args = parser.parse_args()

    # Setup console logging, file logging later
    #logging.basicConfig(level=logging.DEBUG,    format='%(levelname)-8s %(message)s', stream=sys.stderr)
    logging.basicConfig(level=logging.DEBUG,    format='%(message)s', stream=sys.stderr)

    if args.verbosity == 0:
        logging.getLogger('').setLevel(logging.WARN) 
    elif args.verbosity == 1:
        logging.getLogger('').setLevel(logging.INFO) 
    elif args.verbosity == 2:
        logging.getLogger('').setLevel(logging.DEBUG)
    else:
        logging.error('Unrecognized verbosity level {0}'.format(args.verbosity))
        return
    #TODO logging level should/could be chosen in config file    
        #    logging.error('an error message')
        #    logging.info('an info message')
        #    logging.debug('a debug message')
    
    
    
    
    # parameter CONIFG_FILE (configuration file in YAML format)
    CONFIG_FILE = os.path.abspath(args.config_file)    # default is "./config.yaml"
    try:
        with open(CONFIG_FILE) as f:
            conf = yaml.load(f, Loader=yaml.FullLoader)
            f.close()
    except yaml.scanner.ScannerError as e:
        logging.critical('There is a syntax problem with the YAML configuration file {0}.'.format(CONFIG_FILE))
        logging.critical(type(e))
        logging.critical(e)
        return
    except FileNotFoundError as e:        
        logging.critical('No configuration file ({0}) found.'.format(CONFIG_FILE))
        logging.critical('WaveDecActive needs a configuration file. See documentation.')
        # logging.critical(type(e))
        # logging.critical(e)
        return
    except Exception as e:
        logging.critical('There was a problem while loading the YAML configuration file {0}.'.format(CONFIG_FILE))
        logging.critical(type(e))
        logging.critical(e)
        return

    if conf is None: # when we read an empty config file
        conf = dict()
    
    # parameter INPUT
    INPUT_args = args.input # the value of INPUT as parsed from commandline
    INPUT_conf = conf.get('INPUT', None) # the value of INPUT as read from CONFIG_FILE
    if (INPUT_args == None) and (INPUT_conf == None):         INPUT = '.'
    elif (INPUT_args == None) and (INPUT_conf != None):        INPUT = INPUT_conf
    elif (INPUT_args != None) and (INPUT_conf == None):        INPUT = INPUT_args
    else:
        logging.critical('Option clash for \'INPUT\'. Defined both as commandline argument and in config file \'{0}\''.format(CONFIG_FILE))
        return
    INPUT = os.path.abspath(INPUT)
    
    # parameter OUTPUT    
    OUTPUT_args = args.output # the value of OUTPUT as parsed from commandline
    OUTPUT_conf = conf.get('OUTPUT', None) # the value of OUTPUT as read from CONFIG_FILE
    if (OUTPUT_args == None) and (OUTPUT_conf == None):        OUTPUT = '.'
    elif (OUTPUT_args == None) and (OUTPUT_conf != None):    OUTPUT = OUTPUT_conf
    elif (OUTPUT_args != None) and (OUTPUT_conf == None):    OUTPUT = OUTPUT_args
    else:
        logging.critical('Option clash for \'OUTPUT\'. Defined both as commandline argument and in config file \'{0}\''.format(CONFIG_FILE))
        return
    OUTPUT = os.path.abspath(OUTPUT)
    
   
    # create output folder and make sure it is writable
    # TODO ask wether to overwrite output?
    try:
        os.makedirs(OUTPUT)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise    
    
    # Read input, either synthetic wavefield or files from directory
    if os.path.exists(INPUT) == False:
        logging.critical('INPUT folder or file not found ({0})'.format(INPUT))

    # logging to file        
    logFile = logging.FileHandler(filename=os.path.join(OUTPUT, 'WaveDecActive.log'), mode='w')
    logFile.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', "%Y-%m-%d %H:%M:%S")
    logFile.setFormatter(formatter) # tell the handler to use this format
    logging.getLogger('').addHandler(logFile) # add the handler to the root logger
    

    # Read input files
    fileList = None 
    try:
        if os.path.exists(INPUT) == False:
            raise Exception('INPUT folder or file not found ({0})'.format(INPUT))
        if INPUT.endswith('.yaml'):
            (y, ArrayInfo, Ts, _) = readSyntheticWavefield(INPUT)
        elif len(glob.glob(os.path.join(INPUT,'*.sac'))) > 0:        
            (y, ArrayInfo, Ts, fileList) = readSacDir(INPUT)           
        elif len(glob.glob(os.path.join(INPUT,'*.saf'))) > 0:
            raise Exception('Reading of SAF files not yet implemented.') # TODO read also SAF format
        else:
            raise Exception('No suitable input found in \'{0}\''.format(INPUT))
    except Exception as e:
        logging.critical(e)
        return 

    K = shape(y)[0] # number of samples        
    L = shape(y)[1] # number of channels
    Tobservation = Ts*K
    t_axis = linspace(0,Ts*(K-1),K)

    # some quality control on input data
    exitFlag = False
    for ll in range(0, L):
        if sum(abs(y[:, ll])) == 0:
            logging.critical("Channel {0} is empty".format(ll))
            exitFlag = True
            if fileList is not None:
                logging.critical("\t{0}".format(fileList[ll]))
                
    if exitFlag:
        logging.critical("There were some errors in the input data. Exiting.")
        return

    try:
        SourcePosition = conf.get('SourcePosition', None)
        if SourcePosition is None:
            logging.critical('Parameter SourcePosition must be provided in configuration file {0}'.format(CONFIG_FILE))
            return
    except Exception as e:
        logging.critical('There was a problem while parsing the configuration file {0}'.format(CONFIG_FILE))
        logging.critical(type(e))
        logging.critical(e)
        return
    # sensor positions

        # compute resolution limits from array
    pos = array(list(set(tuple(p) for p in ArrayInfo[:,0:3])))
    cmp = ArrayInfo[:,3]
    Nsensors = shape(pos)[0]
    r_dist = unique(cdist(ArrayInfo[:,0:2], array([SourcePosition[0:2]]))[:,0]) # distance of each station from source
    r_min = min(r_dist)         # distance of closest sensor
    r_max = max(r_dist)         # distance of furthermost sensor
    d_max = r_max - r_min       # maximum aperture, offset domain
    d_min = min(diff(r_dist))   # minimum interdistance, offset domain

    
    if (r_min != 0.0) and (r_max != 0.0):
        resolution_Kmax = WDA_C_1 / d_min
        resolution_Kmin = WDA_C_2 / d_max
        KrStep = resolution_Kmin * 0.3 # step size (in wavenumber) for grid search
    else:
        logging.critical("Single sensor analysis, exiting")
        return    

    
    # parsing YAML configuration file    
    try:
        # Several processing parameters
        ShotStarts = conf.get('ShotStarts', None) # this is an array with the start time of each shot
        if ShotStarts is None: ShotStarts = []
        if not isinstance(ShotStarts, list): ShotStarts = [ShotStarts]
        ShotDuration = conf.get('ShotDuration', None) # this is the duration of the shot waveform. It is shared across multiple shots.s
        if ShotDuration is not None: ShotDuration = float(ShotDuration)
        MaxWaves = int(conf.get('MaxWaves', 3)) # maximum number of waves in model
        MaxIterations = int(conf.get('MaxIterations', 10)) # maximum number of iterations when refining estimates
        Ts_tmp = conf.get('Ts', None)
        if Ts_tmp is not None: # this is overriding what found SAC header
            Ts = float(Ts_tmp)
            logging.warning("Overriding sampling interval found in SAC header Ts: {0:.2e}".format(Ts))

        Fmin = float(conf.get('Fmin', 0))
        Fmax = float(conf.get('Fmax',min([20, 0.5/Ts])))
        Fspacing = conf.get('Fspacing','lin')  # 'lin' or 'log'
        Fnum = int(conf.get('Fnum',50))
        Fnum = int(max([Fnum, 1])) # TODO Fnum >= 1 or Fnum = 0 to model all freqs
        Fnyquist = 0.5/Ts 
        if Fmax > Fnyquist:
            Fmax = Fnyquist
            logging.warning("Cannot analyse frequencies higher than Nyquist frequency. Setting Fmax to {0:.2f}.".format(Fmax))
        Estep = float(conf.get('Estep', pi/90.0))
        Vmin = float(conf.get('Vmin', 50))
        Gamma = float(conf.get('Gamma', 1))
        if Gamma < 0:
            logging.critical("Parameter Gamma should be greater than or equal to 0.")
            logging.critical("Gamma: {0}".format(Gamma))
            return
        if Nsensors > 1:
            Kmax = float(conf.get('Kmax', 1.2*resolution_Kmax))
            KrStep = float(conf.get('Kstep', resolution_Kmin)) # step size (in wavenumber) for grid search
            KrStep = float(conf.get('KrStep', resolution_Kmin))
            KiStep = float(conf.get('KiStep', 0.01)) 
        elif Nsensors == 1:
            logging.critical("Single sensor analysis not supported.")
            return
        if ShotDuration is not None and ShotStarts is not None and len(ShotStarts) > 0 and Tobservation < ShotDuration+max(ShotStarts):
            logging.critical("Input signal too short. Or ShotDuration and ShotStarts are improperly set.")
            logging.critical("\tSignal duration {0:.2f} [s] ({1} samples, {2:.2e} [s] sampling interval)".format(Tobservation, K, Ts))
            logging.critical("\tShotDuration: {0:.2f} [s]".format(ShotDuration))
            logging.critical("\tShotStarts: {0} [s]".format(ShotStarts))
            return

        # Which wave models to fit?
        ModelCircularRayleighWaves = db.str2bool(conf.get('ModelCircularRayleighWaves', False))
        ModelCircularDissipativeRayleighWaves = db.str2bool(conf.get('ModelCircularDissipativeRayleighWaves', False))
        ModelNoise = db.str2bool(conf.get('ModelNoise', True))
        ModelCircularVerticalWaves = db.str2bool(conf.get('ModelCircularVerticalWaves', False))
    except Exception as e:
        logging.critical('There was a problem while parsing the configuration file {0}'.format(CONFIG_FILE))
        logging.critical(type(e))
        logging.critical(e)
        return
        
    
    

    
    # Plotting, this may help to find the ShotStarts visually
    if args.plot or ShotDuration is None or len(ShotStarts) == 0:
        logging.info("Plotting traces. You may need these plots to figure out ShotDuration and ShotStarts.")
        plt.ion()
        # Plot array
        plotArray(pos, SourcePosition)

        # Plot entire trace
        fig = plt.figure(figsize=(9,9)) # figsize=(6,9))
        gs = gridspec.GridSpec(3, 1)
        gs.update(wspace=0.05, hspace=0.0, left = 0.15, right = 0.95, bottom = 0.1, top = 0.95)
        ax0 = plt.subplot(gs[0])
        if sum(cmp == UDZ) > 0:
            plt.plot(t_axis, y[:,cmp == UDZ])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.ylabel('UDZ')
        plt.tick_params(labelbottom='off')
        ax1 = plt.subplot(gs[1], sharex=ax0)
        if sum(cmp == EWX) > 0:
            plt.plot(t_axis, y[:,cmp == EWX])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.tick_params(labelbottom='off')
        plt.ylabel('EWX')
        ax2 = plt.subplot(gs[2], sharex=ax0)
        if sum(cmp == NSY) > 0:
            plt.plot(t_axis, y[:,cmp == NSY])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.ylabel('NSY')
        plt.xlabel('[s]')
        #gs.tight_layout(fig)
        plt.show()
        
        plt.ioff()
        # Plot 
        fig = plt.figure(figsize=(9,9))
        gs = gridspec.GridSpec(Nsensors, 1)
        gs.update(wspace=0.05, hspace=0.0, left = 0.15, right = 0.95, bottom = 0.1, top = 0.95)
        for nn in range(0, Nsensors):
            plt.subplot(gs[nn])
            # ndx: this is used to check if the channel is actually present in the input data
            ndx = (cmp == EWX) & (ArrayInfo[:,0] == pos[nn][0]) & ( ArrayInfo[:,1] == pos[nn][1]) & (ArrayInfo[:,2] == pos[nn][2])
            if any(ndx):
                plt.plot(t_axis, y[:, ndx], 'b')
            ndx = (cmp == NSY) & (ArrayInfo[:,0] == pos[nn][0]) & ( ArrayInfo[:,1] == pos[nn][1]) & (ArrayInfo[:,2] == pos[nn][2])
            if any(ndx):
                plt.plot(t_axis, y[:, ndx], 'r')
            ndx = (cmp == UDZ) & (ArrayInfo[:,0] == pos[nn][0]) & ( ArrayInfo[:,1] == pos[nn][1]) & (ArrayInfo[:,2] == pos[nn][2])
            if any(ndx):
                plt.plot(t_axis, y[:, ndx], 'g')
            plt.xlim(t_axis[0], t_axis[-1])
            plt.tick_params(labelbottom='off')
            plt.ylabel(nn)
        plt.tick_params(labelbottom='on')            
        plt.xlim(t_axis[0], t_axis[-1])
        plt.xlabel('[s]')
        plt.show()
    
    
    
    if ShotDuration is None:
        logging.critical('Variable "ShotDuration" is not defined. Check configuration file {0}'.format(CONFIG_FILE))
        return
    if ShotStarts is None:
        logging.critical('Variable "ShotStarts" is not defined. Check configuration file {0}'.format(CONFIG_FILE))
        return 
    if len(ShotStarts) == 0:
        logging.critical('Check definition of variable "ShotStarts". Check configuration file {0}'.format(CONFIG_FILE))
        return
    
    ShotList = []
    for s in ShotStarts:
        ShotList.append({'ShotStart': float(s), 'ShotEnd': float(s)+ShotDuration, 'SourcePosition':SourcePosition})
    Nshots = len(ShotList)
    
  
    if args.plot:
        for ss in range(0, Nshots):
            Shot = ShotList[ss]
            SourcePosition = Shot['SourcePosition']
            ShotStart = Shot['ShotStart']
            ShotEnd = Shot['ShotEnd']
            Kstart = int(round(ShotStart/Ts))
            Kend = int(round(ShotEnd/Ts))
            
            Kw = Kend-Kstart
           
            fig = plt.figure(figsize=(9,9)) # figsize=(6,9))
            gs = gridspec.GridSpec(3, 1)
            gs.update(wspace=0.05, hspace=0.0, left = 0.15, right = 0.95, bottom = 0.1, top = 0.95)
            ax0 = plt.subplot(gs[0])
            plt.title('Shot {0}/{1}'.format(ss+1, Nshots))
            if sum(cmp == UDZ) > 0:
                plt.plot(t_axis[Kstart:Kend], y[Kstart:Kend, cmp == UDZ])
            plt.xlim(t_axis[Kstart], t_axis[Kend-1])
            plt.tick_params(labelbottom='off')
            plt.ylabel('UDZ')
            ax1 = plt.subplot(gs[1], sharex=ax0)
            if sum(cmp == EWX) > 0:
                plt.plot(t_axis[Kstart:Kend], y[Kstart:Kend, cmp == EWX])
            plt.xlim(t_axis[Kstart], t_axis[Kend-1])
            plt.tick_params(labelbottom='off')
            plt.ylabel('EWX')
            ax2 = plt.subplot(gs[2], sharex=ax0)
            if sum(cmp == NSY) > 0:
                plt.plot(t_axis[Kstart:Kend], y[Kstart:Kend, cmp == NSY])
            plt.xlim(t_axis[Kstart], t_axis[Kend-1])
            plt.ylabel('NSY')
            plt.xlabel('[s]')
            plt.show()
        sys.exit()
            
   


    
    # Does not make sense to model Vertical and Rayleigh waves together
    if ModelCircularVerticalWaves and (ModelCircularRayleighWaves or ModelCircularDissipativeRayleighWaves):
        logging.warning("Parameter ModelCircularVerticalWaves is true, setting ModelCircularRayleighWaves and ModelCircularDissipativeRayleighWaves to false")
        ModelCircularRayleighWaves = False
        ModelCircularDissipativeRayleighWaves = False
    
    if ModelCircularRayleighWaves == ModelCircularDissipativeRayleighWaves == True:
        logging.critical("Either ModelCircularRayleighWaves or ModelCircularDissipativeRayleighWaves must be true, not both.")
        sys.exit()
    if ModelCircularVerticalWaves == ModelCircularRayleighWaves == ModelCircularDissipativeRayleighWaves == False: 
        logging.critical("Please choose a wave to model and set it to true.")
        sys.exit()

    # Check whether the needed sensor components are present in input data
    comp=ArrayInfo[:,3]
    if ModelCircularVerticalWaves and any(comp == UDZ) == False:
        logging.critical("ModelCircularVerticalWaves is true, but no vertical components (UDZ) are found in input data")
        sys.exit()
    if ModelCircularRayleighWaves and any(comp == UDZ) == False:
        logging.critical("ModelCircularRayleighWaves is true, but no vertical components (UDZ) are found in input data")
        sys.exit()
    if ModelCircularRayleighWaves and any(comp == EWX) == False:
        logging.critical("ModelCircularRayleighWaves is true, but no east-west (EWX) components are found in input data")
        sys.exit()
    if ModelCircularRayleighWaves and any(comp == NSY) == False:
        logging.critical("ModelCircularRayleighWaves is true, but no north-south (NSY) components are found in input data")
        sys.exit()
    if ModelCircularDissipativeRayleighWaves and any(comp == UDZ) == False:
        logging.critical("ModelCircularDissipativeRayleighWaves is true, but no vertical components (UDZ) are found in input data")
        sys.exit()
    if ModelCircularDissipativeRayleighWaves and any(comp == EWX) == False:
        logging.critical("ModelCircularDissipativeRayleighWaves is true, but no east-west (EWX) components are found in input data")
        sys.exit()
    if ModelCircularDissipativeRayleighWaves and any(comp == NSY) == False:
        logging.critical("ModelCircularDissipativeRayleighWaves is true, but no north-south (NSY) components are found in input data")
        sys.exit()
    
    WavesToModel = {MODEL_CIRCULAR_RAYLEIGH:ModelCircularRayleighWaves, MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH:ModelCircularDissipativeRayleighWaves,
                    MODEL_NOISE:ModelNoise, MODEL_CIRCULAR_VERTICAL:ModelCircularVerticalWaves,
                    MODEL_VERTICAL:False, MODEL_RAYLEIGH:False, MODEL_LOVE:False}
    
    
    # initialize database for storing results
    conn = db.init()
    db.setArrayInfo(conn, ArrayInfo)


    
    Kw = int(ceil(ShotDuration/Ts)) # average samples per shot, used to determine frequencies
    Fvec_fft = fft.fftfreq(Kw, Ts); # array of DFT frequencies
    if (Fspacing.lower() == 'lin') or (Fspacing.lower() == 'linear'):    Fvec_in = linspace(Fmin, Fmax, Fnum)
    elif Fspacing.lower() == 'log':        Fvec_in = logspace(log10(Fmin), log10(Fmax), Fnum)
    else:
        logging.critical("Unrecognized option for frequency spacing ({0})".format(Fspacing))
        return
    
    Fvec=zeros(Fnum) # this will store the frequencies to process
    Fvec_ndx=array([0] * Fnum) # an array of integers, contains the index of Fvec_fft
    
    for ff in range(0,Fnum): # the frequencies actually processed are DFT frequencies
        Freq = Fvec_in[ff]
        ndx = (abs(Fvec_fft-Freq)).argmin()
        Fvec[ff] = Fvec_fft[ndx]
        Fvec_ndx[ff] = ndx
        
    (Fvec, ndx) = unique(Fvec, return_index=True) # make sure Fvec elements are unique    
    Fvec_ndx = Fvec_ndx[ndx]
    ndx = Fvec > 0 # Do not process 0 Hz
    Fvec_ndx = Fvec_ndx[ndx]
    Fvec = Fvec[ndx]
    Fnum = len(Fvec)
    
    cur = conn.cursor()
    for ff in range(0, Fnum):
         cur.execute('''INSERT INTO Frequencies (Fndx, F) VALUES ({0},{1})'''.format(Fvec_ndx[ff], Fvec[ff]))
    conn.commit()
    
   
    # prepare output files
    db.createOutputFiles(conn, OUTPUT, WavesToModel, pos, Fvec, resolution_Kmin, resolution_Kmax, SourcePosition)

    ### Print startup info
    logging.info("")
    logging.info(" _       __                 ____            ___        __  _")
    logging.info("| |     / /___ __   _____  / __ \___  _____/   | _____/ /_(_)   _____")
    logging.info("| | /| / / __ `/ | / / _ \/ / / / _ \/ ___/ /| |/ ___/ __/ / | / / _ \ ")
    logging.info("| |/ |/ / /_/ /| |/ /  __/ /_/ /  __/ /__/ ___ / /__/ /_/ /| |/ /  __/")
    logging.info("|__/|__/\__,_/ |___/\___/_____/\___/\___/_/  |_\___/\__/_/ |___/\___/")
    logging.info("")
    logging.info("© 2017 ETH Zurich, Swiss Seismological Service")
    logging.info("")
    logging.info("")
    logging.info("Input/output info")
    logging.info("\tRead data from")
    logging.info("\t\t{0}".format(INPUT))
    logging.info("\tWith configuration file")
    logging.info("\t\t{0}".format(CONFIG_FILE))
    logging.info("\tOutput will be stored in")
    logging.info("\t\t{0}".format(OUTPUT))
    logging.info("Data info")
    logging.info("\t{0} sensors and {1} channels".format(Nsensors, L))
    logging.info("\tSignal duration {0:.2f} [s] ({1} samples, {2:.2e} [s] sampling interval)".format(Tobservation, K, Ts))
    logging.info("\tSource position (x,y)=({0},{1}) [m]".format(SourcePosition[0], SourcePosition[1]))
    logging.info("\tDistance of closest and farthest sensor: {0:.2f} and {1:.2f} [m], respectively".format(r_min, r_max))
    logging.info("\tSmallest and largest intersensor distance (in offset domain): {0:.2f} and {1:.2f} [m], respectively".format(d_min, d_max))
    logging.info("\tResolution limits, inferred from array geometry,  Kmin={0:.2e} Kmax={1:.2e} [1/m] ".format(resolution_Kmin, resolution_Kmax))
    logging.info("Processing info")
    logging.info("\t{0:.0f} shots in total".format(Nshots))
    logging.info("\tProcessing {0} frequencies between {1:.2f} [Hz] and {2:.2f} [Hz] (with {3} spacing)".format(len(Fvec), Fvec[0], Fvec[-1], Fspacing))
    if Vmin > 0.0: logging.info("\tAnalyzing velocities greater than {0:.2e} [m/s]".format(Vmin))
    logging.info("\tAnalyzing wavenumbers smaller than {0:.2e} [1/m]".format(Kmax))
    logging.info("\tGrid step in wavenumber {0:.2e} [1/m]".format(KrStep))
    if ModelCircularDissipativeRayleighWaves: logging.info("\tGrid step for imaginary wavenumber {0:.2e} [1/m]".format(KiStep))
    logging.info("\tGrid step in ellipticity angle {0:.2e} [rad]".format(Estep))    
    logging.info("Model selection")
    if ModelNoise and Gamma > 0:
        logging.info("\tModeling up to {0} waves for each window and frequency".format(MaxWaves))
    else:
        logging.info("\tModeling exactly {0} waves for each window and frequency".format(MaxWaves))
    if ModelCircularVerticalWaves: logging.info("\tModeling circular vertical waves")
    if ModelCircularRayleighWaves: logging.info("\tModeling circular Rayleigh waves")
    if ModelCircularDissipativeRayleighWaves: logging.info("\tModeling circular dissipative Rayleigh waves")
  
    if Gamma == 0:
        logging.info("\tUsing modified Bayesian information criterion (BIC) with Gamma={0}. This will lead to overfitting.".format(Gamma))
    elif Gamma == 1:
        logging.info("\tUsing Bayesian information criterion (BIC) (Gamma={0})".format(Gamma))
    else:
        logging.info("\tUsing modified Bayesian information criterion (BIC) with Gamma={0:.2f}".format(Gamma))

    logging.info("")
    logging.info("")
    
    

    ### Go with the processing!
    t0 = time.time()
    for ss in range(0, Nshots):
        Shot = ShotList[ss]
        SourcePosition = Shot['SourcePosition']
        ShotStart = Shot['ShotStart']
        ShotEnd = Shot['ShotEnd']
        
        
        logging.info("Processing shot {0} of {1}. ".format(ss+1, Nshots))
        logging.info("\tShot start-end: {0:.3f}-{1:.3f}s".format(ShotStart, ShotEnd))
        if ss > 0:
            elapsedTime = time.time() - t0
            remainingTime = (Nshots - ss)*elapsedTime/ss
            logging.info("\tElapsed time {0:.0f} [s], estimated remaining time {1:.0f} [s]".format(elapsedTime, remainingTime))
        
        Kstart = int(round(ShotStart/Ts))
        Kend = int(round(ShotEnd/Ts))
        WindowId = db.addWindow(conn, Kstart, Kend, Ts)
        
        #        plt.figure()
        #        plt.plot(y[Kstart:Kend,:])
        #        plt.title(ss)
        #        plt.show()
        
        circularDecomposeWavefield(conn, y[Kstart:Kend,:], WindowId, Ts, Fvec_ndx, Kmax, KrStep, KiStep, Estep, Vmin, WavesToModel, MaxWaves, MaxIterations, ArrayInfo, SourcePosition, Gamma)
        
        # save intermediate results
        db.saveResults(conn, OUTPUT, WavesToModel, WindowId)
        


    t1 = time.time()
    elapsedTime = t1 - t0
    logging.info("Processing completed in {0:.0f} seconds. Current time is {1}.".format(elapsedTime, datetime.datetime.fromtimestamp(t1).strftime('%Y-%m-%d %H:%M:%S')))
    logging.info("Output was saved to {0}".format(OUTPUT))

        
    return # end of main()


if  __name__ =='__main__':
    main()


