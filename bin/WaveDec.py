#!/usr/bin/python3
# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
WaveDec a tool for the analysis of seismic waves.
See online documentation.
"""

import glob, os, errno, sys
import yaml
import time, datetime
import logging
import argparse
import numpy, scipy
from numpy import shape, zeros, min, max, any, linspace, logspace, pi, array, fft, floor, ceil, unique, log10
from scipy.spatial.distance import pdist

from EstimationRoutines import decomposeWavefield
from SyntheticWavefield import readSyntheticWavefield
from ReadSAC import readSacDir
import DataUtils as db

from wdSettings import MODEL_NOISE, MODEL_VERTICAL, MODEL_RAYLEIGH, MODEL_LOVE, MODEL_CIRCULAR_VERTICAL, MODEL_CIRCULAR_RAYLEIGH, MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH
from wdSettings import EWX, NSY, UDZ, ROTX, ROTY, ROTZ, WD_C_1, WD_C_2

def main():
    
    
    ### Initialization: commandline args, config file, and other parameters
    parser = argparse.ArgumentParser(description='WaveDec - a tool for seismic wavefield analysis')
    parser.add_argument("--config_file", default='config.yaml', help="path to YAML configuration file [config.yaml]")
    parser.add_argument("--input",  help="path to input folder (or file) [./]")
    parser.add_argument("--output", help="path to output folder [./]")
    parser.add_argument("--verbosity", type = int, default=1, help="increase output verbosity: 0 = only warnings, 1 = info, 2 = debug.  [0]")
    args = parser.parse_args()

    # Setup console logging, file logging later
    #logging.basicConfig(level=logging.DEBUG,    format='%(levelname)-8s %(message)s', stream=sys.stderr)
    logging.basicConfig(level=logging.INFO,    format='%(message)s', stream=sys.stderr)

    if args.verbosity == 0:
        logging.getLogger('').setLevel(logging.WARN) 
    elif args.verbosity == 1:
        logging.getLogger('').setLevel(logging.INFO) 
    elif args.verbosity == 2:
        logging.getLogger('').setLevel(logging.DEBUG)
    else:
        logging.error('Unrecognized verbosity level {0}'.format(args.verbosity))
        return
       

    
    
    
    # parameter CONIFG_FILE (configuration file in YAML format)
    CONFIG_FILE = os.path.abspath(args.config_file)    # default is "./config.yaml"
    try:
        with open(CONFIG_FILE) as f:
            conf = yaml.load(f)
            f.close()
    except yaml.scanner.ScannerError as e:
        logging.critical('There is a syntax problem with the YAML configuration file {0}.'.format(CONFIG_FILE))
        logging.critical(e)
        return
    except FileNotFoundError as e:
        logging.warning('No configuration file ({0}) found. Proceeding with default values.'.format(CONFIG_FILE))
        conf = dict()
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
        if exception.errno != errno.EEXIST:    raise
    
    
    # logging to file        
    logFile = logging.FileHandler(filename=os.path.join(OUTPUT, 'WaveDec.log'), mode='w')
    logFile.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', "%Y-%m-%d %H:%M:%S")
    logFile.setFormatter(formatter) # tell the handler to use this format
    logging.getLogger('').addHandler(logFile) # add the handler to the root logger
    
    logging.debug('Starting WaveDec with Python version:')
    logging.debug(sys.version_info)
    logging.debug('NumPy version: {0}'.format(numpy.__version__))
    logging.debug('SciPy version: {0}'.format(scipy.__version__))

    
    # Read input, either synthetic wavefield or files from directory
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
    
    # some quality control on input data          
    L = shape(y)[1] # number of channels
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
    
    
    # estimate resolution limits from array
    pos = array(list(set(tuple(p) for p in ArrayInfo[:,0:3])))
    Nsensors = shape(pos)[0]
    if Nsensors == 0:
        logging.critical("There is a problem with the array. Appears there are no sensors.")
        return
    elif Nsensors == 1:
        resolution_Kmin = None
        resolution_Kmax = None
        dmin = None
        dmax = None
    elif Nsensors > 1: # normal array setting
        dist = pdist(pos, 'euclidean')
        dmin = min(dist)
        dmax = max(dist)
        if (dmin != 0.0) and (dmax != 0.0):
            resolution_Kmax = WD_C_1 / dmin # As suggested in Asten and Henstridge [1984]
            resolution_Kmin = WD_C_2 / dmax
        
        
    # parsing YAML configuration file    
    try:
        # Several processing parameters
        Twindow = float(conf.get('Twindow', 25.0))
        Tstart = float(conf.get('Tstart', 0.0))
        MaxWindows = int(conf.get('MaxWindows', 0))
        MaxWaves = int(conf.get('MaxWaves', 3)) # maximum number of waves in model
        MaxIterations = int(conf.get('MaxIterations', 10)) # maximum number of iterations when refining estimates
        Ts_tmp = conf.get('Ts', None)
        if Ts_tmp is not None: # this is overriding what found SAC header
            Ts = float(Ts_tmp)
            logging.warning("Overriding sampling interval found in SAC header Ts: {0:.2e}".format(Ts))
        Fmin = float(conf.get('Fmin', 10/Twindow))
        Fmax = float(conf.get('Fmax',min([20, 0.5/Ts])))
        Fspacing = conf.get('Fspacing','lin')  # 'lin' or 'log'
        Fnum = int(conf.get('Fnum',50))
        Fnum = int(max([Fnum, 0])) # Fnum >= 1 or Fnum = 0 to model all freqs
        Estep = float(conf.get('Estep', pi/10.0))
        if Nsensors > 1:
            Kmax = float(conf.get('Kmax', 1.2*resolution_Kmax))
            Kstep = float(conf.get('Kstep', resolution_Kmin)) # step size (in wavenumber) for grid search
        elif Nsensors == 1:
            Kmax = conf.get('Kmax')
            if Kmax == None:
                logging.critical("In the single sensor setting, Kmax must be specified.")
                return
            Kmax = float(Kmax)
            Kstep = float(conf.get('Kstep', Kmax / 100)) # step size (in wavenumber) for grid search
        Vmin = float(conf.get('Vmin', 50))
        Gamma = float(conf.get('Gamma', 1))
        if Gamma < 0:
            logging.critical("Parameter Gamma should be greater than or equal to 0.")
            logging.critical("Gamma: {0}".format(Gamma))
            return
        
        # Which wave models to fit?
        ModelRayleighWaves = db.str2bool(conf.get('ModelRayleighWaves', True))
        ModelLoveWaves = db.str2bool(conf.get('ModelLoveWaves', True))
        ModelNoise = db.str2bool(conf.get('ModelNoise', True))
        ModelVerticalWaves = db.str2bool(conf.get('ModelVerticalWaves', False))
    except Exception as e:
        logging.critical('There was a problem while parsing the configuration file {0}'.format(CONFIG_FILE))
        logging.critical(e)
        raise ValueError('There was a problem while parsing the configuration file {0}'.format(CONFIG_FILE))


    
    # Does not make sense to model Vertical and Rayleigh waves together
    if ModelVerticalWaves and (ModelRayleighWaves == True or ModelLoveWaves == True):
        logging.warning("Parameter ModelVerticalWaves is true, setting ModelRayleighWaves and ModelLoveWaves to false")
        ModelRayleighWaves = ModelLoveWaves = False

    # Check whether the needed sensor components are present in input data
    comp=ArrayInfo[:,3]
    if ModelVerticalWaves and any(comp == UDZ) == False:
        logging.critical("ModelVerticalWaves is true, but no vertical components (UDZ) are found in input data")
        return
    if ModelRayleighWaves and any(comp == UDZ) == False:
        logging.critical("ModelRayleighWaves is true, but no vertical components (UDZ) are found in input data")
        return
    if ModelRayleighWaves and any(comp == EWX) == False:
        logging.critical("ModelRayleighWaves is true, but no east-west (EWX) components are found in input data")
        return
    if ModelRayleighWaves and any(comp == NSY) == False:
        logging.critical("ModelRayleighWaves is true, but no north-south (NSY) components are found in input data")
        return
    if ModelLoveWaves and any(comp == EWX) == False:
        logging.critical("ModelLoveWaves is true, but no east-west (EWX) components are found in input data")
        return
    if ModelLoveWaves and any(comp == NSY) == False:
        logging.critical("ModelLoveWaves is true, but no north-south (NSY) components are found in input data")
        return
    if Nsensors == 1 and (any(comp == ROTX) == False or any(comp == ROTY) == False or any(comp == ROTZ) == False):
        logging.critical("In the single sensor setting, a six components sensor is needed")
        return
        
    
    WavesToModel = {MODEL_RAYLEIGH:ModelRayleighWaves, MODEL_LOVE:ModelLoveWaves,
                    MODEL_NOISE:ModelNoise, MODEL_VERTICAL:ModelVerticalWaves,
                    MODEL_CIRCULAR_RAYLEIGH:False, MODEL_CIRCULAR_VERTICAL:False, MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH:False}
    
    
    # change below for plotting traces
    if False:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        
        K = shape(y)[0] # number of samples        
        L = shape(y)[1] # number of channels
        t_axis = linspace(0,Ts*(K-1),K)
        cmp = ArrayInfo[:,3]
            
        fig = plt.figure(figsize=(9,9)) # figsize=(6,9))
        gs = gridspec.GridSpec(3, 2)
        gs.update(wspace=0.05, hspace=0.0, left = 0.15, right = 0.95, bottom = 0.1, top = 0.95)
        ax0 = plt.subplot(gs[0,0])
        if sum(cmp == UDZ) > 0:
            plt.plot(t_axis, y[:,cmp == UDZ])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.ylabel('UDZ')
        plt.tick_params(labelbottom='off')
        ax1 = plt.subplot(gs[1,0], sharex=ax0)
        if sum(cmp == EWX) > 0:
            plt.plot(t_axis, y[:,cmp == EWX])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.tick_params(labelbottom='off')
        plt.ylabel('EWX')
        ax2 = plt.subplot(gs[2,0], sharex=ax0)
        if sum(cmp == NSY) > 0:
            plt.plot(t_axis, y[:,cmp == NSY])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.ylabel('NSY')
        ax3 = plt.subplot(gs[0,1], sharex=ax0)
        plt.xlabel('[s]')
        if sum(cmp == ROTZ) > 0:
            plt.plot(t_axis, y[:,cmp == ROTZ])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.ylabel('ROTZ')
        plt.tick_params(labelbottom='off')
        ax4 = plt.subplot(gs[1,1], sharex=ax0)
        if sum(cmp == ROTX) > 0:
            plt.plot(t_axis, y[:,cmp == ROTX])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.tick_params(labelbottom='off')
        plt.ylabel('ROTX')
        ax5 = plt.subplot(gs[2,1], sharex=ax0)
        if sum(cmp == ROTY) > 0:
            plt.plot(t_axis, y[:,cmp == ROTY])
        plt.xlim(t_axis[0], t_axis[-1])
        plt.ylabel('ROTY')
        plt.xlabel('[s]')
        #gs.tight_layout(fig)
        #gs.tight_layout(fig)
        plt.show()
        
    
    # initialize database for storing results
    conn = db.init()
    db.setArrayInfo(conn, ArrayInfo)

    if Tstart > 0.0: # ignore intial part of the recording
        y = y[int(ceil(Tstart/Ts)):,:]
    
    K = shape(y)[0] # number of samples
    L = shape(y)[1] # number of channels
    Tobservation = Ts*K
    Nwindows = int(floor(Tobservation/Twindow))
    Nwindows = Nwindows if MaxWindows == 0 else min([Nwindows, MaxWindows])
    Kw = int(ceil(Twindow/Ts)) # samples per window
    

    Fvec_fft = fft.rfftfreq(Kw, Ts); # array of all DFT frequencies
    
    if Fnum == 0: # Fnum = 0 for modeling all DFT frequencies in [Fmin,Fmax]
        Fvec_in = Fvec_fft[(Fvec_fft > 0) & (Fvec_fft >= Fmin) & (Fvec_fft <= Fmax)]
        Fnum = len(Fvec_in)        
        Fspacing = 'lin'
    elif (Fspacing.lower() == 'lin') or (Fspacing.lower() == 'linear'):
        Fvec_in = linspace(Fmin, Fmax, Fnum)
    elif (Fspacing.lower() == 'log') or (Fspacing.lower() == 'logarithmic'):
        Fvec_in = logspace(log10(Fmin), log10(Fmax), Fnum)
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
    db.createOutputFiles(conn, OUTPUT, WavesToModel, pos, Fvec, resolution_Kmin, resolution_Kmax)




    ### Print startup info
    logging.info("")
    logging.info(" _       __                 ____           ")
    logging.info("| |     / /___ __   _____  / __ \___  _____")
    logging.info("| | /| / / __ `/ | / / _ \/ / / / _ \/ ___/")
    logging.info("| |/ |/ / /_/ /| |/ /  __/ /_/ /  __/ /__  ")
    logging.info("|__/|__/\__,_/ |___/\___/_____/\___/\___/  ")
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
    logging.info("\tSignal duration {0:.2f} [s] ({1} samples, {2:.3f} [s] sampling interval)".format(Tobservation, K, Ts))
    if Tstart > 0.0:
        logging.info("\tThe initial {0:.2f} [s] of the recording were ignored".format(Tstart))
    if dmin is not None and dmax is not None and resolution_Kmin is not None and resolution_Kmax is not None:
        logging.info("\tSmallest interstation distance {0:.2f} [m], largest {1:.2f} [m]".format(dmin, dmax))
        logging.info("\tResolution limits, inferred from array geometry,  Kmin={0:.2e} Kmax={1:.2e} [1/m] ".format(resolution_Kmin, resolution_Kmax))
    logging.info("Processing info")
    logging.info("\tWindow duration {0:.2f} [s] ({1:.0f} windows in total)".format(Twindow, Nwindows))
    logging.info("\tProcessing {0} frequencies between {1:.2f} [Hz] and {2:.2f} [Hz] (with {3} spacing)".format(len(Fvec), Fvec[0], Fvec[-1], Fspacing))
    if Vmin > 0.0: logging.info("\tAnalyzing velocities greater than {0:.2e} [m/s]".format(Vmin))
    logging.info("\tAnalyzing wavenumbers smaller than {0:.2e} [1/m]".format(Kmax))
    logging.info("\tGrid step in wavenumber {0:.2e} [1/m]".format(Kstep))
    logging.info("\tGrid step in ellipticity angle {0:.2e} [rad]".format(Estep))

    logging.info("Model selection")
    if ModelVerticalWaves: logging.info("\tModeling vertical waves")
    if ModelLoveWaves: logging.info("\tModeling Love waves")
    if ModelRayleighWaves: logging.info("\tModeling Rayleigh waves")
    if ModelNoise and Gamma > 0:
        logging.info("\tModeling up to {0} wave(s) for each window and frequency".format(MaxWaves))
    else:
        logging.info("\tModeling exactly {0} wave(s) for each window and frequency".format(MaxWaves))

    if Gamma == 0:
        logging.info("\tUsing modified Bayesian information criterion (BIC) with Gamma={0}. This will lead to overfitting.".format(Gamma))
    elif Gamma == 1:
        logging.info("\tUsing Bayesian information criterion (BIC) (Gamma={0})".format(Gamma))
    else:
        logging.info("\tUsing modified Bayesian information criterion (BIC) with Gamma={0:.2f}".format(Gamma))
    if Fnum > 200:
        logging.info('\nWaveDec is about to process {0} frequencies in each time window. It seems a lot.\nComputation time is roughly proportional to the number of frequencies processed.'.format(Fnum))
    logging.info("")
    logging.info("")

    ### Go with the processing!
    t0 = time.time()
    for ww in range(0, Nwindows):
        logging.info("Processing window {0} of {1}".format(ww+1,Nwindows))
        if ww > 0:
            elapsedTime = time.time() - t0
            remainingTime = (Nwindows - ww)*elapsedTime/ww
            logging.info("\tElapsed time {0:.0f} [s], estimated remaining time {1:.0f} [s]".format(elapsedTime, remainingTime))
        
        Kstart = ww*Kw
        Kend = (ww+1)*Kw
        WindowId = db.addWindow(conn, Kstart, Kend, Ts)
        
        decomposeWavefield(conn, y[Kstart:Kend,:], WindowId, Ts, Fvec_ndx, Kmax, Kstep, Estep, Vmin, WavesToModel, MaxWaves, MaxIterations, ArrayInfo, Gamma)
        
        # save intermediate results
        db.saveResults(conn, OUTPUT, WavesToModel, WindowId)
        

    t1 = time.time()
    elapsedTime = t1 - t0
    logging.info("Processing completed in {0:.0f} seconds. Current time is {1}.".format(elapsedTime, datetime.datetime.fromtimestamp(t1).strftime('%Y-%m-%d %H:%M:%S')))
    logging.info("Output was saved to {0}".format(OUTPUT))

        
    return # end of main()


if  __name__ =='__main__':
    main()


