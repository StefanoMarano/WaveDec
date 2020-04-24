# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
Routines for reading SAC files
"""

import glob
import numpy as np
import logging
import struct
from wdSettings import Components


                
                
def readSac(sacFile):
    """Load a single SAC file.

    Parameters
    ----------
    sacFile : string
        Path to the SAC file.

    Returns
    -------
    data : float array
        Data, amplitudes
    t : float array
        Time axis
    sachead : float array
        SAC header

    """

    for ESTR in ["<", ">"]:
        #ESTR=">" # big endian # eg, SESAME format
        #ESTR="<" # little endian # eg, SED format
        #ESTR="@" # same as machine
        SHEAD="%c70f35l5L8s16s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s" % ESTR
        f=open(sacFile, mode="rb")
        sachead_size=struct.calcsize(SHEAD)
        tstr=f.read(sachead_size)
        sachead=struct.unpack(SHEAD, tstr)
        
        nvhdr = sachead[76]
        if nvhdr == 4 | nvhdr == 5:
            # Old sac format. Never tested.
            logging.warning("NVHDR = {0}, file {1} may be from old SAC version.".format(nvhdr, sacFile))
        elif nvhdr != 6:
            # We are reading in the wrong byte order.
            f.close()
        elif nvhdr == 6:
            # Good, we are reading in the propoer byte order.
            break
        else:
            logging.error("NVHDR = {0}, file {1} may be corrupted.".format(nvhdr, sacFile))
            
    dt=sachead[0]
    npts=sachead[79]

    t=np.arange(0, npts*dt, dt)
    dsize=struct.calcsize("%c%df" % (ESTR,npts))
    dat_str=f.read(dsize)
    data=np.array(struct.unpack("%c%df" % (ESTR, npts), dat_str))

    f.close()
    return(data, t, sachead)
    

def readSacDir(sac_dir):
    """Load SAC files from a folder.
 
    Parameters
    ----------
    sac_dir : string
        Path to the folder containing the SAC files.
 
    Returns
    -------
    y : 2d float array
        It is an array of size (K, L). Each column contains the signal at the l-th location.
    info : 2d array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
    Ts : float
        Sampling time in [s].
    
    """
    
    fileList = sorted(glob.glob(sac_dir+'/*.sac'))
    N_files = len(fileList)
    if N_files < 1:
        raise Exception('No files found in {0}'.format(sac_dir))
        
    
    info=np.zeros((N_files,6))
    for ff in range(0,N_files):
        logging.debug("Loading SAC file " + fileList[ff])
        (data, t, sachead)  = readSac(fileList[ff])
        if ff == 0:
            K = len(data)
            Ts = sachead[0] # Sampling interval
            y=np.zeros((K,N_files))
        else:
            if Ts != sachead[0]:
                #print(fileList[ff])
                #print(sachead[0])
                raise Exception('Sampling time Ts is not consistent across input files.')
            if K != len(data):
                #print(fileList[ff])
                #print(len(data))
                # TODO found in SED database a dataset with this problem
                raise Exception('Data length K is not consistent across input files.') 
        y[:,ff] = data
        info[ff,0] = sachead[47] # x coordinate
        info[ff,1] = sachead[48] # y coordinate
        info[ff,2] = sachead[49] # z coordinate
        info[ff,4] = sachead[0] # Sampling interval
        
        cmp = sachead[129].split(b'\x00')[0].strip() # terminate string with null char and remove trailing whitespaces
        try:
            info[ff,3] = Components[cmp] # component
        except KeyError:
            raise Exception('Unknown SAC component {0} found while reading {1}'.format(cmp, fileList[ff]))        

    return(y, info, Ts, fileList)
    
