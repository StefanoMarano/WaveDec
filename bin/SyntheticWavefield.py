# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
Routines for generating syntetic wavefield
"""

import numpy as np
from numpy import size, zeros, vstack, array, shape, linspace, cos, sin, pi, sqrt, random, ones, log, power, angle, abs
from scipy.special import hankel2
import yaml
import logging
import sys
from wdSettings import EWX, NSY, UDZ, ROTX, ROTY, ROTZ, Components



def readSyntheticWavefield(config_file):

    try:
        with open(config_file) as f:
            conf = yaml.load(f)
            f.close()
    except Exception as e:
        logging.critical(e)
        logging.critical('Error: No configuration file "{0}" found.'.format(config_file))
        sys.exit()


    logging.debug('Generating synthetic wavefield')
    
    K = conf.get('NumberSamples')
    Ts = conf.get('SamplingInterval')
    array_info = conf.get('array')
    info=array([])
    
    for aa in array_info:
        position = aa.get('position')
        components = aa.get('components')
        C = size(components)
        info_aa = zeros((C,5))
        for cc in range(0,C):
            info_aa[cc,0] = position[0]
            info_aa[cc,1] = position[1]
            info_aa[cc,2] = 0
            info_aa[cc,3] = Components[components[cc]]
            info_aa[cc,4] = Ts
        if size(info) == 0:
            info = info_aa
        else:
            info = vstack((info,info_aa))
    L=shape(info)[0]
    x=zeros((K,L))
    
    wavefield = conf.get('wavefield')
    for ww in wavefield:
        wave_type = ww.get('wave')
        Amplitude = ww.get('amplitude')
        Phase = ww.get('phase')
        Wavenumber = ww.get('wavenumber')
        Wavenumber_i = ww.get('wavenumber_i',None)
        Azimuth = ww.get('azimuth', None)
        SourcePosition = ww.get('sourceposition', None)
        EllipticityAngle = ww.get('ellipticityangle', None)
        Frequency = ww.get('frequency')
        sigma2 = ww.get('sigma2',None)

        if wave_type.lower() == 'rayleigh':
            logging.debug('\tRayleigh wave')
            logging.debug('\t\tAmplitude: {0:.3f}, Wavenumber: {1:.3f}, Azimuth: {2:.3f}, EllipticityAngle {3:.2f}'.format(Amplitude,Wavenumber, Azimuth, EllipticityAngle))
            logging.debug('\t\tFrequency: {0:.3f}'.format(Frequency))
            x = x + syntheticRayleighWave(Ts, K, info, Amplitude, Phase, Wavenumber, Azimuth, EllipticityAngle, Frequency)
        elif wave_type.lower() == 'circularrayleigh':
            logging.debug('\tCircular Rayleigh wave')
            logging.debug('\t\tAmplitude: {0:.3f}, Wavenumber: {1:.3f}, EllipticityAngle {2:.2f}'.format(Amplitude, Wavenumber, EllipticityAngle))
            logging.debug('\t\tSource Position: ({0:.3f}, {1:.3f})'.format(SourcePosition[0], SourcePosition[0]))
            logging.debug('\t\tFrequency: {0:.3f}'.format(Frequency))
            x = x + syntheticCircularRayleighWave(Ts, K, info, Amplitude, Phase, Wavenumber, EllipticityAngle, Frequency, SourcePosition)
        elif wave_type.lower() == 'circulardissipativerayleigh':
            Wavenumber = Wavenumber + 1j*Wavenumber_i
            logging.debug('\tCircular Dissipative Rayleigh wave')
            logging.debug('\t\tAmplitude: {0:.3f}, Wavenumber: {1:.3f}, Wavenumber_i: {1:.3f}, EllipticityAngle {2:.2f}'.format(Amplitude, Wavenumber, EllipticityAngle))
            logging.debug('\t\tSource Position: ({0:.3f}, {1:.3f})'.format(SourcePosition[0], SourcePosition[0]))
            logging.debug('\t\tFrequency: {0:.3f}'.format(Frequency))
            x = x + syntheticCircularRayleighWave(Ts, K, info, Amplitude, Phase, Wavenumber, EllipticityAngle, Frequency, SourcePosition)

        elif wave_type.lower() == 'love':
            logging.debug('\tLove wave')
            logging.debug('\t\tAmplitude: {0:.3f}, Wavenumber: {1:.3f}, Azimuth: {2:.3f}'.format(Amplitude,Wavenumber, Azimuth))
            logging.debug('\t\tFrequency: {0:.3f}'.format(Frequency))
            x = x + syntheticLoveWave(Ts, K, info, Amplitude, Phase, Wavenumber, Azimuth, Frequency)
        #elif wave_type.lower() == 'circularlove':
        #    logging.debug('\tCircular Love wave')
        #    logging.debug('\t\tAmplitude: {0:.3f}, Wavenumber: {1:.3f}'.format(Amplitude,Wavenumber))
        #    logging.debug('\t\tFrequency: {0:.3f}'.format(Frequency))
        #    x = x + syntheticCircularLoveWave(Ts, K, info, Amplitude, Phase, Wavenumber, Frequency, SourcePosition)
        elif wave_type.lower() == 'noise':
            logging.debug('\tAdditive Gaussian noise')
            logging.debug('\t\tSigma2: {0:.3e}'.format(sigma2))
            z = syntheticAWGN(sigma2, zeros((K,L)))
            
        else:
            logging.warning("Warning: unrecognized wave '{0}'".format(wave_type))
        
    # compute LL under H0 and H1
    y = x + z;
    LL_H0 = logLikelihoodIID(y, zeros((K,L)), sigma2)
    LL_H1 = logLikelihoodIID(y, x, sigma2)
    
    logging.debug('Theoretical log-likelihood values')
    logging.debug('\tH0: {0:.3e} (no signal present)'.format(LL_H0))
    logging.debug('\tH1: {0:.3e} (signal present)'.format(LL_H1))
	
	
	
    return(y, info, Ts, None)


def syntheticRayleighWave(Ts, K, info, Amplitude, Phase, Wavenumber, Azimuth, EllipticityAngle, Frequency):
    """Generate a monochromatic Rayleigh wave.
 
    Parameters
    ----------
    Ts : float
        The sampling time in [s]
    K : int
        Number of samples for each channel. Total duration = Ts*K
    info : 2d array
        Contains information about the sensor positions and channels.
    Amplitude : float
        Wave amplitude in arbitrary units.
    Phase : float
        Wave phase at the origin.
    Wavenumber : float
        Wavenumber in [1/m] (not in [rad/m]!)
    Azimuth : float
        Wave azimuth in [rad]. The wavevector points in the direction of propagation.
        Wavevector = Wavenumber*[cos(Azimuth), sin(Azimuth)]
    EllipticityAngle : float
        Ellipticity angle of the Rayleigh wave in [rad].
        The ellipticity angle is a value in the interval [-pi/2, pi/2]
    Frequency : float
        Temporal frequency of the wave in [Hz].
 
    Returns
    -------
    y : 2d float array
        It is an array of size (K, L). Each column contains the signal at the l-th location.

    """
    
    L = info.shape[0]
    y = zeros((K,L))
    
    
    t_axis = Ts*linspace(0,K-1,K)
    Wavenumber_x = cos(Azimuth)*Wavenumber
    Wavenumber_y = sin(Azimuth)*Wavenumber
    
    for ll in range(0,L):
        Component=info[ll,3]
        pos_x = info[ll,0]
        pos_y =info[ll,1]
        if Component == EWX: 
            y[:,ll] = Amplitude * sin(EllipticityAngle) * cos(Azimuth) * cos(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase)
        elif Component == NSY:
            y[:,ll] = Amplitude * sin(EllipticityAngle) * sin(Azimuth) * cos(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase)
        elif Component == UDZ:
            y[:,ll] = Amplitude * cos(EllipticityAngle) * cos(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase + pi/2)
        elif Component == ROTX:
            y[:,ll] = Amplitude * 2*pi* Wavenumber * cos(EllipticityAngle) * sin(Azimuth) * cos(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase)
        elif Component == ROTY:
            y[:,ll] = - Amplitude * 2*pi* Wavenumber * cos(EllipticityAngle) * cos(Azimuth) * cos(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase)
        elif Component == ROTZ:
            y[:,ll] = zeros((K,))
        else:
            logging.warning("Warning: Unrecognized component. Will produce zero signal.")
        
    return(y)
    
def syntheticLoveWave(Ts, K, info, Amplitude, Phase, Wavenumber, Azimuth, Frequency):
    """Generate a monochromatic Love wave.
 
    Parameters
    ----------
    Ts : float
        The sampling time in [s]
    K : int
        Number of samples for each channel. Total duration = Ts*K
    info : 2d array
        Contains information about the sensor positions and channels.
    Amplitude : float
        Wave amplitude in arbitrary units.
    Phase : float
        Wave phase at the origin.
    Wavenumber : float
        Wavenumber in [1/m] (not in [rad/m]!)
    Azimuth : float
        Wave azimuth in [rad]. The wavevector points in the direction of propagation.
        Wavevector = Wavenumber*[cos(Azimuth), sin(Azimuth)]
    Frequency : float
        Temporal frequency of the wave in [Hz].

    Returns
    -------
    y : 2d float array
        It is an array of size (K, L). Each column contains the signal at the l-th location.

    """
    
    L = info.shape[0]
    y = zeros((K,L))
    
    
    t_axis = Ts*linspace(0,K-1,K)
    Wavenumber_x = cos(Azimuth)*Wavenumber
    Wavenumber_y = sin(Azimuth)*Wavenumber
    
    for ll in range(0,L):
        Component=info[ll,3]
        pos_x = info[ll,0]
        pos_y =info[ll,1]
        if Component == EWX: 
            y[:,ll] = - Amplitude * sin(Azimuth) * cos(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase)
        elif Component == NSY:
            y[:,ll] = Amplitude * cos(Azimuth) * cos(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase)
        elif Component == UDZ:
            y[:,ll] = zeros((K,))
        elif Component == ROTX:
            y[:,ll] = zeros((K,))
        elif Component == ROTY:
            y[:,ll] = zeros((K,))
        elif Component == ROTZ:
            y[:,ll] = Amplitude * pi * Wavenumber * sin(2*pi*(Frequency*t_axis - (Wavenumber_x*pos_x+Wavenumber_y*pos_y) ) + Phase)
        else:
            logging.warning("Warning: Unrecognized component. Will produce zero signal.")
        
    return(y)
    
def syntheticAWGN(sigma2, y):
    """Add white Gaussian noise to a signal.
 
    Parameters
    ----------
    sigma2 : float
        Noise variance.
    y : 2d float array
        The noisless signal.

    Returns
    -------
    y : 2d float array
        The noisy signal.

    """
    
    n = sqrt(sigma2)*random.randn(shape(y)[0], shape(y)[1])
    
    return(y + n)
    
    
    
def logLikelihoodIID(y, x, sigma2):
    """Compute the log-likelihood of the observed signal y
 
    Parameters
    ----------
    y : 2d float array
        Noisy signal
    x : 2d float array
        The noisless signal.
    sigma2 : float
        Noise variance. Can be either a scalar (ie, equal noise variance on all channels)
        or a vector of length L (ie, specifies the different variance on each channel)

    Returns
    -------
    LL : float
        The log-likelihood of the observations y.

    """
    
    # K=shape(y)[0];
    if size(shape(y)) == 1:
        L=1
    else:
        L=shape(y)[1];
    if size(sigma2) == 1:
        sigma2_vec = sigma2 * ones((L,))

    LL=0;
    for ll in range(0,L):
        LL = LL  + sum( -  0.5 * log(2*pi*sigma2_vec[ll]) -  0.5*power(y[:,ll]-x[:,ll],2)/sigma2_vec[ll] )

    return LL


def syntheticCircularRayleighWave(Ts, K, info, Amplitude, Phase, Wavenumber, EllipticityAngle, Frequency, SourcePosition):
    """Generate a monochromatic circular Rayleigh wave.
    
       If Wavenumber_i>0 the wave is dissipative
    """

    if np.imag(Wavenumber) > 0:
        logging.critical("Error: To model a dissipative wave, the imaginary part of the wavenumber should be smaller than zero.")
        sys.exit()
        
        
    L = info.shape[0]
    y = np.zeros((K,L)) +1j*np.zeros((K,L))
    t_axis = Ts*np.linspace(0,K-1,K)
    
    for ll in range(0,L):
        Component=info[ll,3]
        pos_x = info[ll,0]
        pos_y =info[ll,1]
        r = np.sqrt( (SourcePosition[0]-pos_x)**2 + (SourcePosition[1]-pos_y)**2)
        Azimuth = np.arctan2(pos_y-SourcePosition[1], pos_x-SourcePosition[0])

        if Component == EWX:
            mag = abs(-1j*hankel2(1, 2*np.pi*Wavenumber*r))
            ang = angle(-1j*hankel2(1, 2*np.pi*Wavenumber*r))
            y[:,ll] = Amplitude * sin(EllipticityAngle) * cos(Azimuth) * cos(2*pi*Frequency*t_axis + Phase + ang) * mag
        elif Component == NSY:
            mag = abs(-1j*hankel2(1, 2*np.pi*Wavenumber*r))
            ang = angle(-1j*hankel2(1, 2*np.pi*Wavenumber*r))
            y[:,ll] = Amplitude * sin(EllipticityAngle) * sin(Azimuth) * cos(2*pi*Frequency*t_axis + Phase + ang) * mag
        elif Component == UDZ:            
            mag = abs(1j*hankel2(0, 2*np.pi*Wavenumber*r))
            ang = angle(1j*hankel2(0, 2*np.pi*Wavenumber*r))
            y[:,ll] = Amplitude * cos(EllipticityAngle) *                cos(2*pi*Frequency*t_axis + Phase + ang) * mag
        else:
            logging.warning("Warning: Unrecognized component. Will produce zero signal.")
    return np.real(y)



#def syntheticCircularLoveWave(Ts, K, info, Amplitude, Phase, Wavenumber, Frequency, SourcePosition):
#    """Generate a monochromatic circular Love wave.
#       This requires modeling the orientation of the source...        
#    """
      
      