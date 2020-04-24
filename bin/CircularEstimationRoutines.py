# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
Here are estimation routines used in WaveDecActive
"""

from scipy.optimize import minimize, minimize_scalar
import numpy as np
from numpy import shape, fft, ceil, zeros, linspace, concatenate, meshgrid, array
from numpy import isnan, cos, sin, size, matrix, transpose, sum
from numpy import arctan2, log, sqrt, pi, dot, reshape, argmin, arange
import DataUtils as db
import logging
import sys
from scipy.spatial.distance import cdist
from scipy.special import hankel1, hankel2

from EstimationRoutines import fitNoise, estimateNoise, bwMessages
from wdSettings import EWX, NSY, UDZ
from wdSettings import MODEL_NOISE, MODEL_CIRCULAR_VERTICAL, MODEL_CIRCULAR_RAYLEIGH, MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH




def circularDecomposeWavefield(conn, y, WindowId, Ts, Fvec_ndx, Kmax, KrStep, KiStep, Estep, Vmin, WavesToModel, MaxWaves, MaxIterations, ArrayInfo, ShotSource, Gamma):
    """ Fit different wave types at several frequencies.
 
    Parameters
    ----------
    conn : 
        SQL database connection
    y : float array
        input signal
    WindowId : int
        Unique indentifier of the window
    Ts : float
        Sampling time [s]
    Fvec_ndx : int array
        ...
    Kmax : float
        Largest wavenumber [1/m] to analyze
    KrStep : float
        Grid step in the wavenumber plane [1/m]
    KiStep : float
        Grid step in the imaginary part of the wavenumber [1/m]
        Attenuation coefficient alpha must be in [0,1], this set a range for Wavenumber_i [-0.16,0], in practice set it between [-0.1, 0]
    Vmin : float
        Smallest velocity [m/s] to analyze
    WavesToModel : 
        Describe which wave types should be fitted
    MaxWaves : int
        Maximum number of waves to model at each frequency in the time window
    MaxIterations : int
        Maximum number of interations. In each iteration parameters of each wave are re-estimated.
    ArrayInfo : float array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
    Gamma : float
        Controls model complexity. (0 for pure ML, 1 for BIC, other values for intermediate strategies)

    """
    
    K = shape(y)[0]
    L = shape(y)[1]
    Fvec_fft = fft.fftfreq(K, Ts);
    # Fvec = Fvec_fft[Fvec_ndx]
    
    (Sm_bw, Sw_bw, Swm_bw) = bwMessages(y, Ts)
    Sm_fw_all = zeros(shape(Sm_bw))
    
        
    F=shape(Sm_bw)[2]
    
    NumKr = 2*ceil(Kmax/KrStep) # number of points in the wavenumber search grid. Even, so that we do not have 0 in the final search grid
    NumE = ceil(np.pi/Estep) # number of points in the ellipticity search grid Even, so that we have 0 in the final search grid
    
    Wavenumber = linspace(KrStep, Kmax, NumKr)
    EllipticityAngle = linspace(-pi/2, pi/2, NumE, endpoint=False)

    if WavesToModel[MODEL_CIRCULAR_VERTICAL]:
        xx = Wavenumber
        ndx_ok = xx <= Kmax # only wavenumber smaller than Kmax
        xx = xx[ndx_ok];
        X_grid_CV = array([xx])

    if WavesToModel[MODEL_CIRCULAR_RAYLEIGH] or WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH]:
        xx, yy = meshgrid(Wavenumber, EllipticityAngle)
        xx=concatenate(xx); yy=concatenate(yy);
        ndx_ok = xx <= Kmax # only wavenumber smaller than Kmax
        xx = xx[ndx_ok]; yy = yy[ndx_ok];
        X_grid_CR = array([xx, yy]) # Rayleigh waves
    
    
    if WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH]:        
        NumKi = np.max([np.ceil(0.1/KiStep),10])    # minimum 10
        Wavenumber_i = linspace(-0.1, 0, NumKi)     # attenuation coefficient alpha must be in [0,1], this set a range for Wavenumber_i [-0.16,0], in practice set it between [-0.1, 0]
        xx, yy, zz = meshgrid(Wavenumber, EllipticityAngle, Wavenumber_i)
        xx=concatenate(xx); yy=concatenate(yy); zz=concatenate(zz);
        ndx_ok = xx <= Kmax # only wavenumber smaller than Kmax
        xx = xx[ndx_ok]; yy = yy[ndx_ok]; zz = zz[ndx_ok];
        X_grid_CDR = array([xx, yy, zz]) # Rayleigh waves
    

    
    # Fndx
    for ff in Fvec_ndx:
        
        if WavesToModel[MODEL_CIRCULAR_RAYLEIGH] or WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH]:
            X_grid_CR_f = X_grid_CR
        if WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH]:
            X_grid_CDR_f = X_grid_CDR
        if WavesToModel[MODEL_CIRCULAR_VERTICAL]:
            X_grid_CV_f = X_grid_CV

        if Vmin != None and Vmin > 0: # further restrict search grid
            if WavesToModel[MODEL_CIRCULAR_VERTICAL]:
                ndx_ok = X_grid_CV[0,:] <= Fvec_fft[ff]/Vmin
                X_grid_CV_f = X_grid_CV[:,ndx_ok]
            if WavesToModel[MODEL_CIRCULAR_RAYLEIGH]:
                ndx_ok = X_grid_CR[0,:] <= Fvec_fft[ff]/Vmin 
                X_grid_CR_f = X_grid_CR[:,ndx_ok]
            if WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH]:
                ndx_ok = X_grid_CDR[0,:] <= Fvec_fft[ff]/Vmin 
                X_grid_CDR_f = X_grid_CDR[:,ndx_ok]
    
 
        Sm_fw = zeros((2,L,F,MaxWaves)) # shape() = (2,L,F)
        WaveModel = zeros((MaxWaves,))
        WaveAmplitude = zeros((MaxWaves,))
        WavePhase = zeros((MaxWaves,))
        WaveX_ML = [None]* MaxWaves
        
        # gradually increase the number of waves modeled
        NumWaves = 0 # number of waves modeled
        for mm in range(0,MaxWaves):
            logging.debug("Fitting {0}-th model (at {1:.3f} [Hz])".format(mm+1, Fvec_fft[ff]))
            tmpBIC=[]
            tmpModel=[]
            tmpSm_bw = Sm_bw - sum(Sm_fw ,3) + Sm_fw[:,:,:,mm]
            # how about variance messages?
            
            
            # Parameter estimation: attempt to fit different possible models
            if WavesToModel[MODEL_NOISE]:
                (BIC_N, sigma2_ML) = fitNoise(tmpSm_bw, L, K, Gamma)
                tmpBIC.append(BIC_N); tmpModel.append(MODEL_NOISE);               
            if WavesToModel[MODEL_CIRCULAR_VERTICAL] and size(X_grid_CV_f) > 0:
                (BIC_CV, Sm_fw_CV, Amplitude_CV, Phase_CV, X_ML_CV) = fitCircularVerticalWave(X_grid_CV_f, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma)
                tmpBIC.append(BIC_CV); tmpModel.append(MODEL_CIRCULAR_VERTICAL);
            if WavesToModel[MODEL_CIRCULAR_RAYLEIGH] and size(X_grid_CR_f) > 0:
                (BIC_CR, Sm_fw_CR, Amplitude_CR, Phase_CR, X_ML_CR) = fitCircularRayleighWave(X_grid_CR_f, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma)
                tmpBIC.append(BIC_CR); tmpModel.append(MODEL_CIRCULAR_RAYLEIGH);
            if WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH] and size(X_grid_CR_f) > 0 and size(X_grid_CDR_f) > 0:
                (_, _, _, _, X_ML_CR) = fitCircularRayleighWave(X_grid_CR_f, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma)
           
                _WavenumberML = X_ML_CR[0]
                _EllipticityML = X_ML_CR[1]
                _WavenumberRange = 0.025*(np.max(X_grid_CR_f[0,:]) - np.min(X_grid_CR_f[0,:]))
                _EllipticityRange = 0.025*(np.max(X_grid_CR_f[1,:]) - np.min(X_grid_CR_f[1,:]))
                ndx_W = np.logical_and(_WavenumberML-_WavenumberRange <= X_grid_CDR_f[0,:], X_grid_CDR_f[0,:] <= _WavenumberML+_WavenumberRange)
                ndx_E = np.logical_and(_EllipticityML-_EllipticityRange <= X_grid_CDR_f[1,:], X_grid_CDR_f[1,:] <= _EllipticityML+_EllipticityRange )
                ndx = np.logical_and(ndx_W, ndx_E)

                if np.any(ndx==True):
                    X_grid_CDR_f_reduced = X_grid_CDR_f[:,ndx]
                else:
                    xx, yy, zz = meshgrid([_WavenumberML], [_EllipticityML], Wavenumber_i)
                    xx=concatenate(xx); yy=concatenate(yy); zz=concatenate(zz);
                    X_grid_CDR_f_reduced = array([xx, yy, zz]) # Rayleigh waves
                    
                (BIC_CDR, Sm_fw_CDR, Amplitude_CDR, Phase_CDR, X_ML_CDR) = fitCircularDissipativeRayleighWave(X_grid_CDR_f_reduced, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma)
                tmpBIC.append(BIC_CDR); tmpModel.append(MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH);


            # Model selection: choose wave with smallest BIC
            if len(tmpBIC) > 0:
                WaveModel[mm] = tmpModel[np.argmin(tmpBIC)]
                histBIC = zeros((MaxIterations+1,))
                histBIC[0] = np.min(tmpBIC)
            if WaveModel[mm] == MODEL_NOISE:
                break     # when a noise model is chosen no more waves need to be added
            if WaveModel[mm] == 0:
                break     # Nothing was modeled. Perhaps because of stringent Vmin at this frequency
            elif WaveModel[mm] == MODEL_CIRCULAR_VERTICAL:
                Sm_fw[:,:,:,mm] = Sm_fw_CV
                WaveAmplitude[mm] = Amplitude_CV
                WavePhase[mm] = Phase_CV
                WaveX_ML[mm] = X_ML_CV
            elif WaveModel[mm] == MODEL_CIRCULAR_RAYLEIGH:
                Sm_fw[:,:,:,mm] = Sm_fw_CR
                WaveAmplitude[mm] = Amplitude_CR
                WavePhase[mm] = Phase_CR
                WaveX_ML[mm] = X_ML_CR
            elif WaveModel[mm] == MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH:
                Sm_fw[:,:,:,mm] = Sm_fw_CDR
                WaveAmplitude[mm] = Amplitude_CDR
                WavePhase[mm] = Phase_CDR
                WaveX_ML[mm] = X_ML_CDR
            else:
                logging.warning("Warning: unrecognized wave model")
                
                
                
            # refine existing estimates, if needed
            NumWaves = sum( (WaveModel > 0) & (WaveModel != MODEL_NOISE) )
            if (MaxIterations > 0) and (NumWaves  > 1):
                for ii in range(1,MaxIterations):
                    for m in range(0,mm+1):
                        logging.debug("Refining estimates of model {0}/{1}".format(m+1,mm+1))
                        tmpSm_bw = Sm_bw - sum(Sm_fw ,3) + Sm_fw[:,:,:,m]
                        if WaveModel[m] == MODEL_NOISE:
                            continue
                        elif WaveModel[m] == MODEL_CIRCULAR_VERTICAL:
                            X = WaveX_ML[m]
                            (BIC_CV, Sm_fw_CV, Amplitude_CV, Phase_CV, X_ML_CV) = fitCircularVerticalWave(X, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma)
                            WaveX_ML[m] = X_ML_CV
                            WaveAmplitude[m] = Amplitude_CV
                            WavePhase[m] = Phase_CV
                            Sm_fw[:,:,:,m] = Sm_fw_CV
                        elif WaveModel[m] == MODEL_CIRCULAR_RAYLEIGH:
                            X = WaveX_ML[m]
                            (BIC_CR, Sm_fw_CR, Amplitude_CR, Phase_CR, X_ML_CR) = fitCircularRayleighWave(X, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma)
                            WaveX_ML[m] = X_ML_CR
                            WaveAmplitude[m] = Amplitude_CR
                            WavePhase[m] = Phase_CR
                            Sm_fw[:,:,:,m] = Sm_fw_CR
                        elif WaveModel[m] == MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH:
                            X = WaveX_ML[m]
                            (BIC_CR, Sm_fw_CDR, Amplitude_CDR, Phase_CDR, X_ML_CDR) = fitCircularDissipativeRayleighWave(X, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma)
                            WaveX_ML[m] = X_ML_CDR
                            WaveAmplitude[m] = Amplitude_CDR
                            WavePhase[m] = Phase_CDR
                            Sm_fw[:,:,:,m] = Sm_fw_CDR
                        else:
                            logging.warning("Unrecognized wave model {0}".format(WaveModel[m]))
                    tmpSm_bw = Sm_bw - sum(Sm_fw ,3)
                    (BIC_N, sigma2_ML) = fitNoise(tmpSm_bw, L, K, Gamma)
                    histBIC[ii] = BIC_N # TODO is this the right BIC value, multiple waves are they all same?
                    if abs(histBIC[ii]-histBIC[ii-1])/abs(histBIC[ii-1]) < 0.01:
                        break
        # Compute forward messages and estimated paramters
        Sm_fw_all += sum(Sm_fw ,3)
        NumWaves = sum( (WaveModel > 0) & (WaveModel != MODEL_NOISE) )
        for mm in range(0,NumWaves):
            if WaveModel[mm] == MODEL_NOISE:
                continue
            elif WaveModel[mm] == MODEL_CIRCULAR_VERTICAL:
                Wavenumber = WaveX_ML[mm]
                db.addCircularVerticalWave(conn, WindowId, ff, WaveAmplitude[mm], WavePhase[mm], Wavenumber)
            #elif WaveModel[mm] == MODEL_CIRCULAR_LOVE:
            #    Wavenumber = WaveX_ML[mm]
            #    db.addCircularLoveWave(conn, WindowId, ff, WaveAmplitude[mm], WavePhase[mm], Wavenumber)
            elif WaveModel[mm] == MODEL_CIRCULAR_RAYLEIGH:
                Wavenumber = WaveX_ML[mm][0]
                EllipticityAngle = WaveX_ML[mm][1]
                db.addCircularRayleighWave(conn, WindowId, ff, WaveAmplitude[mm], WavePhase[mm], Wavenumber, EllipticityAngle)
            elif WaveModel[mm] == MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH:
                Wavenumber_r = WaveX_ML[mm][0]
                Wavenumber_i = WaveX_ML[mm][2]
                Wavenumber = Wavenumber_r +1j*Wavenumber_i
                EllipticityAngle = WaveX_ML[mm][1]
                db.addCircularDissipativeRayleighWave(conn, WindowId, ff, WaveAmplitude[mm], WavePhase[mm], Wavenumber, EllipticityAngle)
            else:
                logging.warning("Unrecognized wave model {0}".format(WaveModel[mm]))


    # after all freq have been processed
    (BIC_N, sigma2_ML) = fitNoise(Sm_bw - Sm_fw_all, L, K, Gamma)
    db.addNoise(conn, WindowId, sigma2_ML)
    
    # TODO after estimating noise again, do we want to iterate once more on wave params    
        
    return


### Vertical wave



def negLL_CircularVerticalWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    """ Computes the loglikelihood of a circular wave (vertical component only)
 
    Parameters
    ----------
    X_grid : float array
        Array of size (2,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber along the x axes. Second row the wavenumber along the y axes.

    Sw_bw, Swm_bw, SlnGamma_bw :
        sufficient statistics

    ArrayInfo : float array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
 
    Returns
    -------
    negLogLikelihood : float array
        A one dimensional array of N elements with the value of minus the log-likelihood.

    """
    
    N_points=size(X_grid)
    if N_points == 1:
        Wavenumber_vec = np.array([np.asscalar(X_grid)])
    else:
        Wavenumber_vec = X_grid[0,:]
        
    L=shape(Swm_bw)[1]
    negLogLikelihood = zeros((N_points,))
    comp=ArrayInfo[:,3]
    r_source = cdist(ArrayInfo[:,0:2], np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    

    ndx_v = arange(0, L)[comp == UDZ]
    L_v = len(ndx_v) # number of vertical components
    phi = zeros((N_points,))
    alpha = zeros((N_points,))
    Uw = zeros((N_points))
    Uwm = zeros((2,N_points))
    for ll_v in range(0, L_v):
        ll = ndx_v[ll_v]
        hankel_tmp = hankel2(0,2*np.pi*Wavenumber_vec*r_source[ll])
        phi = np.angle(hankel_tmp)
        alpha = np.abs(hankel_tmp)
        Uw[:] += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += + alpha *( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += + alpha *( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
    Um = zeros((2,N_points))
    Um[0,:] = Uwm[0,:] / Uw[:]
    Um[1,:] = Uwm[1,:] / Uw[:]
    negLogLikelihood[:] = -0.5*(Um[0,:] * Uwm[0,:] + Um[1,:] * Uwm[1,:]) + sum(-SlnGamma_bw)
 
    return negLogLikelihood


    
def fwMessages_CircularVerticalWave(X_ML, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    """
    
    RETURN
        amplitude and phase returned refer to outgoing wave. It is assumed there is no ingoing wave.
    """
    
    
    Wavenumber = X_ML
    
    L=shape(ArrayInfo)[0]
    comp=ArrayInfo[:,3]
    
    r_source = cdist(ArrayInfo[:,0:2], np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    
    Sm_fw = zeros((2,L))
    Uw = matrix(zeros((2,2)))
    Uwm = matrix(zeros((2,1)))
    
    # compute the LL
    for ll in range(0,L):
        if comp[ll] == UDZ:
            hankel_tmp = hankel2(0,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        
        Uw = Uw + transpose(H) * matrix(Sw_bw[:,:,ll]) * H
        
        Uwm = Uwm + transpose(H) * transpose(matrix(Swm_bw[:,ll]))
    Um = np.linalg.solve(Uw, Uwm)
    LL = 0.5*transpose(Um) * Uwm + sum(SlnGamma_bw)
    LL = LL[0,0]
    
    # compute the forward messages Sm_fw
    for ll in range(0,L):
        if comp[ll] == UDZ:
            hankel_tmp = hankel2(0,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)], [sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        tmp = H * Um
        Sm_fw[0,ll] = tmp[0]
        Sm_fw[1,ll] = tmp[1]
    
    Amplitude = sqrt(Um[0][0,0]**2 + Um[1][0,0]**2)
    Phase = arctan2(Um[0][0,0], Um[1][0,0])

    return (Sm_fw, Amplitude, Phase, LL)

 
def fitCircularVerticalWave(X_grid, Sm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma):
    """ Computes the loglikelihood of a wave with circular wavefront, vertical component only
 
    Parameters
    ----------
    X_grid : float array
        Array of size (N,). Description of N points to evaluate likelihoodfunction. First row is the wavenumber.

    Sw_bw, Swm_bw, SlnGamma_bw : float array
        sufficient statistics

    ArrayInfo : float array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
 
    Returns
    -------
    negLogLikelihood : float array
        A one dimensional array of N elements with the value of minus the log-likelihood.
"""
    

    
    Sm_fw = zeros(shape(Sm_bw))
    X_grid= X_grid.reshape((1,-1)) # can have grid or single point as input
    

    
    MaxIterations=5
    hist_LL=zeros((MaxIterations,))
    for ii in range(0,MaxIterations):
    
        (sigma2, SlnGamma_bw) = estimateNoise(Sm_bw, Sm_fw)
        
    
        Snw_bw = zeros((2,2,L)); Snwm_bw = zeros((2,L));
        for ll in range(0,L):
            Snw_bw[:,:,ll] = Sw_bw[:,:,ll,ff] / sigma2[ll];
            Snwm_bw[0,ll] = dot( Snw_bw[0,:,ll] , Sm_bw[:,ll,ff] )
            Snwm_bw[1,ll] = dot( Snw_bw[1,:,ll] , Sm_bw[:,ll,ff] )
            
            
        #(LL) = negLL_CircularVerticalWave(np.array([0.02]), Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, d_source)
        #(LL) = negLL_CircularLoveWave(np.array([0.02]), Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource)
        #assert 1 == 2 
        
        if ii == 0:
            # ML estimate with grid search
            (nLL) = negLL_CircularVerticalWave(X_grid, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource)
            
            ndx_min = np.argmin(nLL)
            X_g = X_grid[:,ndx_min];
            


        # Refine ML estimate
        res = minimize_scalar(negLL_CircularVerticalWave, bracket=(np.min(X_grid), X_g[0], np.max(X_grid)), bounds=(0, Kmax), args=(Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource), method='bounded')
        if res.success:
            X_ML = res.x
            LL_ML = -res.fun[0]
        else:
            logging.warning(res)
            X_ML = X_g
            LL_ML = -np.min(nLL)
            

        
        # Compute fw messages
        if isinstance(X_ML, np.ndarray): # Oddly, minimize scalar returns sometimes a scalar, sometimes an array (???)
            X_ML = X_ML[0]
        (Sm_fw_ff, Amplitude, Phase, LL_ML) = fwMessages_CircularVerticalWave(X_ML, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource)
        Sm_fw[:,:,ff] = Sm_fw_ff
        
        
        
        
        hist_LL[ii] = LL_ML
        if (ii > 0) and (abs(hist_LL[ii]-hist_LL[ii-1])/abs(hist_LL[ii-1]) < 0.01):
            break
        
        
    Wavenumber = X_ML
    LogLikelihood = LL_ML
    NumParameters = 3 * L + 3 # Amplitude, phase, wavenumber
    NumPoints = K * L
    BIC = -2* LogLikelihood + Gamma*NumParameters *log(NumPoints)

    
    
    logging.debug('Circular Vertical wave fit')
    logging.debug('\tLL: {0:.3e} BIC: {1:.3e}'.format(LogLikelihood, BIC))
    logging.debug('\tAmplitude: {0:.3f}, Phase: {1:.3f}, Wavenumber: {2:.3f}'.format(Amplitude,Phase,Wavenumber))

    return(BIC, Sm_fw, Amplitude, Phase, Wavenumber)

### Rayleigh wave

def negLL_CircularRayleighWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    """ Computes the loglikelihood of a Rayleigh wave with circular wavefront.
 
    Parameters
    ----------
    X_grid : float array
        Array of size (2,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber. Second row the ellipticity angle.

    Sw_bw, Swm_bw, SlnGamma_bw : float array
        sufficient statistics

    ArrayInfo : float array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
 
    Returns
    -------
    negLogLikelihood : float array
        A one dimensional array of N elements with the value of minus the log-likelihood.
"""


    
    if shape(shape(X_grid))[0] == 1: # single point
        X_grid = reshape(X_grid, (2,-1))
   
    N_points=shape(X_grid)[1]
    if N_points == 1:
        Wavenumber = X_grid[0]
        EllipticityAngle=X_grid[1]
    else:
        Wavenumber = X_grid[0,:]
        EllipticityAngle=X_grid[1,:]

    L=shape(Swm_bw)[1]
    negLogLikelihood = zeros((N_points,))
    pos = ArrayInfo[:,0:2]
    comp = ArrayInfo[:,3]
    r_source = cdist(ArrayInfo[:,0:2], np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    

    Uw = zeros((N_points))
    Uwm = zeros((2,N_points))

    # EWX
    ndx_e = arange(0, L)[comp == EWX]
    L_e = len(ndx_e) # number of EWX components
    alpha = zeros((N_points,))
    phi = zeros((N_points, ))
    for ll_e in range(0, L_e):
        ll = ndx_e[ll_e]
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
        alpha = sin(EllipticityAngle)*cos(Azimuth)*np.abs(hankel_tmp)
        phi = np.angle(hankel_tmp)
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
    
    # NSY
    ndx_n = arange(0, L)[comp == NSY]
    L_n = len(ndx_n) # number of NSY components
    alpha = zeros((N_points,))
    phi = zeros((N_points, ))
    for ll_n in range(0, L_n):
        ll = ndx_n[ll_n]
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
        alpha = sin(EllipticityAngle)*sin(Azimuth)*np.abs(hankel_tmp)
        phi = np.angle(hankel_tmp)
        Uw[:] += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
        
    # UDZ
    ndx_z = arange(0, L)[comp == UDZ]
    L_z = len(ndx_z) # number of NSY components
    alpha = zeros((N_points,))
    phi = zeros((N_points,))
    for ll_z in range(0, L_z):
        ll = ndx_z[ll_z]
        hankel_tmp = 1j*hankel2(0,2*np.pi*Wavenumber*r_source[ll])
        alpha = cos(EllipticityAngle)*np.abs(hankel_tmp)
        phi = np.angle(hankel_tmp)
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )

    
    Um = zeros((2,N_points))
    Um[0,:] = Uwm[0,:] / Uw[:]
    Um[1,:] = Uwm[1,:] / Uw[:]
    negLogLikelihood[:] = -0.5*(Um[0,:] * Uwm[0,:] + Um[1,:] * Uwm[1,:]) + sum(-SlnGamma_bw)
        
    return(negLogLikelihood)
 
    

        
def fwMessages_CircularRayleighWave(X_ML, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    
    
    Wavenumber = X_ML[0]
    EllipticityAngle = X_ML[1]
    
    L=shape(ArrayInfo)[0]
    comp=ArrayInfo[:,3]
    pos = ArrayInfo[:,0:2]
    
    r_source = cdist(pos, np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    
    Sm_fw = zeros((2,L))
    Uw = matrix(zeros((2,2)))
    Uwm = matrix(zeros((2,1)))
    
    # compute the LL
    for ll in range(0,L):
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        if comp[ll] == EWX:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * cos(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif comp[ll] == NSY:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * sin(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif comp[ll] == UDZ:
            hankel_tmp = 1j*hankel2(0,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = cos(EllipticityAngle) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        
        Uw = Uw + transpose(H) * matrix(Sw_bw[:,:,ll]) * H
        
        Uwm = Uwm + transpose(H) * transpose(matrix(Swm_bw[:,ll]))
    Um = np.linalg.solve(Uw, Uwm)
    LL = 0.5*transpose(Um) * Uwm + sum(SlnGamma_bw)
    LL = LL[0,0]
    
    # compute the forward messages Sm_fw
    for ll in range(0,L):
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        if comp[ll] == EWX:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * cos(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)], [sin(phi), cos(phi)]])
        elif  comp[ll] == NSY:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * sin(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == UDZ:
            hankel_tmp = 1j*hankel2(0,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = cos(EllipticityAngle) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        tmp = H * Um
        Sm_fw[0,ll] = tmp[0]
        Sm_fw[1,ll] = tmp[1]
    
    Amplitude = sqrt(Um[0][0,0]**2 + Um[1][0,0]**2)
    Phase = arctan2(Um[0][0,0], Um[1][0,0])

    return (Sm_fw, Amplitude, Phase, LL)
    
    


 
def fitCircularRayleighWave(X_grid, Sm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma):
    
    Sm_fw = zeros(shape(Sm_bw))
    X_grid= X_grid.reshape((2,-1)) # can have grid or single point as input
    
    

    
    MaxIterations=5
    hist_LL=zeros((MaxIterations,))
    for ii in range(0,MaxIterations):
    
        (sigma2, SlnGamma_bw) = estimateNoise(Sm_bw, Sm_fw)
    
        Snw_bw = zeros((2,2,L)); Snwm_bw = zeros((2,L));
        for ll in range(0,L):
            Snw_bw[:,:,ll] = Sw_bw[:,:,ll,ff] / sigma2[ll];
            Snwm_bw[0,ll] = dot( Snw_bw[0,:,ll] , Sm_bw[:,ll,ff] )
            Snwm_bw[1,ll] = dot( Snw_bw[1,:,ll] , Sm_bw[:,ll,ff] )
            
       
        
        if ii == 0:
            # ML estimate with grid search
            (nLL) = negLL_CircularRayleighWave(X_grid, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource)
                        
            ndx_min = argmin(nLL)
            X_g = X_grid[:,ndx_min];
           



        # Refine ML estimate
       
        # TODO 
        # L-BFGS-B often fails giving ABNORMAL_TERMINATION_IN_LNSRCH
        # therefore we use the unconstrained method Nelder-Mead
        # perhaps it would be nice to try L-BFGS-B providing explicit gradient
        
        # res = minimize(negLL_CircularRayleighWave, X_g, (Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource),'L-BFGS-B', bounds=[(0, Kmax), (-pi, pi)])
        # if res.success:
        #     X_ML = res.x
        #     LL_ML = -res.fun[0]
        # else:
        #     logging.warning(res)
        #     X_ML = X_g
        #     LL_ML = -np.min(nLL)
        
        res = minimize(negLL_CircularRayleighWave, X_g, (Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource),'Nelder-Mead')
        if res.success:
            X_ML = res.x
            LL_ML = -res.fun
        else:
            logging.warning(res)
            X_ML = X_g
            LL_ML = -np.min(nLL)
        
                

        
        
        # Compute fw messages
        (Sm_fw_ff, Amplitude, Phase, LL_ML) = fwMessages_CircularRayleighWave(X_ML, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource)

        
        Sm_fw[:,:,ff] = Sm_fw_ff        
        hist_LL[ii] = LL_ML
        if (ii > 0) and (abs(hist_LL[ii]-hist_LL[ii-1])/abs(hist_LL[ii-1]) < 0.01):
            break
        
    Wavenumber = X_ML[0]
    EllipticityAngle = X_ML[1]
    
    LogLikelihood = LL_ML
    NumParameters = 3 * L + 4
    NumPoints = K * L
    BIC = -2* LogLikelihood + Gamma* NumParameters *log(NumPoints)

    logging.debug('Circular Rayleigh wave fit')
    logging.debug('\tLL: {0:.3e} BIC: {1:.3e}'.format(LogLikelihood, BIC))
    logging.debug('\tAmplitude: {0:.3f}, Phase: {1:.3f}, Wavenumber: {2:.3f}, EllipticityAngle: {3:.3f}'.format(Amplitude,Phase,Wavenumber,EllipticityAngle))
    
    return(BIC, Sm_fw, Amplitude, Phase, X_ML)
    


### Dissipative Rayleigh wave

def negLL_CircularDissipativeRayleighWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    """ Computes the loglikelihood of a Dissipative Rayleigh wave with circular wavefront.
 
    Parameters
    ----------
    X_grid : float array
        Array of size (3,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber. Second row the ellipticity angle. Third the imaginary part of the wavenumber.

    Sw_bw, Swm_bw, SlnGamma_bw : float array
        sufficient statistics

    ArrayInfo : float array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
 
    Returns
    -------
    negLogLikelihood : float array
        A one dimensional array of N elements with the value of minus the log-likelihood.
"""


    
    if shape(shape(X_grid))[0] == 1: # single point
        X_grid = reshape(X_grid, (3,-1))
   
    N_points=shape(X_grid)[1]
    if N_points == 1:
        Wavenumber = X_grid[0]
        EllipticityAngle=X_grid[1]
        Wavenumber_i = X_grid[2]
    else:
        Wavenumber = X_grid[0,:]
        EllipticityAngle=X_grid[1,:]
        Wavenumber_i = X_grid[2,:]
        
    

    L=shape(Swm_bw)[1]
    negLogLikelihood = zeros((N_points,))
    Wavenumber = Wavenumber + 1j*Wavenumber_i
    pos = ArrayInfo[:,0:2]
    comp = ArrayInfo[:,3]
    r_source = cdist(ArrayInfo[:,0:2], np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    

    Uw = zeros((N_points))
    Uwm = zeros((2,N_points))

    # EWX
    ndx_e = arange(0, L)[comp == EWX]
    L_e = len(ndx_e) # number of EWX components
    alpha = zeros((N_points,))
    phi = zeros((N_points, ))
    for ll_e in range(0, L_e):
        ll = ndx_e[ll_e]
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
        alpha = sin(EllipticityAngle)*cos(Azimuth)*np.abs(hankel_tmp)
        phi = np.angle(hankel_tmp)
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
    
    # NSY
    ndx_n = arange(0, L)[comp == NSY]
    L_n = len(ndx_n) # number of NSY components
    alpha = zeros((N_points,))
    phi = zeros((N_points, ))
    for ll_n in range(0, L_n):
        ll = ndx_n[ll_n]
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
        alpha = sin(EllipticityAngle)*sin(Azimuth)*np.abs(hankel_tmp)
        phi = np.angle(hankel_tmp)
        Uw[:] += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
        
    # UDZ
    ndx_z = arange(0, L)[comp == UDZ]
    L_z = len(ndx_z) # number of NSY components
    alpha = zeros((N_points,))
    phi = zeros((N_points,))
    for ll_z in range(0, L_z):
        ll = ndx_z[ll_z]
        hankel_tmp = 1j*hankel2(0,2*np.pi*Wavenumber*r_source[ll])
        alpha = cos(EllipticityAngle)*np.abs(hankel_tmp)
        phi = np.angle(hankel_tmp)
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )

    
    Um = zeros((2,N_points))
    Um[0,:] = Uwm[0,:] / Uw[:]
    Um[1,:] = Uwm[1,:] / Uw[:]
    negLogLikelihood[:] = -0.5*(Um[0,:] * Uwm[0,:] + Um[1,:] * Uwm[1,:]) + sum(-SlnGamma_bw)
        
    return(negLogLikelihood)
 
    

        
def fwMessages_CircularDissipativeRayleighWave(X_ML, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    
    
    Wavenumber_r = X_ML[0]
    EllipticityAngle = X_ML[1]
    Wavenumber_i = X_ML[2]
    Wavenumber = Wavenumber_r + 1j*Wavenumber_i
    
    L=shape(ArrayInfo)[0]
    comp=ArrayInfo[:,3]
    pos = ArrayInfo[:,0:2]
    
    r_source = cdist(pos, np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    
    Sm_fw = zeros((2,L))
    Uw = matrix(zeros((2,2)))
    Uwm = matrix(zeros((2,1)))
    
    # compute the LL
    for ll in range(0,L):
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        if comp[ll] == EWX:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * cos(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif comp[ll] == NSY:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * sin(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif comp[ll] == UDZ:
            hankel_tmp = 1j*hankel2(0,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = cos(EllipticityAngle) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        
        Uw = Uw + transpose(H) * matrix(Sw_bw[:,:,ll]) * H
        
        Uwm = Uwm + transpose(H) * transpose(matrix(Swm_bw[:,ll]))
    Um = np.linalg.solve(Uw, Uwm)
    LL = 0.5*transpose(Um) * Uwm + sum(SlnGamma_bw)
    LL = LL[0,0]
    
    # compute the forward messages Sm_fw
    for ll in range(0,L):
        Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        if comp[ll] == EWX:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * cos(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)], [sin(phi), cos(phi)]])
        elif  comp[ll] == NSY:
            hankel_tmp = -1j*hankel2(1,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = sin(EllipticityAngle) * sin(Azimuth) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == UDZ:
            hankel_tmp = 1j*hankel2(0,2*np.pi*Wavenumber*r_source[ll])
            phi = np.angle(hankel_tmp)
            alpha = cos(EllipticityAngle) * np.abs(hankel_tmp)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        tmp = H * Um
        Sm_fw[0,ll] = tmp[0]
        Sm_fw[1,ll] = tmp[1]
    
    Amplitude = sqrt(Um[0][0,0]**2 + Um[1][0,0]**2)
    Phase = arctan2(Um[0][0,0], Um[1][0,0])

    return (Sm_fw, Amplitude, Phase, LL)
    
    


 
def fitCircularDissipativeRayleighWave(X_grid, Sm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, ShotSource, Gamma):
    
    Sm_fw = zeros(shape(Sm_bw))
    X_grid= X_grid.reshape((3,-1)) # can have grid or single point as input
    
    

    
    MaxIterations=5
    hist_LL=zeros((MaxIterations,))
    for ii in range(0,MaxIterations):
    
        (sigma2, SlnGamma_bw) = estimateNoise(Sm_bw, Sm_fw)
    
        Snw_bw = zeros((2,2,L)); Snwm_bw = zeros((2,L));
        for ll in range(0,L):
            Snw_bw[:,:,ll] = Sw_bw[:,:,ll,ff] / sigma2[ll];
            Snwm_bw[0,ll] = dot( Snw_bw[0,:,ll] , Sm_bw[:,ll,ff] )
            Snwm_bw[1,ll] = dot( Snw_bw[1,:,ll] , Sm_bw[:,ll,ff] )
            
       
        
        if ii == 0:
            # ML estimate with grid search
            (nLL) = negLL_CircularDissipativeRayleighWave(X_grid, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource)
                        
            ndx_min = argmin(nLL)
            X_g = X_grid[:,ndx_min];
           



        # Refine ML estimate
       
        # TODO 
        # L-BFGS-B often fails giving ABNORMAL_TERMINATION_IN_LNSRCH
        # therefore we use the unconstrained method Nelder-Mead
        # perhaps it would be nice to try L-BFGS-B providing explicit gradient
        
        # res = minimize(negLL_CircularRayleighWave, X_g, (Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource),'L-BFGS-B', bounds=[(0, Kmax), (-pi, pi)])
        # if res.success:
        #     X_ML = res.x
        #     LL_ML = -res.fun[0]
        # else:
        #     logging.warning(res)
        #     X_ML = X_g
        #     LL_ML = -np.min(nLL)
        
        res = minimize(negLL_CircularDissipativeRayleighWave, X_g, (Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource),'Nelder-Mead')
        if res.success:
            X_ML = res.x
            LL_ML = -res.fun
        else:
            logging.warning(res)
            X_ML = X_g
            LL_ML = -np.min(nLL)
        
                

        
        
        # Compute fw messages
        (Sm_fw_ff, Amplitude, Phase, LL_ML) = fwMessages_CircularDissipativeRayleighWave(X_ML, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo, ShotSource)

        
        Sm_fw[:,:,ff] = Sm_fw_ff        
        hist_LL[ii] = LL_ML
        if (ii > 0) and (abs(hist_LL[ii]-hist_LL[ii-1])/abs(hist_LL[ii-1]) < 0.01):
            break
        
    Wavenumber_r = X_ML[0]
    EllipticityAngle = X_ML[1]
    Wavenumber_i = X_ML[2]
    Wavenumber = Wavenumber_r + 1j*Wavenumber_i
    
    LogLikelihood = LL_ML
    NumParameters = 3 * L + 4
    NumPoints = K * L
    BIC = -2* LogLikelihood + Gamma* NumParameters *log(NumPoints)

    logging.debug('Circular Dissipative Rayleigh wave fit')
    logging.debug('\tLL: {0:.3e} BIC: {1:.3e}'.format(LogLikelihood, BIC))
    logging.debug('\tAmplitude: {0:.3f}, Phase: {1:.3f}, Wavenumber: {2:.3f}, EllipticityAngle: {3:.3f}'.format(Amplitude,Phase,Wavenumber,EllipticityAngle))
    
    return(BIC, Sm_fw, Amplitude, Phase, X_ML)
    




### Old code using a SSM to model Bessel differential equation as presented in
# Stefano Maranò, Donat Fäh, and Hans-Andrea Loeliger, "A state-space approach 
# for the analysis of wave and diffusion fields", in Proceedings IEEE International 
# Conference Acoustics, Speech, and Signal Processing, Brisbane, Australia, Apr. 2015, pp. 2564–2568. 
# doi:10.1109/ICASSP.2015.7178434


def negLL_HelmholtzCylindrical_KrKi(X, S_Wm, S_W, r_x, ReturnState=False):
    """
        X defines the input with real and imaginary wavenumber
        values are in [1/m] and NOT [rad/m]
        
        r is a vector containing the distances from the source
        the state of the first element is returned
    """
    
    X = np.reshape(X,(2,-1))
    N_s = len(r_x)
    N_points = np.shape(X)[1]
    Kr_vec = X[0,:]
    Ki_vec = X[1,:]
    
    assert all(Kr_vec>0)
    assert all(Ki_vec>=0)
    
    LL=np.zeros((N_points,))
    if ReturnState:
        ML_Um=np.zeros((4,N_points))
    
    for kk in range(0, N_points):
        K_r = Kr_vec[kk]
        K_i = Ki_vec[kk]

        Wavenumber =  2*np.pi*(K_r -1j*K_i)
        
        U_wm=np.matrix(np.zeros((4,1)))
        U_w=np.matrix(np.zeros((4,4)))
        
        
        for ss in range(0,N_s):
            
            X_w = np.matrix( np.zeros((4,4)))
            X_wm = np.matrix( np.zeros((4,)))
            # The readout matrix is [1 0 0 0; 0 0 1 0]
            X_wm[0,0] = S_Wm[0,ss]
            X_wm[0,2] = S_Wm[1,ss]
            X_w[0,0] = S_W[0,0,ss]
            X_w[0,2] = S_W[0,1,ss]
            X_w[2,0] = S_W[1,0,ss]
            X_w[2,2] = S_W[1,1,ss]
            
            r0 = np.min(r_x)            # reference receiver
            r1 = r_x[ss]                # current receiver
            #print(r0)
            #print(r1)
            #print(Wavenumber)

            # let Phi be a Fundamental matrix
            # Phi(r_1) = A*Phi(r_0) -> A = Phi(r_1)*inv(Phi(r_0))

#            Phi_r0 = np.matrix([ [yn(0,Wavenumber*r0), jn(0,Wavenumber*r0)],
#                                 [-Wavenumber*yn(1,Wavenumber*r0), -Wavenumber*jn(1,Wavenumber*r0)]])
#            Phi_r1 = np.matrix([ [yn(0,Wavenumber*r1), jn(0,Wavenumber*r1)],
#                                 [-Wavenumber*yn(1,Wavenumber*r1), -Wavenumber*jn(1,Wavenumber*r1)]])
                                 
            Phi_r0 = np.matrix([ [hankel1(0,Wavenumber*r0), hankel2(0,Wavenumber*r0)],
                                 [-Wavenumber*hankel1(1,Wavenumber*r0), -Wavenumber*hankel2(1,Wavenumber*r0)]])
            Phi_r1 = np.matrix([ [hankel1(0,Wavenumber*r1), hankel2(0,Wavenumber*r1)],
                                 [-Wavenumber*hankel1(1,Wavenumber*r1), -Wavenumber*hankel2(1,Wavenumber*r1)]])
                                   
            if any(isnan(Phi_r0)):
                print(r0)
                print(Wavenumber)
                print(Phi_r0)
                
            if any(isnan(Phi_r1)):
                print(r1)
                print(Wavenumber)
                print(Phi_r1)
            
            A = Phi_r1 * Phi_r0.I 
            
            A_real = np.bmat([[A.real, -A.imag], [A.imag, A.real]])
        
            U_w += A_real.T * X_w * A_real
            U_wm += A_real.T * X_wm.T
            
            #print(U_w)
            #print(U_wm)
    
        # with the matrix G, change the state vector space
        G = np.matrix([ [1, 1],
                        [-hankel1(1,Wavenumber*r0)*Wavenumber/hankel1(0,Wavenumber*r0), -hankel2(1,Wavenumber*r0)*Wavenumber/hankel2(0,Wavenumber*r0)]])
        
        G_real = np.bmat([[G.real, -G.imag], [G.imag, G.real]])

        
        U_w = G_real.T * U_w * G_real
        U_wm = G_real.T * U_wm
        U_m = U_w.I * U_wm
        
        # here we force the existence of a single outgoing wave, not two.        
        my_Um = np.copy(U_m)
        my_Um[0] = 0
        my_Um[2] = 0
        
        LL[kk] += -0.5* my_Um.T*U_w*my_Um + my_Um.T * U_wm
        if ReturnState:
            ML_Um[:,kk] = np.reshape(U_m.getA(),(4,)) # save the state at x=0

    if ReturnState:
        RetVal = ML_Um
    else:
        RetVal = -LL
    
    return RetVal
    
    



###########################
###  Old code, not used ###
###########################


    
def negLL_CircularVerticalWave_old(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    """ Computes the loglikelihood of a circular wave (vertical component only)
 
    Parameters
    ----------
    X_grid : float array
        Array of size (2,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber along the x axes. Second row the wavenumber along the y axes.

    Sw_bw, Swm_bw, SlnGamma_bw :
        sufficient statistics

    ArrayInfo : float array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
 
    Returns
    -------
    negLogLikelihood : float array
        A one dimensional array of N elements with the value of minus the log-likelihood.

    """
    
    
    N_points=size(X_grid)
    
    

    if N_points == 1:
        Wavenumber_vec = np.array([np.asscalar(X_grid)])
    else:
        Wavenumber_vec = X_grid[0,:]
        
    L=shape(Swm_bw)[1]
    negLogLikelihood = LogLikelihood = zeros((N_points,))
    comp=ArrayInfo[:,3]
    r_source = cdist(ArrayInfo[:,0:2], np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    

    

    for nn in range(0, N_points):
        Wavenumber =  2*np.pi*Wavenumber_vec[nn]
        
        if Wavenumber == 0:
            LogLikelihood[nn] = - sys.maxsize
        else:

        
            U_wm=np.matrix(np.zeros((4,1)))
            U_w=np.matrix(np.zeros((4,4)))
            
            
            for ll in range(0,L):
                if comp[ll] == UDZ:
                
                    X_w = np.matrix( np.zeros((4,4)))
                    X_wm = np.matrix( np.zeros((4,)))
                    # The readout matrix is [1 0 0 0; 0 0 1 0]
                    X_wm[0,0] = Swm_bw[0,ll]
                    X_wm[0,2] = Swm_bw[1,ll]
                    X_w[0,0] = Sw_bw[0,0,ll]
                    X_w[0,2] = Sw_bw[0,1,ll]
                    X_w[2,0] = Sw_bw[1,0,ll]
                    X_w[2,2] = Sw_bw[1,1,ll]
                    
                    r0 = np.min(r_source)       # reference receiver
                    r1 = r_source[ll]           # current receiver
                    
                           
                    Phi_r0 = np.matrix([ [hankel1(0,Wavenumber*r0), hankel2(0,Wavenumber*r0)],
                                         [-Wavenumber*hankel1(1,Wavenumber*r0), -Wavenumber*hankel2(1,Wavenumber*r0)]])
                    Phi_r1 = np.matrix([ [hankel1(0,Wavenumber*r1), hankel2(0,Wavenumber*r1)],
                                         [-Wavenumber*hankel1(1,Wavenumber*r1), -Wavenumber*hankel2(1,Wavenumber*r1)]])
                                           
                 
                    A = Phi_r1 * Phi_r0.I 
                    
                    A_real = np.bmat([[A.real, -A.imag], [A.imag, A.real]])
                
                    U_w += A_real.T * X_w * A_real
                    U_wm += A_real.T * X_wm.T

                
                #print(U_w)
                #print(U_wm)
        
            # with the matrix G, change the state vector space
            G = np.matrix([ [1, 1],
                            [-hankel1(1,Wavenumber*r0)*Wavenumber/hankel1(0,Wavenumber*r0), -hankel2(1,Wavenumber*r0)*Wavenumber/hankel2(0,Wavenumber*r0)]])
            
            G_real = np.bmat([[G.real, -G.imag], [G.imag, G.real]])
    
            
            U_w = G_real.T * U_w * G_real
            U_wm = G_real.T * U_wm
            U_m = U_w.I * U_wm
            
            # here we force the existence of a single outgoing wave, not two.        
            my_Um = np.copy(U_m)
            my_Um[0] = 0
            my_Um[2] = 0
            
            LogLikelihood[nn] += -0.5* my_Um.T*U_w*my_Um + my_Um.T * U_wm + sum(SlnGamma_bw)
    #if ReturnState:
    #    ML_Um[:,kk] = np.reshape(my_Um.getA(),(4,)) # save the state at x=0

    
    negLogLikelihood = -LogLikelihood 
    return negLogLikelihood
    
    
def fwMessages_CircularVerticalWave_old(X_ML, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    """
    
    RETURN
        amplitude and phase returned refer to outgoing wave. It is assumed there is no ingoing wave.
    """
    
    
    Wavenumber = 2*np.pi*X_ML
    
    L=shape(ArrayInfo)[0]
    comp=ArrayInfo[:,3]
    
    Sm_fw = zeros((2,L))
    
    LogLikelihood = 0
    r_source = cdist(ArrayInfo[:,0:2], np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
        
    # first estimate state vector and force ingoing wave to zero        
    U_wm=np.matrix(np.zeros((4,1)))
    U_w=np.matrix(np.zeros((4,4)))
    C_real = np.matrix([[1, 0, 0, 0], [0, 0, 1, 0]]) # readout matrix

    for ll in range(0,L):
        if comp[ll] == UDZ:
        
            X_w = np.matrix( np.zeros((4,4)))
            X_wm = np.matrix( np.zeros((4,)))
            # The readout matrix is [1 0 0 0; 0 0 1 0]
            X_wm[0,0] = Swm_bw[0,ll]
            X_wm[0,2] = Swm_bw[1,ll]
            X_w[0,0] = Sw_bw[0,0,ll]
            X_w[0,2] = Sw_bw[0,1,ll]
            X_w[2,0] = Sw_bw[1,0,ll]
            X_w[2,2] = Sw_bw[1,1,ll]
            
            r0 = np.min(r_source)            # reference receiver
            r1 = r_source[ll]                # current receiver
                                 
            Phi_r0 = np.matrix([ [hankel1(0,Wavenumber*r0), hankel2(0,Wavenumber*r0)],
                                 [-Wavenumber*hankel1(1,Wavenumber*r0), -Wavenumber*hankel2(1,Wavenumber*r0)]])
            Phi_r1 = np.matrix([ [hankel1(0,Wavenumber*r1), hankel2(0,Wavenumber*r1)],
                                 [-Wavenumber*hankel1(1,Wavenumber*r1), -Wavenumber*hankel2(1,Wavenumber*r1)]])
                                   
         
            A = Phi_r1 * Phi_r0.I 
            
            A_real = np.bmat([[A.real, -A.imag], [A.imag, A.real]])
        
            U_w += A_real.T * X_w * A_real
            U_wm += A_real.T * X_wm.T


    # with the matrix G, change the state vector space
    G = np.matrix([ [1, 1],
                    [-hankel1(1,Wavenumber*r0)*Wavenumber/hankel1(0,Wavenumber*r0), -hankel2(1,Wavenumber*r0)*Wavenumber/hankel2(0,Wavenumber*r0)]])
    G_real = np.bmat([[G.real, -G.imag], [G.imag, G.real]])

    U_w = G_real.T * U_w * G_real
    U_wm = G_real.T * U_wm
    U_m = U_w.I * U_wm
    
    # here we force the existence of a single out going wave, not two.        
    o_Um = np.copy(U_m)
    o_Um[0] = 0
    o_Um[2] = 0
    
    LogLikelihood = -0.5* o_Um.T*U_w*o_Um + o_Um.T * U_wm + sum(SlnGamma_bw)
    
    ML_Um = o_Um # state at position r0, only outgoing wave
    
    Amplitude = np.abs( (ML_Um[1] + 1j*ML_Um[3]) / hankel2(0,Wavenumber*r0) )
    Phase = np.angle( (ML_Um[1] + 1j*ML_Um[3]) / hankel2(0,Wavenumber*r0) )
        
    #print(ML_Um)    
    #print(np.sqrt(ML_Um[1]**2 + ML_Um[3]**2))
    
    # second, compute forward messages Sm_fw
    for ll in range(0,L):
        if comp[ll] == UDZ:
            
            r0 = np.min(r_source)            # reference receiver
            r1 = r_source[ll]                # current receiver
                                 
            Phi_r0 = np.matrix([ [hankel1(0,Wavenumber*r0), hankel2(0,Wavenumber*r0)],
                                 [-Wavenumber*hankel1(1,Wavenumber*r0), -Wavenumber*hankel2(1,Wavenumber*r0)]])
            Phi_r1 = np.matrix([ [hankel1(0,Wavenumber*r1), hankel2(0,Wavenumber*r1)],
                                 [-Wavenumber*hankel1(1,Wavenumber*r1), -Wavenumber*hankel2(1,Wavenumber*r1)]])
            
            A = Phi_r1 * Phi_r0.I
            A_real = np.bmat([[A.real, -A.imag], [A.imag, A.real]])
            
            Sm_fw[:,ll] = np.reshape( C_real * A_real * G_real * ML_Um, (2,) )
            

        
    return (Sm_fw, Amplitude[0], Phase[0], LogLikelihood[0,0])
    
    
def negLL_CircularRayleighWave_old(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo, ShotSource):
    """ Computes the loglikelihood of a Rayleigh wave
 
    Parameters
    ----------
    X_grid : float array
        Array of size (2,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber. Second row the ellipticity angle.

    Sw_bw, Swm_bw, SlnGamma_bw : float array
        sufficient statistics

    ArrayInfo : float array
        An array containing information about sensor location and channels.
        It has L rows, where L is the number of channels.
        Each rows has the following form:
          pos_x, pos_y, pos_z, cmp, Ts
        where the first three fields are the position of the sensor in [m].
        cmp is the component code.
        Ts is the sampling time as read from the SAC file.
 
    Returns
    -------
    negLogLikelihood : float array
        A one dimensional array of N elements with the value of minus the log-likelihood.
"""

    
    if shape(shape(X_grid))[0] == 1: # single point
        X_grid = reshape(X_grid, (2,-1))
   
    N_points=shape(X_grid)[1]
    if N_points == 1:
        Wavenumber_vec = X_grid[0]
        EllipticityAngle_vec=X_grid[1]
    else:
        Wavenumber_vec = X_grid[0,:]
        EllipticityAngle_vec=X_grid[1,:]

    L=shape(Swm_bw)[1]
    LogLikelihood = zeros((N_points,))
    pos = ArrayInfo[:,0:2]
    comp = ArrayInfo[:,3]
    r_source = cdist(ArrayInfo[:,0:2], np.array([ShotSource[0:2]]))[:,0] # distance of each station from source
    r0 = np.min(r_source)            # reference receiver

    

    for nn in range(0, N_points):
        Wavenumber =  2*np.pi*Wavenumber_vec[nn]
        EllipticityAngle = EllipticityAngle_vec[nn]
        
        Phi_r0 = np.array([ [hankel1(0,Wavenumber*r0), hankel2(0,Wavenumber*r0)],
                  [-Wavenumber*hankel1(1,Wavenumber*r0), -Wavenumber*hankel2(1,Wavenumber*r0)]])
                  
        C_real = np.zeros((2,4))

        
        if Wavenumber == 0:
            LogLikelihood[nn] = - sys.maxsize
        else:
        
            U_wm = np.zeros((4,))
            U_w = np.zeros((4,4))
            A = np.zeros((2,2), dtype=complex)
            A_real = np.zeros((4,4))
            G = np.zeros((2,2), dtype=complex)
            G_real = np.zeros((4,4))
            
            for ll in range(0,L):
                r1 = r_source[ll] # distance from source of current receiver
   
                # Naive implementation
                #Phi_r0 = np.array([ [hankel1(0,Wavenumber*r0), hankel2(0,Wavenumber*r0)],
                #                     [-Wavenumber*hankel1(1,Wavenumber*r0), -Wavenumber*hankel2(1,Wavenumber*r0)]])
                #Phi_r1 = np.array([ [hankel1(0,Wavenumber*r1), hankel2(0,Wavenumber*r1)],
                #                     [-Wavenumber*hankel1(1,Wavenumber*r1), -Wavenumber*hankel2(1,Wavenumber*r1)]])
                #A = Phi_r1 * Phi_r0.I
                                       
                A = np.dot(np.array([[hankel1(0,Wavenumber*r1), hankel2(0,Wavenumber*r1)], [-Wavenumber*hankel1(1,Wavenumber*r1), -Wavenumber*hankel2(1,Wavenumber*r1)]]),
                                    np.array([[Phi_r0[1,1], -Phi_r0[0,1]],[-Phi_r0[1,0], Phi_r0[0,0]]]) /(Phi_r0[0,0]*Phi_r0[1,1]-Phi_r0[0,1]*Phi_r0[1,0]))
                
                #A_real = np.bmat([[A.real, -A.imag], [A.imag, A.real]])
                A_real[0:2,0:2] = A_real[2:4,2:4] = A.real
                A_real[0:2,2:4] = -A.imag
                A_real[2:4,0:2] = A.imag
                
                Azimuth = np.arctan2(pos[ll,1]-ShotSource[1], pos[ll,0]-ShotSource[0])
        
                if comp[ll] == EWX:
                    C_real = np.array([[sin(EllipticityAngle)*cos(Azimuth), 0, 0, 0], [0, 0, sin(EllipticityAngle)*cos(Azimuth), 0]])
                elif comp[ll] == NSY:
                    C_real = np.array([[sin(EllipticityAngle)*sin(Azimuth), 0, 0, 0], [0, 0, sin(EllipticityAngle)*sin(Azimuth), 0]])
                elif comp[ll] == UDZ:
                    C_real = np.array([[0, 0, -cos(EllipticityAngle), 0], [cos(EllipticityAngle), 0, 0, 0]]) # pi/2 phase shift

                CA = np.dot(C_real, A_real)
                U_w +=  np.transpose(CA).dot(Sw_bw[:,:,ll]).dot(CA)
                U_wm += np.transpose(CA).dot(Swm_bw[:,ll])


            # with the matrix G, change the state vector space
            G = np.array([[1, 1],
                          [-hankel1(1,Wavenumber*r0)*Wavenumber/hankel1(0,Wavenumber*r0), -hankel2(1,Wavenumber*r0)*Wavenumber/hankel2(0,Wavenumber*r0)]])
            
            # G_real = np.bmat([[G.real, -G.imag], [G.imag, G.real]])
            G_real[0:2,0:2] = G_real[2:4,2:4] = G.real
            G_real[0:2,2:4] = -G.imag
            G_real[2:4,0:2] = G.imag

            U_w = np.transpose(G_real).dot(U_w).dot(G_real)
            U_wm = np.transpose(G_real).dot(U_wm)
            U_m = np.linalg.solve(U_w, U_wm) # faster than U_m = U_w.I * U_wm
            
            # here we force the existence of a single outgoing wave, not two.        
            my_Um = np.copy(U_m)
            my_Um[0] = 0
            my_Um[2] = 0
            
            LogLikelihood[nn] += -0.5* np.transpose(my_Um).dot(U_w).dot(my_Um) + np.transpose(my_Um).dot(U_wm) + np.sum(SlnGamma_bw)

    return -LogLikelihood
 