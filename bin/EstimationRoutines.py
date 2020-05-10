# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
Here are estimation routines used in WaveDec
"""

from scipy.optimize import minimize # , approx_fprime
import numpy as np
from numpy import shape, fft, zeros, real, imag, eye, linalg, ceil, linspace, pi, \
        meshgrid, concatenate, sqrt, array, sum, conjugate, log, mean, size, dot, \
        arctan2, sin, cos, argmin, reshape, matrix, transpose, arange
import DataUtils as db
import logging

from wdSettings import MODEL_NOISE, MODEL_VERTICAL, MODEL_RAYLEIGH, MODEL_LOVE
from wdSettings import EWX, NSY, UDZ, ROTX, ROTY, ROTZ




def decomposeWavefield(conn, y, WindowId, Ts, Fvec_ndx, Kmax, Kstep, Estep, Vmin, WavesToModel, MaxWaves, MaxIterations, ArrayInfo, Gamma):
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
    Kstep : float
        Grid step in the wavenumber plane [1/m]
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
    
    
    NumK = np.int(2*ceil(Kmax/Kstep)) # number of points in the wavenumber search grid. Even, so that we do not have 0 in the final search grid
    NumE = np.int(ceil(np.pi/Estep)) # number of points in the ellipticity search grid Even, so that we have 0 in the final search grid
    Wavenumber_x = linspace(-Kmax, Kmax, NumK)
    Wavenumber_y = linspace(-Kmax, Kmax, NumK)
    EllipticityAngle = linspace(-pi/2, pi/2, NumE, endpoint=False)

    xx, yy, zz = meshgrid(Wavenumber_x,Wavenumber_y,EllipticityAngle)
    xx=concatenate(xx); yy=concatenate(yy); zz=concatenate(zz);
    ndx_ok = sqrt( xx**2 + yy**2 ) <= Kmax # only wavenumber smaller than Kmax
    xx = xx[ndx_ok]; yy = yy[ndx_ok]; zz = zz[ndx_ok];
    X_grid_R = array([xx, yy, zz]) # Rayleigh waves
    
    xx, yy = meshgrid(Wavenumber_x,Wavenumber_y)
    xx=concatenate(xx); yy=concatenate(yy);
    ndx_ok = sqrt( xx**2 + yy**2 ) <= Kmax # only wavenumber smaller than Kmax
    xx = xx[ndx_ok]; yy = yy[ndx_ok];
    X_grid_L = array([xx, yy]) # Love waves
    
    # Fndx
    for ff in Fvec_ndx:
        X_grid_R_f = X_grid_R
        X_grid_L_f = X_grid_L
        X_grid_V_f = X_grid_L

        if Vmin != None and Vmin > 0: # further restrict search grid
            if WavesToModel[MODEL_LOVE] or WavesToModel[MODEL_VERTICAL]:
                ndx_ok = sqrt( X_grid_L[0,:]**2 + X_grid_L[1,:]**2 ) <= Fvec_fft[ff]/Vmin
                if WavesToModel[MODEL_LOVE]:
                    X_grid_L_f = X_grid_L[:,ndx_ok]
                if WavesToModel[MODEL_VERTICAL]:
                    X_grid_V_f = X_grid_L[:,ndx_ok]        
            if WavesToModel[MODEL_RAYLEIGH]:
                ndx_ok = sqrt( X_grid_R[0,:]**2 + X_grid_R[1,:]**2 ) <= Fvec_fft[ff]/Vmin 
                X_grid_R_f = X_grid_R[:,ndx_ok]
 
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
            if WavesToModel[MODEL_VERTICAL] and size(X_grid_V_f) > 0:
                (BIC_V, Sm_fw_V, Amplitude_V, Phase_V, X_ML_V) = fitVerticalWave(X_grid_V_f, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma)
                tmpBIC.append(BIC_V); tmpModel.append(MODEL_VERTICAL);
            if WavesToModel[MODEL_LOVE] and size(X_grid_L_f) > 0:
                (BIC_L, Sm_fw_L, Amplitude_L, Phase_L, X_ML_L) = fitLoveWave(X_grid_L_f, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma)
                tmpBIC.append(BIC_L); tmpModel.append(MODEL_LOVE);
            if WavesToModel[MODEL_RAYLEIGH] and size(X_grid_R_f) > 0:
                (BIC_R, Sm_fw_R, Amplitude_R, Phase_R, X_ML_R) = fitRayleighWave(X_grid_R_f, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma)
                tmpBIC.append(BIC_R); tmpModel.append(MODEL_RAYLEIGH);
            

            # Model selection: choose wave with smallest BIC
            if len(tmpBIC) > 0:
                WaveModel[mm] = tmpModel[np.argmin(tmpBIC)]
                histBIC = zeros((MaxIterations+1,))
                histBIC[0] = np.min(tmpBIC)
            
            if WaveModel[mm] == MODEL_NOISE:
                break     # when a noise model is chosen no more waves need to be added
            if WaveModel[mm] == 0:
                break     # Nothing was modeled. Perhaps because of stringent Vmin at this frequency
            elif WaveModel[mm] == MODEL_VERTICAL:
                Sm_fw[:,:,:,mm] = Sm_fw_V
                WaveAmplitude[mm] = Amplitude_V
                WavePhase[mm] = Phase_V
                WaveX_ML[mm] = X_ML_V
            elif WaveModel[mm] == MODEL_LOVE:
                Sm_fw[:,:,:,mm] = Sm_fw_L
                WaveAmplitude[mm] = Amplitude_L
                WavePhase[mm] = Phase_L
                WaveX_ML[mm] = X_ML_L
            elif WaveModel[mm] == MODEL_RAYLEIGH:
                Sm_fw[:,:,:,mm] = Sm_fw_R
                WaveAmplitude[mm] = Amplitude_R
                WavePhase[mm] = Phase_R
                WaveX_ML[mm] = X_ML_R
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
                        elif WaveModel[m] == MODEL_VERTICAL:
                            X = WaveX_ML[m]
                            (BIC_V, Sm_fw_V, Amplitude_V, Phase_v, X_ML_V) = fitVerticalWave(X, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma)
                            WaveX_ML[m] = X_ML_V
                            WaveAmplitude[m] = Amplitude_V
                            WavePhase[m] = Phase_V
                            Sm_fw[:,:,:,m] = Sm_fw_V
                        elif WaveModel[m] == MODEL_LOVE:
                            X = WaveX_ML[m]
                            (BIC_L, Sm_fw_L, Amplitude_L, Phase_L, X_ML_L) = fitLoveWave(X, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma)
                            WaveX_ML[m] = X_ML_L
                            WaveAmplitude[m] = Amplitude_L
                            WavePhase[m] = Phase_L
                            Sm_fw[:,:,:,m] = Sm_fw_L
                        elif WaveModel[m] == MODEL_RAYLEIGH:
                            X = WaveX_ML[m]
                            (BIC_R, Sm_fw_R, Amplitude_R, Phase_R, X_ML_R) = fitRayleighWave(X, tmpSm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma)
                            WaveX_ML[m] = X_ML_R
                            WaveAmplitude[m] = Amplitude_R
                            WavePhase[m] = Phase_R
                            Sm_fw[:,:,:,m] = Sm_fw_R
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
            elif WaveModel[mm] == MODEL_VERTICAL:
                Kx = WaveX_ML[mm][0]
                Ky = WaveX_ML[mm][1]
                Wavenumber = sqrt( Kx**2 + Ky**2)
                Azimuth = np.mod(arctan2(Ky, Kx), 2*pi) # Note the role reversal, as from documentation
                db.addVerticalWave(conn, WindowId, ff, WaveAmplitude[mm], WavePhase[mm], Wavenumber, Azimuth)
            elif WaveModel[mm] == MODEL_LOVE:
                Kx = WaveX_ML[mm][0]
                Ky = WaveX_ML[mm][1]
                Wavenumber = sqrt( Kx**2 + Ky**2)
                Azimuth = np.mod(arctan2(Ky, Kx), 2*pi) # Note the role reversal, as from documentation
                db.addLoveWave(conn, WindowId, ff, WaveAmplitude[mm], WavePhase[mm], Wavenumber, Azimuth)
            elif WaveModel[mm] == MODEL_RAYLEIGH:
                Kx = WaveX_ML[mm][0]
                Ky = WaveX_ML[mm][1]
                EllipticityAngle = np.mod(WaveX_ML[mm][2] + pi/2, pi)  -pi/2
                Wavenumber = sqrt( Kx**2 + Ky**2)
                Azimuth = np.mod(arctan2(Ky, Kx), 2*pi) # Note the role reversal, as from documentation
                db.addRayleighWave(conn, WindowId, ff, WaveAmplitude[mm], WavePhase[mm], Wavenumber, Azimuth,EllipticityAngle)
            else:
                logging.warning("Unrecognized wave model {0}".format(WaveModel[mm]))


    # after all freq have been processed
    (BIC_N, sigma2_ML) = fitNoise(Sm_bw - Sm_fw_all, L, K, Gamma)
    db.addNoise(conn, WindowId, sigma2_ML)
    
    # TODO after estimating noise again, do we want to iterate once more on wave params  ?  
        
    return

def bwMessages(y, Ts):
    """Compute sufficient statistics, assuming unitary variance
 
    Parameters
    ----------
    y : 2d float array
        It is an array of size (K, L). Each column contains the signal at the l-th location.
    Ts : float
        The sampling time in [s]
 
    Returns
    -------
    Sm_bw : float array
        It is an array of size (2, L). Each column refers to the the signal at the l-th location.
	Contains the mean vector of the message S.
    Sw_bw : float array
        It is an array of size (2, 2, L). Each page refers to the the signal at the l-th location.
	Contains the precision matrix of the message S.
    Swm_bw : float array
        It is an array of size (2, L). Each column refers to the the signal at the l-th location.
	Contains the weighted mean vector of the message S.

    """
    
    K = shape(y)[0]
    L = shape(y)[1]
    #Fvec = fft.fftfreq(K, Ts);
    #Fmin=0
    #Fmax=0.5/Ts
    #Fndx = (Fvec >= Fmin) and (Fvec <= Fmax);
    #Fvec=Fvec[Fndx];
    Fnum=int(K/2+1); # TODO need to make sure K is even
    
    Y_bw= fft.rfft(y, K, 0);

    Swm_bw = zeros((2,L,Fnum));
    Sm_bw= zeros((2,L,Fnum));
    Sw_bw=zeros((2,2,L,Fnum));

    for ff in range(0,Fnum):
        for ll in range(0, L):
            Swm_bw[0,ll,ff] = real(Y_bw[ff,ll]);
            Swm_bw[1,ll,ff] = imag(Y_bw[ff,ll]);
            Sw_bw[:,:,ll,ff] = (K/2)*eye(2);
            Sm_bw[:,ll,ff] = linalg.solve(Sw_bw[:,:,ll,ff], Swm_bw[:,ll,ff]);
    
    return (Sm_bw, Sw_bw, Swm_bw)

### Vertical wave

def negLL_VerticalWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    """ Computes the loglikelihood of a plane wave (vertical component only)
 
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
    
    if shape(shape(X_grid))[0] == 1: # single point
        X_grid = reshape(X_grid, (2,-1))
        
    N_points=shape(X_grid)[1]
    Wavenumber_x=X_grid[0,:]
    Wavenumber_y=X_grid[1,:]

    L=shape(Swm_bw)[1]

    negLogLikelihood = zeros((N_points,))
    
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]

    if False: # TODO check if input matrices are const*identity
        for nn in range(0,N_points):
            Uw = matrix(zeros((2,2)))
            Uwm = matrix(zeros((2,1)))
            for ll in range(0,L):
                if comp[ll] == UDZ:
                    phi = -2*pi*(Wavenumber_x[nn]*pos_x[ll] + Wavenumber_y[nn]*pos_y[ll])
                    H = matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
                else:
                    H = matrix([[0, 0],[0, 0]])
                    
                Uw = Uw + transpose(H) * matrix(Sw_bw[:,:,ll]) * H
                Uwm = Uwm + transpose(H) * transpose(matrix(Swm_bw[:,ll]))
        Um = linalg.solve(Uw, Uwm)
        negLogLikelihood[nn] = -0.5*transpose(Um) * Uwm + sum(-SlnGamma_bw)

    else:    # alternative, faster implementation. Requires Sw_bw to be const*Identity
        ndx_v = arange(0, L)[comp == UDZ]
        L_v = len(ndx_v) # number of vertical components
        phi = zeros((N_points, L_v))
        Uw = zeros((N_points))
        Uwm = zeros((2,N_points))
        for ll_v in range(0, L_v):
            ll = ndx_v[ll_v]
            phi[:,ll_v] = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll])
            Uw[:] = Uw[:] + Sw_bw[0,0,ll]
            Uwm[0,:] += + cos(phi[:,ll_v]) * Swm_bw[0,ll] + sin(phi[:,ll_v]) * Swm_bw[1,ll]
            Uwm[1,:] += - sin(phi[:,ll_v]) * Swm_bw[0,ll] + cos(phi[:,ll_v]) * Swm_bw[1,ll]
        Um = zeros((2,N_points))
        Um[0,:] = Uwm[0,:] / Uw[:]
        Um[1,:] = Uwm[1,:] / Uw[:]
        negLogLikelihood[:] = -0.5*(Um[0,:] * Uwm[0,:] + Um[1,:] * Uwm[1,:]) + sum(-SlnGamma_bw)
        
    return(negLogLikelihood)


def grad_negLL_VerticalWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    """ Computes the gradient of the loglikelihood of a plane wave (vertical component only)
 
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
    grad : float array
        The gradient.

    """
    
        
    
    Wavenumber_x=X_grid[0]
    Wavenumber_y=X_grid[1]

    L=shape(Swm_bw)[1]
    grad = zeros((2,))
    
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]
    
    if False: # TODO check if input matrices are const*identity


        # derivative wrt Wavenumber_x
        Uw = matrix(zeros((2,2)))
        Uwm = matrix(zeros((2,1)))
        dWx_Uw = matrix(zeros((2,2)))
        dWx_Uwm = matrix(zeros((2,1)))
        dWy_Uw = matrix(zeros((2,2)))
        dWy_Uwm = matrix(zeros((2,1)))
        for ll in range(0,L):
            if comp[ll] == UDZ:
                phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
                H = matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
                dWx_H = matrix([[-sin(phi), -cos(phi)],[cos(phi), -sin(phi)]]) * -2*pi*pos_x[ll]
                dWy_H = matrix([[-sin(phi), -cos(phi)],[cos(phi), -sin(phi)]]) * -2*pi*pos_y[ll]
            else:
                H = matrix([[0, 0],[0, 0]])
                dWx_H = matrix([[0, 0],[0, 0]])
                dWy_H = matrix([[0, 0],[0, 0]])
                
            Uw = Uw + transpose(H) * matrix(Sw_bw[:,:,ll]) * H
            Uwm = Uwm + transpose(H) * transpose(matrix(Swm_bw[:,ll]))

            dWx_Uw += transpose(dWx_H) * matrix(Sw_bw[:,:,ll]) * H + transpose(H) * matrix(Sw_bw[:,:,ll]) * dWx_H
            dWx_Uwm += transpose(dWx_H) * transpose(matrix(Swm_bw[:,ll]))
            
            dWy_Uw += transpose(dWy_H) * matrix(Sw_bw[:,:,ll]) * H + transpose(H) * matrix(Sw_bw[:,:,ll]) * dWy_H
            dWy_Uwm += transpose(dWy_H) * transpose(matrix(Swm_bw[:,ll]))

        Um = linalg.solve(Uw, Uwm)
        Uv = linalg.inv(Uw)
                
        grad[0] = -0.5*transpose(-Uv*dWx_Uw*Um + Uv*dWx_Uwm )* Uwm  - 0.5*transpose(Um)*dWx_Uwm
        grad[1] = -0.5*transpose(-Uv*dWy_Uw*Um + Uv*dWy_Uwm )* Uwm  - 0.5*transpose(Um)*dWy_Uwm
        
                
    else:    # alternative, faster implementation. Requires Sw_bw to be const*Identity
        ndx_v = arange(0, L)[comp == UDZ]
        L_v = len(ndx_v) # number of vertical components
        phi = zeros((L_v,))
        Uw = 0 # scalar variance
        Uwm = zeros((2,))
        dWx_Uwm = zeros((2,))
        dWy_Uwm = zeros((2,))
        for ll_v in range(0, L_v):
            ll = ndx_v[ll_v]
            phi[ll_v] = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            Uw = Uw + Sw_bw[0,0,ll]
            Uwm[0] += + cos(phi[ll_v]) * Swm_bw[0,ll] + sin(phi[ll_v]) * Swm_bw[1,ll]
            Uwm[1] += - sin(phi[ll_v]) * Swm_bw[0,ll] + cos(phi[ll_v]) * Swm_bw[1,ll]
            dWx_Uwm[0] += (- sin(phi[ll_v]) * Swm_bw[0,ll] + cos(phi[ll_v]) * Swm_bw[1,ll])*-2*pi*pos_x[ll]
            dWx_Uwm[1] += (- cos(phi[ll_v]) * Swm_bw[0,ll] - sin(phi[ll_v]) * Swm_bw[1,ll])*-2*pi*pos_x[ll]
            dWy_Uwm[0] += (- sin(phi[ll_v]) * Swm_bw[0,ll] + cos(phi[ll_v]) * Swm_bw[1,ll])*-2*pi*pos_y[ll]
            dWy_Uwm[1] += (- cos(phi[ll_v]) * Swm_bw[0,ll] - sin(phi[ll_v]) * Swm_bw[1,ll])*-2*pi*pos_y[ll]
        Um = zeros((2,))
        Um[0] = Uwm[0] / Uw
        Um[1] = Uwm[1] / Uw
        
        grad[0] = -0.5*( dWx_Uwm[0] * Um[0] + dWx_Uwm[1] * Um[1])  - 0.5*(Um[0]*dWx_Uwm[0]  + Um[1]*dWx_Uwm[1])
        grad[1] = -0.5*( dWy_Uwm[0] * Um[0] + dWy_Uwm[1] * Um[1])  - 0.5*(Um[0]*dWy_Uwm[0]  + Um[1]*dWy_Uwm[1])
        
        
    return(grad)


def fwMessages_VerticalWave(X_ML, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    
    Wavenumber_x = X_ML[0]
    Wavenumber_y = X_ML[1]
    
    L=shape(ArrayInfo)[0]
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]
    
    Sm_fw = zeros((2,L))
    
    Uw = matrix(zeros((2,2)))
    Uwm = matrix(zeros((2,1)))
    
    # compute the LL
    for ll in range(0,L):
        if comp[ll] == UDZ:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            H = matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        
        Uw += transpose(H) * matrix(Sw_bw[:,:,ll]) * H
        
        Uwm += transpose(H) * transpose(matrix(Swm_bw[:,ll]))
    Um = linalg.solve(Uw, Uwm)
    LL = 0.5*transpose(Um) * Uwm + sum(SlnGamma_bw)
    LL = LL[0,0]
    
    # compute the forward messages Sm_fw
    for ll in range(0,L):
        if comp[ll] == UDZ:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            H = matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        tmp = H * Um
        Sm_fw[0,ll] = tmp[0]
        Sm_fw[1,ll] = tmp[1]
    
    Amplitude = sqrt(Um[0][0,0]**2 + Um[1][0,0]**2)
    Phase = arctan2(Um[0][0,0], Um[1][0,0])

    return (Sm_fw, Amplitude, Phase, LL)

 
def fitVerticalWave(X_grid, Sm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma):
    
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
            (LL) = negLL_VerticalWave(X_grid, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
            
            ndx_min = argmin(LL)
            X_g = X_grid[:,ndx_min];
            

#            if shape(X_grid)[1] > 3:
#                plt.ion()
#                fig = plt.figure()
#                ax = fig.gca(projection='3d')
#                ax.plot_trisurf(X_grid[0], X_grid[1], -LL, cmap=cm.jet, linewidth=0.2)
#                plt.title("log likelihood " + str(ff) + "Hz")
#                plt.show()

        # Refine ML estimate
        res = minimize(negLL_VerticalWave, X_g, args=(Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo), jac=grad_negLL_VerticalWave, method='TNC', bounds=[(-Kmax, Kmax),(-Kmax, Kmax)])
        
        
        
        X_ML = res.x; LL_ML = -res.fun[0];
        
        # Compute fw messages
        (Sm_fw_ff, Amplitude, Phase, LL_ML) = fwMessages_VerticalWave(X_ML, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
        Sm_fw[:,:,ff] = Sm_fw_ff
        
        hist_LL[ii] = LL_ML
        if (ii > 0) and (abs(hist_LL[ii]-hist_LL[ii-1])/abs(hist_LL[ii-1]) < 0.01):
            break
        
        
    Wavenumber_x = X_ML[0]
    Wavenumber_y = X_ML[1]
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2)
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x)
    
    LogLikelihood = LL_ML
    NumParameters = 3 * L + 4
    NumPoints = K * L
    BIC = -2* LogLikelihood + Gamma *NumParameters *log(NumPoints)

    
    
    logging.debug('Vertical wave fit')
    logging.debug('\tLL: {0:.3e} BIC: {1:.3e}'.format(LogLikelihood, BIC))
    logging.debug('\tAmplitude: {0:.3f}, Phase: {1:.3f}, Wavenumber: {2:.3f}, Azimuth: {3:.3f}'.format(Amplitude,Phase,Wavenumber,Azimuth))

    return(BIC, Sm_fw, Amplitude, Phase, X_ML)
    
### Love wave

def negLL_LoveWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    """ Computes the likelihood of a Love wave
 
    Parameters
    ----------
    X_grid : float array
        Array of size (2,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber along the x axes. Second row the wavenumber along the y axes.

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
    Wavenumber_x=X_grid[0,:]    # wavenumber is measured in [rad / m]
    Wavenumber_y=X_grid[1,:]
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x)
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2) 

#
#    print(Wavenumber_x)    
#    print(Wavenumber_y)       
#    print(Azimuth)
#    print(Wavenumber)

    L=shape(Swm_bw)[1]

    negLogLikelihood = zeros((N_points,))
    
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]
    
    

    # EWX
    ndx_e = arange(0, L)[comp == EWX]
    L_e = len(ndx_e) # number of EWX components
    alpha = -sin(Azimuth)
    phi = zeros((N_points,))
    Uw = zeros((N_points,))
    Uwm = zeros((2,N_points))
    for ll_e in range(0, L_e):
        ll = ndx_e[ll_e]
        
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll])
        Uw[:] += + alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += + alpha *( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += + alpha *( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
    
    # NSY
    ndx_n = arange(0, L)[comp == NSY]
    L_n = len(ndx_n) # number of NSY components
    alpha = cos(Azimuth)
    phi = zeros((N_points,))
    for ll_n in range(0, L_n):
        ll = ndx_n[ll_n]
        
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll])
        Uw[:] += + alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += + alpha *( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += + alpha *( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )

    
    # ROTZ
    ndx_z = arange(0, L)[comp == ROTZ]
    L_z = len(ndx_z) # number of ROTZ components
    alpha = pi*Wavenumber
    phi = zeros((N_points,))
    for ll_z in range(0, L_z):
        ll = ndx_z[ll_z]
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll]) - pi/2
        Uw[:] += + alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += + alpha *( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += + alpha *( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
    
    Um = zeros((2,N_points))
    

    
#        print(Uw)
    Um[0,:] = Uwm[0,:] / Uw[:]
    Um[1,:] = Uwm[1,:] / Uw[:]
    negLogLikelihood[:] = -0.5*(Um[0,:] * Uwm[0,:] + Um[1,:] * Uwm[1,:]) + sum(-SlnGamma_bw)
        
    return(negLogLikelihood)
        

def grad_negLL_LoveWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    """ Computes the gradient of the loglikelihood of a Love wave
 
    Parameters
    ----------
    X_grid : float array
        Array of size (2,1). Description of N points to evaluate likelihoodfunction. First row is the wavenumber along the x axes. Second row the wavenumber along the y axes.

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
    grad : float array
        The gradient of the loglikelihood.
"""

    
       
    Wavenumber_x=X_grid[0]    # wavenumber is measured in [rad / m]
    Wavenumber_y=X_grid[1]
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x) # Note the role reversal
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2) 

    L=shape(Swm_bw)[1]

    grad = zeros((2,))
    
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]
    
    if Wavenumber_x == 0 and Wavenumber_y == 0: # singularity, Azimuth is not defined here
        return(grad)
        
    

    Uw = 0 # scalar variance
    dWx_Uw = 0
    dWy_Uw = 0
    Uwm = zeros((2,))
    dWx_Uwm = zeros((2,))
    dWy_Uwm = zeros((2,))

    # EWX
    ndx_e = arange(0, L)[comp == EWX]
    L_e = len(ndx_e) # number of EWX components
    alpha = -sin(Azimuth)
    dWx_alpha = -cos(Azimuth)* -Wavenumber_y/(Wavenumber_x**2 + Wavenumber_y**2) # TODO here there may be a problem with division by zero            
    dWy_alpha = -cos(Azimuth)* Wavenumber_x/(Wavenumber_x**2 + Wavenumber_y**2)
    for ll_e in range(0, L_e):
        ll = ndx_e[ll_e]
        
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]
    
    # NSY
    ndx_n = arange(0, L)[comp == NSY]
    L_n = len(ndx_n) # number of NSY components
    alpha = cos(Azimuth)
    dWx_alpha = -sin(Azimuth)* -Wavenumber_y/(Wavenumber_x**2 + Wavenumber_y**2)
    dWy_alpha = -sin(Azimuth)* Wavenumber_x/(Wavenumber_x**2 + Wavenumber_y**2)

    for ll_n in range(0, L_n):
        ll = ndx_n[ll_n]
        
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]            
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]

    # ROTZ
    ndx_z = arange(0, L)[comp == ROTZ]
    L_z = len(ndx_z) # number of ROTZ components
    alpha = pi*Wavenumber
    dWx_alpha = pi*Wavenumber_x/Wavenumber
    dWy_alpha = pi*Wavenumber_y/Wavenumber
    for ll_z in range(0, L_z):
        ll = ndx_z[ll_z]
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll]) - pi/2
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]
  
    
    Um = zeros((2,))
    Um[0] = Uwm[0] / Uw
    Um[1] = Uwm[1] / Uw
    
    
    grad[0] =  -0.5*( -(Um[0] * dWx_Uw * Um[0] + Um[1] *dWx_Uw * Um[1]) + dWx_Uwm[0] * Um[0] + dWx_Uwm[1] * Um[1])  - 0.5*(Um[0]*dWx_Uwm[0]  + Um[1]*dWx_Uwm[1])
    grad[1] =  -0.5*( -(Um[0] * dWy_Uw * Um[0] + Um[1] *dWy_Uw * Um[1]) + dWy_Uwm[0] * Um[0] + dWy_Uwm[1] * Um[1])  - 0.5*(Um[0]*dWy_Uwm[0]  + Um[1]*dWy_Uwm[1])
        
    return(grad)      
        
def fwMessages_LoveWave(X_ML, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    """ Computes the estimated wavefield of a Love wave.
 
    Parameters
    ----------
    X_ML : float array
        Array of size (2,1). The Maximum Likelihood estimate of the wave parameters.
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
    Sm_fw : float array
        Wavefield
    Amplitude : float
        Wave amplitude
    Phase : float
        Wave phase
    LL : float
        Model loglikelihood

    """
    
    Wavenumber_x = X_ML[0]
    Wavenumber_y = X_ML[1]
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x)
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2) 
    
    L=shape(ArrayInfo)[0]
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]

    
    Sm_fw = zeros((2,L))
    
    Uw = matrix(zeros((2,2)))
    Uwm = matrix(zeros((2,1)))
    for ll in range(0,L):    # compute the LL
        if comp[ll] == EWX:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = -sin(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == NSY:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = cos(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == ROTZ:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll]) -pi/2
            alpha = pi*Wavenumber
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        
        Uw = Uw + transpose(H) * matrix(Sw_bw[:,:,ll]) * H
        
        Uwm = Uwm + transpose(H) * transpose(matrix(Swm_bw[:,ll]))
    Um = linalg.solve(Uw, Uwm)
    LL = 0.5*transpose(Um) * Uwm + sum(SlnGamma_bw)
    LL = LL[0,0]
    
    for ll in range(0,L):    # compute the LL
        if comp[ll] == EWX:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = -sin(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == NSY:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = cos(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == ROTZ:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll]) -pi/2
            alpha = pi*Wavenumber
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        tmp = H * Um
        Sm_fw[0,ll] = tmp[0]
        Sm_fw[1,ll] = tmp[1]
    
    Amplitude = sqrt(Um[0][0,0]**2 + Um[1][0,0]**2)
    Phase = arctan2(Um[0][0,0], Um[1][0,0])

    return (Sm_fw, Amplitude, Phase, LL)

 
def fitLoveWave(X_grid, Sm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma):
    """ Finds the maximum likelihood fit a monochromatic Love wave
 
    Parameters
    ----------
    X_grid : float array
        Array of size (2,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber along the x axes. Second row the wavenumber along the y axes.

    Sw_bw, Swm_bw : float array
        sufficient statistics

    ff : int
        index of the frequency of interest

    L : int
        Number of channels

    K : int
        number of samples

    Kmax : float
        Largest wavenumber [1/m] to analyze

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
    BIC : float
        Value of the Bayesian information criterion
    Sm_fw : float array
        Wavefield generated by the estimated wave
    Amplitude : float
        Wave amplitude
    Phase : float
        Wave phase
    X_ML : float array
        Array of size (2,1). The Maximum Likelihood estimate of the wave parameters.

"""
    
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
            (LL) = negLL_LoveWave(X_grid, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
            
            ndx_min = argmin(LL)
            X_g = X_grid[:,ndx_min];

            ### test gradient
            #print("Love wave likelihood: testing gradient")
            #print(X_grid[:,0])
            #grad1 = approx_fprime(X_grid[:,0], negLL_LoveWave, sqrt(np.finfo(float).eps), Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
            #print(grad1)
            #grad2 = grad_negLL_LoveWave(X_grid[:,0], Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
            #print(grad2)
            #
            #if shape(X_grid)[1] > 3:
            #    from mpl_toolkits.mplot3d import Axes3D # this line is needed for projection='3d'
            #    from matplotlib import cm
            #    import matplotlib.pyplot as plt
            #    plt.ion()
            #    fig = plt.figure()
            #    ax = fig.gca(projection='3d')
            #    ax.view_init(90, 270)
            #    ax.plot_trisurf(X_grid[0], X_grid[1], -LL, cmap=cm.jet, linewidth=0.0)
            #    plt.title("Love wave log likelihood - " + str(ff))
            #    plt.show()
            #    plt.pause(5)
                

        # Refine ML estimate
        res = minimize(negLL_LoveWave, X_g, (Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo),'TNC', bounds=[(-Kmax, Kmax),(-Kmax, Kmax)], jac=grad_negLL_LoveWave)
        X_ML = res.x; LL_ML = -res.fun[0];
        
        
        #Kx = X_ML[0]
        #Ky = X_ML[1]
        #Wavenumber = sqrt( Kx**2 + Ky**2)
        #Azimuth = np.mod(arctan2(Ky, Kx), 2*pi)
        #print("Love wave ML estimate: {:.2e} {:.2e}".format(Wavenumber, Azimuth))

        
        # Compute fw messages
        (Sm_fw_ff, Amplitude, Phase, LL_ML) = fwMessages_LoveWave(X_ML, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
        Sm_fw[:,:,ff] = Sm_fw_ff
        
        hist_LL[ii] = LL_ML
        if (ii > 0) and (abs(hist_LL[ii]-hist_LL[ii-1])/abs(hist_LL[ii-1]) < 0.01):
            break
        
        
    Wavenumber_x = X_ML[0]
    Wavenumber_y = X_ML[1]
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2)
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x)
    
    LogLikelihood = LL_ML
    NumParameters = 3 * L + 4
    NumPoints = K * L
    BIC = -2* LogLikelihood + Gamma *NumParameters *log(NumPoints)

    
    
    logging.debug('Love wave fit')
    logging.debug('\tLL: {0:.3e} BIC: {1:.3e}'.format(LogLikelihood, BIC))
    logging.debug('\tAmplitude: {0:.3f}, Phase: {1:.3f}, Wavenumber: {2:.3f}, Azimuth: {3:.3f}'.format(Amplitude,Phase,Wavenumber,Azimuth))

    return(BIC, Sm_fw, Amplitude, Phase, X_ML)


### Rayleigh wave

def negLL_RayleighWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    """ Computes the loglikelihood of a Rayleigh wave
 
    Parameters
    ----------
    X_grid : float array
        Array of size (3,N). Description of N points to evaluate likelihoodfunction. First row is the wavenumber along the x axes. Second row the wavenumber along the y axes. Third row the ellipticity angle.

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
    Wavenumber_x=X_grid[0,:]    # wavenumber is measured in [rad / m]
    Wavenumber_y=X_grid[1,:]
    EllipticityAngle=X_grid[2,:]
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x)
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2) 

    L=shape(Swm_bw)[1]

    negLogLikelihood = zeros((N_points,))
    
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]
    
    Uw = zeros((N_points))
    Uwm = zeros((2,N_points))

    # EWX
    ndx_e = arange(0, L)[comp == EWX]
    L_e = len(ndx_e) # number of EWX components
    alpha = sin(EllipticityAngle)*cos(Azimuth)
    phi = zeros((N_points, ))
    for ll_e in range(0, L_e):
        ll = ndx_e[ll_e]
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll])
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
    
    # NSY
    ndx_n = arange(0, L)[comp == NSY]
    L_n = len(ndx_n) # number of NSY components
    alpha = sin(EllipticityAngle)*sin(Azimuth)
    phi = zeros((N_points, ))
    for ll_n in range(0, L_n):
        ll = ndx_n[ll_n]
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll])
        Uw[:] += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
        
    # UDZ
    ndx_z = arange(0, L)[comp == UDZ]
    L_z = len(ndx_z) # number of NSY components
    alpha = cos(EllipticityAngle)
    phi = zeros((N_points, ))
    for ll_z in range(0, L_z):
        ll = ndx_z[ll_z]
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll]) + pi/2
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )

    # ROTX
    ndx_x = arange(0, L)[comp == ROTX]
    L_x = len(ndx_x) # number of ROTX components
    alpha = 2*pi*Wavenumber*sin(Azimuth)*cos(EllipticityAngle)
    phi = zeros((N_points, ))
    for ll_x in range(0, L_x):
        ll = ndx_x[ll_x]
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll])
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )

    # ROTY
    ndx_y = arange(0, L)[comp == ROTY]
    L_y = len(ndx_y) # number of ROTY components
    alpha = -2*pi*Wavenumber*cos(Azimuth)*cos(EllipticityAngle)
    phi = zeros((N_points, ))
    for ll_y in range(0, L_y):
        ll = ndx_y[ll_y]
        phi = -2*pi*(Wavenumber_x[:]*pos_x[ll] + Wavenumber_y[:]*pos_y[ll])
        Uw[:] += alpha**2*Sw_bw[0,0,ll]
        Uwm[0,:] += alpha * ( + cos(phi) * Swm_bw[0,ll] + sin(phi) * Swm_bw[1,ll] )
        Uwm[1,:] += alpha * ( - sin(phi) * Swm_bw[0,ll] + cos(phi) * Swm_bw[1,ll] )
    
    Um = zeros((2,N_points))
    Um[0,:] = Uwm[0,:] / Uw[:]
    Um[1,:] = Uwm[1,:] / Uw[:]
    negLogLikelihood[:] = -0.5*(Um[0,:] * Uwm[0,:] + Um[1,:] * Uwm[1,:]) + sum(-SlnGamma_bw)
        
    return(negLogLikelihood)
 

 
def grad_negLL_RayleighWave(X_grid, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    """ Computes the gradient of the loglikelihood of a Rayleigh wave
 
    Parameters
    ----------
    X_grid : float array
        Array of size (3,1). Description of N points to evaluate likelihoodfunction. First row is the wavenumber along the x axes. Second row the wavenumber along the y axes. Third row the ellipticity angle.

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
    grad : float array
        The gradient of the loglikelihood.
"""

    
    Wavenumber_x=X_grid[0]    # wavenumber is measured in [rad / m]
    Wavenumber_y=X_grid[1]
    EllipticityAngle=X_grid[2]
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x) # Note the role reversal
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2) 

    L=shape(Swm_bw)[1]

    grad = zeros((3,))
    
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]
    
    if Wavenumber_x == 0 and Wavenumber_y == 0: # singularity, Azimuth is not defined here
        return(grad)
    
    Uw = 0 # scalar variance
    dWx_Uw = 0
    dWy_Uw = 0
    dEl_Uw = 0
    Uwm = zeros((2,))
    dWx_Uwm = zeros((2,))
    dWy_Uwm = zeros((2,))
    dEl_Uwm = zeros((2,))

    # EWX
    ndx_e = arange(0, L)[comp == EWX]
    L_e = len(ndx_e) # number of EWX components
    alpha = sin(EllipticityAngle)*cos(Azimuth)
    dWx_alpha = sin(EllipticityAngle)*-sin(Azimuth)* -Wavenumber_y/(Wavenumber_x**2 + Wavenumber_y**2)
    dWy_alpha = sin(EllipticityAngle)*-sin(Azimuth)* Wavenumber_x/(Wavenumber_x**2 + Wavenumber_y**2)
    dEl_alpha = cos(EllipticityAngle)*cos(Azimuth)
    dEl_phi = 0
    for ll_e in range(0, L_e):
        ll = ndx_e[ll_e]
        
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dEl_Uwm[0] += dEl_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dEl_phi
        dEl_Uwm[1] += dEl_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dEl_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]
        dEl_Uw += 2*alpha*dEl_alpha*Sw_bw[0,0,ll]
    
    # NSY
    ndx_n = arange(0, L)[comp == NSY]
    L_n = len(ndx_n) # number of NSY components
    alpha = sin(EllipticityAngle)*sin(Azimuth)
    dWx_alpha = sin(EllipticityAngle)*cos(Azimuth)* -Wavenumber_y/(Wavenumber_x**2 + Wavenumber_y**2)
    dWy_alpha = sin(EllipticityAngle)*cos(Azimuth)* Wavenumber_x/(Wavenumber_x**2 + Wavenumber_y**2)
    dEl_alpha = cos(EllipticityAngle)*sin(Azimuth)
    dEl_phi = 0

    for ll_n in range(0, L_n):
        ll = ndx_n[ll_n]
        
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]            
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dEl_Uwm[0] += dEl_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dEl_phi
        dEl_Uwm[1] += dEl_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dEl_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]
        dEl_Uw += 2*alpha*dEl_alpha*Sw_bw[0,0,ll]
        
    # UDZ
    ndx_z = arange(0, L)[comp == UDZ]
    L_z = len(ndx_z) # number of UDZ components
    alpha = cos(EllipticityAngle)
    dWx_alpha = 0
    dWy_alpha = 0
    dEl_alpha = -sin(EllipticityAngle)
    dEl_phi = 0
    for ll_z in range(0, L_z):
        ll = ndx_z[ll_z]
        
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll]) + pi/2
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dEl_Uwm[0] += dEl_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dEl_phi
        dEl_Uwm[1] += dEl_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dEl_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]
        dEl_Uw += 2*alpha*dEl_alpha*Sw_bw[0,0,ll]

    # ROTX
    ndx_x = arange(0, L)[comp == ROTX]
    L_x = len(ndx_x) # number of NSY components
    alpha = 2*pi*Wavenumber*sin(Azimuth)*cos(EllipticityAngle)
    dWx_alpha = 2*pi*cos(EllipticityAngle)*( Wavenumber_x/Wavenumber*sin(Azimuth)+ Wavenumber*cos(Azimuth)*-Wavenumber_y/(Wavenumber_x**2 + Wavenumber_y**2) )
    dWy_alpha = 2*pi*cos(EllipticityAngle)*( Wavenumber_y/Wavenumber*sin(Azimuth)+ Wavenumber*cos(Azimuth)*Wavenumber_x/(Wavenumber_x**2 + Wavenumber_y**2) )
    dEl_alpha = 2*pi*Wavenumber*sin(Azimuth)*-sin(EllipticityAngle)
    dEl_phi = 0
    for ll_x in range(0, L_x):
        ll = ndx_x[ll_x]
        
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dEl_Uwm[0] += dEl_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dEl_phi
        dEl_Uwm[1] += dEl_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dEl_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]
        dEl_Uw += 2*alpha*dEl_alpha*Sw_bw[0,0,ll]

    # ROTY
    ndx_y = arange(0, L)[comp == ROTY]
    L_y = len(ndx_y) # number of NSY components
    alpha = -2*pi*Wavenumber*cos(Azimuth)*cos(EllipticityAngle)
    dWx_alpha = -2*pi*cos(EllipticityAngle)*( Wavenumber_x/Wavenumber*cos(Azimuth)+ Wavenumber*-sin(Azimuth)*-Wavenumber_y/(Wavenumber_x**2 + Wavenumber_y**2) )
    dWy_alpha = -2*pi*cos(EllipticityAngle)*( Wavenumber_y/Wavenumber*cos(Azimuth)+ Wavenumber*-sin(Azimuth)*Wavenumber_x/(Wavenumber_x**2 + Wavenumber_y**2) )
    dEl_alpha = -2*pi*Wavenumber*cos(Azimuth)*-sin(EllipticityAngle)
    dEl_phi = 0
    for ll_y in range(0, L_y):
        ll = ndx_y[ll_y]
        phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        dWx_phi = -2*pi*pos_x[ll]
        dWy_phi = -2*pi*pos_y[ll]
        
        Uw += alpha**2*Sw_bw[0,0,ll]
        Uwm[0] += alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )
        Uwm[1] += alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )
        
        dWx_Uwm[0] += dWx_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWx_phi
        dWx_Uwm[1] += dWx_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWx_phi
        dWy_Uwm[0] += dWy_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dWy_phi
        dWy_Uwm[1] += dWy_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dWy_phi
        dEl_Uwm[0] += dEl_alpha * ( + cosPhi * Swm_bw[0,ll] + sinPhi * Swm_bw[1,ll] )+ alpha * ( -sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] )*dEl_phi
        dEl_Uwm[1] += dEl_alpha * ( - sinPhi * Swm_bw[0,ll] + cosPhi * Swm_bw[1,ll] ) + alpha * ( -cosPhi * Swm_bw[0,ll] - sinPhi * Swm_bw[1,ll] )*dEl_phi
        dWx_Uw += 2*alpha*dWx_alpha*Sw_bw[0,0,ll]
        dWy_Uw += 2*alpha*dWy_alpha*Sw_bw[0,0,ll]
        dEl_Uw += 2*alpha*dEl_alpha*Sw_bw[0,0,ll]
    
    Um = zeros((2,))
    Um[0] = Uwm[0] / Uw
    Um[1] = Uwm[1] / Uw
    
    
    grad[0] =  -0.5*( -(Um[0] * dWx_Uw * Um[0] + Um[1] *dWx_Uw * Um[1]) + dWx_Uwm[0] * Um[0] + dWx_Uwm[1] * Um[1])  - 0.5*(Um[0]*dWx_Uwm[0]  + Um[1]*dWx_Uwm[1])
    grad[1] =  -0.5*( -(Um[0] * dWy_Uw * Um[0] + Um[1] *dWy_Uw * Um[1]) + dWy_Uwm[0] * Um[0] + dWy_Uwm[1] * Um[1])  - 0.5*(Um[0]*dWy_Uwm[0]  + Um[1]*dWy_Uwm[1])
    grad[2] =  -0.5*( -(Um[0] * dEl_Uw * Um[0] + Um[1] *dEl_Uw * Um[1]) + dEl_Uwm[0] * Um[0] + dEl_Uwm[1] * Um[1])  - 0.5*(Um[0]*dEl_Uwm[0]  + Um[1]*dEl_Uwm[1])
        
    return(grad)
        
def fwMessages_RayleighWave(X_ML, Sw_bw, Swm_bw, SlnGamma_bw, ArrayInfo):
    
    Wavenumber_x = X_ML[0]
    Wavenumber_y = X_ML[1]
    EllipticityAngle = X_ML[2]
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x)
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2) 
    
    L=shape(ArrayInfo)[0]
    pos_x=ArrayInfo[:,0]
    pos_y=ArrayInfo[:,1]
    comp=ArrayInfo[:,3]
    
    Sm_fw = zeros((2,L))
    
    Uw = matrix(zeros((2,2)))
    Uwm = matrix(zeros((2,1)))
    for ll in range(0,L):    # compute the LL
        if comp[ll] == EWX:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = sin(EllipticityAngle)*cos(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == NSY:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = sin(EllipticityAngle)*sin(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == UDZ:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll]) +pi/2
            alpha = cos(EllipticityAngle)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == ROTX:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha =2*pi*Wavenumber*sin(Azimuth)*cos(EllipticityAngle)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == ROTY:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = -2*pi*Wavenumber*cos(Azimuth)*cos(EllipticityAngle)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        
        Uw = Uw + transpose(H) * matrix(Sw_bw[:,:,ll]) * H
        
        Uwm = Uwm + transpose(H) * transpose(matrix(Swm_bw[:,ll]))
    Um = linalg.solve(Uw, Uwm)
    LL = 0.5*transpose(Um) * Uwm + sum(SlnGamma_bw)
    LL = LL[0,0]
    
    for ll in range(0,L):    # compute the forward messages Sm_fw
        if comp[ll] == EWX:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = sin(EllipticityAngle)*cos(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == NSY:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = sin(EllipticityAngle)*sin(Azimuth)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == UDZ:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll]) +pi/2
            alpha = cos(EllipticityAngle)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == ROTX:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha =2*pi*Wavenumber*sin(Azimuth)*cos(EllipticityAngle)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        elif  comp[ll] == ROTY:
            phi = -2*pi*(Wavenumber_x*pos_x[ll] + Wavenumber_y*pos_y[ll])
            alpha = -2*pi*Wavenumber*cos(Azimuth)*cos(EllipticityAngle)
            H = alpha * matrix([[cos(phi), -sin(phi)],[sin(phi), cos(phi)]])
        else:
            H = matrix([[0, 0],[0, 0]])
        tmp = H * Um
        Sm_fw[0,ll] = tmp[0]
        Sm_fw[1,ll] = tmp[1]
    
    Amplitude = sqrt(Um[0][0,0]**2 + Um[1][0,0]**2)
    Phase = arctan2(Um[0][0,0], Um[1][0,0])

    return (Sm_fw, Amplitude, Phase, LL)

 
def fitRayleighWave(X_grid, Sm_bw, Sw_bw, ff, L, K, Kmax, ArrayInfo, Gamma):
    
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
            (LL) = negLL_RayleighWave(X_grid, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
            
            ndx_min = argmin(LL)
            X_g = X_grid[:,ndx_min];
            

            ### test gradient
            #print("Rayleigh wave likelihood: testing gradient")
            #print(X_grid[:,0])
            #grad1 = approx_fprime(X_grid[:,0], negLL_RayleighWave, sqrt(np.finfo(float).eps), Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
            #print(grad1)
            #grad2 = grad_negLL_RayleighWave(X_grid[:,0], Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
            #print(grad2)

            #
            #if shape(X_grid)[1] > 3:
            #    from mpl_toolkits.mplot3d import Axes3D # this line is needed for projection='3d'
            #    from matplotlib import cm
            #    import matplotlib.pyplot as plt
            #    plt.ion()
            #    fig = plt.figure()
            #    ax = fig.gca(projection='3d')
            #    ax.view_init(90, 270)
            #    ax.plot_trisurf(X_grid[0], X_grid[1], -LL, cmap=cm.jet, linewidth=0.0)
            #    plt.title("Love wave log likelihood - " + str(ff))
            #    plt.show()
            #    plt.pause(5)


#            if shape(X_grid)[1] > 3:
#                plt.ion()
#                fig = plt.figure()
#                ax = fig.gca(projection='3d')
#                ax.plot_trisurf(X_grid[0], X_grid[1], -LL, cmap=cm.jet, linewidth=0.2)
#                plt.title("log likelihood " + str(ff) + "Hz")
#                plt.show()

        # Refine ML estimate
        res = minimize(negLL_RayleighWave, X_g, (Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo),'TNC', bounds=[(-Kmax, Kmax),(-Kmax, Kmax), (-pi, pi)], jac=grad_negLL_RayleighWave)
        # numerical optimization with 'L-BFGS-B', showed some issues. Estimates were clustered around grid search points.
        X_ML = res.x; LL_ML = -res.fun[0];
        
        #Kx = X_ML[0]
        #Ky = X_ML[1]
        #Ellipticity = X_ML[2]
        #Wavenumber = sqrt( Kx**2 + Ky**2)
        #Azimuth = np.mod(arctan2(Ky, Kx), 2*pi)
        #print("Rayleigh wave ML estimate: {:.2e} {:.2e} {:.2e}".format(Wavenumber, Azimuth, Ellipticity))

        
        # Compute fw messages
        (Sm_fw_ff, Amplitude, Phase, LL_ML) = fwMessages_RayleighWave(X_ML, Snw_bw, Snwm_bw, SlnGamma_bw, ArrayInfo)
        Sm_fw[:,:,ff] = Sm_fw_ff
        
        hist_LL[ii] = LL_ML
        if (ii > 0) and (abs(hist_LL[ii]-hist_LL[ii-1])/abs(hist_LL[ii-1]) < 0.01):
            break
        
        
    Wavenumber_x = X_ML[0]
    Wavenumber_y = X_ML[1]
    EllipticityAngle = X_ML[2]
    Wavenumber = sqrt(Wavenumber_x**2 + Wavenumber_y**2)
    Azimuth = arctan2(Wavenumber_y, Wavenumber_x)
    
    LogLikelihood = LL_ML
    NumParameters = 3 * L + 5
    NumPoints = K * L
    BIC = -2* LogLikelihood + Gamma *NumParameters *log(NumPoints)

    
    
    logging.debug('Rayleigh wave fit')
    logging.debug('\tLL: {0:.3e} BIC: {1:.3e}'.format(LogLikelihood, BIC))
    logging.debug('\tAmplitude: {0:.3f}, Phase: {1:.3f}, Wavenumber: {2:.3f}, Azimuth: {3:.3f}, EllipticityAngle: {4:.3f}'.format(Amplitude,Phase,Wavenumber,Azimuth,EllipticityAngle))

    return(BIC, Sm_fw, Amplitude, Phase, X_ML)
    

### Noise model       
        
def estimateNoise(Sm_bw, Sm_fw):
    """ Estimate noise parameters
 
    Parameters
    ----------
    Sm_bw : float array
        Sufficient statistics from measurements
    Sm_fw : float array
        Wavefield explained by modeled waves

    Returns
    -------
    sigma2 : float array
        Noise variance estimated from residual signal. One dimensional array of length L.
    SlnGamma_bw : float array
	Array of gamma scale factors computed with the estimated noise variance
"""
    
    L=shape(Sm_bw)[1]
    K=int(2*(shape(Sm_bw)[2] -1))
    SlnGamma_bw=zeros((L,))
    sigma2=zeros((L,))
    for ll in range(0,L):
        Y_bw=zeros((K,))+1j*zeros((K,))
        Y_fw=zeros((K,))+1j*zeros((K,))

        Kd2p1 = int(K/2+1)
        Kd2   = int(K/2)
        Y_bw[0:Kd2p1] = Sm_bw[0,ll,:] + 1j*Sm_bw[1,ll,:]
        Y_bw[Kd2p1:K] = conjugate(Y_bw[1:Kd2])
        Y_fw[0:Kd2p1] = Sm_fw[0,ll,:] + 1j*Sm_fw[1,ll,:]
        Y_fw[Kd2p1:K] = conjugate(Y_fw[1:Kd2])
        
        sigma2[ll] = sum( real(conjugate(Y_bw)*Y_bw) -2*real(conjugate(Y_bw)*Y_fw) + real(conjugate(Y_fw)*Y_fw)) /4
        # sigma2_b[ll] = sum((y_fw[:,ll] - y_bw[:,ll])**2)/K # time domain
        
        SlnGamma_bw[ll] = - (K/2)* log(2*pi*sigma2[ll]) -K*sum(real(conjugate(Y_bw)*Y_bw))/sigma2[ll]/8
    
    #if EQUAL_SIGMA2:
    #    sigma2 = np.mean(sigma2)*np.ones((L,))
    #    SlnGamma_bw = - (K/2)* log(2*pi*sigma2[0])
    
    return (sigma2, SlnGamma_bw)
    

def fitNoise(Sm_bw, L, K, Gamma):
    """ Estimate noise parameters
 
    Parameters
    ----------
    Sm_bw : float array
        Sufficient statistics from measurements
    Sm_fw : float array
        Wavefield explained by modeled waves

    Returns
    -------
    BIC : float
        Value of the Bayesian information criterion
    sigma2_ML : float array
        Noise variance estimated from residual signal. One dimensional array of length L.
"""

    
    Sm_fw = zeros(shape(Sm_bw))
    
    (sigma2, SlnGamma_bw) = estimateNoise(Sm_bw, Sm_fw)
    
    LogLikelihood = sum(SlnGamma_bw)
    NumParameters = 3 * L
    NumPoints = K * L
    BIC = -2* LogLikelihood + Gamma *NumParameters *log(NumPoints)
    sigma2_ML = sigma2
    
    logging.debug('Additive Gaussian noise fit')
    logging.debug('\tLL: {0:.3e} BIC: {1:.3e}'.format(LogLikelihood, BIC))
    logging.debug('\tsigma2: {0:.2e} ... {1:.2e} ...  {2:.2e}'.format(min(sigma2), mean(sigma2), max(sigma2)))

    return(BIC, sigma2_ML)
    
