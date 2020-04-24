# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
from numpy.matlib import repmat
import CircStat as cs
from wdSettings import DEFAULT_SCALE_WAVENUMBER, DEFAULT_SCALE_ELLIPTICITY, DEFAULT_SCALE_AZIMUTH, DEFAULT_SCALE_VELOCITY
from wdSettings import DEFAULT_USETEX, DEFAULT_PLOT_ESTIMATE_DOT
from wdSettings import DEFAULT_PLOT_FREQUENCY_LOG, DEFAULT_PLOT_VELOCITY_LOG




def density_estimation(X_vec, Y_vec, X_data, Y_data, scale=1, modularity=None):
    
    assert len(X_data) == len(Y_data)    
    
    NumX = len(X_vec) # eg, frequency
    NumY = len(Y_vec) # eg, wavenumber, ellipticity
    density = np.zeros((NumY, NumX))
    
    for xx in range(0,NumX):
        ndx_x = np.where(X_data == X_vec[xx])[0]
        Y_x = Y_data[ndx_x]
        #ndx_Xnz = np.where(Y_x > 0) # remove zero elements
        #Y_x = Y_x[ndx_Xnz]
        
        if len(Y_x) > 1 and sum(abs(Y_x)) > 0: # TODO. gaussian_kde needs more points than dimension
            kernel = cs.CircKde(Y_x, bw_method=scale, modularity=modularity)
            density[:,xx] = kernel(Y_vec)
            if np.max(density[:,xx]) > 0:
                density[:,xx] /= np.max(density[:,xx]) # for visualization, better /max() than /sum()
    
    return (density)
    
    
    
def plotWavenumber(F, K, WaveLabel=None, ax=None, NumK=300, scale=DEFAULT_SCALE_WAVENUMBER, xlim=None, ylim=None):
    
    assert len(F) == len(K)
    
    Fvec = np.unique(F) # WARNING this vector is not necessarily equally spaced! (eg, log-spacing, approximation)
    NumF = len(Fvec)
    
    Kvec = np.linspace(0,np.max(K)*1.05,NumK)
    
    density_K = np.zeros((NumK,NumF+1)) # the +1 is intended to trick pcolor() which ignores the last column and row
    density_K[:,0:-1] = density_estimation(Fvec, Kvec, F, K, scale)
    
    if NumF > 1:
        Fvec_pcolor = 0.5*( np.concatenate((np.array([2*Fvec[0] -Fvec[1]]), Fvec )) + np.concatenate((Fvec, np.array([2*Fvec[-1] -Fvec[-2]]))) )            
        Flim = np.array([Fvec[0], Fvec[-1]]) + np.array([-0.5, +0.5])*(Fvec[1]-Fvec[0])
    else:
        Fvec_pcolor = np.array([Fvec[0], Fvec[0]+0.5])
        Flim = [Fvec_pcolor.min(), Fvec_pcolor.max()]
    
    # mesh grid for pcolor()
    K_mx, K_my = np.meshgrid(Fvec_pcolor,Kvec)
    
    # Plotting
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    else:
        fig = ax.get_figure()
    
    ax.set_xlim([Fvec_pcolor.min(), Fvec_pcolor.max()])    
    ax.pcolor(K_mx, K_my, density_K, cmap='hot_r') # TODO this line prints an error to the console
    if DEFAULT_PLOT_ESTIMATE_DOT: ax.plot(F, K,'r.', markersize=3)
    
    if DEFAULT_PLOT_FREQUENCY_LOG: 
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xticks(np.kron(np.power(10, np.arange(np.floor(np.log10(Flim[0])), np.ceil(np.log10(Flim[1])))), np.array([1, 2, 5])))
    ax.set_xlim(Flim) if xlim is None else ax.set_xlim(xlim)
    Kedge = (Kvec.max() - Kvec.min())*0.05
    ax.set_ylim([Kvec.min()-Kedge, Kvec.max()+Kedge]) if ylim is None else ax.set_ylim(ylim)    
    ax.set_ylabel('Wavenumber [1/m]')
    ax.set_xlabel('Frequency [Hz]')
    
    if WaveLabel is not None: ax.set_title('{0} wave dispersion'.format(WaveLabel))
    
    return ax
    
def plotVelocity(F, V, WaveLabel=None, ax=None, NumV=150, scale=DEFAULT_SCALE_VELOCITY, xlim=None, ylim=None):
    
    Fvec = np.unique(F) # WARNING this vector is not necessarily equally spaced! (eg, log-spacing, approximation)
    NumF = len(Fvec)
    Vmin = np.max(np.array([np.min(V), 50]))
    Vmax = np.min(np.array([np.max(V), 5000]))
    Vvec = np.logspace(np.log10(Vmin),np.log10(Vmax),NumV)
    
    density_V = np.zeros((NumV,NumF+1)) # the +1 is intended to trick pcolor() which ignores the last column and row
    density_V[:,0:-1] = density_estimation(Fvec, Vvec, F, V, scale)
    
    if NumF > 1:
        Fvec_pcolor = 0.5*( np.concatenate((np.array([2*Fvec[0] -Fvec[1]]), Fvec )) + np.concatenate((Fvec, np.array([2*Fvec[-1] -Fvec[-2]]))) )            
        Flim = np.array([Fvec[0], Fvec[-1]]) + np.array([-0.5, +0.5])*(Fvec[1]-Fvec[0])
    else:
        Fvec_pcolor = np.array([Fvec[0], Fvec[0]+0.5])
        Flim = [Fvec_pcolor.min(), Fvec_pcolor.max()]
    
    # mesh grid for pcolor()
    V_mx, V_my = np.meshgrid(Fvec_pcolor, Vvec)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    else:
        fig = ax.get_figure()
    
       
    ax.set_xlim([F.min(), F.max()])
    
    ax.pcolor(V_mx, V_my, density_V, cmap='hot_r') # TODO this line prints an error to the console
    if DEFAULT_PLOT_ESTIMATE_DOT: ax.plot(F, V, 'r.', markersize=3)
    if DEFAULT_PLOT_FREQUENCY_LOG:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xticks(np.kron(np.power(10, np.arange(np.floor(np.log10(Flim[0])), np.ceil(np.log10(Flim[1])))), np.array([1, 2, 5])))
    if DEFAULT_PLOT_VELOCITY_LOG:
        ax.set_yscale("log")        
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.set_yticks(np.kron(np.power(10, np.arange(np.floor(np.log10(Flim[0])), np.ceil(np.log10(Flim[1])))), np.array([1, 2, 5])))
    ax.set_xlim(Flim) if xlim is None else ax.set_xlim(xlim)
    ax.set_ylim([Vmin, Vmax]) if ylim == None else ax.set_ylim(ylim)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Velocity [m/s]')
    
    if WaveLabel is not None: ax.set_title('{0} wave dispersion'.format(WaveLabel))
    plt.show()
    
    return ax
    

def plotEllipticity(F, E, WaveLabel=None, ax=None, NumE=300, scale=DEFAULT_SCALE_ELLIPTICITY, xlim=None, ylim=None):
    """
        generates and plots histogram for ellipticity angle
    """
    
    Fvec = np.unique(F) # WARNING this vector is not necessarily equally spaced! (eg, log-spacing, approximation)
    NumF = len(Fvec)
    
    Evec = np.linspace(-np.pi/2,np.pi/2,NumE)
    
    
    density_E = np.zeros((NumE,NumF+1))  # the +1 is intended to trick pcolor() which ignores the last column and row
        
    density_E[:,0:-1] = density_estimation(Fvec, Evec, F, E, scale, np.pi)
    
    if NumF > 1:
        Fvec_pcolor = 0.5*( np.concatenate((np.array([2*Fvec[0] -Fvec[1]]), Fvec )) + np.concatenate((Fvec, np.array([2*Fvec[-1] -Fvec[-2]]))) )            
        Flim = np.array([Fvec[0], Fvec[-1]]) + np.array([-0.5, +0.5])*(Fvec[1]-Fvec[0])
    else:
        Fvec_pcolor = np.array([Fvec[0], Fvec[0]+0.5])
        Flim = [Fvec_pcolor.min(), Fvec_pcolor.max()]
    
    
    # mesh grid for pcolor()
    E_mx, E_my = np.meshgrid(Fvec_pcolor,Evec)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    else:
        fig = ax.get_figure()
    
    # Plotting            
    ax.pcolor(E_mx, E_my, density_E, cmap='hot_r')
    ax.plot([0, 100], [0,0], 'k:', linewidth=1)
    if DEFAULT_PLOT_ESTIMATE_DOT: ax.plot(F, E, 'r.', markersize=3)
    if DEFAULT_PLOT_FREQUENCY_LOG:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xticks(np.kron(np.power(10, np.arange(np.floor(np.log10(Flim[0])), np.ceil(np.log10(Flim[1])))), np.array([1, 2, 5])))
    ax.set_yticks([-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2])
    if DEFAULT_USETEX:
        ax.set_yticklabels(['$\\texttt{-}\\frac{\pi}{2}$', '$\\texttt{-}\\frac{\pi}{4}$', '$0$', '$\\texttt{+}\\frac{\pi}{4}$', '$\\texttt{+}\\frac{\pi}{2}$'])
    else:
        ax.set_yticklabels(['-pi/2', '-pi/4', '0', '+pi/4', '+pi/2'])
    ax.set_xlim(Flim) if xlim is None else ax.set_xlim(xlim)
    ax.set_ylim([-np.pi/2, np.pi/2]) if ylim is None else ax.set_ylim(ylim)   
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Ellipticity Angle [rad]')
    if WaveLabel is not None: ax.set_title('{0} ellipticity curve'.format(WaveLabel))
            
            
    return fig
    
def plotAzimuth(F, A, WaveLabel=None, ax=None, NumA=300, scale=DEFAULT_SCALE_AZIMUTH, xlim=None, ylim=None):
    
    Fvec = np.unique(F) # WARNING this vector is not necessarily equally spaced! (eg, log-spacing, approximation)
    NumF = len(Fvec)
    Avec = np.linspace(0,2*np.pi,NumA)
    
    density_A = np.zeros((NumA,NumF+1))  # the +1 is intended to trick pcolor() which ignores the last column and row
    density_A[:,0:-1] = density_estimation(Fvec, Avec, F, A, scale, 2*np.pi)    
        
    # mesh grid for pcolor()
    if NumF > 1:
        Fvec_pcolor = 0.5*( np.concatenate((np.array([2*Fvec[0] -Fvec[1]]), Fvec )) + np.concatenate((Fvec, np.array([2*Fvec[-1] -Fvec[-2]]))) )            
        Flim = np.array([Fvec[0], Fvec[-1]]) + np.array([-0.5, +0.5])*(Fvec[1]-Fvec[0])
    else:
        Fvec_pcolor = np.array([Fvec[0], Fvec[0]+0.5])
        Flim = [Fvec_pcolor.min(), Fvec_pcolor.max()]
    A_mx, A_my = np.meshgrid(Fvec_pcolor,Avec)
    
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    else:
        fig = ax.get_figure()
    
    ax_right = ax.twinx() 
    # Plotting            
    ax.pcolor(A_mx, A_my, density_A, cmap='hot_r')
    if DEFAULT_PLOT_ESTIMATE_DOT: ax.plot(F, A, 'r.', markersize=3)
    if DEFAULT_PLOT_FREQUENCY_LOG:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xticks(np.kron(np.power(10, np.arange(np.floor(np.log10(Flim[0])), np.ceil(np.log10(Flim[1])))), np.array([1, 2, 5])))
    ax.set_xlim(Flim) if xlim is None else ax.set_xlim(xlim)
    ax.set_ylim([0, 2*np.pi]) if ylim == None else ax.set_ylim(ylim)    
    ax.set_yticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax_right.set_ylim([0, 2*np.pi]) if ylim == None else ax.set_ylim(ylim)    
    ax_right.set_yticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    if DEFAULT_USETEX:
        ax.set_yticklabels(['$0$', '$\\frac{\pi}{2}$', '$\pi$', '$\\frac{3\pi}{2}$', '$2 \pi$'])
    else:
        ax.set_yticklabels(['0', 'pi/2', 'pi', '3pi/2', '2pi'])
        
    ax_right.set_yticklabels(['W', 'S', 'E', 'N', 'W'])
    ax_right.set_ylabel('Direction of arrival')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Azimuth [rad]')
            
    return fig



def plotEllipticityAngle(F, E, E_std=None, color='r', marker='.', linestyle='solid', linewidth=1, alpha=1, label=None):
    """
    This function to plots nicely the ellipticity angle, avoiding unecessary jumps
    """
    assert len(F) == len(E)
    E = np.mod(E+np.pi/2, np.pi) - np.pi/2
    
    for ff in range(0, len(F)-1):
        if np.abs(E[ff] - E[ff+1]) > np.pi/2: # if there is a large jump, plot differently
            if E[ff] > 0:
                plt.plot(F[ff:ff+2], [E[ff], E[ff+1]+np.pi], color=color, marker=marker, linestyle=linestyle, linewidth=linewidth, alpha=alpha, label=label)
                plt.plot(F[ff:ff+2], [E[ff]-np.pi, E[ff+1]], color=color, marker=marker, linestyle=linestyle, linewidth=linewidth, alpha=alpha, label=label)
                if E_std is not None:
                    plt.errorbar(F[ff:ff+2], [E[ff], E[ff+1]+np.pi], yerr=E_std[ff:ff+2], fmt='none', color=color, marker=marker, ecolor=color, zorder=10, label=label)
                    plt.errorbar(F[ff:ff+2], [E[ff]-np.pi, E[ff+1]], yerr=E_std[ff:ff+2], fmt='none', color=color, marker=marker, ecolor=color, zorder=10, label=label)
            elif E[ff] < 0:
                plt.plot(F[ff:ff+2], [E[ff], E[ff+1]-np.pi], color=color, marker=marker, linestyle=linestyle, linewidth=linewidth, alpha=alpha, label=label)
                plt.plot(F[ff:ff+2], [E[ff]+np.pi, E[ff+1]], color=color, marker=marker, linestyle=linestyle, linewidth=linewidth, alpha=alpha, label=label)
                if E_std is not None:
                    plt.errorbar(F[ff:ff+2], [E[ff], E[ff+1]-np.pi], yerr=E_std[ff:ff+2], fmt='none', color=color, marker=marker, ecolor=color, zorder=10, label=label)
                    plt.errorbar(F[ff:ff+2], [E[ff]+np.pi, E[ff+1]], yerr=E_std[ff:ff+2], fmt='none', color=color, marker=marker, ecolor=color, zorder=10, label=label)
            else:
                plt.plot(F[ff:ff+2], E[ff:ff+2], color=color, marker=marker, linewidth=linewidth,  linestyle=linestyle, alpha=alpha, label=label)
                if E_std is not None:
                    plt.errorbar(F[ff:ff+2], E[ff:ff+2], yerr=E_std[ff:ff+2], fmt='none', marker=marker, color=color, ecolor=color, zorder=10, label=label)
        else:
            plt.plot(F[ff:ff+2], E[ff:ff+2], color=color, linewidth=linewidth, marker=marker, linestyle=linestyle, alpha=alpha, label=label)
            if E_std is not None:
                plt.errorbar(F[ff:ff+2], E[ff:ff+2], yerr=E_std[ff:ff+2], fmt='none', marker=marker, color=color, ecolor=color, zorder=10, label=label)
    return
    
    
    

def plotBounds(ub_x, ub_y, lb_x, lb_y, ax=None):
    """
    This function to plots the upper and lower bounds used in picking
    """
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    ax.plot(ub_x, ub_y, 'bv--', linewidth=3, markersize=10)
    ax.plot(lb_x, lb_y, 'b^--', linewidth=3, markersize=10)

def plotArray(array_info, SourcePosition=None):
    """
    Plots array layout
    """
    
    pos_x = array_info[:,0]
    pos_y = array_info[:,1]
    
    xMin = np.min(pos_x)
    xMax = np.max(pos_x)
    yMin = np.min(pos_y)
    yMax = np.max(pos_y)
    if SourcePosition is not None:
        xMin = np.min([xMin, SourcePosition[0]])
        xMax = np.max([xMax, SourcePosition[0]])
        yMin = np.min([yMin, SourcePosition[1]])
        yMax = np.max([yMax, SourcePosition[1]])
        
    xAperture = np.max([xMax - xMin, 0.8*(yMax - yMin)]) # avoid really stretched aspect ratios
    yAperture = np.max([yMax - yMin, 0.8*(xMax - xMin)])
   
    fig = plt.figure()
    ax = fig.gca()
    ax.set_aspect('equal')
    ax.grid()
    ax.plot(pos_x, pos_y, 'o', ms=10, lw=2, alpha=0.7, mfc='orange')
    if SourcePosition is not None:
        ax.plot(SourcePosition[0], SourcePosition[1], '*', ms=18, lw=2, alpha=0.8, mfc='red')
    ax.set_xlim(xMin-0.05*xAperture, xMax+0.05*xAperture)
    ax.set_ylim(yMin-0.05*yAperture, yMax+0.05*yAperture)
    ax.set_xlabel('Easting [m]')
    ax.set_ylabel('Northing [m]')
    plt.show()
    
    return fig
    
    
def plotArrayResponse(array_info, Kmax, Knum=100, normalized=True):
    """
        plot the squared modulo of the array response H^2
        
        pos the sensor positions
        Kmax the largest wavenumber of interest [1/m]
        Knum grid size
        normalized if True, sets the largest value to one
    """
    
    Nsensors = np.shape(array_info)[0]
    pos_x = array_info[:,0]
    pos_y = array_info[:,1]
    K_vec = np.linspace(-Kmax, Kmax, num=Knum)
    K_plane = repmat(K_vec, Knum, 1)
    H = np.zeros((Knum,Knum))  # the squared modulos is real
    
    # loop solely on the elements above diagonal
    for n1 in range(0,Nsensors):
        for n2 in range(n1+1, Nsensors):
            H += np.cos(2*np.pi*(K_plane*(pos_x[n1]-pos_x[n2]) + np.transpose(K_plane*(pos_y[n1]-pos_y[n2]))))
    
    H *= 2        # elements below diagonal
    H += Nsensors # elements on the diagonal (all equal to 1)

    if normalized:
        H = H / Nsensors**2
        
        
       
    fig = plt.figure()
    plt.pcolor(K_vec, K_vec, H, cmap='hot_r')
    cbar = plt.colorbar()

    if normalized:
        plt.clim(0,1)
        cbar.ax.set_ylabel('Normalized array response')
    else:
        plt.clim(0,Nsensors**2)
        cbar.ax.set_ylabel("Array response")
    plt.xlim(-Kmax, Kmax)
    plt.ylim(-Kmax, Kmax)
    
    plt.axes().set_aspect('equal')
    plt.xlabel('Wavenumber x [1/m]')
    plt.ylabel('Wavenumber y [1/m]')
    plt.show()
    return fig
    
def plotArrayResponseCuts(array_info, Kmax, Knum=200, Thetanum=60, half=False, normalized=True):
    """
        plot sections the squared modulo of the array response H^2
        
        pos the sensor positions
        Kmax the largest wavenumber of interest [1/m]
        Knum grid size
        Thetanum number of sections to be displayed
        normalized if True, sets the largest value to one
    """
    Nsensors = np.shape(array_info)[0]
    pos_x = array_info[:,0]
    pos_y = array_info[:,1]
    if half:
        K_vec = np.linspace(0, Kmax, num=Knum)
    else:
        K_vec = np.linspace(-Kmax, Kmax, num=Knum)
    Theta_vec = np.linspace(0, np.pi, num=Thetanum, endpoint=False)
    H = np.zeros((Knum,Thetanum))  # the squared modulos is real
    
    for tt in range(0,Thetanum):
        Kx_vec = K_vec * np.cos(Theta_vec[tt])
        Ky_vec = K_vec * np.sin(Theta_vec[tt])
        
        # loop solely on the elements above diagonal
        for n1 in range(0,Nsensors):
            for n2 in range(n1+1, Nsensors):
                H[:,tt] += np.cos(2*np.pi*(Kx_vec*(pos_x[n1]-pos_x[n2]) + np.transpose(Ky_vec*(pos_y[n1]-pos_y[n2]))))
    
    H *= 2        # elements below diagonal
    H += Nsensors # elements on the diagonal (all equal to 1)

    if normalized:
        H = H / Nsensors**2
        
       
    fig = plt.figure()
    plt.plot(K_vec, H, color='red', linewidth=1.0, alpha=0.5)
    if normalized:
        plt.ylim(0,1)
        plt.ylabel('Normalized array response')
    else:
        plt.ylim(0,Nsensors**2)
        plt.ylabel('Array response')
    plt.xlim(np.min(K_vec), np.max(K_vec))
    
    plt.xlabel('Wavenumber [1/m]')
    plt.show()
    return fig
