# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
This class provides functionalities to pick dispersion and ellipticty curves
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import CircStat as cs
from PlotUtils import plotVelocity, density_estimation
from numpy import pi, array

from wdSettings import DEFAULT_SCALE_WAVENUMBER, DEFAULT_SCALE_VELOCITY, DEFAULT_SCALE_ELLIPTICITY, DEFAULT_PLOT_ESTIMATE_DOT, DEFAULT_PLOT_FREQUENCY_LOG


class Picker:
    """
    This class provide the graphical interface for selecting points from 
    dispersion curve and ellipticity curve.
    """
    def __init__(self, x, y, mode, arrayResolutionLimits=None):
        """
        Initialization
        """
        
        self.dispersion = True if mode.lower() == 'dispersion' else False
        self.ellipticity = True if mode.lower() == 'ellipticity' else False
        assert self.dispersion ^ self.ellipticity # dispersion XOR ellipticity
        
        self.F = x
        self.Fvec = np.unique(self.F) # WARNING this vector is not necessarily equally spaced! (eg, log-spacing, approximation)
        self.NumF = len(self.Fvec)
        self.Flim = np.array([np.min(self.Fvec), np.max(self.Fvec)])

        self.scale_wavenumber = DEFAULT_SCALE_WAVENUMBER
        self.scale_velocity = DEFAULT_SCALE_VELOCITY
        self.scale_ellipticity = DEFAULT_SCALE_ELLIPTICITY
        if self.dispersion:            
            self.K = y
            self.V = x/y
            self.NumK=200
            self.Kvec = np.linspace(0,np.max(self.K),self.NumK)
        if self.ellipticity:
            self.E = y
            self.F3 = np.concatenate((self.F,self.F,self.F))
            self.E3 = np.concatenate((self.E,self.E-pi,self.E+pi))
            self.NumE=200
            self.Evec = np.linspace(-pi/2, pi/2, self.NumE)
            self.E3vec = np.linspace(-pi/2-pi, pi/2+pi, 3*self.NumE)
            
            
            
        if arrayResolutionLimits is None:
            self.plotArrayResolution = False
        else:
            self.plotArrayResolution = True
            (self.Fres, self.Kmin, self.Kmax, self.Vmin, self.Vmax) = arrayResolutionLimits
            
        self.lb_x = list() # these lists contain the x,y upper and lower bound for coordiantes.
        self.lb_y = list()
        self.ub_x = list()
        self.ub_y = list()

        self.pickingUpperBound = True
        self.pickingLowerBound = False
        
        
    def plotAndSelect(self):
        """
        Plot and get selection bounds from user
        """
        
        self.fig = plt.figure(figsize=(8,12))
        self.cid1 = self.fig.canvas.mpl_connect('button_press_event', self.onMouseClick)
        self.cid2 = self.fig.canvas.mpl_connect('key_press_event', self.onKeyPress)
        self.cid3 = self.fig.canvas.mpl_connect('scroll_event', self.onMouseScroll)
        
        
        if self.dispersion:
            self.ax1 = self.fig.add_subplot(2,1,1)
            
            density_K = np.zeros((self.NumK,self.NumF+1)) # the +1 is intended to trick pcolor() which ignores the last column and row
            density_K[:,0:-1] = density_estimation(self.Fvec, self.Kvec, self.F, self.K, self.scale_wavenumber)
            
            Fvec_pcolor = 0.5*( np.concatenate((array([2*self.Fvec[0] -self.Fvec[1]]), self.Fvec )) + np.concatenate((self.Fvec, array([2*self.Fvec[-1] -self.Fvec[-2]]))) )
            self.Flim = [Fvec_pcolor.min(), Fvec_pcolor.max()]
            
            # mesh grid for pcolor()
            K_mx, K_my = np.meshgrid(Fvec_pcolor,self.Kvec)
    
            self.ax1.set_xlim(self.Flim)    
            self.pColor = self.ax1.pcolor(K_mx, K_my, density_K, cmap='hot_r') # TODO this line prints an error to the console
            if DEFAULT_PLOT_FREQUENCY_LOG:
                self.ax1.set_xscale("log")        
                self.ax1.xaxis.set_major_formatter(ScalarFormatter())
                self.ax1.set_xticks(np.kron(np.power(10, np.arange(np.floor(np.log10(self.Flim[0])), np.ceil(np.log10(self.Flim[1])))), np.array([1, 2, 5])))
 
            if DEFAULT_PLOT_ESTIMATE_DOT: self.ax1.plot(self.F, self.K,'r.', markersize=3)
            self.ax1.set_xlim(self.Flim)
            Kedge = (self.Kvec.max() - self.Kvec.min())*0.05
            self.ax1.set_ylim([self.Kvec.min()-Kedge, self.Kvec.max()+Kedge])
            self.ax1.set_ylabel('Wavenumber [1/m]')
            self.ax1.set_xlabel('Frequency [Hz]')
            
            if self.plotArrayResolution: self.ax1.plot(self.Fres, self.Kmin,'y-',linewidth=2),self.ax1.plot(self.Fres, self.Kmax,'y-',linewidth=2)
                            
            self.ax2 = self.fig.add_subplot(2,1,2, sharex=self.ax1)
            plotVelocity(self.F, self.V, WaveLabel=None, ax=self.ax2, scale=self.scale_velocity)
            if self.plotArrayResolution: self.ax2.plot(self.Fres, self.Vmin,'y-',linewidth=2),self.ax2.plot(self.Fres, self.Vmax,'y-',linewidth=2)
                
            self.fig.axes[0].set_title("Please read instructions in console.")        
            # xlim = self.fig.axes[0].get_xlim()
            # ylim = self.fig.axes[0].get_ylim()
            
            a = self.fig.axes[0] # frequency-wavenumber plot
            self.ub, = a.plot(a.axis()[0], a.axis()[2], 'bv--', linewidth=3, markersize=10)
            self.lb, = a.plot(a.axis()[0], a.axis()[2], 'b^--', linewidth=3, markersize=10)

            a = self.fig.axes[1] # frequency-velocity plot
            self.ub_v_points, = a.plot(a.axis()[0], a.axis()[2], 'bv', linewidth=3, markersize=10)
            self.lb_v_points, = a.plot(a.axis()[0], a.axis()[2], 'b^', linewidth=3, markersize=10)
            self.ub_v_lines, = a.plot(a.axis()[0], a.axis()[2], 'b--', linewidth=3, markersize=10)
            self.lb_v_lines, = a.plot(a.axis()[0], a.axis()[2], 'b--', linewidth=3, markersize=10)
            
        elif self.ellipticity:
            self.ax1 = self.fig.add_subplot(1,1,1)

            density_E = np.zeros((3*self.NumE, self.NumF+1)) # the +1 is intended to trick pcolor() which ignores the last column and row
            density_E[0:self.NumE,0:-1] = density_estimation(self.Fvec, self.Evec, self.F, self.E, self.scale_ellipticity, np.pi)
            density_E[self.NumE:2*self.NumE,0:-1] = density_E[0:self.NumE,0:-1]
            density_E[2*self.NumE:,0:-1] = density_E[0:self.NumE,0:-1]
            
            Fvec_pcolor = 0.5*( np.concatenate((array([2*self.Fvec[0] -self.Fvec[1]]), self.Fvec )) + np.concatenate((self.Fvec, array([2*self.Fvec[-1] -self.Fvec[-2]]))) )
            self.Flim = [Fvec_pcolor.min(), Fvec_pcolor.max()]
            
            E_mx, E_my = np.meshgrid(Fvec_pcolor, self.E3vec) # mesh grid for pcolor()
    
            self.ax1.set_xlim(self.Flim)
            self.pColor = self.ax1.pcolor(E_mx, E_my, density_E, cmap='hot_r') # TODO this line prints an error to the console
            if DEFAULT_PLOT_FREQUENCY_LOG:
                self.ax1.set_xscale("log")        
                self.ax1.xaxis.set_major_formatter(ScalarFormatter())
                self.ax1.set_xticks(np.kron(np.power(10, np.arange(np.floor(np.log10(self.Flim[0])), np.ceil(np.log10(self.Flim[1])))), np.array([1, 2, 5])))
            if DEFAULT_PLOT_ESTIMATE_DOT: self.ax1.plot(self.F, self.E,'r.', markersize=3)
            self.ax1.plot(self.Flim, [0, 0],'b:', markersize=2)
            self.ax1.plot(self.Flim, [-pi/2, -pi/2],'b:', markersize=2)
            self.ax1.plot(self.Flim, [+pi/2, +pi/2],'b:', markersize=2)
            self.ax1.set_xlim(self.Flim)
            self.ax1.set_ylim([self.E3vec.min(), self.E3vec.max()])
            self.ax1.set_ylabel('Ellipticity angle [rad]')
            self.ax1.set_xlabel('Frequency [Hz]')
            self.ax1.set_yticks([-3*np.pi/2, -np.pi,-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2, np.pi, 3*np.pi/2])
            self.ax1.set_yticklabels(['pi/2','0','-pi/2', '-pi/4', '0', '+pi/4', '+pi/2', '0', 'pi/2'])
            
            self.fig.axes[0].set_title("Please read instructions in console.")        
            # xlim = self.fig.axes[0].get_xlim()
            # ylim = self.fig.axes[0].get_ylim()
            
            a = self.fig.axes[0]
            self.ub, = a.plot(a.axis()[0], a.axis()[2], 'bv--', linewidth=3, markersize=10)
            self.lb, = a.plot(a.axis()[0], a.axis()[2], 'b^--', linewidth=3, markersize=10)

        
        print("\tSelect UPPER bound using the mouse. When finished press the SPACE KEY.")
        return
            

    def redrawDensity(self):
        """
        This function redraws the density
        """
        if self.dispersion:
            
            # TODO testare se anche il wavenumber si puo plottare chiamando plotwavenumber come plotvelocity
            
            density_K = np.zeros((self.NumK,self.NumF+1)) # the +1 is intended to trick pcolor() which ignores the last column and row
            density_K[:,0:-1] = density_estimation(self.Fvec, self.Kvec, self.F, self.K, self.scale_wavenumber)
            self.pColor.set_array(density_K[0:-1,0:-1].ravel(order='C')) # pcolor() throws away last row and last column
            
            plotVelocity(self.F, self.V, WaveLabel=None, ax=self.ax2, scale=self.scale_velocity)
            
        elif self.ellipticity:
            density_E = np.zeros((3*self.NumE, self.NumF+1)) # the +1 is intended to trick pcolor() which ignores the last column and row
            density_E[0:self.NumE,0:-1] = density_estimation(self.Fvec, self.Evec, self.F, self.E, self.scale_velocity, np.pi)
            density_E[self.NumE:2*self.NumE,0:-1] = density_E[0:self.NumE,0:-1]
            density_E[2*self.NumE:,0:-1] = density_E[0:self.NumE,0:-1]
            
            self.pColor.set_array(density_E[0:-1,0:-1].ravel(order='C')) # pcolor() throws away last row and last column

        plt.draw()
        return

    def onMouseClick(self, event):
        # print("Mouse click")
        # print(event)
    
        # here are all the picking possibilities      
        if self.dispersion and event.inaxes==self.fig.axes[0]: # freq-wavenumber plane
            xdata = event.xdata
            ydata = np.max([event.ydata, 1e-9]) # make sure we cannot pick negative wavenumbers
        elif self.dispersion and event.inaxes==self.fig.axes[1]: # freq-velocity plane
            xdata = event.xdata # frequency
            ydata = event.xdata / event.ydata # convert velocity to wavenumber
        elif self.ellipticity: # ellipticity plane
            xdata = event.xdata
            ydata = event.ydata
        else:
            # do nothing
            return
        
        
        if self.pickingUpperBound and not self.pickingLowerBound: # picking upper bound
            if event.button == 1: # adding a point
                if len(self.ub_x) == 0 or self.ub_x[-1] < xdata: # accept picks only to the right of last point
                    self.ub_x.append(xdata)
                    self.ub_y.append(ydata)
            elif event.button == 3 and len(self.ub_x) > 0: # remove last point
                self.ub_x.pop()
                self.ub_y.pop()
                
            self.ub.set_data(self.ub_x, self.ub_y)
            if self.dispersion:
                self.ub_v_points.set_data(self.ub_x, np.array(self.ub_x)/np.array(self.ub_y))
                v_line_x=[]
                v_line_y=[]
                for ii in range(0,len(self.ub_x)-1):
                    v_line_x += list(np.linspace(self.ub_x[ii], self.ub_x[ii+1], 10))
                    v_line_y += list(np.array(self.ub_y[ii]) + (np.array(self.ub_y[ii+1])-np.array(self.ub_y[ii]))*np.linspace(0,1,10))
                self.ub_v_lines.set_data(v_line_x, np.array(v_line_x)/np.array(v_line_y))
  
        elif not self.pickingUpperBound and self.pickingLowerBound: # picking lower bound
            if event.button == 1:
                if len(self.lb_x) == 0 or self.lb_x[-1] < xdata: # accept picks only to the right of last point
                    self.lb_x.append(xdata)
                    self.lb_y.append(ydata)
            elif event.button == 3 and len(self.lb_x) > 0: # right click, remove data point
                self.lb_x.pop()
                self.lb_y.pop()
        
            self.lb.set_data(self.lb_x, self.lb_y)            
            if self.dispersion:
                self.lb_v_points.set_data(self.lb_x, np.array(self.lb_x)/np.array(self.lb_y))
                v_line_x=[]
                v_line_y=[]
                for ii in range(0,len(self.lb_x)-1):
                    v_line_x += list(np.linspace(self.lb_x[ii], self.lb_x[ii+1], 10))
                    v_line_y += list(np.array(self.lb_y[ii]) + (np.array(self.lb_y[ii+1])-np.array(self.lb_y[ii]))*np.linspace(0,1,10))
                self.lb_v_lines.set_data(v_line_x, np.array(v_line_x)/np.array(v_line_y))
                
        self.lb.figure.canvas.draw() # redraw figure

    def onKeyPress(self, event):
        
        if event.key == ' ':
            if self.pickingUpperBound and not self.pickingLowerBound:
                self.pickingUpperBound = False
                self.pickingLowerBound = True
                print("\tSelect LOWER bound using the mouse. When finished press the SPACE KEY.")
            elif not self.pickingUpperBound and self.pickingLowerBound:
                plt.close('all')
                
        elif event.key == '+':
            self.scale_velocity *= 1.5
            self.scale_wavenumber *= 1.5
            self.scale_ellipticity *= 1.5
            self.redrawDensity()
        elif event.key == '-':
            self.scale_velocity /= 1.5
            self.scale_wavenumber /= 1.5
            self.scale_ellipticity /= 1.5
            self.redrawDensity()
        #else:
        #    print(event)
        #    print("Key pressed: {0}".format(event.key))
        return
            
            
            
            
    def onMouseScroll(self, event):
        """
        Scrolling the mouse wheel it is possible to control the variance of the density
        """

        if event.button == 'up':
            self.scale_velocity *= 1.5
            self.scale_wavenumber *= 1.5
            self.scale_ellipticity *= 1.5
            self.redrawDensity()
        elif event.button == 'down':
            self.scale_velocity /= 1.5
            self.scale_wavenumber /= 1.5
            self.scale_ellipticity /= 1.5
            self.redrawDensity()
        
        # TODO would be nice to draw only once in case of multiple scroll
        return

            
    def filterPoints(self, x, y, data=None):
        """
        x   is the frequency axes
        y   is wavenumber / velocity / ellipticity
            the lower and upper bound are referred to this
        data   table of points to filter accordingly
        """
        
        assert len(self.lb_x) == len(self.lb_y)
        assert len(self.ub_x) == len(self.ub_y)
        assert len(x) == len(y) and (data is None or len(x) == np.size(data,0))
        
        if self.ellipticity:
            x  = np.concatenate((x, x, x))
            y  = np.concatenate((y, y-pi, y+pi))
        
        N = len(x)
        if len(self.lb_x) > 1 and len(self.ub_x) > 1: # a meaningful selection was done
            Nub = len(self.ub_x)
            Nlb = len(self.lb_x)
            
            ndxUB = np.array([False] * N)
            for ii in range(1, Nub):
                # find elements in right frequency range:
                ndx_tmp = (np.array(x)>self.ub_x[ii-1]) & (np.array(x)<self.ub_x[ii])
                # find elements also below the line:
                ndx_tmp = ndx_tmp & np.array(y < (x - self.ub_x[ii-1])*(self.ub_y[ii] - self.ub_y[ii-1])/(self.ub_x[ii] -self.ub_x[ii-1]) + self.ub_y[ii-1])
                ndxUB = ndxUB | ndx_tmp
                
            ndxLB = np.array([False] * N)
            for ii in range(1, Nlb):
                # find elements in right frequency range:
                ndx_tmp = (np.array(x)>self.lb_x[ii-1]) & (np.array(x)<self.lb_x[ii])
                # find elements also above the line:
                ndx_tmp = ndx_tmp & np.array(y > (x - self.lb_x[ii-1])*(self.lb_y[ii] - self.lb_y[ii-1])/(self.lb_x[ii] -self.lb_x[ii-1]) + self.lb_y[ii-1])
                ndxLB = ndxLB | ndx_tmp
            
            ndx_selection = ndxUB & ndxLB
            
        else: # no meaningful selection was done, perhaps no selection. So we keep all the points.
            ndx_selection = np.array([True] * N)

        if self.ellipticity:
            assert np.mod(N,3) == 0
            ndx_selection = ndx_selection[0:int(N/3)] | ndx_selection[int(N/3):int(2*N/3)] | ndx_selection[int(2*N/3):]

        
        if data is None:
            ret = (x[ndx_selection], y[ndx_selection])
        else:
            ret = (x[ndx_selection], y[ndx_selection], data[ndx_selection,:])
        return ret


    def pickPoints(self, x, y, mode='dispersion'):
        """
        select median curve from points
        """        
        
        assert len(x) == len(y)
        
        x_vec = np.unique(x)
        Nx = len(x_vec)
        pX = np.zeros((Nx,))
        pY = np.zeros((Nx,))
        pY_std = np.zeros((Nx,))
        if mode == 'dispersion':
            for xx in range(0,Nx):
                ndx = x == x_vec[xx]
                pX[xx] = x_vec[xx]
                pY[xx] = np.median(y[ndx])
                pY_std[xx] = np.std(y[ndx])
        elif mode == 'ellipticity':
            # ellipticity angle ranges between -pi/2 and pi/2
            y = np.mod(y +pi/2, pi) - pi/2
            for xx in range(0,Nx):
                ndx = x == x_vec[xx]
                pX[xx] = x_vec[xx]
                pY[xx] = np.mod( cs.circ_median( 2*(y[ndx]) ) /2 +pi/2, pi) -pi/2
                pY_std[xx] = cs.circ_std( 2*y[ndx] ) /2
                
        return (pX, pY, pY_std)
        
    def getScaleEllipticity(self):
        return self.scale_ellipticity
    
    def getScaleVelocity(self):
        return self.scale_velocity

    def getScaleWavenumber(self):
        return self.scale_wavenumber
        
    def getUpperBound(self):
        return (self.ub_x, self.ub_y)
        
    def getLowerBound(self):
        return (self.lb_x, self.lb_y)