#!/usr/bin/python3
# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
"""
Tool for picking and plotting WaveDec or WaveDecActive output files
See online documentation.
"""



# TODO
# Bounds should be normalized within limits otherwise they are plotted outsides

# plot array and array response, need to save array coordinates
# perhaps "Waves" could be avoided in file names
# make nice velocity plots, especially for publication




import numpy as np
import os, errno, csv, glob, argparse, shutil, inspect
import matplotlib
import ntpath
import matplotlib.pyplot as plt
from Picker import Picker
from PlotUtils import plotWavenumber, plotEllipticity, plotAzimuth, plotVelocity, \
                     plotArray, plotBounds, plotEllipticityAngle, plotArrayResponse, plotArrayResponseCuts

from wdSettings import MODEL_NOISE, MODEL_VERTICAL, MODEL_RAYLEIGH, MODEL_LOVE, MODEL_CIRCULAR_VERTICAL, MODEL_CIRCULAR_RAYLEIGH, MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH
from wdSettings import DEFAULT_IMAGE_FORMAT, DEFAULT_DPI
from wdSettings import DEFAULT_DIR_PLOTS, DEFAULT_DIR_PICKED, DEFAULT_DIR_FILTERED, DEFAULT_DIR_FILTERED_PLOTS
from wdSettings import DEFAULT_SAVE_PLOT, DEFAULT_SAVE_CSV
from wdSettings import DEFAULT_USETEX, DEFAULT_FONT_SIZE_1, DEFAULT_FONT_SIZE_2
from wdSettings import DEFAULT_SWITCH_ELLIPTICITY
    




imageFormat = DEFAULT_IMAGE_FORMAT
matplotlib.rcParams['text.usetex'] = DEFAULT_USETEX
matplotlib.rcParams['text.latex.unicode'] = True
plt.rc('font', family='serif', size=DEFAULT_FONT_SIZE_2)
plt.rc('xtick', labelsize=DEFAULT_FONT_SIZE_1)
plt.rc('ytick', labelsize=DEFAULT_FONT_SIZE_1)
plt.rc('legend', fontsize=DEFAULT_FONT_SIZE_2)

def main():
    
    
    ### Initialization: commandline args
    parser = argparse.ArgumentParser(description='wdPicker - a tool for the analysis of  WaveDec output files')
    parser.add_argument("--input",  default='./', help="Path to input folder [./]")
    parser.add_argument("--output",  default='./', help="Path to output folder [./]")
    
    parser.set_defaults(savePlotsToFile=DEFAULT_SAVE_PLOT)
    if DEFAULT_SAVE_PLOT:
         parser.add_argument('--dontsave_plots', dest='savePlotsToFile', action='store_false',
                             help="Do not save pictures to disk. By default plots are saved to disk.")
    else:
        parser.add_argument('--save_plots', dest='savePlotsToFile', action='store_true',
                            help="Save pictures to disk. By default plots are not saved to disk.")
        
    parser.set_defaults(saveCSVToFile=DEFAULT_SAVE_CSV)
    if DEFAULT_SAVE_CSV:
        parser.add_argument('--dontsave_csv', dest='saveCSVToFile', action='store_false',
                            help="Do not save CSV files to disk. By default CSV are saved to disk.")
    else:
        parser.add_argument('--save_csv', dest='saveCSVToFile', action='store_true',
                            help="Save CSV files to disk. By default CSV are not saved to disk.")
    
    parser.add_argument("--format", default=DEFAULT_IMAGE_FORMAT,
                        help="Image format [{0}]. Possible formats include 'png', 'pdf', 'eps'".format(DEFAULT_IMAGE_FORMAT))
    parser.add_argument("--dpi", default=DEFAULT_DPI, help="Image resolution in dpi [{0}]".format(DEFAULT_DPI))
    args = parser.parse_args()
    
    plt.ion()

    # Look for input files
    fileVertical            = glob.glob(os.path.join(args.input, 'VerticalWaves*.csv'))
    fileLove                = glob.glob(os.path.join(args.input, 'LoveWaves*.csv'))
    fileRayleigh            = glob.glob(os.path.join(args.input, 'RayleighWaves*.csv'))
    fileCircularVertical    = glob.glob(os.path.join(args.input, 'CircularVerticalWaves*.csv'))
    fileCircularRayleigh    = glob.glob(os.path.join(args.input, 'CircularRayleighWaves*.csv'))
    fileCircularDissipativeRayleigh    = glob.glob(os.path.join(args.input, 'CircularDissipativeRayleighWaves*.csv'))
    fileArrayResolution     = os.path.join(args.input, 'ArrayResolutionLimits.csv')
    fileArrayLayout         = os.path.join(args.input, 'ArrayLayout.csv')
    fileSourcePosition      = os.path.join(args.input, 'SourcePosition.csv')


    fileTheoreticalDispersion = None
    fileTheoreticalEllipticity = None
    
    ### M2.1
    ## DEFAULT_SWITCH_ELLIPTICITY set it to True in wdSettings.py
#    colorTheoretical = 'lime'
#    fileTheoreticalDispersion = '/home/marra/Dropbox/ETH/EllipticityPaper/doc/pictures/tikz/M2.1_RayleighVelocity_0.out'
#    fileTheoreticalEllipticity = '/home/marra/Dropbox/ETH/EllipticityPaper/doc/pictures/tikz/M2.1_RayleighEllipticity_0.out'
#    fileTheoreticalDispersion = '/home/marra/Dropbox/ETH/EllipticityPaper/doc/pictures/tikz/M2.1_RayleighVelocity_1.out'
#    fileTheoreticalEllipticity = '/home/marra/Dropbox/ETH/EllipticityPaper/doc/pictures/tikz/M2.1_RayleighEllipticity_1.out'
#    fileTheoreticalDispersion = '/home/marra/Dropbox/ETH/EllipticityPaper/doc/pictures/tikz/M2.1_LoveVelocity_0.out'
    
    
    
    savePlotsToFile = args.savePlotsToFile
    saveCSVToFile = args.saveCSVToFile
    
    imageFormat = args.format
    imageDPI = args.dpi
    
    outputDirPlot = os.path.join(args.output, DEFAULT_DIR_PLOTS)
    outputDirPicked = os.path.join(args.output, DEFAULT_DIR_PICKED)
    outputDirFiltered = os.path.join(args.output, DEFAULT_DIR_FILTERED)
    outputDirFilteredPlot = os.path.join(args.output, DEFAULT_DIR_FILTERED_PLOTS)


    
    waveFiles = fileVertical + fileLove + fileRayleigh \
                + fileCircularVertical + fileCircularRayleigh
    waveTypes = [MODEL_VERTICAL]*len(fileVertical) + [MODEL_LOVE]*len(fileLove) + [MODEL_RAYLEIGH]*len(fileRayleigh) \
                + [MODEL_CIRCULAR_VERTICAL]*len(fileCircularVertical) + [MODEL_CIRCULAR_RAYLEIGH]*len(fileCircularRayleigh)
                
    waveLabels = {MODEL_VERTICAL:'Vertical', MODEL_LOVE:'Love', MODEL_RAYLEIGH:'Rayleigh',
                  MODEL_CIRCULAR_VERTICAL:'Circular Vertical', MODEL_CIRCULAR_RAYLEIGH:'Circular Rayleigh', MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH:'Circular Dissipative Rayleigh',
                  MODEL_NOISE:'Noise'}
    
    if len(waveFiles) < 1:
        path = os.path.abspath(args.input)
        print("No WaveDec output files found in {0}".format(path))
        return
    else:
        print("Reading the following WaveDec output files:\n\t{0}".format(waveFiles))
        
    if savePlotsToFile:
        print("Will save figures to file ({0})".format(imageFormat))
        if DEFAULT_USETEX: print("\tUsing LaTex text rendering")
    if saveCSVToFile: print("Will save CSV to file")
        
        
    # Columns index
    COL_Frequency = 0
    COL_Amplitude = 1
    COL_Wavenumber = 2
    COL_Velocity = 3
    COL_Azimuth = 4     # this column is missing for circular waves
    COL_EllipticityAngle = 5
    COL_Wavenumber_i = 6
    COL_Attenuation = 7
    
    
    if savePlotsToFile:
        try: os.makedirs(outputDirPlot)
        except OSError as exception:
            if exception.errno != errno.EEXIST: raise
            

    
    
   
    # if available, read array resolution limits
    plotArrayResolution = False # plot array resolution limits. True only if the file is present!
    if os.path.isfile(fileArrayResolution):
        ArrayResolution = np.atleast_2d(np.loadtxt(fileArrayResolution,comments='#', delimiter='\t'))
        if np.size(ArrayResolution) > 0:
            Fres = ArrayResolution[:,COL_Frequency]
            Kmin = ArrayResolution[:,1]
            Kmax = ArrayResolution[:,2]
            Vmin = ArrayResolution[:,3]
            Vmax = ArrayResolution[:,4]
            plotArrayResolution = True # Yes, we will plot resolution limits
        else:
            print("Array resolution limit file ({0}) is empty.".format(fileArrayResolution))
            
    if os.path.isfile(fileArrayLayout):
        ArrayLayout = np.atleast_2d(np.loadtxt(fileArrayLayout,comments='#', delimiter='\t'))
        if np.size(ArrayLayout) > 0:
            pos = ArrayLayout[:,0:2]
            if os.path.isfile(fileSourcePosition):
                SourcePosition = np.atleast_1d(np.loadtxt(fileSourcePosition, comments='#', delimiter='\t'))
            else:
                SourcePosition = None
            plotArray(pos, SourcePosition)

            if savePlotsToFile: plt.savefig('{0}/ArrayLayout.{1}'.format(outputDirPlot, imageFormat), dpi=imageDPI)
                
            if plotArrayResolution and not os.path.isfile(fileSourcePosition): # for WaveDecActive, do not plot array response
                plotArrayResponse(pos, 1.5*Kmax[0])
                if plotArrayResolution:
                    Thetavec = np.linspace(0, 2*np.pi, 180)
                    plt.plot(Kmin[0]*np.cos(Thetavec), Kmin[0]*np.sin(Thetavec),'y-',linewidth=2)
                    plt.plot(Kmax[0]*np.cos(Thetavec), Kmax[0]*np.sin(Thetavec),'y-',linewidth=2)
                if savePlotsToFile: plt.savefig('{0}/ArrayResponse.{1}'.format(outputDirPlot, imageFormat), dpi=imageDPI)
                    
                
                plotArrayResponseCuts(pos, 1.5*Kmax[0])
                if plotArrayResolution:
                    plt.plot([Kmin[0], Kmin[0]], [0, 1],'y-',linewidth=2)
                    plt.plot([-Kmin[0], -Kmin[0]], [0, 1],'y-',linewidth=2)
                    plt.plot([Kmax[0], Kmax[0]], [0, 1],'y-',linewidth=2)
                    plt.plot([-Kmax[0], -Kmax[0]], [0, 1],'y-',linewidth=2)
                if savePlotsToFile: plt.savefig('{0}/ArrayResponseCuts.{1}'.format(outputDirPlot, imageFormat), dpi=imageDPI)
            
        else:
            print("Array layout file ({0}) is empty.".format(fileArrayResolution))
            
    # load files with syntetic dispersion curves
    if fileTheoreticalDispersion is not None:
        tData = np.atleast_2d(np.loadtxt(fileTheoreticalDispersion, comments='#', delimiter=' '))
        tDispFreq = tData[:,0]
        tDispSlowness = tData[:,1]
        tDispVelocity = 1/ tDispSlowness
        tDispWavenumber = tDispFreq * tDispSlowness
        
    if fileTheoreticalEllipticity is not None:
        tData = np.atleast_2d(np.loadtxt(fileTheoreticalEllipticity, comments='#', delimiter=' '))
        tEllipticityFreq = tData[:,0]
        tEllipticityHV = tData[:,1]
        tEllipticityAngle = np.arctan(tEllipticityHV)

            
            
    for waveFile in waveFiles:
        ww = waveFiles.index(waveFile)
        data = np.atleast_2d(np.loadtxt(waveFile, comments='#', delimiter='\t'))

        if np.size(data) > 0:
            waveType = waveTypes[ww]
            waveLabel = waveLabels[waveType]
            fileBasename = ntpath.basename(waveFile).split('.')[0]
            F = data[:,COL_Frequency]
            K = data[:,COL_Wavenumber]
            V = data[:,COL_Velocity]
            
            
            plotWavenumber(F, K)
            if fileTheoreticalDispersion is not None: plt.plot(tDispFreq, tDispWavenumber, colorTheoretical, linewidth=2)
            if plotArrayResolution: plt.plot(Fres, Kmin,'y-',linewidth=2),plt.plot(Fres, Kmax,'y-',linewidth=2)
            if savePlotsToFile: plt.savefig('{0}/{1}_Wavenumber.{2}'.format(outputDirPlot, fileBasename, imageFormat), dpi=imageDPI)
            plt.title(waveLabel)
            

            
            plotVelocity(F, V)
            if fileTheoreticalDispersion is not None: plt.plot(tDispFreq, tDispVelocity, colorTheoretical, linewidth=2)
            if plotArrayResolution:
                plt.plot(Fres, Vmin,'y-',linewidth=2),plt.plot(Fres, Vmax,'y-',linewidth=2)
                Fvec = np.unique(F) # WARNING this vector is not necessarily equally spaced! (eg, log-spacing, approximation)
                if len(Fvec) > 1: Fvec_pcolor = 0.5*( np.concatenate((np.array([2*Fvec[0] -Fvec[1]]), Fvec )) + np.concatenate((Fvec, np.array([2*Fvec[-1] -Fvec[-2]]))) )            
                else: Fvec_pcolor = np.array([Fvec[0], Fvec[0]+0.5])
                plt.xlim([Fvec_pcolor.min(), Fvec_pcolor.max()])
            if savePlotsToFile: plt.savefig('{0}/{1}_Velocity.{2}'.format(outputDirPlot, fileBasename, imageFormat), dpi=imageDPI)
            plt.title(waveLabel)
 
            # plotting azimuth
            if waveType == MODEL_RAYLEIGH or waveType == MODEL_LOVE or waveType == MODEL_VERTICAL:
                A = np.mod(data[:,COL_Azimuth], 2*np.pi)
                plotAzimuth(F, A)
                if savePlotsToFile: plt.savefig('{0}/{1}_Azimuth.{2}'.format(outputDirPlot, fileBasename, imageFormat), dpi=imageDPI)
                plt.title(waveLabel)
                
            # plotting ellipticity    
            if waveType == MODEL_RAYLEIGH or waveType == MODEL_CIRCULAR_RAYLEIGH:
                if waveType == MODEL_RAYLEIGH:
                    E = np.mod(data[:,COL_EllipticityAngle] + np.pi/2, np.pi) - np.pi/2
                elif waveType == MODEL_CIRCULAR_RAYLEIGH:
                    E = np.mod(data[:,COL_EllipticityAngle - 1] + np.pi/2, np.pi) - np.pi/2 # For circular waves there is no azimuth colum. This explain the -1
                else:
                    return
                    
                if DEFAULT_SWITCH_ELLIPTICITY:
                    print("\n*** Switching ellipticity sign! ***")
                    print("\t this may be necessary if the UDZ component points downwards\n")
                    E = -E
                
                plotEllipticity(F, E)
                if fileTheoreticalEllipticity is not None: plotEllipticityAngle(tEllipticityFreq, tEllipticityAngle, E_std=None, color=colorTheoretical, linewidth=2, alpha=1)
                if savePlotsToFile: plt.savefig('{0}/{1}_EllipticityAngle.{2}'.format(outputDirPlot, fileBasename, imageFormat), dpi=imageDPI)
                plt.title(waveLabel)
    
                # Dispersion and Ellipticity angle on same plot
                fig, axes = plt.subplots(2,1, sharex=True)
                plotWavenumber(F, K, ax=axes[0])
                if plotArrayResolution: axes[0].plot(Fres, Kmin,'y-',linewidth=2),axes[0].plot(Fres, Kmax,'y-',linewidth=2)
                axes[0].set_title('')
                axes[0].set_xlabel('')
                plotEllipticity(F, E, ax=axes[1])
                axes[1].set_title('')
                #fig.tight_layout()
                plt.show()
                if savePlotsToFile: plt.savefig('{0}/{1}_Wavenumber_EllipticityAngle.{2}'.format(outputDirPlot, fileBasename, imageFormat), dpi=imageDPI)
                    
          
        else: # if any file is empty, we remove it from the list
            print("File '{0}' is empty.".format(waveFile))
            del waveFiles[ww]
            del waveTypes[ww]
   

        
        
    ### and now manually select dispersion boundaries

    
    while True: # A while loop for selection/filter/pick.

        while True: # Select which wave to pick
            print("\nWould you like to manually select and filter ML estimates?")
            print("What would you like to filter?")
            print("0 - Nothing. Exit.")
                
            for ww in range(0,len(waveFiles)):
                if waveTypes[ww] == MODEL_VERTICAL:
                    print("{0} - Vertical Wave          \t({1})".format(ww+1, waveFiles[ww]))
                elif waveTypes[ww] == MODEL_LOVE:
                    print("{0} - Love Wave              \t({1})".format(ww+1, waveFiles[ww]))
                elif waveTypes[ww] == MODEL_RAYLEIGH:
                    print("{0} - Rayleigh Wave          \t({1})".format(ww+1, waveFiles[ww]))
                elif waveTypes[ww] == MODEL_CIRCULAR_VERTICAL:
                    print("{0} - Circular Vertical Wave \t({1})".format(ww+1, waveFiles[ww]))
                elif waveTypes[ww] == MODEL_CIRCULAR_RAYLEIGH:
                    print("{0} - Circular Rayleigh Wave \t({1})".format(ww+1, waveFiles[ww]))
            
            try: cc = int(input(":"))
            except: cc = -1

            if 0 <= cc <= len(waveFiles):
                plt.close('all')
                break
        
        if cc == 0:
            print("Exiting.")
            break
        else:
            cc = cc - 1
            if savePlotsToFile or saveCSVToFile:
                try: os.makedirs(outputDirFiltered)
                except OSError as exception:
                    if exception.errno != errno.EEXIST: raise
                try: os.makedirs(outputDirPicked)
                except OSError as exception:
                    if exception.errno != errno.EEXIST: raise
            if savePlotsToFile:
                try: os.makedirs(outputDirFilteredPlot)
                except OSError as exception:
                    if exception.errno != errno.EEXIST: raise            
                
                if plotArrayResolution:
                    try: shutil.copyfile(fileArrayResolution, os.path.join(outputDirFiltered, ntpath.basename(fileArrayResolution)))
                    except OSError as exception:
                        if exception.errno != errno.EEXIST: raise 
            
        waveType = waveTypes[cc]
        waveFile = waveFiles[cc]
        waveLabel = waveLabels[waveType]
        
        if waveType == MODEL_VERTICAL:
            print("\nPicking Vertical Wave        ({0})".format(waveFile))
            print("Insert label for the mode. (Eg: V0, V1, ...)")
        elif waveType == MODEL_LOVE:
            print("\nPicking Love Wave              ({0})".format(waveFile))
            print("Insert label for the mode. (Eg: L0, L1, ...)")
        elif waveType == MODEL_RAYLEIGH:
            print("\nPicking Rayleigh Wave          ({0})".format(waveFile))
            print("Insert label for the mode. (Eg: R0, R1, ...)")
        elif waveType == MODEL_CIRCULAR_VERTICAL:
            print("Picking Circular Vertical Wave ({0})".format(waveFile))
            print("Insert label for the mode. (Eg: V0, V1, ...)")
        elif waveType == MODEL_CIRCULAR_RAYLEIGH:
            print("Picking Circular Rayleigh Wave ({0})".format(waveFile))
            print("Insert label for the mode. (Eg: R0, R1, ...)")
        else:
            print("Unrecognized option. Exiting.")
            break

        ModeLabel = str(input(":"))
        # check the input or Latex will complain later
        while not all(s.isnumeric() or s.isalpha() for s in ModeLabel):
            print("Label may contain only letters and numbers (a-zA-Z0-9)")
            ModeLabel = str(input(":"))

        
        data = np.atleast_2d(np.loadtxt(waveFile, comments='#', delimiter='\t'))
        fileBasename = ntpath.basename(waveFile).split('.')[0]
        F = data[:,COL_Frequency]
        K = data[:,COL_Wavenumber]
        V = data[:,COL_Velocity]        


        if waveType == MODEL_RAYLEIGH:
            E = np.mod(data[:,COL_EllipticityAngle] + np.pi/2, np.pi) - np.pi/2
        elif waveType == MODEL_CIRCULAR_RAYLEIGH:
            E = np.mod(data[:,COL_EllipticityAngle - 1] + np.pi/2, np.pi) - np.pi/2 # For circular waves there is no azimuth colum. This explain the -1
            
        if DEFAULT_SWITCH_ELLIPTICITY:
            E = -E
            
        # axes limits used in all plots
        Flim = np.array([np.min(F), np.max(F)])
        Klim = np.array([np.min(K), np.max(K)*1.05])

        print("Step one, selecting points in the wavenumber-frequency plane")

        if plotArrayResolution:
          picker_K = Picker(F, K, 'dispersion', (Fres, Kmin, Kmax, Vmin, Vmax) )
        else:
          picker_K = Picker(F, K, 'dispersion')
        picker_K.plotAndSelect()
        plt.ioff()
        plt.show()            
        (fF_1, fK_1, fData_1) = picker_K.filterPoints(F, K, data)                        
        plt.ion()
        
        if waveType == MODEL_RAYLEIGH or waveType == MODEL_CIRCULAR_RAYLEIGH:
            print("Step two, selecting points in the ellipticity angle-frequency plane")
            if waveType == MODEL_RAYLEIGH:
                fE_1 = np.mod(fData_1[:,COL_EllipticityAngle] + np.pi/2, np.pi) - np.pi/2
            elif  waveType == MODEL_CIRCULAR_RAYLEIGH:
                fE_1 = np.mod(fData_1[:,COL_EllipticityAngle-1] + np.pi/2, np.pi) - np.pi/2
                
            if DEFAULT_SWITCH_ELLIPTICITY:
                fE_1 = -fE_1
                    
            picker_E = Picker(fF_1, fE_1, 'ellipticity')
            picker_E.plotAndSelect()
            plt.ioff()
            plt.show()
            (fF_2, fE_2, fData_2) = picker_E.filterPoints(fF_1, fE_1, fData_1)
            fK_2 = fData_2[:,COL_Wavenumber]
            plt.ion()
        else:
            fF_2 = fF_1
            fK_2 = fK_1
            fData_2 = fData_1
            
        fV_2 =  fData_2[:,COL_Velocity] 
        # Pick curves
        (pF, pK, pK_std) = picker_K.pickPoints(fF_2, fK_2, 'dispersion')
        pV = pF/pK
        if waveType == MODEL_RAYLEIGH or waveType == MODEL_CIRCULAR_RAYLEIGH:
            (pF, pE, pE_std) = picker_E.pickPoints(fF_2, fE_2, 'ellipticity')
        

        
        # plotting wavenumber
        fig = plotWavenumber(F, K, scale=picker_K.getScaleWavenumber(), xlim=Flim, ylim=Klim)
        if fileTheoreticalDispersion is not None: plt.plot(tDispFreq, tDispWavenumber, colorTheoretical,linewidth=2)
        plt.plot(pF, pK, 'b-', linewidth=2)
        plt.errorbar(pF, pK, yerr=pK_std, fmt='none', color='b', zorder=10)
        if plotArrayResolution: plt.plot(Fres, Kmin,'y-',linewidth=2),plt.plot(Fres, Kmax,'y-',linewidth=2)
        if savePlotsToFile: plt.savefig('{0}/{1}_{2}_Wavenumber.{3}'.format(outputDirFilteredPlot, fileBasename, ModeLabel, imageFormat), dpi=imageDPI)
        (ub_x, ub_y) = picker_K.getUpperBound(); (lb_x, lb_y) = picker_K.getLowerBound()
        plotBounds(ub_x, ub_y, lb_x, lb_y, fig.axes)
        if savePlotsToFile: plt.savefig('{0}/{1}_{2}_Wavenumber_B.{3}'.format(outputDirFilteredPlot, fileBasename, ModeLabel, imageFormat), dpi=imageDPI)
        plt.title('{0} {1}'.format(waveLabel, ModeLabel))
            
        # plotting velocity
        fig = plotVelocity(F, V, '{0} {1}'.format(waveLabel, ModeLabel), scale=picker_K.getScaleVelocity(), xlim=Flim)
        if fileTheoreticalDispersion is not None: plt.plot(tDispFreq, tDispVelocity, colorTheoretical,linewidth=2)
        if plotArrayResolution: plt.plot(Fres, Vmin,'y-',linewidth=2),plt.plot(Fres, Vmax,'y-',linewidth=2)
        plt.plot(pF, pV, 'b-', linewidth=2)
        if savePlotsToFile: plt.savefig('{0}/{1}_{2}_Velocity.{3}'.format(outputDirFilteredPlot, fileBasename, ModeLabel, imageFormat), dpi=imageDPI)
        
        fig = plotVelocity(fF_2, fV_2, '{0} {1}'.format(waveLabel, ModeLabel), xlim=Flim)
        if plotArrayResolution: plt.plot(Fres, Vmin,'y-',linewidth=2),plt.plot(Fres, Vmax,'y-',linewidth=2)
        plt.plot(pF, pV, 'b-', linewidth=2)
        if savePlotsToFile: plt.savefig('{0}/{1}_{2}_Velocity_2.{3}'.format(outputDirFilteredPlot, fileBasename, ModeLabel, imageFormat), dpi=imageDPI)


        # plotting ellipticity angle
        if waveType == MODEL_RAYLEIGH or waveType == MODEL_CIRCULAR_RAYLEIGH:
            
            fig = plotEllipticity(fF_1, fE_1, scale=picker_E.getScaleEllipticity(), xlim=Flim)
            if fileTheoreticalEllipticity is not None: plotEllipticityAngle(tEllipticityFreq, tEllipticityAngle, E_std=None, color=colorTheoretical, linewidth=2, alpha=1)
            plotEllipticityAngle(pF, pE, E_std=pE_std, color='b', linewidth=2, alpha=1)
            if savePlotsToFile: plt.savefig('{0}/{1}_{2}_EllipticityAngle.{3}'.format(outputDirFilteredPlot, fileBasename, ModeLabel, imageFormat), dpi=imageDPI)
            (ub_x, ub_y) = picker_E.getUpperBound(); (lb_x, lb_y) = picker_E.getLowerBound()
            plotBounds(ub_x, ub_y, lb_x, lb_y, fig.axes[0])
            if savePlotsToFile: plt.savefig('{0}/{1}_{2}_EllipticityAngle_B.{3}'.format(outputDirFilteredPlot, fileBasename, ModeLabel, imageFormat), dpi=imageDPI)
            plt.title('{0} {1}'.format(waveLabel, ModeLabel))
                
        # plotting azimuth
        if waveType == MODEL_RAYLEIGH or waveType == MODEL_LOVE or waveType == MODEL_VERTICAL:
            fA_2 = np.mod(fData_2[:,COL_Azimuth], 2*np.pi)
            plotAzimuth(fF_2, fA_2, 'Azimuth {0}'.format(ModeLabel))
            if savePlotsToFile: plt.savefig('{0}/{1}_{2}_Azimuth.{3}'.format(outputDirFilteredPlot, fileBasename, ModeLabel, imageFormat), dpi=imageDPI)
            plt.title('{0} {1}'.format(waveLabel, ModeLabel))
                   
                
        

        ### save to CSV file filtered estimates and picked points
        if saveCSVToFile:
            if waveType == MODEL_VERTICAL:
                # save CSV with filtered points
                with open(os.path.join(outputDirFiltered, '{0}_{1}.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec output file for vertical-component waves'])
                    writer.writerow(['# MODEL_VERTICAL'])
                    writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'Azimuth'])
                    writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]'])
                    for row in fData_2:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                # save CSV with picked curves               
                with open(os.path.join(outputDirPicked, '{0}_{1}_Dispersion_picked.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec picked dispersion curve for Vertical wave'])
                    writer.writerow(['# MODEL_VERTICAL'])
                    writer.writerow(['# Frequency', 'Wavenumber' ,'Std(Wavenumber)', 'Velocity'])
                    writer.writerow(['# [Hz]', '[1/m]', '[1/m]', '[m/s]'])                
                    for row in np.vstack((pF.T, pK.T, pK_std.T, pV.T)).T:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])

            elif waveType == MODEL_LOVE:
                with open(os.path.join(outputDirFiltered, '{0}_{1}.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec output file for Love waves'])
                    writer.writerow(['# MODEL_LOVE'])
                    writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'Azimuth'])
                    writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]'])
                    for row in fData_2:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                # save CSV with picked curves               
                with open(os.path.join(outputDirPicked, '{0}_{1}_Dispersion_picked.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec picked dispersion curve for Love wave'])
                    writer.writerow(['# MODEL_LOVE'])
                    writer.writerow(['# Frequency', 'Wavenumber' ,'Std(Wavenumber)', 'Velocity'])
                    writer.writerow(['# [Hz]', '[1/m]', '[1/m]', '[m/s]'])                
                    for row in np.vstack((pF.T, pK.T, pK_std.T, pV.T)).T:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])

            elif waveType == MODEL_RAYLEIGH:
                with open(os.path.join(outputDirFiltered, '{0}_{1}.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec output file for Rayleigh waves'])
                    writer.writerow(['# MODEL_RAYLEIGH'])
                    writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'Azimuth', 'EllipticityAngle'])
                    writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]', '[rad]'])
                    for row in fData_2:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                with open(os.path.join(outputDirPicked, '{0}_{1}_Dispersion_picked.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec picked dispersion curve for Rayleigh wave'])
                    writer.writerow(['# MODEL_RAYLEIGH'])
                    writer.writerow(['# Frequency', 'Wavenumber' ,'Std(Wavenumber)', 'Velocity'])
                    writer.writerow(['# [Hz]', '[1/m]', '[1/m]', '[m/s]'])                
                    for row in np.vstack((pF.T, pK.T, pK_std.T, pV.T)).T:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                with open(os.path.join(outputDirPicked, '{0}_{1}_EllipticityAngle_picked.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec picked ellipticity angle curve for Rayleigh wave'])
                    writer.writerow(['# MODEL_RAYLEIGH'])
                    writer.writerow(['# Frequency', 'EllipticityAngle', 'Std(EllipticityAngle)'])
                    writer.writerow(['# [Hz]', '[rad]', '[rad]'])
                    for row in np.vstack((pF.T,pE.T,pE_std.T)).T:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                        
            elif waveType == MODEL_CIRCULAR_VERTICAL:
                with open(os.path.join(outputDirFiltered, '{0}_{1}.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec output file for circular vertical-component waves'])
                    writer.writerow(['# MODEL_CIRCULAR_VERTICAL'])
                    writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity'])
                    writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]'])
                    for row in fData_2:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                with open(os.path.join(outputDirPicked, '{0}_{1}_Dispersion_picked.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec picked dispersion curve for Circular Vertical wave'])
                    writer.writerow(['# MODEL_CIRCULAR_VERTICAL'])
                    writer.writerow(['# Frequency', 'Wavenumber' ,'Std(Wavenumber)', 'Velocity'])
                    writer.writerow(['# [Hz]', '[1/m]', '[1/m]', '[m/s]'])                
                    for row in np.vstack((pF.T, pK.T, pK_std.T, pV.T)).T:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
            elif waveType == MODEL_CIRCULAR_RAYLEIGH:
                with open(os.path.join(outputDirFiltered, '{0}_{1}.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec output file for Circular Rayleigh waves'])
                    writer.writerow(['# MODEL_CIRCULAR_RAYLEIGH'])
                    writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'EllipticityAngle'])
                    writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]'])
                    for row in fData_2:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                with open(os.path.join(outputDirPicked, '{0}_{1}_Dispersion_picked.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec picked dispersion curve for Circular Rayleigh wave'])
                    writer.writerow(['# MODEL_CIRCULAR_RAYLEIGH'])
                    writer.writerow(['# Frequency', 'Wavenumber' ,'Std(Wavenumber)', 'Velocity'])
                    writer.writerow(['# [Hz]', '[1/m]', '[1/m]', '[m/s]'])                
                    for row in np.vstack((pF.T, pK.T, pK_std.T, pV.T)).T:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
                with open(os.path.join(outputDirPicked, '{0}_{1}_EllipticityAngle_picked.csv'.format(fileBasename, ModeLabel)), 'w', encoding='UTF8') as ff:
                    writer = csv.writer(ff, delimiter='\t')
                    writer.writerow(['# WaveDec picked ellipticity angle curve for Circular Rayleigh wave'])
                    writer.writerow(['# MODEL_CIRCULAR_RAYLEIGH'])
                    writer.writerow(['# Frequency', 'EllipticityAngle', 'Std(EllipticityAngle)'])
                    writer.writerow(['# [Hz]', '[rad]', '[rad]'])
                    for row in np.vstack((pF.T,pE.T,pE_std.T)).T:
                        writer.writerow(['{:6.3e}'.format(c) for c in row])
            else:
                print("Unrecognized option. Exiting.")
                break


if  __name__ =='__main__':
    main()

