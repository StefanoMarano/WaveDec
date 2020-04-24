# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################


import sqlite3
import numpy as np
import os
import csv
import logging
from wdSettings import MODEL_VERTICAL, MODEL_RAYLEIGH, MODEL_LOVE, MODEL_CIRCULAR_VERTICAL, MODEL_CIRCULAR_RAYLEIGH, MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH


def init():
    
    conn = sqlite3.connect(':memory:')
    cur = conn.cursor()
    # cur.execute('SELECT SQLITE_VERSION()') 
    # data = cur.fetchone()
    # print("SQLite version: {0}".format(data))
    cur.execute('''CREATE TABLE Windows (
        WindowId INTEGER PRIMARY KEY,
        Kstart INTEGER, Kend INTEGER,
        Ts FLOAT)''')
    cur.execute('''CREATE TABLE Frequencies (
        Fndx INTEGER PRIMARY KEY,
        F FLOAT)''')
    cur.execute('''CREATE TABLE RayleighWaves (
        WaveId INTEGER PRIMARY KEY,
        WindowId INTEGER,
        Fndx INTEGER,
        Amplitude FLOAT, Phase FLOAT, Wavenumber FLOAT, Azimuth FLOAT, EllipticityAngle FLOAT,
        FOREIGN KEY(Fndx) REFERENCES Frequencies(Fndx),
        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
        )''')
    cur.execute('''CREATE TABLE CircularRayleighWaves (
        WaveId INTEGER PRIMARY KEY,
        WindowId INTEGER,
        Fndx INTEGER,
        Amplitude FLOAT, Phase FLOAT, Wavenumber FLOAT, EllipticityAngle FLOAT,
        FOREIGN KEY(Fndx) REFERENCES Frequencies(Fndx),
        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
        )''')
    cur.execute('''CREATE TABLE CircularDissipativeRayleighWaves (
        WaveId INTEGER PRIMARY KEY,
        WindowId INTEGER,
        Fndx INTEGER,
        Amplitude FLOAT, Phase FLOAT, Wavenumber FLOAT, Wavenumber_i FLOAT, EllipticityAngle FLOAT,
        FOREIGN KEY(Fndx) REFERENCES Frequencies(Fndx),
        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
        )''')
    cur.execute('''CREATE TABLE VerticalWaves (
        WaveId INTEGER PRIMARY KEY,
        WindowId INTEGER,
        Fndx INTEGER,
        Amplitude FLOAT, Phase FLOAT, Wavenumber FLOAT, Azimuth FLOAT,
        FOREIGN KEY(Fndx) REFERENCES Frequencies(Fndx),
        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
        )''')
    cur.execute('''CREATE TABLE CircularVerticalWaves (
        WaveId INTEGER PRIMARY KEY,
        WindowId INTEGER,
        Fndx INTEGER,
        Amplitude FLOAT, Phase FLOAT, Wavenumber FLOAT,
        FOREIGN KEY(Fndx) REFERENCES Frequencies(Fndx),
        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
        )''')
    cur.execute('''CREATE TABLE LoveWaves (
        WaveId INTEGER PRIMARY KEY,
        WindowId INTEGER,
        Fndx INTEGER,
        Amplitude FLOAT, Phase FLOAT, Wavenumber FLOAT, Azimuth FLOAT,
        FOREIGN KEY(Fndx) REFERENCES Frequencies(Fndx),
        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
        )''')
#    cur.execute('''CREATE TABLE CircularLoveWaves (
#        WaveId INTEGER PRIMARY KEY,
#        WindowId INTEGER,
#        Fndx INTEGER,
#        Amplitude FLOAT, Phase FLOAT, Wavenumber FLOAT,
#        FOREIGN KEY(Fndx) REFERENCES Frequencies(Fndx),
#        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
#        )''')
    cur.execute('''CREATE TABLE Noise (
        WindowId INTEGER,
        ChannelId INTEGER,
        sigma2 FLOAT,
        FOREIGN KEY(WindowId) REFERENCES Windows(WindowId)
        )''')
    cur.execute('''CREATE TABLE Sensors (
        SensorId INTEGER PRIMARY KEY,
        x FLOAT,
        y FLOAT,
        z FLOAT)''')
    cur.execute('''CREATE TABLE Channels (
        ChannelId INTEGER PRIMARY KEY,
        SensorId INTEGER,
        type TEXT    ,
        FOREIGN KEY(SensorId) REFERENCES Sensors(SensorId))''')

    conn.commit()
    cur.close()
    
    return conn

def createOutputFiles(conn, OUTPUT, WavesToModel, pos, Fvec, resolution_Kmin, resolution_Kmax, SourcePosition=None):
    
    if resolution_Kmin is not None and resolution_Kmax is not None:
        with open(os.path.join(OUTPUT, 'ArrayResolutionLimits.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDec output file with array resolution limits'])
            writer.writerow(['# Frequency', 'Kmin', 'Kmax', 'Vmin', 'Vmax'])
            writer.writerow(['# [Hz]', '[1/m]', '[1/m]', '[m/s]', '[m/s]'])
            for F in Fvec:
                Vmin = F / resolution_Kmax
                Vmax = F / resolution_Kmin
                row = [F, resolution_Kmin, resolution_Kmax, Vmin, Vmax]
                writer.writerow(['{:6.3e}'.format(c) for c in row])
            
    with open(os.path.join(OUTPUT, 'ArrayLayout.csv'), 'w', encoding='UTF8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['# WaveDec output file with array layout'])
        writer.writerow(['# Easting', 'Northing', 'Elevation'])
        writer.writerow(['# [m]', '[m]', '[m]'])
        cur = conn.cursor()
        cur.execute('''SELECT x, y, z FROM Sensors''')
        data = cur.fetchall()
        for row in data:
            writer.writerow(['{:.3f}'.format(c) for c in row])
    
    if SourcePosition is not None:
        with open(os.path.join(OUTPUT, 'SourcePosition.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDecActive output file with source position'])
            writer.writerow(['# Easting', 'Northing'])
            writer.writerow(['# [m]', '[m]'])
            writer.writerow(['{:.3f}'.format(SourcePosition[0]), '{:.3f}'.format(SourcePosition[1])])
            
    
    if MODEL_VERTICAL in WavesToModel and WavesToModel[MODEL_VERTICAL]:
        with open(os.path.join(OUTPUT, 'VerticalWaves.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDec output file for vertical-component waves'])
            writer.writerow(['# MODEL_VERTICAL'])
            writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'Azimuth'])
            writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]'])
            
    if MODEL_LOVE in WavesToModel and WavesToModel[MODEL_LOVE]:
        with open(os.path.join(OUTPUT, 'LoveWaves.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDec output file for Love waves'])
            writer.writerow(['# MODEL_LOVE'])
            writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'Azimuth'])
            writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]'])
            
    if MODEL_RAYLEIGH in WavesToModel and WavesToModel[MODEL_RAYLEIGH]:
        with open(os.path.join(OUTPUT, 'RayleighWaves.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDec output file for Rayleigh waves'])
            writer.writerow(['# MODEL_RAYLEIGH'])
            writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'Azimuth', 'EllipticityAngle'])
            writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]', '[rad]'])
            
    if MODEL_CIRCULAR_VERTICAL in WavesToModel and WavesToModel[MODEL_CIRCULAR_VERTICAL]:
        with open(os.path.join(OUTPUT, 'CircularVerticalWaves.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDecActive output file for circular vertical-component waves'])
            writer.writerow(['# MODEL_CIRCULAR_VERTICAL'])
            writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity'])
            writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]'])
            
    if MODEL_CIRCULAR_RAYLEIGH in WavesToModel and WavesToModel[MODEL_CIRCULAR_RAYLEIGH]:
        with open(os.path.join(OUTPUT, 'CircularRayleighWaves.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDecActive output file for Circular Rayleigh waves'])
            writer.writerow(['# MODEL_CIRCULAR_RAYLEIGH'])
            writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber', 'Velocity', 'EllipticityAngle'])
            writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]'])
    if MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH in WavesToModel and WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH]:
        with open(os.path.join(OUTPUT, 'CircularDissipativeRayleighWaves.csv'), 'w', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# WaveDecActive output file for Circular Dissipative Rayleigh waves'])
            writer.writerow(['# MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH'])
            writer.writerow(['# Frequency', 'Amplitude', 'Wavenumber_r', 'Velocity', 'EllipticityAngle', 'Wavenumber_i', 'Alpha'])
            writer.writerow(['# [Hz]', '[a.u.]', '[1/m]', '[m/s]', '[rad]', '[1/m]', '[rad/m]'])
#    if WavesToModel[MODEL_NOISE]:        
#        with open(os.path.join(OUTPUT, 'Noise.csv'), 'w',encoding='UTF8') as f:
#            writer = csv.writer(f, delimiter='\t')
#            writer.writerow(['# WaveDec output file'])
#            writer.writerow(['# ChannelId', 'Sigma2', 'Tstart', 'Tend'])
#            writer.writerow(['#', '[a.u.^2]', '[s]', '[s]'])
    return
    
def saveResults(conn, OUTPUT, WavesToModel, WindowId):
    cur = conn.cursor()
    if WavesToModel[MODEL_VERTICAL]:
        cur.execute('''SELECT Frequencies.F, VerticalWaves.Amplitude, VerticalWaves.Wavenumber, Frequencies.F/VerticalWaves.Wavenumber, VerticalWaves.Azimuth 
                FROM VerticalWaves JOIN Frequencies ON VerticalWaves.Fndx = Frequencies.Fndx
                WHERE VerticalWaves.WindowId = {0}'''.format(WindowId))
        data = cur.fetchall()
        with open(os.path.join(OUTPUT, 'VerticalWaves.csv'), 'a', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            for row in data:
                # TODO there may be a problem with division by zero, in case wavenumber=0, the velocity is None
                writer.writerow(['{:6.3e}'.format(c) for c in row])
    if WavesToModel[MODEL_LOVE]:  
        cur.execute('''SELECT Frequencies.F, LoveWaves.Amplitude, LoveWaves.Wavenumber, Frequencies.F/LoveWaves.Wavenumber, LoveWaves.Azimuth 
                FROM LoveWaves JOIN Frequencies ON LoveWaves.Fndx = Frequencies.Fndx
                WHERE LoveWaves.WindowId = {0}'''.format(WindowId))
        data = cur.fetchall()
        with open(os.path.join(OUTPUT, 'LoveWaves.csv'), 'a', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')            
            for row in data:
                writer.writerow(['{:6.3e}'.format(c) for c in row])
                
    if WavesToModel[MODEL_RAYLEIGH]:                    
        cur.execute('''SELECT Frequencies.F, RayleighWaves.Amplitude, RayleighWaves.Wavenumber, Frequencies.F/RayleighWaves.Wavenumber, RayleighWaves.Azimuth, RayleighWaves.EllipticityAngle
                FROM RayleighWaves JOIN Frequencies ON RayleighWaves.Fndx = Frequencies.Fndx
                WHERE RayleighWaves.WindowId = {0}'''.format(WindowId))
        data = cur.fetchall()
        with open(os.path.join(OUTPUT, 'RayleighWaves.csv'), 'a', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            for row in data:
                writer.writerow(['{:6.3e}'.format(c) for c in row])
                
    if WavesToModel[MODEL_CIRCULAR_VERTICAL]:
        cur.execute('''SELECT Frequencies.F, CircularVerticalWaves.Amplitude, CircularVerticalWaves.Wavenumber, Frequencies.F/CircularVerticalWaves.Wavenumber
                FROM CircularVerticalWaves JOIN Frequencies ON CircularVerticalWaves.Fndx = Frequencies.Fndx
                WHERE CircularVerticalWaves.WindowId = {0}'''.format(WindowId))
        data = cur.fetchall()
        with open(os.path.join(OUTPUT, 'CircularVerticalWaves.csv'), 'a', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            for row in data:
                # TODO there may be a problem with division by zero, in case wavenumber=0, the velocity is None
                writer.writerow(['{:6.3e}'.format(c) for c in row])

    if WavesToModel[MODEL_CIRCULAR_RAYLEIGH]:                    
        cur.execute('''SELECT Frequencies.F, CircularRayleighWaves.Amplitude, CircularRayleighWaves.Wavenumber, Frequencies.F/CircularRayleighWaves.Wavenumber, CircularRayleighWaves.EllipticityAngle
                FROM CircularRayleighWaves JOIN Frequencies ON CircularRayleighWaves.Fndx = Frequencies.Fndx
                WHERE CircularRayleighWaves.WindowId = {0}'''.format(WindowId))
        data = cur.fetchall()
        with open(os.path.join(OUTPUT, 'CircularRayleighWaves.csv'), 'a', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            for row in data:
                writer.writerow(['{:6.3e}'.format(c) for c in row])
                
    if WavesToModel[MODEL_CIRCULAR_DISSIPATIVE_RAYLEIGH]:                    
        cur.execute('''SELECT Frequencies.F, CircularDissipativeRayleighWaves.Amplitude, CircularDissipativeRayleighWaves.Wavenumber, Frequencies.F/CircularDissipativeRayleighWaves.Wavenumber, CircularDissipativeRayleighWaves.EllipticityAngle, CircularDissipativeRayleighWaves.Wavenumber_i, -6.28318530718*CircularDissipativeRayleighWaves.Wavenumber_i 
                FROM CircularDissipativeRayleighWaves JOIN Frequencies ON CircularDissipativeRayleighWaves.Fndx = Frequencies.Fndx
                WHERE CircularDissipativeRayleighWaves.WindowId = {0}'''.format(WindowId))
        data = cur.fetchall()
        with open(os.path.join(OUTPUT, 'CircularDissipativeRayleighWaves.csv'), 'a', encoding='UTF8') as f:
            writer = csv.writer(f, delimiter='\t')
            for row in data:
                writer.writerow(['{:6.3e}'.format(c) for c in row])               
#    if WavesToModel[MODEL_NOISE]:                    
#        cur.execute('''SELECT Noise.ChannelId, Noise.sigma2, Windows.Kstart*Windows.Ts AS Tstart, Windows.Kend*Windows.Ts AS Tend 
#                FROM Noise JOIN Windows ON Noise.WindowId = Windows.WindowId
#                WHERE Noise.WindowId = {0}'''.format(WindowId))
#        data = cur.fetchall()
#        with open(os.path.join(OUTPUT, 'Noise.csv'), 'a',encoding='UTF8') as f:
#            writer = csv.writer(f, delimiter='\t')
#            for row in data:
#                writer.writerow(['{:6.3e}'.format(c) for c in row])
    conn.commit()
    cur.close()
    return

def setArrayInfo(conn,array_info):
    

    Nchannels = np.shape(array_info)[0]
    cur = conn.cursor()
    
    for cc in range(0,Nchannels):
        ChannelId = cc
        pos_x = array_info[cc,0]
        pos_y = array_info[cc,1]
        pos_z = array_info[cc,2]
        component = array_info[cc,3]
        # Ts = array_info[cc,4]

        cur.execute('''SELECT count(*) FROM Sensors WHERE x={0} AND y={1} AND z={2}'''.format(pos_x, pos_y, pos_z))
        data = cur.fetchall()
        if data[0][0] == 0:
            cur.execute('''INSERT INTO Sensors (x, y, z) 
                VALUES ({0},{1},{2})'''.format(pos_x, pos_y, pos_z))
            SensorId = cur.lastrowid
        cur.execute('''INSERT INTO Channels (ChannelId, SensorId, type) 
                VALUES ({0},{1},{2})'''.format(ChannelId, SensorId, component))
            
    conn.commit()
    cur.close()
    return
    
def addWindow(conn, Kstart, Kend, Ts):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO Windows (Kstart, Kend, Ts) 
                VALUES ({0},{1},{2})'''.format(Kstart,Kend,Ts))
    WindowId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return WindowId
    
def addVerticalWave(conn, WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO VerticalWaves (WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth) 
                VALUES ({0},{1},{2},{3},{4},{5})'''.format(WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth))
    VerticalWaveId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return VerticalWaveId
    
def addCircularVerticalWave(conn, WindowId, Fndx, Amplitude, Phase, Wavenumber):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO CircularVerticalWaves (WindowId, Fndx, Amplitude, Phase, Wavenumber) 
                VALUES ({0},{1},{2},{3},{4})'''.format(WindowId, Fndx, Amplitude, Phase, Wavenumber))
    CircularVerticalWaveId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return CircularVerticalWaveId

def addLoveWave(conn, WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO LoveWaves (WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth) 
                VALUES ({0},{1},{2},{3},{4},{5})'''.format(WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth))
    LoveWaveId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return LoveWaveId
    
def addCircularLoveWave(conn, WindowId, Fndx, Amplitude, Phase, Wavenumber):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO CircularLoveWaves (WindowId, Fndx, Amplitude, Phase, Wavenumber) 
                VALUES ({0},{1},{2},{3},{4})'''.format(WindowId, Fndx, Amplitude, Phase, Wavenumber))
    CircularLoveWaveId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return CircularLoveWaveId

def addRayleighWave(conn, WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth, EllipticityAngle):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO RayleighWaves (WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth, EllipticityAngle) 
                VALUES ({0},{1},{2},{3},{4},{5},{6})'''.format(WindowId, Fndx, Amplitude, Phase, Wavenumber, Azimuth, EllipticityAngle))
    RayleighWaveId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return RayleighWaveId
    
def addCircularRayleighWave(conn, WindowId, Fndx, Amplitude, Phase, Wavenumber, EllipticityAngle):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO CircularRayleighWaves (WindowId, Fndx, Amplitude, Phase, Wavenumber, EllipticityAngle) 
                VALUES ({0},{1},{2},{3},{4},{5})'''.format(WindowId, Fndx, Amplitude, Phase, Wavenumber, EllipticityAngle))
    CircularRayleighWaveId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return CircularRayleighWaveId
    
def addCircularDissipativeRayleighWave(conn, WindowId, Fndx, Amplitude, Phase, Wavenumber, EllipticityAngle):
    
    cur = conn.cursor()
    cur.execute('''INSERT INTO CircularDissipativeRayleighWaves (WindowId, Fndx, Amplitude, Phase, Wavenumber, Wavenumber_i, EllipticityAngle) 
                VALUES ({0},{1},{2},{3},{4},{5},{6})'''.format(WindowId, Fndx, Amplitude, Phase, np.real(Wavenumber), np.imag(Wavenumber), EllipticityAngle))
    CircularDissipativeRayleighWaveId = cur.lastrowid
    conn.commit()
    cur.close()
    
    return CircularDissipativeRayleighWaveId
    
def addNoise(conn, WindowId, sigma2):
    
    L = len(sigma2)
    cur = conn.cursor()
    for ll in range(0, L):
        cur.execute('''INSERT INTO Noise (WindowId, ChannelId, sigma2) 
                VALUES ({0},{1},{2})'''.format(WindowId, ll, sigma2[ll]))
    conn.commit()
    cur.close()
    
    return
    
    
def str2bool(s):
    
    if type(s) is bool:
        return s

    if s.lower() in ['true', '1', 't', 'y', 'yes']:
        return True
    elif s.lower() in ['false', '0', 'f', 'n', 'no']:
        return False
    else:
        logging.critical("Cannot covert \"{}\" to a bool.".format(s))
        raise ValueError("Cannot covert \"{}\" to a bool.".format(s))