import os, sys
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt



def run(N, state='Dicke', inte='RK4', pol=0.1, hthr=0.3, gamma=0.8, ti=1e-3, tf=10.0, meas=50):
    '''
    '''
    main_dir = os.getcwd()
    exe = main_dir + "/OpenNu"

    options = ["--N", str(N), "--thr", str(hthr), "--sthr", "0.00001",
               "--ti", str(ti), "--tf", str(tf), "--meas", str(meas),
               "--gamma", str(gamma), "--state", state, "--mtime", '1',
               "--integrator", inte, "--pol", str(pol)]
    command = [exe] + options
    subprocess.run(command)

def driverN():
    for n in [20, 40, 60, 80, 100]:
        run(n)

def driverStep():
    N = 100
    for thr in [0.01, 0.05, 0.1, 0.5, 1]:
        run(N, hthr=thr)

def driverGamma():
    N = 100
    for g in np.linspace(0.95,0.99,5):
        run(N, gamma=g, state='Ground', inte='Euler', tf=120.0, hthr=0.4, meas=100)
        #run(N, gamma=g, state='Ground', tf=120.0, hthr=0.1, meas=100)

def driverSxN():
    for n in [10, 20, 30, 40, 50, 60, 70]:
        run(n, state='Ground', tf=100.0, gamma=0.8, hthr=0.1)

def driverPol():

    li0 = np.array([1e5, 1e-4, 1e-3, 1e-2, 1e-1])
    li1 = np.linspace(0.1, 0.9, 9)
    li2 = np.linspace(0.91, 0.99, 9)

    pols = np.concatenate([li0, li1, li2])
    #pols = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.9]
    #pols = [0.5]

    for n in [100,]:
#        for g in np.arange(0.95, 0.999, 0.001):
        for p in pols:
            run(n, state='Thermal', inte='Euler', pol=p, tf=100.0, gamma=0.99, hthr=0.5)



if __name__ == "__main__":
    #driverN()
    #driverStep()
    #driverGamma()
    #driverSxN()
    driverPol()
