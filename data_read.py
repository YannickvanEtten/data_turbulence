'''
Filename: question.py

Purpose:
    
Date:
    15 August 2023
Author:
    Yannick van Etten 2688877  
'''
###########################################################
### Imports
import numpy.random as rnd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts
import scipy.integrate as ing
from scipy.integrate import quad
import pandas as pd
import cdsapi
import xarray as xr
import pygrib
c = cdsapi.Client()
###########################################################

###########################################################

def main():
    grib = "C://Users//yanni//OneDrive//Documenten//Universiteit//Dataproject//ERA5_global_raw_[date].grib" 
    #ds = xr.open_dataset(file_path, engine='cfgrib')
    #print(ds)
    #grib = 'cams_aod.grib' # Set the file name of your input GRIB file
    grbs = pygrib.open(grib)
    grbs = pygrib.open(grib)
    grb = grbs.select()[0]
    data = grb.values
    #print('test')

###########################################################
### call main
if __name__ == "__main__":
    main()