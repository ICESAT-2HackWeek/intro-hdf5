# -*- coding: utf-8 -*-
"""
Utility functions for the ICESat2-HackWeek 2019.

This utils.py module is part of the following tutorials:

    - Intro to HDF5 and ICESat-2 files (by Fernando)
    - Gridding ICESat-2 data (by Johan)

Jet Propulsion Laboratory
California Institute of Technology

Fernando Paolo (paolofer@jplnasa.govl)
Johan Nilsson (johan.nilsson@jpl.nasa.gov)

"""
#import re
#import os
#import sys
#import h5py
#import pyproj
#import numpy as np
from astropy.time import Time
#from scipy.ndimage import map_coordinates

'''

try:
    from gdalconst import *
    from osgeo import gdal, osr
except:
    print('Proceeding without GDAL!')

# Hide anoying warnings
import warnings
warnings.filterwarnings('ignore')


def seg_diff_filter(dh_fit_dx, h_li, tol=2):
    """ Filter segments by consecutive difference. """
    dAT=20.
    
    if h_li.shape[0] < 3:
        mask = np.ones_like(h_li, dtype=bool)
        return mask

    EPplus  = h_li + dAT * dh_fit_dx
    EPminus = h_li - dAT * dh_fit_dx

    segDiff       = np.zeros_like(h_li)
    segDiff[0:-1] = np.abs(EPplus[0:-1] - h_li[1:])
    segDiff[1:]   = np.maximum(segDiff[1:], np.abs(h_li[0:-1] - EPminus[1:]))
    
    mask = segDiff < tol

    return mask

'''

def gps2dyr(time):
    """ Converte GPS time to decimal years. """
    return Time(time, format='gps').decimalyear

'''

def list_files(path, endswith='.h5'):
    """ List files in dir recursively."""
    return [os.path.join(dpath, f)
            for dpath, dnames, fnames in os.walk(path)
            for f in fnames if f.endswith(endswith)]


def transform_coord(proj1, proj2, x, y):
    """
    Transform coordinates from proj1 to proj2 (EPSG num).

    Example EPSG projs:
        Geodetic (lon/lat): 4326
        Polar Stereo AnIS (x/y): 3031
        Polar Stereo GrIS (x/y): 3413
    """
    # Set full EPSG projection strings
    proj1 = pyproj.Proj("+init=EPSG:"+str(proj1))
    proj2 = pyproj.Proj("+init=EPSG:"+str(proj2))
    return pyproj.transform(proj1, proj2, x, y)  # convert

'''

def track_type(time, lat, tmax=1):
    """
    Separate tracks into ascending and descending.

    Defines tracks as segments with time breaks > tmax,
    and tests whether lat increases or decreases w/time.
    """
    tracks = np.zeros(lat.shape)  # generate track segment
    tracks[0:np.argmax(np.abs(lat))] = 1  # set values for segment
    i_asc = np.zeros(tracks.shape, dtype=bool)  # output index array
    
    # Loop trough individual secments
    for track in np.unique(tracks):
        
        i_track, = np.where(track == tracks)  # get all pts from seg 
        
        if len(i_track) < 2: continue
        
        # Test if lat increases (asc) or decreases (des) w/time
        i_min = time[i_track].argmin()
        i_max = time[i_track].argmax()
        lat_diff = lat[i_track][i_max] - lat[i_track][i_min]
        
        # Determine track type
        if lat_diff > 0:  i_asc[i_track] = True

    return i_asc, np.invert(i_asc)  # index vectors

'''

def h5read(ifile, vnames):
    """ Simple HDF5 reader. """
    with h5py.File(ifile, 'r') as f:
        return [f[v][:] for v in vnames]


def ncread(ifile, vnames):
    """ Simple NetCDF4 reader. """
    ds = Dataset(ifile, "r")   # NetCDF4
    d = ds.variables
    return [d[v][:] for v in vnames]


def tifread(ifile, metaData):
    """ Simple GeoTIFF reader. """
    file = gdal.Open(ifile, GA_ReadOnly)
    projection = file.GetProjection()
    src = osr.SpatialReference()
    src.ImportFromWkt(projection)
    proj = src.ExportToWkt()
    Nx = file.RasterXSize
    Ny = file.RasterYSize
    trans = file.GetGeoTransform()
    dx = trans[1]
    dy = trans[5]
    if metaData == "A":
        xp = np.arange(Nx)
        yp = np.arange(Ny)
        (Xp, Yp) = np.meshgrid(xp,yp)
        X = trans[0] + (Xp+0.5)*trans[1] + (Yp+0.5)*trans[2]  #FIXME: bottleneck!
        Y = trans[3] + (Xp+0.5)*trans[4] + (Yp+0.5)*trans[5]
    if metaData == "P":
        xp = np.arange(Nx)
        yp = np.arange(Ny)
        (Xp, Yp) = np.meshgrid(xp,yp)
        X = trans[0] + Xp*trans[1] + Yp*trans[2]  #FIXME: bottleneck!
        Y = trans[3] + Xp*trans[4] + Yp*trans[5]
    band = file.GetRasterBand(1)
    Z = band.ReadAsArray()
    dx = np.abs(dx)
    dy = np.abs(dy)
    #return X, Y, Z, dx, dy, proj
    return X, Y, Z


def interp2d(xd, yd, data, xq, yq, **kwargs):
    """ Bilinear interpolation from grid. """
    
    xd = np.flipud(xd)
    yd = np.flipud(yd)
    data = np.flipud(data)
    xd = xd[0,:]
    yd = yd[:,0]
    
    nx, ny = xd.size, yd.size
    (x_step, y_step) = (xd[1]-xd[0]), (yd[1]-yd[0])
    
    assert (ny, nx) == data.shape
    assert (xd[-1] > xd[0]) and (yd[-1] > yd[0])
    
    if np.size(xq) == 1 and np.size(yq) > 1:
        xq = xq*ones(yq.size)
    elif np.size(yq) == 1 and np.size(xq) > 1:
        yq = yq*ones(xq.size)
    
    xp = (xq-xd[0])*(nx-1)/(xd[-1]-xd[0])
    yp = (yq-yd[0])*(ny-1)/(yd[-1]-yd[0])
    
    coord = np.vstack([yp,xp])
    zq = map_coordinates(data, coord, **kwargs)
    
    return zq

'''
