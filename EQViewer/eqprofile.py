import pygmt
import types
import random
from scipy import interpolate
import os
import string
import pandas as pd
import numpy as np
from . import utils as ut
import datetime as dt
import geopandas as gpd
from obspy.geodetics.base import gps2dist_azimuth
pygmt.config(FORMAT_GEO_MAP="ddd.xx")

def map(region,
        catalogs=None,
        stations=None,
        profiles=None,
        fms = None,
        wells=None,
        shapes_before_catalog=None,
        shapes_after_catalog=None,
        fig = None,
        coast_dict = {"projection":"M?",
                    "shorelines":"True",
                    "land":"gray",
                    "water":'lightblue',
                    "rivers":"['2/blue','2/blue']",
                    "borders":'1/1p,black',
                    "frame":["afg","WNse"]},
        grdimage_dict={},
        map_scale_loc = {"xloc":None,
                        "yloc":None,
                        "distance":5},
        legend_loc = {"xloc":None, 
                        "yloc":None},
        legend=True
        ):
    """
    Parameters:
    -----------
    region: list
        [lonw , lone, lats, latn]
    catalogs: Catalogs
        Seismological catalogs
    stations: Stations
        Seismological stations
    profiles: Profiles
        Seismic Profiles
    fms: FocalMechanisms
        Seismic fochal mechanisms
    wells: Wells
        Wells trajectory
    shapes_before_catalog: Shapes
        Shapes before to plot the catalogs
    shapes_after_catalog: Shapes
        Shapes after to plot the catalogs
    fig : pygmt.Figure
        Previous pygmt figure
    coast_dict : dict
        See pygmt.coast to know the items  of the dictionary.
    grdimage_dict : dict
        See pygmt.grdimage to know the items  of the dictionary.
    map_scale_loc: dict
        keys values       represents
        ----------------------------
        xloc  float       longitude in degrees
        yloc  float       latitude in degrees
        distance float    distance in km  
    legend_loc: dict
        keys values       represents
        ----------------------------
        xloc  float       longitude in degrees
        yloc  float       latitude in degrees
    legend: bool
        True to show the legend
    """

    if fig == None:
        fig = pygmt.Figure()  
    
    if grdimage_dict:
        fig.grdimage(**grdimage_dict)
 
    coast_dict.pop('region', None)
    fig.coast(region=region,
            **coast_dict)

    if shapes_before_catalog != None:
        