# conda forge pygmt obspy geopandas openpyxl

import sys
import os
repository_path = r"/home/emmanuel/EDCT/EQviewer"  ##change this path where is located the main directory
rep_data = os.path.join(repository_path,"data")
rep_out = os.path.join(repository_path,"example")
sys.path.insert(0,repository_path)

import pygmt
import numpy as np
import pandas as pd
import datetime as dt
import geopandas as gpd
from EQViewer import Catalog,Station,Shape,Well,Wells,Cbar
from EQViewer import seismic_profile,xy_seismic_profile
import EQViewer.utils as equt

# data_path = os.path.join(rep_data,"well","well_example.csv")
data_path = "/home/emmanuel/EDCT/EQviewer/data/well/survey_proc2.xlsx"
data = pd.read_excel(data_path)
data = data.rename(columns={"MD (ft)":"MD",
                            "Lon (°)":"lon",
                            "Lat (°)":"lat",
                            "Z (m)": "z",
                            "TVD (ft)":"TVD"})
# print(data)
injection = pd.read_csv("/home/emmanuel/EDCT/EQviewer/data/well/CASTILLA_ft_ini.dat")
injection = injection[injection["name"]=="CAD01"]
injection = injection[injection["name"]=="CAD01"]
print(injection.info())
# exit()

well = Well(data,"PAD",injection=injection,cmap=True)

wells = Wells(wells=[well],
                injection_cmap=None)
                # injection_cmap=Cbar(
                #                 color_target="water_flow",
                #                 cmap='cool',
                #                 transparency=80,
                #                 label=["bbl/day"],
                #                 overrule_bg=True,
                #                 reverse=True,
                #                 series=[0,1e5,1e4])
                # )

lats = 3.78;latn = 3.90;lonw = -73.72;lone = -73.60 #castilla
reg = [lonw , lone, lats, latn ]
fig,new_catalog = seismic_profile.map(reg,catalogs=None,
                shapes_before_catalog=None,
                wells=wells,
                cmap_args = Cbar(cmap='rainbow', #name of the colorbar. see gmt colorbars
                            transparency=0,
                            color_target="depth", #Here, colorbar refers to the depth column of the catalog. 
                                                    #You can use origin_time column as well
                            label="Depth (m)", #label of the colorbar
                            overrule_bg=True, #overrule True to color white and black for outliers.
                            reverse=True, # reverse the colorbar
                            series=[0, 3e3,250] ), #remove 200 if you want continous colorbar.
                                                    # [] empty to select the minimum and maximum values of the catalog
                map_scale_args = {"xloc":-70, #xloc,yloc are the x-y map coordinates to locate the scale bar. 
                                "yloc":-5.1, #if they are None, put bottom left as default.
                                "distance":500}, #distance en km
                rivers=None,
                legend=False
                )
fig.show()