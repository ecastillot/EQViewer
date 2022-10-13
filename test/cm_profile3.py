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
from EQViewer import Catalog,Station,Shape,Profile
from EQViewer import seismic_profile,xy_seismic_profile,utils
import EQViewer.utils as equt

# 407580
# 292766

lats = -5;latn = 15;lonw = -85;lone = -65 #castilla
# cat_csv = os.path.join(rep_data,"earthquakes","CM_20160101T20220901_ok.csv")
cat_csv = os.path.join(rep_data,"earthquakes","events_20160101T20220901_NLLOC_ok.csv")
# cat_csv = os.path.join("/media/emmanuel/TOSHIBA EXT/ColSeismicity/events_20160101T20220901.csv")
tecto_shp = os.path.join(rep_data,"shapes","Tectónica.shp")


df = pd.read_csv(cat_csv)
df = df.drop_duplicates(subset="id",ignore_index=True)
events = equt.transform_to_fmt_catalog(cat_csv,
        columns={"time_event":"origin_time"})
events = events.sort_values(by="depth",ascending=True)
events.loc[(events["depth"]<0),"depth"] = 0
events.loc[(events["depth"]>200),"depth"] = 200

reg = [lonw , lone, lats, latn ]

catalogs = [Catalog(data=events, color = "lightblue",
                    style="c0.1c",
                    size=None,
                    cmap = True,
                    pen = "black"
            )
            ]


all_data = gpd.read_file(tecto_shp)
all_data_l = all_data[all_data["NombreEstr"] == "Falla de Jordan"]
all_data_r = all_data[all_data["NombreEstr"] != "Falla de Jordan"]
shape_r = Shape(data=all_data_r,color="white", pen=["0.02c,black,-"],
        connection="rr",
        style="f1c/0.1c+r+t")
shape_l = Shape(data=all_data_l,color="white", pen=["0.02c,black,-"],
        connection="rr",
        style="f1c/0.1c+l+t")

shapes = [shape_r,shape_l]


profile1 = Profile(name=("B","B'"),      
        # coords=((-74.868,11.356),(-70.272,8.588)), 
        coords=((-78.517,5.354),(-70.352,3.966)), 
        width=(-150,150), #left and right width
        colorline="magenta", #only for map figure
        color= "blue", # only for profile figure.
                        #color of the events in the profile. Only if cmap is False
        # size=lambda x= 0.1 * np.sqrt(1.5 ** (x*1.5)), # only for profile figure.
        grid=[50,50],
        size=None, # only for profile figure.
        style ="c0.1c",# only for profile figure.
        pen="black", # only for profile figure. 
        cmap=True, # only for profile figure. # cmap controlled by cbar_profile_args
        cbar_profile_args = {"cmap":'viridis', 
                                "color_target":"depth",
                                "label":"Depth (km)",
                                "reverse":False,
                                "series":[0, 200] } # only for profile figure. # cbar profile args
        )
profile2 = Profile(name=("A","A'"),      
        coords=((-75.583,11.494),(-68.714,8.799)), 
        width=(-150,150), #left and right width
        colorline="magenta", #only for map figure
        color= "blue", # only for profile figure.
                        #color of the events in the profile. Only if cmap is False
        # size=lambda x= 0.1 * np.sqrt(1.5 ** (x*1.5)), # only for profile figure.
        grid=[50,50],
        size=None, # only for profile figure.
        style ="c0.1c",# only for profile figure.
        pen="black", # only for profile figure. 
        cmap=True, # only for profile figure. # cmap controlled by cbar_profile_args
        cbar_profile_args = {"cmap":'viridis', 
                                "color_target":"depth",
                                "label":"Depth (km)",
                                "reverse":False,
                                "series":[0, 200] } # only for profile figure. # cbar profile args
        )

profile3 = Profile(name=("C","C'"),      
        coords=((-77.075,5.403),(-72.161,4.526)), 
        width=(-150,150), #left and right width
        colorline="magenta", #only for map figure
        color= "blue", # only for profile figure.
                        #color of the events in the profile. Only if cmap is False
        # size=lambda x= 0.1 * np.sqrt(1.5 ** (x*1.5)), # only for profile figure.
        grid=[50,50],
        size=None, # only for profile figure.
        style ="c0.1c",# only for profile figure.
        pen="black", # only for profile figure. 
        cmap=True, # only for profile figure. # cmap controlled by cbar_profile_args
        cbar_profile_args = {"cmap":'viridis', 
                                "color_target":"depth",
                                "label":"Depth (km)",
                                "reverse":False,
                                "series":[0, 200] } # only for profile figure. # cbar profile args
        )

profiles = [profile1,profile2]
# profiles = [profile1,profile2,profile3]


fig,new_catalog = seismic_profile.map(reg,catalogs=catalogs,
                shapes_before_catalog=shapes,
                profiles=profiles,
                cmap_args = {"cmap":'rainbow', #name of the colorbar. see gmt colorbars
                            "color_target":"depth", #Here, colorbar refers to the depth column of the catalog. 
                                                    #You can use origin_time column as well
                            "label":"Depth (km)", #label of the colorbar
                            "overrule_bg":True, #overrule True to color white and black for outliers.
                            "reverse":True, # reverse the colorbar
                            "series":[0, 200,25] }, #remove 200 if you want continous colorbar.
                                                    # [] empty to select the minimum and maximum values of the catalog
                map_scale_args = {"xloc":-70, #xloc,yloc are the x-y map coordinates to locate the scale bar. 
                                "yloc":-5.1, #if they are None, put bottom left as default.
                                "distance":500}, #distance en km
                rivers=None,
                legend=False
                )

proj = 'G-70/0/1i'
fig.shift_origin(xshift="0i",yshift="0i")
fig.grdimage(
            '@earth_relief_10m',
            region='g',
            projection=proj,
            cmap='globe',
            shading=True,
          )
fig.coast(
            region='g',
            projection=proj,
            shorelines=True,
            water='white',
            borders='1/1p,black',
            land='grey',
            frame=True,
        )
x_reg = [reg[0], reg[1], reg[1], reg[0], reg[0]]
y_reg = [reg[2], reg[2], reg[3], reg[3], reg[2]]
fig.plot(x=x_reg,y=y_reg,
        pen="2p,red")



fig.show()
outpath = os.path.join(rep_out,"col_map_lon.pdf")
fig.savefig(outpath,dpi=300)

prof_fig = seismic_profile.profile_plots(region=reg,
                catalogs=catalogs,
                profiles= profiles,
                depth=[0,200],
                subsize = ("14c", "14c"),
                depth_unit="km",
                save=True)
prof_fig.show()
outpath = os.path.join(rep_out,"col_prof_lon.pdf")
prof_fig.savefig(outpath,dpi=300)
