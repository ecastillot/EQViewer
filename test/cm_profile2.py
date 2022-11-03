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


profile1 = Profile(name=("A","A'"),      
        coords=((-80,10),(-70,10)), 
        width=(-300,300), #left and right width
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
# profiles = [profile1,profile1]
profiles_tmp = profile1
profile1 = utils.get_cross_profiles(((-74.0875,5.2868),
                                     (-73.3667,6.7886)),
                                7,90,profiles =profiles_tmp)
profile2 = utils.get_cross_profiles(((-73.3138,6.7503),
                                     (-73.2623,6.9633)),
                                2,90,number_of_profile=len(profile1),
                                profiles =profiles_tmp)
profile3 = utils.get_cross_profiles(((-73.2651,7.3209),
                                     (-73.1207,8.8896)),
                                7,90,
                                number_of_profile=len(profile1)+len(profile2),
                                profiles =profiles_tmp)

profiles = profile1+profile2+profile3

w = 0.8
points = [x.coords for x in profiles ]
reg = np.array(points)
reg = reg.reshape(reg.shape[0]*reg.shape[1],-1)
latm = min(reg[:,1])
latM = max(reg[:,1])
lonm = min(reg[:,0])
lonM = max(reg[:,0])
reg = [lonm-w, lonM+w , latm-w, latM+w]


# fig = pygmt.Figure()


# fig.show()


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
                map_scale_args = {"xloc":-71.4, #xloc,yloc are the x-y map coordinates to locate the scale bar. 
                                "yloc":3.5, #if they are None, put bottom left as default.
                                "distance":200}, #distance en km
                rivers=None,
                legend=False,
                )
fig.shift_origin(xshift="0.08c",yshift="0.08c")

lats = -15;latn = 35;lonw = -100;lone = -60 #castilla
col_reg = [lonw , lone, lats, latn ]

with pygmt.config(MAP_FRAME_TYPE="plain"):
        fig.coast(
                region=col_reg,
                projection='M1i',
                shorelines=True,
                water='white',
                borders='1/1p,black',
                land='grey',
                    frame="f",
                )

with pygmt.config(FONT_TITLE=4):
    fig.basemap(rose="jBL+w1c+lO,E,S,N+o-0.1c/0.5c+p,red", 
        # map_scale="jBL+w500k+o0.5c/0.5c+f+lkm+at"
        ) 
x_reg = [reg[0], reg[1], reg[1], reg[0], reg[0]]
y_reg = [reg[2], reg[2], reg[3], reg[3], reg[2]]
fig.plot(x=x_reg,y=y_reg,
        pen="2p,red")




fig.show()
outpath = os.path.join(rep_out,"vmm_map.pdf")
fig.savefig(outpath,dpi=300)

prof_fig = seismic_profile.profile_plots(region=reg,
                catalogs=catalogs,
                profiles= profiles,
                depth=[0,200],
                subsize = ("14c", "14c"),
                depth_unit="km")
prof_fig.show()
outpath = os.path.join(rep_out,"vmm_prof.pdf")
prof_fig.savefig(outpath,dpi=300)
