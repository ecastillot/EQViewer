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
from EQViewer import Catalog,Station,Shape
from EQViewer import seismic_profile,xy_seismic_profile
import EQViewer.utils as equt


lats = -5;latn = 15;lonw = -85;lone = -65 #castilla
cat_csv = os.path.join(rep_data,"earthquakes","events_20160101T20220901_NLLOC_ok.csv")
# cat_csv = os.path.join(rep_data,"earthquakes","CM_20160101T20220901_ok.csv")
tecto_shp = os.path.join(rep_data,"shapes","Tectónica.shp")
faults_shp = os.path.join(rep_data,"shapes","Fallas.shp")


df = pd.read_csv(cat_csv)
df = df.drop_duplicates(subset="id",ignore_index=True)
events = equt.transform_to_fmt_catalog(cat_csv,
        columns={"time_event":"origin_time"})
events = events.sort_values(by="depth",ascending=False)
events.loc[(events["depth"]<0),"depth"] = 0
events.loc[(events["depth"]>200),"depth"] = 200
events_sup = events[events["depth"]<=40]
events_dep= events[events["depth"]>40]


reg = [lonw , lone, lats, latn ]

catalogs = [
            Catalog(data=events_dep, color = "lightblue",
                    style="c0.07c",
                    size=None,
                    cmap = True,
                    pen = None,
            ),
        Catalog(data=events_sup, color = "lightblue",
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

all_data = gpd.read_file(faults_shp)
shape_f = Shape(data=all_data,color="white", pen=["0.01c,black,-"],connection="rr")

shapes_before_catalog = [shape_r,shape_l]
shapes_after_catalog = [shape_f]

fig,new_catalog = seismic_profile.map(reg,catalogs=catalogs,
                shapes_before_catalog=shapes_before_catalog,
                # shapes_after_catalog = shapes_after_catalog,
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
outpath = os.path.join(rep_out,"cm_sup.png")
fig.savefig(outpath,dpi=300)