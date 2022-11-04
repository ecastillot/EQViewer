
import sys
import os
import numpy as np
repository_path = r"/home/emmanuel/EDCT/EQviewer"  ##change this path where is located the main directory
rep_data = os.path.join(repository_path,"data")
rep_out = os.path.join(repository_path,"example")
sys.path.insert(0,repository_path)
import geopandas as gpd
from EQViewer.objects import *
import EQViewer.utils as equt
import matplotlib.pyplot as plt
lats = -5;latn = 15;lonw = -85;lone = -65 #castilla
cat_csv = os.path.join(rep_data,"earthquakes","events_20160101T20220901_NLLOC_ok.csv")
tecto_shp = os.path.join(rep_data,"shapes","Tectónica.shp")


df = pd.read_csv(cat_csv)
df = df.drop_duplicates(subset="id",ignore_index=True)
events = equt.transform_to_fmt_catalog(cat_csv,
        columns={"time_event":"origin_time"})
reg = [lonw , lone, lats, latn ]

catalog = Catalog(data=events, 
                    color = "lightblue",
                    style="cc",
                    size=lambda x: 0.11 * np.sqrt(1.2 ** (x*1.4)),
                    # style="c0.1c",
                    # size=None,
                    apply_cbar = True,
                    pen = None)
# fig = catalog.plot()
# fig = catalog.mplot()
# fig.show()
# plt.show()

# stations_path = "/home/emmanuel/EDCT/EQviewer/data/stations/castilla.csv"
# stations = pd.read_csv(stations_path)
# stations = Station(data = stations)
# stations.plot()
# plt.show()
# print(stations)

# tecto_shp = os.path.join(rep_data,"shapes","Tectónica.shp")
# all_data = gpd.read_file(tecto_shp)
# all_data_l = all_data[all_data["NombreEstr"] == "Falla de Jordan"]
# all_data_r = all_data[all_data["NombreEstr"] != "Falla de Jordan"]
# shape_r = Shape(data=all_data_r,color="white", pen=["0.02c,black,-"],
#         connection="rr",
#         style="f1c/0.1c+r+t")
# print(shape_r.data.area)
# plt.show()


# # data_path = os.path.join(rep_data,"well","well_example.csv")
# data_path = "/home/emmanuel/EDCT/EQviewer/data/well/survey_proc2.xlsx"
# data = pd.read_excel(data_path)
# data = data.rename(columns={"MD (ft)":"MD",
#                             "Lon (°)":"longitude",
#                             "Lat (°)":"latitude",
#                             "Z (m)": "z",
#                             "TVD (ft)":"TVD"})
# # print(data)
# injection = pd.read_csv("/home/emmanuel/EDCT/EQviewer/data/well/CASTILLA_ft_ini.dat")
# injection = injection[injection["name"]=="CAD01"]
# injection = injection[injection["name"]=="CAD01"]
# print(injection.info())
# # exit()

# well = Well(data,"PAD",injection=injection,cbar=True)
# well.plot()
# plt.show()
# # wells = Wells(wells=[well],
# #                 injection_cmap=None)
    