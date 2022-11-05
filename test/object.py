
import sys
import os
import numpy as np
import datetime as dt
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

cat1 = catalog.filter_datetime(starttime=dt.datetime(2022,1,1))
x = cat1.project((-80,10),(-70,10),(-501,501),verbose=False)

print(cat1)

seismicity = Seismicity([cat1,cat1])
x = seismicity.project((-80,10),(-70,10),(-501,501),verbose=False)

# filter_domain=[-76,-75,5,6]
# polygon = [(filter_domain[0],filter_domain[2]),
#                 (filter_domain[0],filter_domain[3]),
#                 (filter_domain[1],filter_domain[3]),
#                 (filter_domain[1],filter_domain[2]),
#                 (filter_domain[0],filter_domain[2])
#                 ]
# cat1 = cat1.filter_region(polygon)
# print(cat1)
# print(cat1.__str__(True))
# exit()
# print(cat1.__str__(True))

# a = Catalog()
# a = a.sort_values(by="origin_time")
# print(a)
# # cat1.color = "gray"
# fig = a.plot()
# # # fig = a.plot()
# fig.show()
# # print(cat1)
# # cat1.append(cat1.data)
# # print(cat1)
# exit()
# fig = cat1.matplot()
# plt.show()
# cat1.apply_cbar = False
# cat2 = catalog.filter_datetime(starttime=dt.datetime(2022,1,1))    
# # fig = cat2.plot()
# # fig.show()
# catalogs = Catalogs([cat1,cat2])  
# fig = catalogs.plot()
# fig.show()
# catalogs = catalogs.sort_values(by="origin_time")
# catalogs = catalogs.filter_datetime(endtime=dt.datetime(2022,5,1))
# catalogs = catalogs.filter_region(polygon)
# for catalog in catalogs:
#         print(catalog.__str__(True))



# # print(catalogs) 
# print(catalogs.__str__(True)) 
# fig = catalogs.plot(cbar=Cbar("depth","catalogs",cmap="oleron",
#                                 series=[0,200]))
# fig.show()
# catalogs 
# region = catalogs.get_region()           
# print(region)
# print(catalogs)
# print(cat2)

# fig = catalog.plot()
# fig.show()
# plt.show()
# print(catalog.__dict__)
# catalog.color = "red"
# print(catalog.__dict__)
# catalog = Catalogs([catalog,catalog,catalog])
# fig = catalog.mplot()
# fig.show()
# plt.show()

stations_path = "/home/emmanuel/EDCT/EQviewer/data/stations/castilla.csv"
station = pd.read_csv(stations_path)
station = Station(data = station,color="red")
sta1 = station.copy()
network = Network([station,sta1])
print(network)
station.remove({"station":["CA01","CA04"]})
print(network)
exit()
station2 = station.copy()

filter_domain=[-73.67,-73.61,3.79,3.87]
polygon = [(filter_domain[0],filter_domain[2]),
                (filter_domain[0],filter_domain[3]),
                (filter_domain[1],filter_domain[3]),
                (filter_domain[1],filter_domain[2]),
                (filter_domain[0],filter_domain[2])
                ]
station2 = station2.filter_region(polygon )
station2.color = "yellow"
network = Network([station,station2])
network.append(station2.copy())
print(network)
# print(networks)
# print(network.__str__(True))
# fig = station.plot()
# fig.show()
# fig = network.plot()
# fig.show()

# fig = station.matplot()
# plt.show()
# print(station.__str__(True))
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
    