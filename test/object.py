
import sys
import os
import numpy as np
import datetime as dt
repository_path = r"/home/emmanuel/EQviewer"  ##change this path where is located the main directory
rep_data = os.path.join(repository_path,"data")
rep_out = os.path.join(repository_path,"example")
sys.path.insert(0,repository_path)
import geopandas as gpd
from EQViewer.new_objects import *
import EQViewer.utils as equt
import matplotlib.pyplot as plt
import datetime as t
import random




# # data_path = os.path.join(rep_data,"well","well_example.csv")
# data_path = "/home/emmanuel/EQviewer/data/well/survey_proc2.xlsx"
# data = pd.read_excel(data_path)
# data = data.rename(columns={"MD (ft)":"MD",
#                             "Lon (°)":"longitude",
#                             "Lat (°)":"latitude",
#                             "Z (m)": "depth",
#                             "TVD (ft)":"TVD"})
# # print(data)
# injection = pd.read_csv("/home/emmanuel/EQviewer/data/well/CASTILLA_ft_ini.dat")
# injection = injection.rename(columns={"min":"min_depth",
#                                     "max":"max_depth",
#                                     "water_flow":"measurement"})
# injection = injection[injection["name"]=="CAD01"]
# injection = injection[injection["name"]=="CAD01"]
# # print(injection.info())
# # exit()
# print(data)
# print(injection)
# injection = Injection(injection,depth_type="MD",
#                     # baseplot=BasePlot(cmap=True,
#                     #                 style="g0.3")
#                                     )

# zmin = data.depth.min()
#                 zmax = data.depth.max()
#                 survey_cpt = CPT(color_target="depth",
#                             label="depth",
#                             cmap="rainbow",
#                             series=[zmin,zmax],
#                             reverse=True,
#                             overrule_bg=True)
# well = Well(data,"PAD",
#             survey_baseplot = BasePlot(
#                         # size=None,
#                         # style="g0.05",
#                         cmap=False,
#                         color=None,
#                         label=None,
#                         transparency=None,
#                         pen="1p"
#                         ),
#             injection=injection)
# # print(well)
# # print(well.data)
# fig = well.plot()
# fig.show()
# # plt.show()
# exit()
# wells = Wells(wells=[well],
#                 injection_cmap=None)


lats = -5;latn = 15;lonw = -85;lone = -65 #castilla
cat_csv = os.path.join(rep_data,"earthquakes","events_20160101T20220901_NLLOC_ok.csv")
# tecto_shp = os.path.join(rep_data,"shapes","Tectónica.shp")





# fm_csv = os.path.join(rep_data,"mf","MF_quifa_201301_202203.csv")
# events = pd.read_csv(fm_csv)
# events = events.rename(columns={"Origin time":"origin_time",
#                         "Latitude (deg)":"latitude",
#                         "Longitude (deg)":"longitude",
#                         "Depth (m)":"depth",
#                         "Mag. (Mw)":"magnitude",
#                         "Dip n1 (deg)":"dip",
#                         "Rake n1 (deg)":"rake",
#                         "Strike n1 (deg)":"strike",
#                         "Dip n2 (deg)":"dip_n2",
#                         "Rake n2 (deg)":"rake_n2",
#                         "Strike n2 (deg)":"strike_n2",
#                         }
#                         )
# events["depth"]= events["depth"]/1e3
# ev2 = events.copy()
# fm = FM(events,BaseMeca(cmap=False))
# region = fm.get_region()
# events["plot_longitude"] = pd.Series(np.random.uniform(region[0], region[1], size=len(events)))
# events["plot_latitude"] = pd.Series(np.random.uniform(region[2], region[3], size=len(events)))
# fm1 = FM(events,BaseMeca(cmap=True))
# # print(events)
# fm2 = FM(ev2.iloc[0:6],BaseMeca(cmap=False))
# fm = FMs([fm1,fm2])
# fig = fm.plot()
# fig.show()
# exit()





df = pd.read_csv(cat_csv)
df = df.drop_duplicates(subset="id",ignore_index=True)
events = equt.transform_to_fmt_catalog(cat_csv,
        columns={"time_event":"origin_time"})

df1 = events[events["origin_time"]<=dt.datetime(2022,1,1)]
df2 = events[events["origin_time"]>=dt.datetime(2022,1,1)]
reg = [lonw , lone, lats, latn ]
baseplot2 = BasePlot(color = "lightblue",
                    style="cc",
                    size=lambda x: 0.11 * np.sqrt(1.2 ** (x.magnitude*1.4)),
                #     style="c0.1c",
                #     size=None,
                    cmap = True,
                    pen = "black")
baseplot1 = BasePlot(color = "gray",
                    style="cc",
                    size=lambda x: 0.11 * np.sqrt(1.2 ** (x.magnitude*1.4)),
                #     style="c0.1c",
                #     size=None,
                    cmap = False,
                    pen = "black")
catalog = Catalog(df1,baseplot1 )
catalog2 = Catalog(df2,baseplot2 )
profile = Profile(name=("A","A'"),      
        coords=((-80,10),(-70,10)), 
        width=(-300,300)
            )
baseprofile = BaseProfile(projection="X10/-10",
                        depth_lims=[0,200],output_unit="km")

# fig = catalog.plot_profile(profile,baseprofile,
#                             depth_unit="km",
#                             verbose=False)
# fig = catalog2.plot_profile(profile,baseprofile,
#                             depth_unit="km",fig=fig,
#                             verbose=False)
# fig.show()

multicatalog = MulCatalog(catalogs=[catalog,catalog2],
                            show_cpt=True)
# print(multicatalog.__str__(True))
fig = multicatalog.plot_profile(profile,baseprofile,
                            depth_unit="km",
                            )
fig.show()
exit()
# mulcatalog = MulCatalog([catalog])
# fig = pygmt.Figure()
# with fig.subplot(nrows=2, ncols=3, 
#             subsize = ("12c", "12c"),
#            margins=["3c", "-3.5c"],frame="lrtb"):
#     for i in range(2):  # row number starting from 0
#         for j in range(3):  # column number starting from 0
#             index = i * 3 + j  # index number starting from 0
#             with fig.set_panel(panel=index):  # sets the current panel
#                 # fig.text(
#                 #     position="MC",
#                 #     text=f"index: {index}; row: {i}, col: {j}",
#                 #     region=[0, 1, 0, 1],
#                 # )
#                 # p_o = aa.mul_objects[0]["p(o)"]
#                 # fig = p_o.plot(fig=fig)
#                 fig = catalog.plot_profile(profile,baseprofile,depth_unit="km",fig=fig,)
#                 break
# fig.show()
# exit()               
aa = Profile(name=("A","A'"),      
        coords=((-80,10),(-70,10)), 
        width=(-300,300)
            )
aa.add(catalog,depth_unit="km")



p_o = aa.objects[0]["p(o)"]
fig = p_o.plot()
fig.show()
# print(aa.objects)
exit()
# cat1 = catalog.copy()
# cat1 = cat1.filter_datetime(starttime=dt.datetime(2022,1,1))
# cat2 = catalog.copy()
# baseplot = BasePlot(color = "gray",
#                     style="c0.1c",
#                     size=None,
#                     cmap = False,
#                     pen = "black")
# cat2.baseplot = baseplot
# cats = Catalogs([cat2,cat1])
# fig = cats.plot()
# fig.show()
# cat1.matplot()
# plt.show()
# exit()

# print(cat2.info2pygmt)
# fig = cat2.plot()
# # fig.show()
# # exit()
# # x = cat1.project((-80,10),(-70,10),(-501,501),verbose=False)
# # fig = cat1.plot()
# # fig.show()

# fig = s.plot()
# fig.show()
# print(cat1.info2pygmt)
# exit()
# print(cat1)

# seismicity = Seismicity([cat1,cat1])
# x = seismicity.project((-80,10),(-70,10),(-501,501),verbose=False)

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

# stations_path = "/home/emmanuel/EQviewer/data/stations/castilla.csv"
# station = pd.read_csv(stations_path)
# baseplot = BasePlot(color="black",
#                 label="stations",
#                 transparency = 0,
#                 style="i0.3c",
#                 pen="black")
# basetext = BaseText(font="10p,Helvetica,black",
#                     fill=None,
#                     offset="-0.05c/0.15c")

# station = Station(data = station)
# sta2 = station.copy()
# sta2.baseplot.color = "red"
# print(station.baseplot.__dict__)
# print(sta2.baseplot.__dict__)
# station = station.remove_data({"station":["CA03"]})
# stations = Stations([sta2,station])
# fig = stations.plot()
# fig.show()
# # station.matplot()
# # plt.show()
# exit()
# sta1 = station.copy()
# network = Network([station,sta1])
# fig = network.plot()
# fig.show()
# print(network)
# station.remove({"station":["CA01","CA04"]})
# print(network)
# exit()
# station2 = station.copy()

# filter_domain=[-73.67,-73.61,3.79,3.87]
# polygon = [(filter_domain[0],filter_domain[2]),
#                 (filter_domain[0],filter_domain[3]),
#                 (filter_domain[1],filter_domain[3]),
#                 (filter_domain[1],filter_domain[2]),
#                 (filter_domain[0],filter_domain[2])
#                 ]
# station2 = station2.filter_region(polygon )
# station2.color = "yellow"
# network = Network([station,station2])
# network.append(station2.copy())
# print(network)
# print(networks)
# print(network.__str__(True))

# fig = network.plot()
# fig.show()

# fig = station.matplot()
# plt.show()
# print(station.__str__(True))
# stations.plot()
# plt.show()
# print(stations)

tecto_shp = os.path.join(rep_data,"shapes","Tectónica.shp")
all_data = gpd.read_file(tecto_shp)
all_data_l = all_data[all_data["NombreEstr"] == "Falla de Jordan"]
all_data_r = all_data[all_data["NombreEstr"] != "Falla de Jordan"]

df = all_data["geometry"].apply(lambda x: x.bounds)
df = pd.DataFrame(df.tolist(),columns=["min_x","min_y","max_x","max_y"])                                  
# print(df)
# print(df)

# print(all_data[["min_x","min_y","max_x","max_y"]])
# print(all_data_l.geometry.iloc[0].bounds)
shape_r = Shape(data=all_data_r,
        projection="EPSG:4326")
s = shape_r.copy()
# print(s.data)
s = s.select_data({"OBJECTID":[2]})
s.baseplot.pen = ["0.02c,red,-"]
# print(s.__dict__)
# fig = s.plot()
# fig.show()
# print(s.data)
ss = Shapes([shape_r,s])
fig = ss.plot()
fig.show()
# fig = shape_r.plot()
# fig.show()
# print(shape_r)
# print(shape_r.__str__(True))
# print(shape_r.get_region())
# plt.show()



    