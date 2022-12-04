
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

cat_csv = os.path.join(rep_data,"earthquakes",
                        "CAS_20190515-20220831.csv")
events = pd.read_csv(cat_csv,sep=";")
events = events.rename(columns={"Origin time":"origin_time",
                                "Depth (m)":"depth",
                                "Mag. (Mw)":"magnitude",
                                "Latitude (deg)":"latitude",
                                "Longitude (deg)":"longitude"})
events["origin_time"] = pd.to_datetime(events['origin_time']).dt.tz_localize(None)

# cat_csv = os.path.join(rep_data,"earthquakes",
#                         "events_20160101T20220901.csv")
# events = pd.read_csv(cat_csv)
# events["origin_time"] = pd.to_datetime(events['origin_time'])

df1 = events[events["origin_time"]<=dt.datetime(2022,1,1)]
df2 = events[events["origin_time"]>=dt.datetime(2022,1,1)]

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

cat1 = Catalog(df1,baseplot1 )
cat2 = Catalog(df2,baseplot2 )

mulcatalog = MulCatalog([cat1,cat2])

data_path = "/home/emmanuel/EQviewer/data/well/survey_proc2.xlsx"
data = pd.read_excel(data_path)
data = data.rename(columns={"MD (ft)":"MD",
                            "Lon (°)":"longitude",
                            "Lat (°)":"latitude",
                            "Z (m)": "depth",
                            "TVD (ft)":"TVD"})
# data_path = "/home/emmanuel/EQviewer/data/well/well_example.csv"
# data = pd.read_csv(data_path)
# data = data.rename(columns={"MD":"MD",
#                             "lon":"longitude",
#                             "lat":"latitude",
#                             "z": "depth",
#                             "TVD":"TVD"})
injection = pd.read_csv("/home/emmanuel/EQviewer/data/well/CASTILLA_ft_ini.dat")
injection = injection.rename(columns={"min":"min_depth",
                                    "max":"max_depth",
                                    "water_flow":"measurement"})
injection = injection[injection["name"]=="CAD01"]
print(injection)
injection = Injection(injection,depth_type="MD",
                    baseplot=BasePlot(cmap=True,
                                color="blue",
                                    style="g0.3")
                                    )
well = Well(data,"PAD",
            baseplot = BasePlot(
                        # size=None,
                        style="c0.5",
                        cmap=True,
                        color=None,
                        label=None,
                        transparency=None,
                        # pen="1p"
                        ),
            injection=injection
                        
                )
# p = well.project((-73.67,3.82),(-73.66,3.81),(-0.1,0.1),with_injection=True)
# print(p)
wellfig = well.plot_map()
wellfig.show()

exit()
mulwell = MulWell([well])
baseprofile = BaseProfile(projection="X10/-10",
                        depth_lims=[0,3e3],output_unit="m")
profile = Profile(name=("A","A'"),      
        coords=((-73.67,3.82),(-73.66,3.81)), 
        # coords=((-73.690175,3.869714),(-73.666203,3.895990)), 
        # coords=((-80,10),(-70,10)), 
        width=(-0.1,0.1),
        baseprofile=baseprofile
            )

# wellfig = mulwell.plot_map()
# wellfig = mulcatalog.plot_map(fig = wellfig )
# wellfig = profile.plot_in_map(wellfig )
# wellfig.show()

# w.show()
profile.add_mulobject(mulcatalog,depth_unit="m",
                    verbose=True)
profile.add_mulobject(mulwell,depth_unit="m",
                    verbose=True)
# x = profile.plot_map(mapfig)
fig = profile.plot_profile()
fig.show()
# print(profile.projections)
# profile.add_mulobject()
# prof_fig = profile.plot_map()
# prof_fig.show()

# fig = catalog.plot_map(profiles=[profile])
# fig.show()
# prof_fig = catalog.plot_profile(profile,
#                             depth_unit="km",
#                             verbose=False)
# prof_fig.show()
# # fig = catalog2.plot_profile(profile,baseprofile,
#                             depth_unit="km",fig=fig,
#                             verbose=False)
# fig.show()

# multicatalog = MulCatalog(catalogs=[catalog,catalog2],
#                             show_cpt=True)
# print(multicatalog.__str__(True))
# fig = multicatalog.plot()
# fig.show()
# prof_fig = multicatalog.plot_profile(profile,
#                             depth_unit="km",
#                             verbose=False)
# prof_fig.show()