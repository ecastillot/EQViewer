
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
                        "events_20160101T20220901.csv")
events = pd.read_csv(cat_csv)
events["origin_time"] = pd.to_datetime(events['origin_time'])

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


baseprofile = BaseProfile(projection="X10/-10",
                        depth_lims=[0,200],output_unit="km")
profile = Profile(name=("A","A'"),      
        coords=((-80,10),(-70,10)), 
        width=(-300,300),
        baseprofile=baseprofile
            )

profile.add_mulobject(mulcatalog,depth_unit="km")
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