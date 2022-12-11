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

fm_csv = os.path.join(rep_data,"mf","MF_quifa_201301_202203.csv")
events = pd.read_csv(fm_csv)
events = events.rename(columns={"Origin time":"origin_time",
                        "Latitude (deg)":"latitude",
                        "Longitude (deg)":"longitude",
                        "Depth (m)":"depth",
                        "Mag. (Mw)":"magnitude",
                        "Dip n1 (deg)":"dip",
                        "Rake n1 (deg)":"rake",
                        "Strike n1 (deg)":"strike",
                        "Dip n2 (deg)":"dip_n2",
                        "Rake n2 (deg)":"rake_n2",
                        "Strike n2 (deg)":"strike_n2",
                        }
                        )
fm = FM(events,BaseMeca(cmap=False))
# fig = fm.plot()
# fig.show()
# exit()
mulfm = MulFM([fm])
# x = fm.project((-73.67,3.82),(-73.66,3.81),(-10,10))
# print(x)

baseprofile = BaseProfile(
                        projection="x0.002c/-0.002",
                        # projection="X10/-10",
                        depth_lims=[0,3e3],output_unit="m")
profile = Profile(name=("A","A'"),      
        coords=((-73.67,3.82),(-73.66,3.81)), 
        # coords=((-73.690175,3.869714),(-73.666203,3.895990)), 
        # coords=((-80,10),(-70,10)), 
        # width=(-0.1,0.1),
        width=(-10,10),
        baseprofile=baseprofile
            )
profile.add_mulobject(mulfm,depth_unit="m",
                    verbose=True)
fig = profile.plot_profile()
fig.show()