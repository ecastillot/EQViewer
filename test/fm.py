import sys
import os
import numpy as np
import datetime as dt
repository_path = r"/home/emmanuel/EQviewer"  ##change this path where is located the main directory
rep_data = os.path.join(repository_path,"data")
rep_out = os.path.join(repository_path,"example")
sys.path.insert(0,repository_path)
import geopandas as gpd
from EQViewer.eqviewer import *
import EQViewer.utils as equt
import matplotlib.pyplot as plt
import datetime as t
import random
import pandas as pd

df = pd.read_csv("/home/emmanuel/EQviewer/data/mf/mf.csv")

fm = FM(df,BaseMeca(cmap=False))
# fig= fm.plot_map()
# fig.show()
mulfm = MulFM([fm])

baseprofile = BaseProfile(
                        projection="x0.002c/-0.002",
                        # projection="X10/-10",
                        depth_lims=[0,3e3],output_unit="m")
profile = Profile(name=("A","A'"),      
        coords=((-122.86865,38.86238),(-122.70583,38.74421)), 
        # coords=((-73.690175,3.869714),(-73.666203,3.895990)), 
        # coords=((-80,10),(-70,10)), 
        # width=(-0.1,0.1),
        width=(-2,2),
        baseprofile=baseprofile
            )
profile.add_mulobject(mulfm,depth_unit="m",
                    verbose=True)
fig = profile.plot()
fig.show()