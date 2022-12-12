import pandas as pd
# df = pd.read_csv("/home/emmanuel/EQviewer/data/earthquakes/catalog1.csv",delim_whitespace=True)
# df["origin_time"]  = df['Date'] + ' ' + df['Time']
# df["origin_time"] = pd.to_datetime(df["origin_time"])
# df = df.rename(columns={"Lat":"latitude","Lon":"longitude","Depth":"depth","Mag":"magnitude"})
# df = df[["origin_time","latitude","longitude","depth","magnitude"]]
# df.to_csv("/home/emmanuel/EQviewer/data/earthquakes/catalog.csv",index=False)
# print(df)

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

# https://ncedc.org/egs/catalog-search.html
df = pd.read_csv("/home/emmanuel/EQviewer/data/earthquakes/catalog.csv",
                parse_dates=["origin_time"])
baseplot1 = BasePlot(color = "gray",
                    style="cc",
                    size=lambda x: 0.11 * np.sqrt(1.2 ** (x.magnitude*1.4)),
                #     style="c0.1c",
                #     size=None,
                    cmap = False,
                    pen = "black")
cat1 = Catalog(df,baseplot1 )

data = pd.read_csv("/home/emmanuel/EQviewer/data/wells/survey/well_7.csv")
well = Well(data,"PAD",
            baseplot = BasePlot(
                        # size=None,
                        # style="c0.2",
                        cmap=False,
                        color=None,
                        label=None,
                        transparency=None,
                        pen="1p,red"
                        )
                )
print(well.__str__(True))
# print(cat1)
# print(cat1.__str__(True))
fig = cat1.plot_map()
fig = well.plot_map(fig=fig)
fig.show()

