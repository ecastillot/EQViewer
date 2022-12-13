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

# myshp = ""
# all_data = gpd.read_file(myshp)