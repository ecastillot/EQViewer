import pandas as pd

df = pd.read_csv("/home/emmanuel/EDCT/EQviewer/data/well/well_example.csv")
df["lon"] = df["lon"]+5
df.to_csv("/home/emmanuel/EDCT/EQviewer/data/well/well_example.csv",index=False)