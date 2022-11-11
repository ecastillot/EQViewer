import pandas as pd
index = pd.Index([10, 20, 30, 40, 50, 60, 70])
print("Pandas Index\n",index)
print("\nGet the indexes\n",index.get_indexer([30, 25, 58, 50, 69], method="nearest"))