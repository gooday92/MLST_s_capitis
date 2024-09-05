import pandas as pd
import os

df = pd.read_excel("stastic.xlsx", index_col=0)
print(df)
print(df.describe())