import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import re

df = pd.read_excel(r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\Task 24_AuburnSamples_JVSummary_SampleLevel_SabinSamples_forDDD.xlsx")

print(df)