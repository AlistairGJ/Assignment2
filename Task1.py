#Task 1: Data Retrieving
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import urllib2
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00342/Data_Cortex_Nuclear.xls"
RequestURL = urllib2.Request(url)
AllProtein_filename = urllib2.urlopen(RequestURL)
AllProtein = pd.read_excel(AllProtein_filename, headers=0)