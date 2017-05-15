#Task 1: Data Retrieving
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import urllib2
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00342/Data_Cortex_Nuclear.xls"
RequestURL = urllib2.Request(url)
AllProtein_filename = urllib2.urlopen(RequestURL)
AllProtein = pd.read_excel(AllProtein_filename, headers=0)

#Task 2: Data Exploration and Visualisation
#Once split into separate scripts should start with
#%run Task1.py
AllProtein.describe()

#Makes a unique mouse column
MakeMouseID = AllProtein['MouseID'].str.split('_').apply(pd.Series, 1)
MakeMouseID.name = 'MouseIDavg'
AllProteinID = AllProtein.join(MakeMouseID)
AllProteinID.rename(columns = {0:'MouseIDavg'}, inplace = True)
AllProteinID.columns
AllProteinID.head(10)

del AllProteinID[1]
del AllProtein

#Summary of data 1
AllProteinID.dtypes
AllProteinID.describe()
AllProteinID['Genotype'].value_counts()
AllProteinID['Treatment'].value_counts()
AllProteinID['Behavior'].value_counts()
AllProteinID['class'].value_counts()

#Cuts data to include only 11 proteins we are going to use
Protein11 = AllProteinID[['MouseID', 'BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Genotype', 'Treatment', 'Behavior', 'class', 'MouseIDavg']]

#Summary of data 2
Protein11.describe()
Protein11['Genotype'].value_counts()
Protein11['Treatment'].value_counts()
Protein11['Behavior'].value_counts()
Protein11['class'].value_counts()

#Data checking
#Null
mask_BRAF_null = pd.isnull(Protein11['BRAF_N'])
Protein11.loc[mask_BRAF_null, 'MouseID']
mask_pERK_null = pd.isnull(Protein11['pERK_N'])
Protein11.loc[mask_pERK_null, 'MouseID']
mask_pNUMB_null = pd.isnull(Protein11['pNUMB_N'])
Protein11.loc[mask_pNUMB_null, 'MouseID']
mask_DYRK1A_null = pd.isnull(Protein11['DYRK1A_N'])
Protein11.loc[mask_DYRK1A_null, 'MouseID']
mask_ITSN1_null = pd.isnull(Protein11['ITSN1_N'])
Protein11.loc[mask_ITSN1_null, 'MouseID']
mask_SOD1_null = pd.isnull(Protein11['SOD1_N'])
Protein11.loc[mask_SOD1_null, 'MouseID']

#Outliers
mask_BRAF_outlier = Protein11['BRAF_N'] > 2
print Protein11.loc[mask_BRAF_outlier, 'MouseID']
mask_pERK_outlier = Protein11['pERK_N'] > 3
print Protein11.loc[mask_pERK_outlier, 'MouseID']
mask_DYRK1A_outlier = Protein11['DYRK1A_N'] > 2
print Protein11.loc[mask_DYRK1A_outlier, 'MouseID']
mask_ITSN1_outlier = Protein11['ITSN1_N'] > 2
print Protein11.loc[mask_ITSN1_outlier, 'MouseID']

#Initial graphs per variable and scatter matrix
#Histograms
ProteinList = ['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']
colours = ['#a6cee3', '#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99']
for index, data in enumerate(ProteinList):
    Protein11[data].plot(kind="hist", bins=20)
    plt.title("Distribution of " + data)
    plt.xlabel("Mice")
    #plt.show()

#Scatter plots
##for index, data in enumerate(ProteinList):
##    Protein11.plot(kind='scatter', figsize=(6,6), x=16, y=index, color=colours[index], label=data)
##    plt.grid()
##    plt.legend(loc='upper left')
##    plt.xlim(-1,73)
##    plt.ylim(0,3)
##    plt.ylabel(data + ' Expression Level')
##    plt.xlabel('Mouse')
##    plt.title(data + " Expression")
##    #plt.show()
##    #plt.savefig(data + "_scatter.png") #need to get this working

from pandas.tools.plotting import scatter_matrix
scatter_matrix(Protein11, alpha=0.2,figsize=(16,16),diagonal='hist')
#plt.show()

#Creates one line per mouse by averaging variables
ProteinMeans = Protein11.groupby(['MouseIDavg'], group_keys=True).mean()
ProteinMeans['MouseIDavg'] = ProteinMeans.index
ProteinObjects = AllProteinID[['MouseIDavg', 'Genotype', 'Treatment', 'Behavior', 'class']]
ProteinObjectsGrp = ProteinObjects.groupby(['MouseIDavg']).first()
ProteinObjectsGrp['MouseIDavg'] = ProteinObjectsGrp.index
ProteinData = pd.merge(ProteinMeans, ProteinObjectsGrp, on='MouseIDavg', how='left')

#Summary of data 3
ProteinData.describe()
ProteinData['Genotype'].value_counts()
ProteinData['Treatment'].value_counts()
ProteinData['Behavior'].value_counts()
ProteinData['class'].value_counts()

#Clean up old names
del AllProtein_filename
del AllProteinID
del ProteinMeans
del ProteinObjects
del ProteinObjectsGrp

#Setting up plotting basics
font = {'family' : 'monospace', 'weight' : 'regular', 'size' : '10'}
plt.rc('font', **font)

#Graphs to check for outliers
colours3 = ['#a6cee3', '#1f78b4','#ff7f00','#cab2d6']
ProteinList3 = ['BRAF_N', 'pERK_N', 'DYRK1A_N', 'ITSN1_N']

for index, data in enumerate(ProteinList3):
    ProteinData.plot(kind='box', y=index, color=colours3[index], label=data)
    plt.grid()
    plt.title(data + " Expression")
    #plt.show()

#Graphs
from pandas.tools.plotting import scatter_matrix
scatter_matrix(ProteinData, alpha=0.2,figsize=(16,16),diagonal='hist')
#plt.show()

#Getting a numeric value for mouse to use in scatter plots
ProteinData['MouseNo'] = ProteinData.index
ProteinData.columns

#Scatter plot showing all proteins on the one plot
colours2 = ['#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99']
ProteinList2 = ['pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']
ax = ProteinData.plot(kind='scatter', figsize=(12,12), x=16, y=0, color='#a6cee3', label='BRAF_N');
for index, data in enumerate(ProteinList2):
    ProteinData.plot(kind='scatter', x=16, y=index+1, color=colours2[index], label=data, ax=ax);
plt.grid()
plt.legend(loc='upper left')
plt.xlim(-1,73)
plt.ylim(0,3)
plt.ylabel('Protein Expression Level')
plt.xlabel('Mouse')
#plt.show()

colours = ['#a6cee3', '#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99']
ProteinList = ['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']

for index, data in enumerate(ProteinList):
    ProteinData.plot(kind='scatter', figsize=(6,6), x=16, y=index, color=colours[index], label=data)
    plt.grid()
    plt.legend(loc='upper left')
    plt.xlim(-1,73)
    plt.ylim(0,3)
    plt.ylabel(data + ' Expression Level')
    plt.xlabel('Mouse')
    plt.title(data + " Expression")
    #plt.show()
 
#Task 3: Data Modelling (Classification)
#Once split into separate scripts should start with
#%run Task1.py

#Change object variables into boolean for analysis
ProteinData['Genotype'] = (ProteinData['Genotype']!='Control').astype(int) # Makes Control=0 and TS65DN=1
ProteinData['Treatment'] = (ProteinData['Treatment']!='Saline').astype(int) # Makes Saline=0 and Memantine=1
ProteinData['Behavior'] = (ProteinData['Behavior']!='S/C').astype(int) # Makes S/C=0 and C/S=1
ProteinData.dtypes

#Classification 1: K Nearest Neighbor
from sklearn.cross_validation import train_test_split
ProteinData_kneighbors = ProteinData[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Genotype', 'Treatment', 'Behavior']]
ProteinData_kneighbors.dtypes
ProteinData_kneighbors.describe()

X_train, X_test, y_train, y_test = train_test_split(ProteinData_kneighbors, ProteinData_kneighbors['Behavior'], test_size=0.4)

X_train.shape
y_train.shape

#Analysis: 3 neighbours
from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(3)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "3 neighbours test: "
print predicted
predicted.shape
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "3 neighbours predicted: "
print cm

#Analysis: 8 neighbours
#from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(8)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "8 neighbours test: "
print predicted
predicted.shape
#from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "8 neighbours predicted: "
print cm

#Analysis: 10 neighbours
#from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(10)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "10 neighbours test: "
print predicted
predicted.shape
#from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "10 neighbours predicted: "
print cm

#Analysis: 5 neighbours
#from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(5)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "5 neighbours test: "
print predicted
predicted.shape
#from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm

#Analysis: 5 neighbours, weighted distance
#from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(5, weights='distance')
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "5 neighbours test: "
print predicted
predicted.shape
#from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm

#Analysis: 5 neighbours, p=1
#from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(5, p=1)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "5 neighbours test: "
print predicted
predicted.shape
#from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm

#Analysis: 5 neighbours, p=3
#from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(5, p=3)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "5 neighbours test: "
print predicted
predicted.shape
#from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm

X_train, X_test, y_train, y_test = train_test_split(ProteinData_kneighbors, ProteinData_kneighbors['Genotype'], test_size=0.4)

X_train.shape
y_train.shape

#Analysis: 5 neighbours
from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(5, weights='distance')
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
print "5 neighbours test: "
print predicted
predicted.shape
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm


