import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import urllib2
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00342/Data_Cortex_Nuclear.xls"
requestURL = urllib2.Request(url)
allProtein_filename = urllib2.urlopen(requestURL)
allProtein = pd.read_excel(allProtein_filename, headers=0)
allProtein.rename(columns={'class': 'Class'}, inplace=True)

proteinNames = ['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']
miceIDs = allProtein['MouseID'].str.split('_').apply(pd.Series, 1)[0].unique()

#print AllProtein.describe()

protein11 = allProtein[['MouseID', 'BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Genotype', 'Treatment', 'Behavior', 'Class']]

#print Protein11.describe()

# Prevents stupid error messages
pd.options.mode.chained_assignment = None

# Renaming Genotype Column
# Control = 1
# Ts65Dn = 0

protein11['Genotype'].replace('Control', '1', inplace=True)
protein11['Genotype'].replace('Ts65Dn', '0', inplace=True)
protein11['Genotype'] = protein11['Genotype'].astype(int)

# Renaming Behavior Column
# C/S = 1
# S/C = 0

protein11['Behavior'].replace('C/S', '1', inplace=True)
protein11['Behavior'].replace('S/C', '0', inplace=True)
protein11['Behavior'] = protein11['Behavior'].astype(int)

# Renaming Treatment Column
# Saline = 1
# Memantine = 0

protein11['Treatment'].replace('Saline', '1', inplace=True)
protein11['Treatment'].replace('Memantine', '0', inplace=True)
protein11['Treatment'] = protein11['Treatment'].astype(int)

# Class = Genotype,Behaviour,Treatment into binary and then into decimal
def change_class(row):
    row['Class'] = row['Genotype'] * 4 + row['Behavior'] * 2 + row['Treatment']
    return row

protein11 = protein11.apply(change_class, axis=1)

#print protein11


def makeTables():
    proteinRows = []
    outlierMiceRows = []

    for proteinName in proteinNames:

        count = protein11[proteinName].count()
        mean = protein11[proteinName].mean()
        sd = protein11[proteinName].std()
        minusThreeSD = mean - (3 * sd)
        minusTwoSD = mean - (2 * sd)
        twoSD = mean + (2 * sd)
        threeSD = mean + (3 * sd)
        outliers = protein11.query(proteinName + ' < ' + str(minusThreeSD) + ' | ' + proteinName + ' > ' + str(threeSD))

        row = {'Protein': proteinName,
               'Count': count,
               'Mean': mean,
               'SD': sd,
               '-3SD': minusThreeSD,
               '-2SD': minusTwoSD,
               '+2SD': twoSD,
               '+3SD': threeSD,
               'Outliers': outliers[proteinName].count()}
        proteinRows.append(row)

        if outliers.empty:
            row = {'Protein': proteinName,
                   'MouseID': '-',
                   '# Instances': '-',
                   'Genotype': '-',
                   'Treatment': '-',
                   'Behavior': '-',
                   'Class': '-'}
            outlierMiceRows.append(row)
        else:
            for mouseID in miceIDs:
                mouseOutlierRows = outliers[outliers['MouseID'].str.contains(mouseID)]
                if not mouseOutlierRows.empty:
                    row = {'Protein': proteinName,
                           'MouseID': mouseID,
                           '# Instances': len(mouseOutlierRows),
                           'Genotype': mouseOutlierRows['Genotype'].iloc[0],
                           'Treatment': mouseOutlierRows['Treatment'].iloc[0],
                           'Behavior': mouseOutlierRows['Behavior'].iloc[0],
                           'Class': mouseOutlierRows['Class'].iloc[0]}
                    outlierMiceRows.append(row)

    nnPctRangeDF = pd.DataFrame(proteinRows, index=proteinNames,
                                columns=['Count', 'Mean', 'SD', '-3SD', '-2SD', '+2SD', '+3SD', 'Outliers'])
    print nnPctRangeDF

    outliersDF = pd.DataFrame(outlierMiceRows,
                              columns=['Protein', 'MouseID', '# Instances', 'Genotype', 'Treatment', 'Behavior', 'Class'])
    print outliersDF


    # Average for each class for each protein
    # Rows: Each of the 11 Proteins
    # Columns: Average for All, Std Dev All, Average for Class 0, Std Dev 0, ... Average for Class 7, Std Dev 7

    classAveragesRows = []

    for proteinName in proteinNames:
        row = {}
        row['Protein'] = proteinName
        row['Avg All'] = protein11[proteinName].mean()
        row['SD All'] = protein11[proteinName].std()

        for c in range(8):
            row['Avg Class ' + str(c)] = protein11[protein11.Class == c][proteinName].mean()
            row['SD Class ' + str(c)] = protein11[protein11.Class == c][proteinName].std()

        classAveragesRows.append(row)

    outliersDF = pd.DataFrame(classAveragesRows, index=proteinNames,
                              columns=['Avg All', 'SD All',
                                       'Avg Class 0', 'SD Class 0',
                                       'Avg Class 1', 'SD Class 1',
                                       'Avg Class 2', 'SD Class 2',
                                       'Avg Class 3', 'SD Class 3',
                                       'Avg Class 4', 'SD Class 4',
                                       'Avg Class 5', 'SD Class 5',
                                       'Avg Class 6', 'SD Class 6',
                                       'Avg Class 7', 'SD Class 7'])

    print outliersDF


print "Main tables with outliers included"

makeTables()

print "BEFORE"

for proteinName in proteinNames:
    print protein11[protein11[proteinName].isnull()]

def make_nans(row):
    for proteinName in proteinNames:
        mean = protein11[proteinName].mean()
        sd = protein11[proteinName].std()
        minusThreeSD = mean - (3 * sd)
        threeSD = mean + (3 * sd)

        if row[proteinName] < minusThreeSD or row[proteinName] > threeSD:
            row[proteinName] = None
    return row

protein11 = protein11.apply(make_nans, axis=1)

print "AFTER"

for proteinName in proteinNames:
    print protein11[protein11[proteinName].isnull()]

# Removes mouse 3484_n
protein11 = protein11[~protein11['MouseID'].str.contains('3484')]
indexOfMouse = np.where(miceIDs=='3484')[0]
miceIDs = np.delete(miceIDs, indexOfMouse)


print "Main tables after conversion of outliers to NaN and removal of 3484_n"

makeTables()

def make_averages(row):
    for proteinName in proteinNames:
        if np.isnan(row[proteinName]):
            average = protein11[protein11.Class == row['Class']][proteinName].mean()
            row[proteinName] = average
    return row

protein11 = protein11.apply(make_averages, axis=1)

print "NANS CONVERTED TO AVERAGES"

for proteinName in proteinNames:
    print protein11[protein11[proteinName].isnull()]


print "Main tables after conversion of NAN to Average"

makeTables()

protein11.to_csv("finalData.csv")

mouseRows = []
for mouseID in miceIDs:
    row = {}
    row['MouseID'] = mouseID
    for proteinName in proteinNames:
        row[proteinName] = protein11[protein11['MouseID'].str.contains(mouseID)][proteinName].mean()

    row['Genotype'] = protein11[protein11['MouseID'].str.contains(mouseID)]['Genotype'].iloc[0]
    row['Treatment'] = protein11[protein11['MouseID'].str.contains(mouseID)]['Treatment'].iloc[0]
    row['Behavior'] = protein11[protein11['MouseID'].str.contains(mouseID)]['Behavior'].iloc[0]
    row['Class'] = protein11[protein11['MouseID'].str.contains(mouseID)]['Class'].iloc[0]
    mouseRows.append(row)

cols = proteinNames + ['Genotype', 'Treatment', 'Behavior', 'Class']

miceDF = pd.DataFrame(mouseRows, index=miceIDs, columns=cols)

print miceDF

miceDF.to_csv('miceAverages.csv')

from pandas.tools.plotting import scatter_matrix
scatterProteins = protein11[['MouseID', 'BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Class']]
scatter_matrix(scatterProteins, alpha=0.2,figsize=(16,16),diagonal='hist')
plt.show()

