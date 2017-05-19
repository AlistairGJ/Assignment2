import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import urllib2
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00342/Data_Cortex_Nuclear.xls"
requestURL = urllib2.Request(url)
allProtein_filename = urllib2.urlopen(requestURL)
allProtein = pd.read_excel(allProtein_filename, headers=0)


proteinNames = ['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']

proteinRows = []
outlierMiceRows = []

miceIDs = allProtein['MouseID'].str.split('_').apply(pd.Series, 1)[0].unique()

#print AllProtein.describe()

protein11 = allProtein[['MouseID', 'BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Genotype', 'Treatment', 'Behavior', 'class']]

#print Protein11.describe()

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
                       'Class': mouseOutlierRows['class'].iloc[0]}
                outlierMiceRows.append(row)

nnPctRangeDF = pd.DataFrame(proteinRows,
                            columns=['Protein', 'Count', 'Mean', 'SD', '-3SD', '-2SD', '+2SD', '+3SD', 'Outliers'])
print nnPctRangeDF

outliersDF = pd.DataFrame(outlierMiceRows,
                          columns=['Protein', 'MouseID', '# Instances', 'Genotype', 'Treatment', 'Behavior', 'Class'])
print outliersDF
