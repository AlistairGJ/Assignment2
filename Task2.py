#Task 2: Data Exploration and Visualisation
%run Task1.py
#Cuts data to include only 11 proteins we are going to use
Protein11 = AllProtein[['MouseID', 'BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Genotype', 'Treatment', 'Behavior', 'class']]

#Summary of data
Protein11.dtypes
Protein11.shape
Protein11.columns
Protein11['Genotype'].value_counts()
Protein11['Treatment'].value_counts()
Protein11['Behavior'].value_counts()
Protein11['class'].value_counts()
Protein11.describe()

#Data checking
#Null - this will print out all the mice will null values
proteinNames = ['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']
miceIDs = Protein11['MouseID'].str.split('_').apply(pd.Series, 1)[0].unique()

# Prevents stupid error messages
pd.options.mode.chained_assignment = None

# Renaming Genotype Column
Protein11['Genotype'].replace('Control', '1', inplace=True)
Protein11['Genotype'].replace('Ts65Dn', '0', inplace=True)
Protein11['Genotype'] = Protein11['Genotype'].astype(int)

# Renaming Behavior Column
Protein11['Behavior'].replace('C/S', '1', inplace=True)
Protein11['Behavior'].replace('S/C', '0', inplace=True)
Protein11['Behavior'] = Protein11['Behavior'].astype(int)

# Renaming Treatment Column
Protein11['Treatment'].replace('Saline', '1', inplace=True)
Protein11['Treatment'].replace('Memantine', '0', inplace=True)
Protein11['Treatment'] = Protein11['Treatment'].astype(int)

# Class = Genotype,Behaviour,Treatment into binary and then into decimal
def change_class(row):
    row['Class'] = row['Genotype'] * 4 + row['Behavior'] * 2 + row['Treatment']
    return row

Protein11 = Protein11.apply(change_class, axis=1)

#Outlier checking
def makeTables():
    proteinRows = []
    outlierMiceRows = []

    for proteinName in proteinNames:
        count = Protein11[proteinName].count()
        mean = Protein11[proteinName].mean()
        sd = Protein11[proteinName].std()
        minusThreeSD = mean - (3 * sd)
        minusTwoSD = mean - (2 * sd)
        twoSD = mean + (2 * sd)
        threeSD = mean + (3 * sd)
        outliers = Protein11.query(proteinName + ' < ' + str(minusThreeSD) + ' | ' + proteinName + ' > ' + str(threeSD))
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
        row['Avg All'] = Protein11[proteinName].mean()
        row['SD All'] = Protein11[proteinName].std()

        for c in range(8):
            row['Avg Class ' + str(c)] = Protein11[Protein11.Class == c][proteinName].mean()
            row['SD Class ' + str(c)] = Protein11[Protein11.Class == c][proteinName].std()

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
    print Protein11[Protein11[proteinName].isnull()]

def make_nans(row):
    for proteinName in proteinNames:
        mean = Protein11[proteinName].mean()
        sd = Protein11[proteinName].std()
        minusThreeSD = mean - (3 * sd)
        threeSD = mean + (3 * sd)

        if row[proteinName] < minusThreeSD or row[proteinName] > threeSD:
            row[proteinName] = None
    return row

Protein11 = Protein11.apply(make_nans, axis=1)

print "AFTER"

for proteinName in proteinNames:
    print Protein11[Protein11[proteinName].isnull()]

# Removes mouse 3484_n
Protein11 = Protein11[~Protein11['MouseID'].str.contains('3484')]
indexOfMouse = np.where(miceIDs=='3484')[0]
miceIDs = np.delete(miceIDs, indexOfMouse)

print "Main tables after conversion of outliers to NaN and removal of 3484_n"

makeTables()

def make_averages(row):
    for proteinName in proteinNames:
        if np.isnan(row[proteinName]):
            average = Protein11[Protein11.Class == row['Class']][proteinName].mean()
            row[proteinName] = average
    return row

Protein11 = Protein11.apply(make_averages, axis=1)

print "NANS CONVERTED TO AVERAGES"

for proteinName in proteinNames:
    print Protein11[Protein11[proteinName].isnull()]

print "Main tables after conversion of NAN to Average"

makeTables()
Protein11.to_csv("finalData.csv")

mouseRows = []
for mouseID in miceIDs:
    row = {}
    row['MouseID'] = mouseID
    for proteinName in proteinNames:
        row[proteinName] = Protein11[Protein11['MouseID'].str.contains(mouseID)][proteinName].mean()
    row['Genotype'] = Protein11[Protein11['MouseID'].str.contains(mouseID)]['Genotype'].iloc[0]
    row['Treatment'] = Protein11[Protein11['MouseID'].str.contains(mouseID)]['Treatment'].iloc[0]
    row['Behavior'] = Protein11[Protein11['MouseID'].str.contains(mouseID)]['Behavior'].iloc[0]
    row['Class'] = Protein11[Protein11['MouseID'].str.contains(mouseID)]['Class'].iloc[0]
    mouseRows.append(row)

cols = proteinNames + ['Genotype', 'Treatment', 'Behavior', 'Class']

miceDF = pd.DataFrame(mouseRows, index=miceIDs, columns=cols)

#Summary of mouse average data
miceDF.describe()
miceDF['Genotype'].value_counts()
miceDF['Treatment'].value_counts()
miceDF['Behavior'].value_counts()
miceDF['Class'].value_counts()

#print miceDF
miceDF.to_csv('miceAverages.csv')

#This makes a unique mouse field
MakeMouseID = Protein11['MouseID'].str.split('_').apply(pd.Series, 1)
MakeMouseID.name = 'MouseIDavg'
Protein11 = Protein11.join(MakeMouseID)
Protein11.rename(columns = {0:'MouseIDavg'}, inplace = True)
Protein11.columns
del Protein11[1]
Protein11.head(10)

#Setting up plotting basics
font = {'family' : 'monospace', 'weight' : 'regular', 'size' : '10'}
plt.rc('font', **font)

#Initial graphs per variable and scatter matrix - use all values for all mice
from pandas.tools.plotting import scatter_matrix
scatterProteins = Protein11[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']]
scatter_matrix(scatterProteins, alpha=0.2,figsize=(16,16),diagonal='hist')
plt.suptitle('Protein Distribution Scatter Matrix')
plt.show()

#Histograms
ProteinList = ['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']
colours = ['#a6cee3', '#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99']
for index, data in enumerate(ProteinList):
    Protein11[data].plot(kind="hist", bins=20, color=colours[index])
    plt.title("Distribution of " + data)
    plt.show()

#Graphs using avg mice data
miceDF_scatter = miceDF[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']]
scatter_matrix(miceDF_scatter, alpha=0.2,figsize=(16,16),diagonal='hist')
plt.suptitle('Protein Distribution - One Value per Mouse')
plt.show()

for index, data in enumerate(ProteinList):
    miceDF[data].plot(kind="hist", bins=20, color=colours[index])
    plt.title("Distribution of " + data)
    plt.show()

#Getting a numeric value for mouse to use in scatter plots
Protein11['MouseNo'] = Protein11.index
Protein11.columns

#Scatter plot showing all proteins on the one plot
colours2 = ['#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99']
ProteinList2 = ['pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']
ax = Protein11.plot(kind='scatter', figsize=(12,12), x=18, y=1, color='#a6cee3', label='BRAF_N',s=50);
for index, data in enumerate(ProteinList2):
    Protein11.plot(kind='scatter', x=18, y=index+1, color=colours2[index], label=data, ax=ax, s=50);
plt.grid()
plt.legend(loc='upper left')
plt.ylim(0,3)
plt.ylabel('Protein Expression Level')
plt.xlabel('Mouse')
plt.show()

for index, data in enumerate(ProteinList):
    Protein11.plot(kind='scatter', figsize=(6,6), x=18, y=index+1, color=colours[index], label=data, s=30, facecolor='0.5')
    plt.grid()
    plt.legend(loc='upper left')
    plt.ylabel(data + ' Expression Level')
    plt.xlabel('Mouse')
    plt.title(data + " Expression")
    plt.show()

#Scatter matrix for 4 highly correlated variables
Protein_correlated = Protein11[['MouseID', 'BRAF_N', 'pERK_N', 'DYRK1A_N','ITSN1_N']]
scatter_matrix(Protein_correlated, alpha=0.2,figsize=(16,16),diagonal='hist')
plt.suptitle("Scatter Matrix for Highly Correlated Proteins")
plt.show()

#Scatter plot per protein by genotype
import matplotlib.patches as mpatches
Protein11['Genotype'].value_counts()
colour_palette = {1 :'#ed2939', 0:'#2cc2e8'}
colors = [colour_palette[c] for c in Protein11['Genotype']]
colours = ['#ed2939', '#2cc2e8']
for index, data in enumerate(ProteinList):
    Protein11.plot(kind='scatter', x=18, y=index+1, s=50, c=colors)
    plt.title(str(data) + ' Expression by Mouse Genotype')
    plt.xlabel('Mouse')
    plt.ylabel(data)
    plt.grid(True, which='major', color='#131313', linestyle='-')
    plt.minorticks_on()
    recs = []
    labels = 'Control', 'Ts65Dn'
    for i in range(0, len(colours)):
        recs.append(mpatches.Rectangle((0,0),1,1,fc=colours[i]))
    plt.legend(recs, labels, loc=2)
    plt.show()

#Scatter plot per protein by treatment
Protein11['Treatment'].value_counts()
colour_palette = {0 :'#ade799', 1:'#d399e7'}
colors = [colour_palette[c] for c in Protein11['Treatment']]
colours = ['#ade799', '#d399e7']
for index, data in enumerate(ProteinList):
    Protein11.plot(kind='scatter', x=18, y=index+1, s=50, c=colors)
    plt.title(str(data) + ' Expression by Treatment')
    plt.xlabel('Mouse')
    plt.ylabel(data)
    plt.grid(True, which='major', color='#131313', linestyle='-')
    plt.minorticks_on()
    recs = []
    labels = 'Memantine', 'Saline'
    for i in range(0, len(colours)):
        recs.append(mpatches.Rectangle((0,0),1,1,fc=colours[i]))
    plt.legend(recs, labels, loc=2)
    plt.show()

#Scatter plot per protein by behaviour
Protein11['Behavior'].value_counts()
colour_palette = {0 :'#ffdd0f', 1:'#f22f74'}
colors = [colour_palette[c] for c in Protein11['Behavior']]
colours = ['#ffdd0f', '#f22f74']
for index, data in enumerate(ProteinList):
    Protein11.plot(kind='scatter', x=18, y=index+1, s=50, c=colors)
    plt.title(str(data) + ' Expression by Mouse Behaviour')
    plt.xlabel('Mouse')
    plt.ylabel(data)
    plt.grid(True, which='major', color='#131313', linestyle='-')
    plt.minorticks_on()
    recs = []
    labels = 'S/C', 'C/S'
    for i in range(0, len(colours)):
        recs.append(mpatches.Rectangle((0,0),1,1,fc=colours[i]))
    plt.legend(recs, labels, loc=2)
    plt.show()

#Scatter plot per protein by class
Protein11['class'].value_counts()
colour_palette = {'c-SC-m' :'#e41a1c', 'c-CS-m':'#377eb8', 'c-CS-s' :'#4daf4a', 't-CS-m':'#984ea3','t-SC-s' :'#ff7f00', 't-SC-m':'#ffff33', 'c-SC-s' :'#a65628', 't-CS-s':'#f781bf'}
colors = [colour_palette[c] for c in Protein11['class']]
colours = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf']
for index, data in enumerate(ProteinList):
    Protein11.plot(kind='scatter', x=18, y=index+1, s=30, c=colors)
    plt.title(str(data) + ' Expression by Mouse Class')
    plt.xlabel('Mouse')
    plt.ylabel(data)
    plt.grid(True, which='major', color='#131313', linestyle='-')
    plt.minorticks_on()
    recs = []
    labels = 'c-SC-m','c-CS-m','c-CS-s','t-CS-m','t-SC-s','t-SC-m','c-SC-s','t-CS-s'
    for i in range(0, len(colours)):
        recs.append(mpatches.Rectangle((0,0),1,1,fc=colours[i]))
    plt.legend(recs, labels, loc=2)
    plt.show()

