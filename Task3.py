#Task 3: Data Modelling (Classification)
%run Task2.py

#Classification 1: K Nearest Neighbor
#Import correct modules
from sklearn.cross_validation import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix

#Using SOD1, Genotype & Treatment to predict Behavior
#Function for classification
Protein11_SOD1 = Protein11[['SOD1_N', 'Genotype', 'Treatment']]
Protein11_SOD1.dtypes
Protein11_SOD1.describe()
X_train, X_test, y_train, y_test = train_test_split(Protein11_SOD1, Protein11['Behavior'], test_size=0.4)
X_train.shape
y_train.shape
#Analysis: 5 neighbours
clf = KNeighborsClassifier(5)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm
#Analysis: 2 neighbours
clf = KNeighborsClassifier(2)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "2 neighbours predicted: "
print cm
#Analysis: 8 neighbours
clf = KNeighborsClassifier(8)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "8 neighbours predicted: "
print cm

#Using ITSN1, Genotype & Treatment to predict Behavior
Protein11_ITSN1 = Protein11[['ITSN1_N', 'Genotype', 'Treatment']]
print Protein11_ITSN1.dtypes
Protein11_ITSN1.describe()
X_train, X_test, y_train, y_test = train_test_split(Protein11_ITSN1, Protein11['Behavior'], test_size=0.4)
X_train.shape
y_train.shape
#Analysis: 5 neighbours
clf = KNeighborsClassifier(5)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm
#Analysis: 2 neighbours
clf = KNeighborsClassifier(2)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "2 neighbours predicted: "
print cm
#Analysis: 8 neighbours
clf = KNeighborsClassifier(8)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "8 neighbours predicted: "
print cm
#Analysis: 5 neighbours, distance weighted
clf = KNeighborsClassifier(5, weights='distance')
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "5 neighbours, distance weighted predicted: "
print cm
#Analysis: 8 neighbours, distance weighted
clf = KNeighborsClassifier(8, weights='distance')
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "8 neighbours, distance weighted predicted: "
print cm
#Analysis: 5 neighbours, distance weighted, p=1
clf = KNeighborsClassifier(5, weights='distance', p=1)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "5 neighbours, distance weighted predicted, p=1: "
print cm
#Analysis: 8 neighbours, distance weighted, p=1
clf = KNeighborsClassifier(8, weights='distance', p=1)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "8 neighbours, distance weighted predicted, p=1: "
print cm

#Use all protein expression levels to clasify class
Protein11_class = Protein11[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N']]
print Protein11_class.dtypes
print Protein11_class.describe()
X_train, X_test, y_train, y_test = train_test_split(Protein11_class, Protein11['class'], test_size=0.4)
print X_train.shape
print y_train.shape

#Analysis: 5 neighbours
clf = KNeighborsClassifier(5)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "5 neighbours predicted: "
print cm
#Analysis: 2 neighbours
clf = KNeighborsClassifier(2)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "2 neighbours predicted: "
print cm
#Analysis: 8 neighbours
clf = KNeighborsClassifier(8)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "8 neighbours predicted: "
print cm
#Analysis: 5 neighbours, distance weighted
clf = KNeighborsClassifier(5, weights='distance')
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "5 neighbours, distance weighted predicted: "
print cm
#Analysis: 8 neighbours, distance weighted
clf = KNeighborsClassifier(8, weights='distance')
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "8 neighbours, distance weighted predicted: "
print cm
#Analysis: 5 neighbours, distance weighted, p=1
clf = KNeighborsClassifier(5, weights='distance', p=1)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "5 neighbours, distance weighted predicted, p=1: "
print cm
#Analysis: 8 neighbours, distance weighted, p=1
clf = KNeighborsClassifier(8, weights='distance', p=1)
fit = clf.fit(X_train, y_train)
predicted = fit.predict(X_test)
cm = confusion_matrix(y_test, predicted)
print "8 neighbours, distance weighted predicted, p=1: "
print cm