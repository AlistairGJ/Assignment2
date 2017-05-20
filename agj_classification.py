from sklearn.tree import DecisionTreeClassifier
from sklearn.cross_validation import train_test_split
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn import tree
import os

import pandas as pd

df = pd.read_csv("finalData.csv")
df_class = df[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Class']]

y = df_class.pop('Class')
X = df_class

X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.2, random_state = 0)

#X_train, X_test, y_train, y_test = train_test_split(proteins, proteins, random_state = 0)
clf = DecisionTreeClassifier()
fit = clf.fit(X_train, y_train)

y_pre = fit.predict(X_test)
cm = confusion_matrix(y_test, y_pre)

print "Class"
print cm
print classification_report(y_test, y_pre)
with open("class_tree.dot", 'w') as f:
    f = tree.export_graphviz(clf, out_file=f,
                             feature_names=['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N',
                                            'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N'], class_names="01234567",
                             filled=True, rounded=True, special_characters=True)
os.system("dot class_tree.dot -o class_tree.png -Tpng")

####

df_genotype = df[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Genotype']]

y = df_genotype.pop('Genotype')
X = df_genotype

X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.2, random_state = 0)

#X_train, X_test, y_train, y_test = train_test_split(proteins, proteins, random_state = 0)
clf = DecisionTreeClassifier()
fit = clf.fit(X_train, y_train)

y_pre = fit.predict(X_test)
cm = confusion_matrix(y_test, y_pre)

print "Genotype"
print cm
print classification_report(y_test, y_pre)
with open("genotype_tree.dot", 'w') as f:
    f = tree.export_graphviz(clf, out_file=f,
                             feature_names=['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N',
                                            'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N'], class_names="01",
                             filled=True, rounded=True, special_characters=True)
os.system("dot genotype_tree.dot -o genotype_tree.png -Tpng")


####

df_behavior = df[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Behavior']]

y = df_behavior.pop('Behavior')
X = df_behavior

X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.2, random_state = 0)

#X_train, X_test, y_train, y_test = train_test_split(proteins, proteins, random_state = 0)
clf = DecisionTreeClassifier()
fit = clf.fit(X_train, y_train)

y_pre = fit.predict(X_test)
cm = confusion_matrix(y_test, y_pre)

print "Behavior"
print cm
print classification_report(y_test, y_pre)
with open("behavior_tree.dot", 'w') as f:
    f = tree.export_graphviz(clf, out_file=f,
                             feature_names=['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N',
                                            'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N'], class_names="01",
                             filled=True, rounded=True, special_characters=True)
os.system("dot behavior_tree.dot -o behavior_tree.png -Tpng")


####

df_treatment = df[['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N', 'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N', 'Treatment']]

y = df_treatment.pop('Treatment')
X = df_treatment

X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.2, random_state = 0)

#X_train, X_test, y_train, y_test = train_test_split(proteins, proteins, random_state = 0)
clf = DecisionTreeClassifier()
fit = clf.fit(X_train, y_train)

y_pre = fit.predict(X_test)
cm = confusion_matrix(y_test, y_pre)

print "Treatment"
print cm
print classification_report(y_test, y_pre)
with open("treatment_tree.dot", 'w') as f:
    f = tree.export_graphviz(clf, out_file=f,
                             feature_names=['BRAF_N', 'pERK_N', 'S6_N', 'pGSK3B_N', 'CaNA_N', 'CDK5_N', 'pNUMB_N',
                                            'DYRK1A_N', 'ITSN1_N', 'SOD1_N', 'GFAP_N'],
                             class_names=["Mamantine","Saline"],
                             filled=True, rounded=True, special_characters=True)
os.system("dot treatment_tree.dot -o treatment_tree.png -Tpng")

####

df_class_braf = df[['BRAF_N','Class']]

y = df_class_braf.pop('Class')
X = df_class_braf

X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.2, random_state = 0)

#X_train, X_test, y_train, y_test = train_test_split(proteins, proteins, random_state = 0)
clf = DecisionTreeClassifier()
fit = clf.fit(X_train, y_train)

y_pre = fit.predict(X_test)
cm = confusion_matrix(y_test, y_pre)

print "Class by BRAF"
print cm
print classification_report(y_test, y_pre)
with open("classbybraf_tree.dot", 'w') as f:
    f = tree.export_graphviz(clf, out_file=f,
                             feature_names=['BRAF_N'], class_names="01234567",
                             filled=True, rounded=True, special_characters=True)
os.system("dot classbybraf_tree.dot -o classbybraf_tree.png -Tpng")

