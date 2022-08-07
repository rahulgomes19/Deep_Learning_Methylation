#!/usr/bin/env python
# coding: utf-8

# In[1]:


import random
random.seed(0)


# In[2]:


from platform import python_version
print(python_version())


# In[3]:


import pandas as pd

df = pd.read_csv("data/BetaData_27K_SimpleImpute_Mean_1.csv")


# In[4]:


df = df.iloc[: , 1:]  #deleting 1st column 
df = df.drop(['Donor_Sample'], axis=1) #deleting this column
df = df.fillna(0)


# In[5]:


#target = df['is_tumor']
df.is_tumor = df.is_tumor.astype(int) #converted is_tumor value float to int


# In[6]:


print("Number of markers/ features", len(df.columns)-1)
print("Number of samples", len(df))
print(df['is_tumor'].value_counts())
#print(df)


# Splitting data into features and target

# In[7]:


#split dataset into features and target
X = df.iloc[:,0:len(df.columns)-1].values
print(type(X))
print("X = ", X.shape)
Y = df.iloc[:,len(df.columns)-1].values
print("Y = ",Y.shape)


# Feature scaling

# In[20]:


from sklearn.preprocessing import StandardScaler
sc = StandardScaler()
X = sc.fit_transform(X)
#print(X)
print(X.shape)
print(Y.shape)
Y_original = Y


# Handling oversampling with SMOTE

# In[9]:


from imblearn.over_sampling import SMOTE

oversample = SMOTE()
X, Y = oversample.fit_resample(X, Y)
print(X.shape)
print(Y.shape)
#print(target.shape)


# In[10]:


target = pd.DataFrame(Y)  #This is the target (Y) after oversampling


# Reduced features using ANOVA and Random forest

# In[30]:


from sklearn.feature_selection import f_classif
import numpy as np
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression


print("Total number of features:", X.shape[1])
a = round(0.85*(X.shape[1]))
print("85% of total feature", a)
X_new = SelectKBest(f_classif, k=a).fit_transform(X, Y)
#df=pd.DataFrame(X_new)
fvals, pvals = f_classif(X_new, Y)
to_remove = pvals >= (0.05/X_new.shape[1])
X_anova = np.delete(X_new, obj=to_remove, axis=1)
print("Number of selected features using ANOVA F-test:", X_anova.shape[1])
#print(X_anova)

print("\n")
print("Only random forest")
sel = SelectFromModel(RandomForestClassifier(n_estimators = 100))
sel.fit(X, Y)
rf_sel_features=sel.get_support()
X_rf = np.delete(X, obj=~np.array(rf_sel_features), axis=1)
print("Number of selected features using RF:", X_rf.shape[1])

print("\n")
print("Both ANOVA and random forest")
sel = SelectFromModel(RandomForestClassifier(n_estimators = 100))
sel.fit(X_anova, Y)
anova_rf_sel_features=sel.get_support()
X_anova_rf = np.delete(X_anova, obj=~np.array(anova_rf_sel_features), axis=1)
print("Number of selected features using anova and RF:", X_anova_rf.shape[1])


# Train-test-validation split on total features (before feature selection)

# In[12]:


#import imblearn
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score, recall_score, precision_score, cohen_kappa_score, f1_score, roc_auc_score



print("train and test split on total dataset", len(X))
X_train0, X_test, y_train0, y_test = train_test_split(X, Y, test_size=0.2, random_state = 0)
print("train and validation split on training dataset", len(X_train0))
X_train, X_val, y_train, y_val = train_test_split(X_train0, y_train0, test_size=0.2, random_state = 0)
print("Training data:",X_train.shape)
print("Validation data", X_val.shape)
print("test data",X_test.shape)


# Performance of Random forest model on total features

# In[13]:


RF = RandomForestClassifier()
RF.fit(X_train, y_train)
print("on validation dataset")
y_pred = RF.predict(X_val)
print("tn, fp, fn, tp=", confusion_matrix(y_val, y_pred).ravel())
tn, fp, fn, tp = confusion_matrix(y_val, y_pred).ravel()
print("Accuracy: " + str(round((tp + tn) / (tn + fp + tp + fn), 2)))
print("Recall/Sensitivity : " + str(round(tp / (tp + fn), 2)))
print("Specificity : " + str(round(tn / (tn + fp), 2)))
print("F1 score : " + str(round(f1_score(y_val, y_pred))))
print("Cohens kappa : " + str(round(tn / (tn + fp), 2)))
print("ROC AUC : " + str(round(roc_auc_score(y_val, y_pred))))

print("\n")
print("on test dataset")
y_pred = RF.predict(X_test)
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
print("tn, fp, fn, tp=", confusion_matrix(y_test, y_pred).ravel())
print("Accuracy: " + str(round((tp + tn) / (tn + fp + tp + fn), 2)))
print("Recall/Sensitivity : " + str(round(tp / (tp + fn), 2)))
print("Specificity : " + str(round(tn / (tn + fp), 2)))
print("F1 score : " + str(round(f1_score(y_test, y_pred))))
print("Cohens kappa : " + str(round(tn / (tn + fp), 2)))
print("ROC AUC : " + str(round(roc_auc_score(y_test, y_pred))))


# Train-test-validation split on selected features

# In[14]:


print("train and test split on total dataset", len(X))
X_train0, X_test, y_train0, y_test = train_test_split(X_anova_rf, Y, test_size=0.2, random_state = 0)
print("train and validation split on training dataset", len(X_train0))
X_train, X_val, y_train, y_val = train_test_split(X_train0, y_train0, test_size=0.2, random_state = 0)
print("Training data:",X_train.shape)
print("Validation data", X_val.shape)
print("test data",X_test.shape)


# Performance of Random forest model on reduced features

# In[15]:


from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score, recall_score, precision_score, cohen_kappa_score, f1_score, roc_auc_score

RF = RandomForestClassifier()
RF.fit(X_train, y_train)
print("on validation dataset")
y_pred = RF.predict(X_val)
print("tn, fp, fn, tp=", confusion_matrix(y_val, y_pred).ravel())
tn, fp, fn, tp = confusion_matrix(y_val, y_pred).ravel()
print("Accuracy: " + str(round((tp + tn) / (tn + fp + tp + fn), 2)))
print("Recall/Sensitivity : " + str(round(tp / (tp + fn), 2)))
print("Specificity : " + str(round(tn / (tn + fp), 2)))
print("F1 score : " + str(round(f1_score(y_val, y_pred))))
print("Cohens kappa : " + str(round(tn / (tn + fp), 2)))
print("ROC AUC : " + str(round(roc_auc_score(y_val, y_pred))))

print("\n")
print("on test dataset")
y_pred = RF.predict(X_test)
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
print("tn, fp, fn, tp=", confusion_matrix(y_test, y_pred).ravel())
print("Accuracy: " + str(round((tp + tn) / (tn + fp + tp + fn), 2)))
print("Recall/Sensitivity : " + str(round(tp / (tp + fn), 2)))
print("Specificity : " + str(round(tn / (tn + fp), 2)))
print("F1 score : " + str(round(f1_score(y_test, y_pred))))
print("Cohens kappa : " + str(round(tn / (tn + fp), 2)))
print("ROC AUC : " + str(round(roc_auc_score(y_test, y_pred))))


# Array of selected features to CSV for further analysis

# In[16]:


anovaRan = pd.DataFrame(X_anova_rf)
anovaRan['is_tumor'] = target
#print(anovaRan)
#anovaRan.to_csv("450_sim_anova_RF_woSMOTE.csv")


# Common features using ANOVA and Random forest

# In[29]:


df1=df.drop(["is_tumor"], axis=1)
print("Features using ANOVA")
bestfeatures = SelectKBest(score_func=f_classif, k=X_anova.shape[1])
#df1 = pd.DataFrame(X)
anova = bestfeatures.fit(df1, Y_original)
dfscores = pd.DataFrame(anova.scores_)
dfcolumns = pd.DataFrame(df1.columns)
# concat two dataframes for better visualization
featureScores = pd.concat([dfcolumns, dfscores], axis=1)
featureScores.columns = ['fFeature', 'fScore']  # naming the dataframe columns
best_features = featureScores.nlargest(X_anova.shape[1], 'fScore')
#print(best_features)
anova_bestfeatures = best_features['fFeature'].tolist()
print(len(anova_bestfeatures))

print("Features using Random forest")
RF = RandomForestClassifier()
RF.fit(X, Y)
# feature importance
rf_importance = RF.feature_importances_

dfscores = pd.DataFrame(RF.feature_importances_)
dfcolumns = pd.DataFrame(df1.columns)
# concat two dataframes for better visualization
fs = pd.concat([dfcolumns, dfscores], axis=1)
fs.columns = ['rfFeature', 'rfScore']  # naming the dataframe columns
rf_best_feature = fs.nlargest(X_rf.shape[1], 'rfScore')
#print(rf_best_feature)
rf_bestfeatures = rf_best_feature['rfFeature'].tolist()
print(len(rf_bestfeatures))


print("Common features of ANOVA and Random forest")
#df_anova = 
RF = RandomForestClassifier()
RF.fit(X_anova, Y)
# feature importance
importance = RF.feature_importances_
dfscores = pd.DataFrame(importance)
dfcolumns = pd.DataFrame(df1.columns)
# concat two dataframes for better visualization
fs = pd.concat([dfcolumns, dfscores], axis=1)
fs.columns = ['rfFeature', 'rfScore']  # naming the dataframe columns
anovaRF_best_features = fs.nlargest(X_anova_rf.shape[1], 'rfScore')
anova_rf_bestfeatures = anovaRF_best_features['rfFeature'].tolist()
print(len(anova_rf_bestfeatures))
anova_rf_markers = pd.DataFrame (anova_rf_bestfeatures,columns=['Reduced_Markers_AnovaRF'])
#anova_rf_markers.to_csv("450_sim_AnovaRF_markers_woSMOTE.csv")

