import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, cohen_kappa_score, roc_auc_score, confusion_matrix

### Adjustable Paths and Parameters

curr_dir = sys.path[0]

csv_path = curr_dir + '/BetaData_27K_SimpleImpute_Mean_1.csv'

OUTPUT_NAME = "27K_SimpleImpute_Mean_1"

# Testing Split
testing_split = 0.3

# Filter Sizes
FILTERS = [10, 20, 30, 20]


### Importing Data and Organizing DataFrames
# Read Full CSV
df_pos = pd.read_csv(csv_path)

# Sort DataFrame by is_tumor
df_pos = df_pos.sort_values(by=['is_tumor'])

# Split DataFrame into pos and neg DataFrames
tumor_counts = df_pos['is_tumor'].value_counts()
df_neg = df_pos.iloc[:tumor_counts[0],:]
df_pos = df_pos.iloc[tumor_counts[0]:,:]

# Randomly Shuffled Split DataFrames
df_neg = df_neg.sample(frac=1)
df_pos = df_pos.sample(frac=1)

# Remove extras to equilize number and pos and neg in both dataFrames
if tumor_counts[0] > tumor_counts[1]:
	df_extra = df_neg.iloc[tumor_counts[1]:,:]
	df_neg = df_neg.iloc[:tumor_counts[1],:]
else:
	df_extra = df_pos.iloc[tumor_counts[0]:,:]
	df_pos = df_pos.iloc[:tumor_counts[0],:]

tumor_counts = df_pos['is_tumor'].value_counts()

# Recombine into Training and Testing Dataframes
num_testing = int(tumor_counts[1] * testing_split)

df_testing = pd.concat([df_neg.iloc[:num_testing],df_pos.iloc[:num_testing],df_extra])
df_testing = df_testing.sample(frac=1)
df_testing = df_testing.reset_index()

df_training = pd.concat([df_neg.iloc[num_testing:], df_pos.iloc[num_testing:]])
df_training = df_training.sample(frac=1)
df_training = df_training.reset_index()

print("--- Training is_tumor Count ---")
print(df_training['is_tumor'].value_counts())
print("--- Testing is_tumor Count ---")
print(df_testing['is_tumor'].value_counts())

# Delete Temporary DataFrames
del df_neg
del df_pos



### Splitting Outcome from Features

df_trainingY = df_training['is_tumor'].copy()
try:
     df_trainingX = df_training.drop(columns=['index', 'Donor_Sample', 'Unnamed: 0', 'is_tumor'])
except KeyError:
     df_trainingX = df_training.drop(columns=['index', 'Unnamed: 0', 'is_tumor'])

df_testingY = df_testing['is_tumor'].copy()
try:
     df_testingX = df_testing.drop(columns=['index', 'Donor_Sample', 'Unnamed: 0', 'is_tumor'])
except KeyError:
     df_testingX = df_testing.drop(columns=['index', 'Unnamed: 0', 'is_tumor'])

print("-------")
print(df_trainingY)
print("-----")
print(df_trainingX)
print("-----")
print(df_testingY)
print("-----")
print(df_testingX)
print("-------")


# Number of Features
input_size = len(df_trainingX.columns)


### Building Sequential Model 

model = Sequential()

#Hidden Layer 1 (Input)
model.add(Dense(FILTERS[0], input_dim=input_size, activation="relu"))
model.add(Dropout(0.25))

#Hidden Layer 2
model.add(Dense(FILTERS[1], activation="relu"))
model.add(Dropout(0.25))

#Hidden Layer 3
model.add(Dense(FILTERS[2], activation="relu"))
model.add(Dropout(0.25))

#Hidden Layer 4
model.add(Dense(FILTERS[3], activation="relu"))
model.add(Dropout(0.25))

#Output Layer
model.add(Dense(1, activation="sigmoid"))
opt = Adam(learning_rate=0.001)

model.compile(loss="binary_crossentropy", optimizer=opt, metrics=['accuracy'])
print(model.summary())



# Early Stopping and Learning Rate Reduction

earlystop = EarlyStopping(patience=20)
learning_rate_reduction = ReduceLROnPlateau(monitor="val_accuracy",
                                            patience=5,
                                            verbose=1,
                                            factor=0.5,
                                            min_lr=0.00001)
callback = [earlystop, learning_rate_reduction]



### Fitting Model

epochs = 30
history = model.fit(x=df_trainingX,
                    y=df_trainingY,
                    epochs=epochs,
                    validation_split=0.3,
                    callbacks=callback)

model.save_weights(OUTPUT_NAME + sys.argv[1] + ".h5")



### Saving History Data as a CSV for graphing in Excel

history_file = "History_" + sys.argv[1] + ".csv"

df_history = pd.DataFrame()
df_history['loss'] = history.history['loss']
df_history['val_loss'] = history.history['val_loss']
df_history['accuracy'] = history.history['accuracy']
df_history['val_accuracy'] = history.history['val_accuracy']
df_history = df_history.transpose()
df_history.to_csv(history_file)


### Testing Model

predict = model.predict(df_testingX)
predict_classes = model.predict_classes(df_testingX)

predicted_prob = predict[:, 0]
predicted_prob_class = predict_classes[:, 0]

### Saving statistical data as CSV for graphing in Excel

statistics_file = "Statistics_" + sys.argv[1] + ".csv"

accuracy = accuracy_score(df_testingY, predicted_prob_class)
precision = precision_score(df_testingY, predicted_prob_class)
recall = recall_score(df_testingY, predicted_prob_class)
f1 = f1_score(df_testingY, predicted_prob_class)
cohen_kappa = cohen_kappa_score(df_testingY, predicted_prob_class)
roc_auc = roc_auc_score(df_testingY, predicted_prob_class)

df_stats = pd.DataFrame()
df_stats['Accuracy'] = [accuracy]
df_stats['Precision'] = [precision]
df_stats['Recall'] = [recall]
df_stats['F1 Score'] = [f1]
df_stats['Cohen Kappa'] = [cohen_kappa]
df_stats['ROC AUC'] = [roc_auc]
df_stats.to_csv(statistics_file)


prediction_file = "Prediction_Summary_" + sys.argv[1] + ".txt"


with open(prediction_file, 'w') as file:

	sys.stdout = file

	print("--- Prediction Summary ---")

   	 # How many predictions (Pos and Neg) were correct.
	print("Accuracy: %f" %accuracy)

   	 # Precision = T(Pos) / (T(Pos) + F(Pos))
   	 # Out of all the Predicted Positives, How many are actually positive

	print("Precision: %f" %precision)

   	 # Recall = T(Pos) / (T(Pos) + F(Neg))
   	 # Out of all the Actual Positives, How many were predicted positive

	print("Recall: %f" %recall)

    	 # Balance between Accuracy and Recall

	print("F1 Score: %f" %f1)

	print("Cohen Kappa: %f" %cohen_kappa)

	print("ROC AUC: %f" %roc_auc)

	matrix = confusion_matrix(df_testingY, predicted_prob_class)
	print(matrix)
