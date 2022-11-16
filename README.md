# Deep_Learning_Methylation
This research explores machine learning and deep learning approach for analyzing DNA methylation data

Please cite the following paper associated with this reseacrch: 

Gomes, Rahul, Nijhum Paul, Nichol He, Aaron Francis Huber, and Rick J. Jansen. 2022. "Application of Feature Selection and Deep Learning for Cancer Prediction Using DNA Methylation Markers" Genes 13, no. 9: 1557. https://doi.org/10.3390/genes13091557

Author: Nijhum Paul

-------final_feature_analysis.py------

ANOVA and Random forest methods on DNA methylation dataset to get a reduced set of features/CPG markers. 

Python version 3.9.7

Input: CSV file containing CpG markers. Columns represent markers and each row represents each sample. is_tumor column is the last column that represents if a sample is tumor(1) or normal(0).
BetaData_27K_SimpleImpute_Mean_1.csv and BetaData_450K_SimpleImpute_Zero_1.

Output: 4 sets of reduced features: 1) 27_sim_AnovaRF_markers_woSMOTE, 2) 27_sim_AnovaRF_markers_wSMOTE, 3) 450K_sim_AnovaRF_markers_woSMOTE, 4) 450K_sim_AnovaRF_markers_wSMOTE.

Description: First, a csv file with methylation value is imported as a dataframe (df) using pandas library. The first column 'Donor_Sample' does not add any meaning to the model, so I dropped it. 
The df is split into features (X) and target arrays (Y) for analysis. SMOTE is applied to features for balancing the data using python's imblearn package. Then the value of target column 
is updated (target) as number of samples are different now after SMOTE application. 
ANOVA and Random forest models are applied both separately (only ANOVA, only RF) and together (ANOVA_RF) to get the reduced setS of features. This process was done by using Scikit-learn
package and it is done on both the balanced and imbalanced dataset.
Dataset is divided into 3 sets (train, validation and test) to verify the importance of reduced features on Random forest model. Both arrays of total features and reduced features are 
applied to Random forest model to see if reduced features perform well in classification. 
I also calculated number of common features in the list of only ANOVA and only RF. 
Finally, array of of reduced features (methylation value) and list of reduced features (CPG markers) are converted to a csv file for further analysis.

------annotation_final.py-----

In this file, I compared the closest genes associated with reduced markers with the cancer related genes collected from several resources (COSMIC, TSGene database, etc).

Python version 3.9.7

Input: CSV file containing sets of genes in column 'Gene_Symbol'.

Output: Common genes between the resources and our reduced sets of genes.

Description: Genes from COSMIC, TSGene, etc were collected and python's intersection method was used to find common genes among different sets.


------EA.R-------

Performed GSEA analysis on the closest genes associated with overall markers before feature reduction and reduced markers.

R version 4.2.0 

Input: CSV file containing sets of genes in column 'Gene_Symbol'.

Output: GSEA figures for all 6 sets of genes.

Description: R's TCGAbiolinks package is used for GSEA analysis. TCGAanalyze_EAcomplete function is used on sets of genes to identify classes of genes or proteins that are over-represented using annotations for that gene set. Then 
TCGAvisualize_EAbarplot was used to show canonical pathways significantly overrepresented (enriched) by the DEGs (differentially expressed genes) with the number of genes for the main categories of three ontologies 
(GO:biological process, GO:cellular component, and GO:molecular function, respectively).


#########################################################################################################################

Author: NICHOL Y. HE
Purpose: Algorithm for the general preprocessing of methylation data, building and training of a simple deep learning model, 
		running a prediction using the built model on a subset of methylation data, and reporting statistics and results.

Files:
methylation_dl_model.py
- Input: .csv file containing methylation data. 
	Eg. BetaData_27K_SimpleImpute_Mean_1.csv
- Output: 
	Histoy.csv - Contains training history data including accuracy, loss, validation accuracy, validation with 1 entry/epoch
	Statistics.csv - Contains statistical results in the form of accuracy, precision, Recall, F1 SCore, Cohen Kappa, and ROC AUC
	Prediction Summary.txt - Contains statistical results in text form along with Confusion Matrix 
	Output.h5 - Trained model weights. Name is modifiable 

- Process:
	The preprocessing revolves around obtaining 70% of the samples from the limiting outcome and matching an equal number of 
	  samples of the excess outcome to generate a training subset containing equal numbers of postive and negative tumor outcomes.
	The remaining 30% of the samples of the limiting outcome is combined with all the rest of the excess outcome to generate a 
	  testing subset for prediction.

	This process is to ensure the training is performed on a balanced dataset. Retaining 30% of the limiting outcome for prediction
	ensures there is some amount of each outcome to predict on. 

	Eg. 10 samples have tumors. 13 samples have no tumors. 
		7 samples w/ tumors is randomly selected and 7 samples w/o tumors is randomly selected to form a training subset.
		Reamining 3 samples w/ tumors is combined with the remaining 6 samples w/o tumors to form a testing subset. 

	The model is a simple sequential model built using the Tensorflow.Keras library. There are four fully-connected layers, each 
		followed by a 25% dropout. The final layer generates a binary output. 

	The model history is saved. Prediction is ran on the training dataset to generate predictions and those outputs are used 
		to calculate prediction statistics, which are then saved to files.

run_methylation_trials.sh
- SLURM Script to run methylation_dl_model.py
	Utilizes array=0-4 to run the algorithm 5 times, generating 5 independant trials.
	Outputs can be summarized manually. 

wide_to_long.ipynb
- Jupyter Notebook used to convert wide formatted data to long formatted data
- Input: Wide format Marker Metadata in .csv format
- Output: Long format Marker Metadata in .csv format


	


