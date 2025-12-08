# vEDA_AhR_Elbe_sediment

This file contains two pipelines: (A) Machine learning models and (B) Suspect screening analysis.

(A)	Machine learning models
A binary classification model was trained to distinguish AhR agonists from inactive molecular structures, and a regression model was trained to predict the effect concentration (EC10) of these agonists.
The workflow includes:

(1)	Data curation for training dataset and unseen dataset from our in-house data

(2)	Molecular descriptor calculation

(3)	Classification model training

(4)	Regression model training

(5)	Prediction for molecular structures

The inputs and outputs for the pipeline calculation are available in the Zenodo database under DOI: 10.5281/zenodo.17791247.

Dependencies:
Python version (3.12.12)
numpy, 1.26.4
pandas 2.2.3
rdkit, 2025.9.6
mordred 1.2.0
xgboost 3.0.2
sklearn 1.6.1
pubchempy 1.0.5

(B)	Suspect screening analysis
A suspect screening analysis workflow was established to identify AhR agonists from GC-HRMS data.
The workflow includes:

(1)	Peak picking and deconvolution by MS-DIAL software

(2)	Pre-filter GC-HRMS features

(3)	Mass spectra reference library matching

(4)	AhR activation prediction and EC10 value prediction for molecular structure candidates

(5)	Identification of filtered features 

(6)	Alignment


The inputs and outputs for the pipeline calculation are available in the Zenodo database under DOI: 10.5281/zenodo.17791247.

Dependencies:
R version (4.4.3)
mssearchr 0.2.0
tidyr 1.3.1
dplyr 1.1.4
ggplot2 3.5.2
UpsetR 1.4.0

Here we provided an example of using the trained machine learning (ML) models. Please see the Example.ipynb in ‘Example of using ML models’ folder.

•	Retrieve 3D SDF from the PubChem database using the CID number.

•	Calculate the molecular descriptor

•	Predict their AhR activity (1: active; 0: inactive)

•	Predict their log (EC10 values (μM)) values

