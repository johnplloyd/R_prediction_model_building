# R_prediction_models

This repository includes scripts for running 4 algorithms to build prediction models in R: logistic regression, regularized regression, and both classification (binary) and regression random forest.

* The scripts assume the class and feature data have been pre-processed and are saved as RDS files.
* The model-building is done with a valdiation set that is held out for each model. The validation set is used to report model performance.
* The scripts make use of a bagging approach where multiple models are developed while selecting a random subset of training instances. The final performance for validation instances are reported as the mean of the model repetitions.

Each script takes the same 8 inputs:
1) RDS file with class data matrix
2) Index of class data matrix to use
3) RDS file with feature data matrix
4) RDS file with vector of feature IDs to include in prediction model
5) RDS file with vector of instance IDs to hold out as validation set
6) Proportion of training instances (i.e. non-validation) to include in prediction model: [0, 1]
7) \# of repetitions of prediction models to build, each with a random subset of the training data: Integer
8) Prefix character string for output files
