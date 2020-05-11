# R_prediction_models

This repository includes scripts for running 4 algorithms to build prediction models in R: logistic regression, regularized regression, and both classification (binary) and regression random forest.

* The scripts assume the class and feature data have been pre-processed and are saved as RDS files.
* The model-building is done with a valdiation set that is held out for each model. The validation set is used to report model performance.
* The scripts make use of a bagging approach where multiple models are developed while selecting a random subset of training instances. The final prediction score for a validation instance is reported as the mean of model repetitions.
* Each script performs a parameter sweep to identify optimal parameters using cross-validation within the training set. Regularized regression tunes lamba (beta penalization) and alpha (LASSO vs. ridge regression vs. elastic nets) parameters. Logistic and both random forest algorithms utilize LASSO as a feature selection tool prior to model-building. The lambda beta penalization for feature selection is tuned.

Each script takes the same 8 inputs:
1) Location of RDS file with class data matrix (scripts were developed for data with multiple response values per instance)
2) Colunm index of class data matrix to use: Integer
3) Location of RDS file with feature data matrix
4) Location of RDS file with vector of feature IDs to include in prediction model
5) Location of RDS file with vector of instance IDs to hold out as validation set
6) Proportion of training instances (i.e. non-validation) to include in prediction model: [0, 1]
7) \# of repetitions of prediction models to build, each with a random subset of the training data: Integer
8) Prefix character string for output files

The scripts produce four output files:
1) [prefix].prediction.RDS: Object with prediction on the validation set
2) [prefix].featSelWeightsRDS: Object with feature selection weights (penalized betas)
3) [prefix].parameters.RDS: Object with the parameters utilized to train the model, for reporting purposes.
4) [prefix].trainTest.RDS: Object with training and testing instances for each repetition of model-building, for reporting purposes.
