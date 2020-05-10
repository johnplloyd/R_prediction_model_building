
args = commandArgs(TRUE)

print("Reading arguments")

FILE.m.class = args[1]
class_i = as.numeric(args[2])
FILE.m.feat = args[3]
FILE.features = args[4]
FILE.IDs.val = args[5]
train_prop = as.numeric(args[6])
n_iter = as.numeric(args[7])
out_prefix = args[8]

#FILE.m.class = "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/Klijn2015/class.IC50.Klijn.RDS"
#class_i = 3
#FILE.m.feat = "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/Klijn2015/feat.ALL.Klijn.RDS"
#FILE.features = "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/IDs.all_feat.RDS"
#FILE.IDs.val = "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/IDs.CL.overlap.IC50.RDS"
#train_prop = 0.7
#n_iter = 3
#out_prefix = "test_REGU"

m.class <- readRDS(file = FILE.m.class)
log_class <- 'y'

m.feat <- readRDS(file = FILE.m.feat)
features <- readRDS(file = FILE.features)

IDs.val <- readRDS(file = FILE.IDs.val)

#n_train <- 125
#n_test <- 50
#n_iter <- 25

argsObj <- list( class_matrix_file = FILE.m.class, class_i = class_i,
                 feat_matrix_file = FILE.m.feat, use_feat_file = FILE.features,
                 val_IDs_file = FILE.IDs.val, train_prop = train_prop, n_iter = n_iter )

IDs.exclude <- "n"

## Pre-process class data frame
print("Processing class data frame")

m.class <- m.class[class_i]
if( tolower(log_class) == "y" ){
  m.class <- log(m.class)
}

if(anyNA(m.class)){
  NA_ind <- which( is.na(m.class[,1]) )
  m.class <- m.class[ -NA_ind, , drop = F ]
}

IDs.class_in_feat <- intersect(x = row.names(m.class), y = row.names(m.feat))
m.class <- m.class[ IDs.class_in_feat, , drop = F]

## Pre-process features
print("Processing features data frame")

m.feat <- m.feat[features]

IDs.class <- row.names(m.class)
IDs.noClass <- row.names(m.feat) [which( ! row.names(m.feat) %in% IDs.class )]

m.feat.class <- m.feat[ IDs.class, ]
m.feat.noClass <- m.feat[ IDs.noClass, ]

dim(m.feat.class)
dim(m.feat.noClass)

dim(m.feat)
nrow(m.feat.class)+nrow(m.feat.noClass)

## Pre-process train, test, validation IDs
print("Processing train, test, validation IDs")

subset_train_test <- function(vec, train_prop, n_iter){
  l_train <- list()
  l_test <- list()
  
  n_train <- floor( length(vec) * train_prop )
  
  for(i in 1:n_iter){
    train_samp <- sample( x = vec, size = n_train )
    test_samp <- setdiff(vec, train_samp)
    
    #vec.noTrain <- setdiff( x = vec, train_samp )
    #test_samp <- sample( x = vec.noTrain, size = n_test )
    
    l_train[[i]] <- train_samp
    l_test[[i]] <- test_samp
  }
  return( list( train = l_train, test = l_test ) )
}

IDs.val <- intersect( IDs.class, IDs.val )
IDs.trainTest <- setdiff(IDs.class, IDs.val)

ID_obj.trainTest <- subset_train_test(vec = IDs.trainTest, train_prop = train_prop, n_iter = n_iter)

## Grid search 
print("performing grid search")

alphas <- seq(from = 0, to = 1, by = 0.1)
lambdas <- c( 1e-3, 1e-2, 1e-1, 1e0, 1e1 )

sort_LASSO_weights <- function (fit_coef){
  nonzero_ind <- which(fit_coef != 0)
  if(length(nonzero_ind) == 1){
    return(NA)
  }else{
    fit_coef.weights <- fit_coef[ nonzero_ind,  ]
    intercept_ind <- which(names(fit_coef.weights) == "(Intercept)")
    fit_coef.weights <- fit_coef.weights[-intercept_ind]
    pos_weights.names <- labels(which( fit_coef.weights > 0))
    neg_weights.names <- labels(which( fit_coef.weights < 0))
    fit_coef.weights.abs.sort <- sort(abs(fit_coef.weights), decreasing = TRUE)
    
    df.weights <- data.frame(rank = 1:length(fit_coef.weights.abs.sort), weight=fit_coef.weights.abs.sort, sign=rep(NA, length(fit_coef.weights.abs.sort)), row.names = labels(fit_coef.weights.abs.sort))
    df.weights[pos_weights.names, "sign"] = "+"
    df.weights[neg_weights.names, "sign"] = "-"
    #fit_coef.weights.order_names <- names(fit_coef.weights.abs)
    return (df.weights)
  }
}

glmnet_grid_search.1 <- function(df_class, df_feat, ID_obj, alphas, lambdas){
  library(glmnet)
  
  dim(df_class)
  dim(df_feat)
  
  predictions_object <- list()
  weights_object <- list()
  for( i in 1:length(ID_obj[[1]]) ){
    print(i)
    IDs.train <- ID_obj[[1]][[i]]
    IDs.test <- ID_obj[[2]][[i]]
    #intersect(IDs.train, IDs.test)
    
    df_class.train <- df_class[ IDs.train, , drop = F ]
    df_class.test <- df_class[ IDs.test, , drop = F ]
    df_feat.train <- df_feat[ IDs.train, , drop = F ]
    df_feat.test <- df_feat[ IDs.test, , drop = F ]
    
    y = df_class.train[,1]
    x = as.matrix(df_feat.train)
    
    newx <- as.matrix(df_feat.test)
    
    for( j in 1:length(alphas) ){
      alpha <- alphas[j]
      for(k in 1:length(lambdas) ){
        lambda <- lambdas[k]
        
        model <- glmnet(x, y, alpha = alpha, lambda = lambda)
        model.coef  <- predict(model, type = 'coefficients', s = lambda)
        weights <- sort_LASSO_weights(model.coef)
        
        fit <- predict(model, s = lambda, newx = newx)[,1]
        
        df_pred <- data.frame( actual = df_class.test[,1], pred = fit, row.names = row.names(df_class.test) )
        
        grid_position <- paste( i, alpha, lambda, sep = "_" )
        predictions_object[[ length(predictions_object)+1 ]] <- df_pred
        names(predictions_object)[ length(predictions_object) ] <- grid_position
        
        weights_object[[ length(weights_object)+1 ]] <- weights
        names(weights_object)[ length(weights_object) ] <- grid_position
      }
    }
  }
  
  returnObj <- list( predObj = predictions_object, weightObj = weights_object )
  
  return(returnObj)
}

glmnet_grid_predictions_weights <- glmnet_grid_search.1(df_class = m.class, df_feat = m.feat.class, ID_obj = ID_obj.trainTest, alphas = alphas, lambdas = lambdas)
glmnet_grid_predictions <- glmnet_grid_predictions_weights[[1]]
glmnet_grid_weights <- glmnet_grid_predictions_weights[[2]]

outNm.parameters <- paste(out_prefix, ".parameters.RDS", sep = "")
outNm.trainTest <- paste(out_prefix, ".trainTest.RDS", sep = "")
outNm.predictions <- paste(out_prefix, ".predictions.RDS", sep = "")
outNm.weights <- paste(out_prefix, ".featSelWeights.RDS", sep = "")

saveRDS( object = argsObj, file = outNm.parameters )
saveRDS( object = ID_obj.trainTest, file = outNm.trainTest )
saveRDS( object = glmnet_grid_predictions, file = outNm.predictions )
saveRDS( object = glmnet_grid_weights, file = outNm.weights )
