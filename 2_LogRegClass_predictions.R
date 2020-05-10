
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

#FILE.m.class <- "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/Klijn2015/class.IC50.Klijn.RDS"
#class_i <- 3
#FILE.m.feat <- "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/Klijn2015/feat.ALL.Klijn.RDS"
#FILE.features <- "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/IDs.all_feat.RDS"
#FILE.IDs.val <- "C:/Users/johnp/OneDrive/_Projects/1_Z_pancancer_biomarkers/data/IDs.CL.overlap.IC50.RDS"
#train_prop <- 0.75
#n_iter <- 3
#out_prefix <- "test"

m.class <- readRDS(file = FILE.m.class)
log_class <- 'y'

threshold <- 1
log_thresh <- 'y'
if(log_thresh == "y"){
  threshold <- log(threshold)
}

m.feat <- readRDS(file = FILE.m.feat)
features <- readRDS(file = FILE.features)

IDs.val <- readRDS(file = FILE.IDs.val)

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

#head(m.class)
pos_ind <- which(m.class <= threshold)
neg_ind <- which(m.class > threshold)
#m.class[ pos_ind, ]
#m.class[ neg_ind, ]
m.class[ pos_ind, ] <- 1
m.class[ neg_ind, ] <- 0

## Pre-process features
print("Processing features data frame")

m.feat <- m.feat[ features]

IDs.class <- row.names(m.class)
IDs.noClass <- row.names(m.feat) [which( ! row.names(m.feat) %in% IDs.class )]

m.feat.class <- m.feat[ IDs.class, ]
m.feat.noClass <- m.feat[ IDs.noClass, ]

dim(m.feat.class)
dim(m.feat.noClass)

dim(m.feat)
nrow(m.feat.class)+nrow(m.feat.noClass)

dim(m.class)
dim(m.feat.class)
all.equal( row.names(m.class), row.names(m.feat.class) )

## Pre-process train, test, validation IDs
print("Processing train, test, validation IDs")

subset_train_test <- function(vec, df_class, train_prop, n_iter){
  
  nms.pos <- row.names( df_class ) [ which(df_class == 1) ]
  nms.neg <- row.names( df_class ) [ which(df_class == 0) ]
  
  vec.pos <- intersect( vec, nms.pos )
  vec.neg <- intersect( vec, nms.neg )
  
  min_class <- min( c( length(vec.pos), length(vec.neg) ) )
  n_train <- floor( min_class * train_prop )
  n_test <- min_class-n_train-1
  
  l_train <- list()
  l_test <- list()
  for(i in 1:n_iter){
    train_samp.pos <- sample( x = vec.pos, size = n_train )
    vec.pos.noTrain <- setdiff( x = vec.pos, y = train_samp.pos )
    test_samp.pos <- sample( x = vec.pos.noTrain, size = n_test )
    
    intersect(train_samp.pos, vec.pos.noTrain)
    intersect(train_samp.pos, test_samp.pos)
    intersect(vec.pos.noTrain, test_samp.pos)
    
    train_samp.neg <- sample( x = vec.neg, size = n_train )
    vec.neg.noTrain <- setdiff( x = vec.neg, y = train_samp.neg )
    test_samp.neg <- sample( x = vec.neg.noTrain, size = n_test )
    
    intersect(train_samp.neg, vec.neg.noTrain)
    intersect(train_samp.neg, test_samp.neg)
    intersect(vec.neg.noTrain, test_samp.neg)
    
    train_samp <- c( train_samp.pos, train_samp.neg )
    test_samp <- c( test_samp.pos, test_samp.neg )
    
    intersect(train_samp, test_samp)
    
    df_class[ train_samp, ]
    df_class[ test_samp, ]
    table( df_class[ train_samp, ] )
    table( df_class[ test_samp, ] )
    
    l_train[[i]] <- train_samp
    l_test[[i]] <- test_samp
  }
  return( list( train = l_train, test = l_test ) )
}

IDs.val <- intersect( IDs.class, IDs.val )
IDs.trainTest <- setdiff(IDs.class, IDs.val)

ID_obj.trainTest <- subset_train_test(vec = IDs.trainTest, df_class = m.class, train_prop = train_prop, n_iter = n_iter)

## Grid search 
print("performing grid search")

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

LogRegClass_grid_search.1 <- function(df_class, df_feat, ID_obj, lambdas){
  library(glmnet)
  library(randomForest)
  
  dim(df_class)
  dim(df_feat)
  head(df_class)
  
  predictions_object <- list()
  weights_object <- list()
  for( i in 1:length(ID_obj[[1]]) ){
    print(i)
    IDs.train <- ID_obj[[1]][[i]]
    IDs.test <- ID_obj[[2]][[i]]
    intersect(IDs.train, IDs.test)
    length(IDs.train)
    length(IDs.test)
    
    df_class.train <- df_class[ IDs.train, , drop = F ]
    df_class.test <- df_class[ IDs.test, , drop = F ]
    df_feat.train <- df_feat[ IDs.train, , drop = F ]
    df_feat.test <- df_feat[ IDs.test, , drop = F ]
    
    y <- df_class.train[,1]
    #y <- as.factor(df_class.train[,1])
    #y.glmnet <- as.numeric(as.character(y))
    
    x <- as.matrix(df_feat.train) # for feature selection
    newx <- df_feat.test
    
    data <- cbind(y, df_feat.train)
    
    for( j in 1:length(lambdas) ){
      if(j == 0){
        # NO FEATURE SELECTION
        lambda <- "noFeatSel"
        print(lambda)
        weights <- NA
      }else{
        lambda <- lambdas[j]
        print(lambda)
        
        model <- glmnet(x, y, alpha = 1, lambda = lambda)
        model.coef  <- predict(model, type = 'coefficients', s = lambda)
        weights <- sort_LASSO_weights(model.coef)
        
        head(weights)
        tail(weights)
        SNP_ind <- grep( pattern = "_SNP", x = row.names(weights) )
        CNV_ind <- grep( pattern = "_CNV", x = row.names(weights) )
        head(SNP_ind)
        head(CNV_ind)
        
        df_feat.train.sel <- df_feat.train[ row.names(weights) ]
        df_feat.test.sel <- df_feat.test[ row.names(weights) ]
        
        data <- cbind(y, df_feat.train.sel)
        newx <- df_feat.test.sel
      }
      
      logit_model <- glm( formula = y ~ ., data = data, family = binomial(link="logit") )
      fit <- predict(object = logit_model, newdata = newx, type = "response" )
      df_pred <- data.frame( actual = df_class.test[,1], pred = fit, row.names = row.names(df_class.test) )
      #boxplot( df_pred[,2] ~ df_pred[,1] )
      
      grid_position <- paste( i, lambda, sep = "_" )
      
      predictions_object[[ length(predictions_object)+1 ]] <- df_pred
      names(predictions_object)[ length(predictions_object) ] <- grid_position
      
      weights_object[[ length(weights_object)+1 ]] <- weights
      names(weights_object)[ length(weights_object) ] <- grid_position
    }
  }
  returnObj <- list( predObj = predictions_object, weightObj = weights_object )
  return(returnObj)
}


lambdas <- c(  1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)

LogReg_grid_predictions_weights <- LogRegClass_grid_search.1(df_class = m.class, df_feat = m.feat.class, ID_obj = ID_obj.trainTest, lambdas = lambdas)
LogReg_grid_predictions <- LogReg_grid_predictions_weights[[1]]
LogReg_grid_weights <- LogReg_grid_predictions_weights[[2]]

outNm.parameters <- paste(out_prefix, ".parameters.RDS", sep = "")
outNm.trainTest <- paste(out_prefix, ".trainTest.RDS", sep = "")
outNm.predictions <- paste(out_prefix, ".predictions.RDS", sep = "")
outNm.weights <- paste(out_prefix, ".featSelWeights.RDS", sep = "")

saveRDS( object = argsObj, file = outNm.parameters )
saveRDS( object = ID_obj.trainTest, file = outNm.trainTest )
saveRDS( object = LogReg_grid_predictions, file = outNm.predictions )
saveRDS( object = LogReg_grid_weights, file = outNm.weights )
