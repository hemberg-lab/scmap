#' Run support vector machines (\code{SVM}) prediction
#'
#' Train an \code{SVM} classifier on a training dataset (\code{train}) and then
#' classify a study dataset (\code{study}) using the classifier.
#'
#' @param train training data.frame with colnames, corresponding to training labels
#' @param study study data.frame
#' @param kern kernel to be used with SVM
#' @return classification of the study dataset
#'
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @export
# support_vector_machines <- function(train, study, kern = "linear") {
#     train <- t(train)
#     labs <- factor(rownames(train))
#     rownames(train) <- NULL
#     model <- tryCatch(e1071::svm(train, labs, kernel = kern), error = function(cond) return(NA))
#     pred <- stats::predict(model, t(study))
#     return(pred = pred)
# }

#' Run random forest prediction
#'
#' Create a random forest calssifier based on a training dataset (\code{train}) and then
#' classify a study dataset (\code{study}) using the classifier.
#'
#' @param train training data.frame with colnames, corresponding to training labels
#' @param study study data.frame
#' @param ntree number of trees to be used in random forest
#' @return classification of the study dataset
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export
# random_forest <- function(train, study, ntree = 50) {
#     train <- t(train)
#     train <- as.data.frame(train)
#     y <- as.factor(rownames(train))
#     study <- t(study)
#     study <- as.data.frame(study)
#     rownames(train) <- NULL
#     rownames(study) <- NULL
#     train_rf <- randomForest::randomForest(x = train, y = y, ntree = ntree)
#     Prediction <- stats::predict(train_rf, study)
#     return(Prediction)
# }
