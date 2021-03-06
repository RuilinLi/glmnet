% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assess.glmnet.R, R/confusion.glmnet.R,
%   R/roc.glmnet.R
\name{assess.glmnet}
\alias{assess.glmnet}
\alias{confusion.glmnet}
\alias{roc.glmnet}
\title{assess performance of a 'glmnet' object using test data.}
\usage{
assess.glmnet(
  object,
  newx = NULL,
  newy,
  weights = NULL,
  family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
  ...
)

confusion.glmnet(
  object,
  newx = NULL,
  newy,
  family = c("binomial", "multinomial"),
  ...
)

roc.glmnet(object, newx = NULL, newy, ...)
}
\arguments{
\item{object}{Fitted \code{"glmnet"} or \code{"cv.glmnet"}, \code{"relaxed"}
or \code{"cv.relaxed"} object, or a matrix of predictions (for
\code{roc.glmnet} or \code{assess.glmnet}). For \code{roc.glmnet} the model
must be a 'binomial', and for \code{confusion.glmnet} must be either
'binomial' or 'multinomial'}

\item{newx}{If predictions are to made, these are the 'x' values. Required
for \code{confusion.glmnet}}

\item{newy}{required argument for all functions; the new response values}

\item{weights}{For observation weights for the test observations}

\item{family}{The family of the model, in case predictions are passed in as
'object'}

\item{...}{additional arguments to \code{predict.glmnet} when "object" is a
"glmnet" fit, and predictions must be made to produce the statistics.}
}
\value{
\code{assess.glmnet} produces a list of vectors of measures.
\code{roc.glmnet} a list of 'roc' two-column matrices, and
\code{confusion.glmnet} a list of tables. If a single prediction is
provided, or predictions are made from a CV object, the latter two drop the
list status and produce a single matrix or table.
}
\description{
Given a test set, produce summary performance measures for the glmnet
model(s)
}
\details{
\code{assess.glmnet} produces all the different performance measures
provided by \code{cv.glmnet} for each of the families. A single vector, or a
matrix of predictions can be provided, or fitted model objects or CV
objects. In the case when the predictions are still to be made, the
\code{...} arguments allow, for example, 'offsets' and other prediction
parameters such as values for 'gamma' for 'relaxed' fits.  \code{roc.glmnet}
produces for a single vector a two column matrix with columns TPR and FPR
(true positive rate and false positive rate). This object can be plotted to
produce an ROC curve. If more than one predictions are called for, then a
list of such matrices is produced.  \code{confusion.glmnet} produces a
confusion matrix tabulating the classification results. Again, a single
table or a list, with a print method.
}
\examples{

data(QuickStartExample)
set.seed(11)
train = sample(seq(length(y)),70,replace=FALSE)
fit1 = glmnet(x[train,], y[train])
assess.glmnet(fit1, newx = x[-train,], newy = y[-train])
preds = predict(fit1, newx = x[-train, ], s = c(1, 0.25))
assess.glmnet(preds, newy = y[-train], family = "gaussian")
fit1c = cv.glmnet(x, y, keep = TRUE)
fit1a = assess.glmnet(fit1c$fit.preval, newy=y,family="gaussian")
plot(fit1c$lambda, log="x",fit1a$mae,xlab="Log Lambda",ylab="Mean Absolute Error")
abline(v=fit1c$lambda.min, lty=2, col="red")
data(BinomialExample)
fit2 = glmnet(x[train,], y[train], family = "binomial")
assess.glmnet(fit2,newx = x[-train,], newy=y[-train], s=0.1)
plot(roc.glmnet(fit2, newx = x[-train,], newy=y[-train])[[10]])
fit2c = cv.glmnet(x, y, family = "binomial", keep=TRUE)
idmin = match(fit2c$lambda.min, fit2c$lambda)
plot(roc.glmnet(fit2c$fit.preval, newy = y)[[idmin]])
data(MultinomialExample)
set.seed(103)
train = sample(seq(length(y)),100,replace=FALSE)
fit3 = glmnet(x[train,], y[train], family = "multinomial")
confusion.glmnet(fit3, newx = x[-train, ], newy = y[-train], s = 0.01)
fit3c = cv.glmnet(x, y, family = "multinomial", type.measure="class", keep=TRUE)
idmin = match(fit3c$lambda.min, fit3c$lambda)
confusion.glmnet(fit3c$fit.preval, newy = y, family="multinomial")[[idmin]]

}
\seealso{
\code{cv.glmnet}, \code{glmnet.measures} and \code{vignette("relax",package="glmnet")}
}
\author{
Trevor Hastie and Rob Tibshirani\cr Maintainer: Trevor Hastie
<hastie@stanford.edu>
}
\keyword{classification}
\keyword{models}
