\name{pgbart_predict}
\title{Make Predictions Using Bayesian Additive Regression Trees}
\alias{pgbart_predict}
\description{
Make predictions for a new test data set after building a model using trainng data by function \code{\link{pgbart_train}}.}
\usage{
pgbart_predict(x.test, model)}

\arguments{
\item{x.test}{Explanatory variables for test (out of sample) data.\cr
   Should have same structure as \code{x.train} in \code{pgbart_train}.}

\item{model}{The path to save the model file as specified in \code{pgbart_train}.}
}
\details{
   PGBART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case
   and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is a row from the test data (x.test).
}
\value{
   \code{pgbart_predict} returns a list assigned class \sQuote{pgbart}.
   In the numeric \eqn{y} case, the list has components:

   \item{yhat.test}{
   A matrix with (ndpost/keepevery) rows and nrow(x.test) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.test.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.test. Burn-in is dropped.
   }
   \item{yhat.test.mean}{Test data fits = mean of yhat.test columns. Only exists when \eqn{y} is not binary.}

   In the binary \eqn{y} case, the returned list has the components
   yhat.test and binaryOffset.

   Note that in the binary \eqn{y} case, yhat.test is
   \eqn{f(x)} + binaryOffset.  If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
}
\references{
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298.

Lakshminarayanan B, Roy D, Teh Y W. (2015)
   Particle Gibbs for Bayesian Additive Regression Trees
   \emph{Artificial Intelligence and Statistics}, 553-561.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.
   \emph{Advances in Neural Information Processing Systems} \bold{19},
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

Friedman, J.H. (1991)
   Multivariate Adaptive Regression Splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
        
Breiman, L. (1996)
   Bias, Variance, and Arcing Classifiers. 
   \emph{Tech. Rep.} \bold{460}, 
   Statistics Department, University of California, Berkeley, CA, USA.        
}

\seealso{
\code{\link{pgbart_train}}, \code{\link{pdpgbart}}
}
\examples{
##Example 1: simulated continuous outcome data (example from section 4.3 of Friedman's MARS paper)
f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x = matrix(runif(n*10), n, 10) 
Ey = f(x)
y = Ey+sigma*rnorm(n)
model_path = file.path(tempdir(),'pgbart.model')
pgbartFit = pgbart_train(x[1:(n*.75),], y[1:(n*.75)], 
                         model=model_path,
                         ndpost=200, ntree=5, usepg=TRUE)

pgbartPredict = pgbart_predict(x[(n*.75+1):n,], model=model_path)

cor(pgbartPredict$yhat.test.mean, y[(n*.75+1):n])

##Example 2: simulated binary outcome data (two normal example from Breiman) 
f <- function (n, d = 20) 
{
  x <- matrix(0, nrow = n, ncol = d)
  c1 <- sample.int(n, n/2)
  c2 <- (1:n)[-c1]
  a <- 2/sqrt(d)
  x[c1, ] <- matrix(rnorm(n = d * length(c1), mean = -a), ncol = d)
  x[c2, ] <- matrix(rnorm(n = d * length(c2), mean = a), ncol = d)
  
  x.train <- x
  y.train <- rep(0, n)
  y.train[c2] <- 1
  list(x.train=x.train, y.train=as.factor(y.train))
}


set.seed(99)
n <- 200
train <- f(n)
model_path = file.path(tempdir(),'pgbart.model')
pgbartFit = pgbart_train(train$x.train[1:(n*.75),], train$y.train[1:(n*.75)], 
                         model=model_path, ndpost=200, ntree=5, usepg=TRUE)
pgbartPredict = pgbart_predict(train$x.train[(n*.75+1):n,], model=model_path)
class.pred = ifelse(colMeans(apply(pgbartPredict$yhat.test, 2, pnorm)) <= 0.5, 0, 1)
table(class.pred, train$y.train[(n*.75+1):n])
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
