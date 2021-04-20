import sys, os, re, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import rpy2.robjects as ro
import xgboost as xgb
from rpy2.robjects.vectors import IntVector, FloatVector
from RTassignment import getRT


def loess():
    rstring = """
    loess.as = function(x, y, degree = 1, criterion="aicc", family="gaussian", user.span=NULL, plot=FALSE, ...) {

        criterion <- match.arg(criterion)
        family <- match.arg(family)
        x <- as.matrix(x)

        if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
        if (!is.numeric(x)) stop("argument 'x' must be numeric!")
        if (!is.numeric(y)) stop("argument 'y' must be numeric!")
        if (any(is.na(x))) stop("'x' contains missing values!")
        if (any(is.na(y))) stop("'y' contains missing values!")
        if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) 
            stop("argument 'user.span' must be a numerical number!")
        if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
        if(length(y) < 3) stop("not enough observations!")

        data.bind <- data.frame(x=x, y=y)
        if (ncol(x) == 1) {
            names(data.bind) <- c("x", "y")
        } else { names(data.bind) <- c("x1", "x2", "y") }

        opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
            as.crit <- function (x) {
                span <- x$pars$span
                traceL <- x$trace.hat
                sigma2 <- sum(x$residuals^2 ) / (x$n-1)
                aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
                gcv <- x$n*sigma2 / (x$n-traceL)^2
                result <- list(span=span, aicc=aicc, gcv=gcv)
                return(result)
            }
            criterion <- match.arg(criterion)
            fn <- function(span) {
                mod <- update(model, span=span)
                as.crit(mod)[[criterion]]
            }
            result <- optimize(fn, span.range)
            return(list(span=result$minimum, criterion=result$objective))
        }

        control = loess.control(surface = "direct")
        if (ncol(x)==1) {
            if (is.null(user.span)) {
                fit0 <- loess(y ~ x, degree=degree, family=family, data=data.bind, control=control, ...)
                span1 <- opt.span(fit0, criterion=criterion)$span
            } else {
                span1 <- user.span
            }		
            fit <- loess(y ~ x, degree=degree, span=span1, family=family, data=data.bind, control=control, ...)
        } else {
            if (is.null(user.span)) {
                fit0 <- loess(y ~ x1 + x2, degree=degree,family=family, data.bind, control=control, ...)
                span1 <- opt.span(fit0, criterion=criterion)$span
            } else {
                span1 <- user.span
            }		
            fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family=family, data=data.bind, control=control...)
        }
        return(fit)
    }
    """
    return ro.r(rstring)


# inputDf = pd.read_pickle("df_key_run_RT.pickle")
mzXML = [r"../Data/FTLD_Batch2_F50.mzXML", r"../Data/FTLD_Batch2_F51.mzXML", r"../Data/FTLD_Batch2_F52.mzXML"]
idTxt = "../Data/ID.txt"
# df = getRT(mzXML, idTxt)
df = pd.read_csv("rtAlignment.txt", sep="\t")

# Find the keys fully identified across the runs
    # To-do: almost fully identified keys may be used

# LOESS
res = df.copy()
rLoess = loess()
rPredict = ro.r("predict")
colN = [c for c in df.columns if c.endswith("nPSMs")]
colRt = df.columns.drop(colN)
colRt = colRt.drop("key")
fullRt = df[df[colRt].isna().sum(axis=1) == 0][colRt]
fullN = df[df[colRt].isna().sum(axis=1) == 0][colN]
y = (fullRt * fullN.values).sum(axis=1) / fullN.sum(axis=1) # Reference as the weighted mean of RTs
for col in colRt:
    x = fullRt[col]
    mod = rLoess(FloatVector(x), FloatVector(y))
    res[col] = rPredict(mod, FloatVector(df[col]))
print()


# mod = rLoess(FloatVector(compRt), FloatVector(refRt))
# compRtHat = rPredict(mod, FloatVector(compRt))
# plt.plot(refRt, compRt, '.')
# plt.plot(refRt, compRtHat, '+')
# plt.grid()
# plt.show()

