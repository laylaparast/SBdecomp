\name{sbdecomp}
\alias{sbdecomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection Bias Decomposition
}
\description{
This function decomposes the estimated selection bias to quantify what proportion of the estimated selection bias is explained by each observed confounder used in the propensity score model when estimating propensity score weighted treatment effects. The function offers two approaches - confounder inclusion or removal, and offers two estimation approaches - parametric or nonparametric.
}
\usage{
sbdecomp(outcome, treatment, confounders, data=NULL, type = "inclusion", estimation 
= "parametric", Bonly = T, balance = T, n.trees = 20000, interaction.depth = 4, 
shrinkage = 0.005, verbose = FALSE, stop.method = c("es.max"), cv.folds = 0,
standard.error = F, boot.rep=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outcome}{
Continuous outcome; numeric vector or name of the variable in supplied data
}
  \item{treatment}{
Binary treatment indicator, 0 or 1; numeric vector or name of the variable in supplied data
}
  \item{confounders}{
a data frame of confounders with names or a vector containing the names of the confounders in supplied data (each confounder must be either numeric or factor)
}
  \item{data}{
data containing all variables; required if names provided in outcome, treatment, and/or confounders arguments above
}
  \item{type}{
The type of method for selection bias deocomposition; options are inclusion or removal, default is inclusion
}
  \item{estimation}{
The type of estimation that should be used; options are parametric (based on logistic regression) or nonparametric (based on generalized boosted models [GBM]), default is parametric
}
  \item{Bonly}{
TRUE or FALSE; whether the user only wants the quantities B which are the proportions of the selection bias explained by each confounder, default is TRUE
}
  \item{balance}{
%%     ~~Describe \code{balance} here~~
TRUE or FALSE; whether the user wants balance information, default is TRUE
}
  \item{n.trees}{
%%     ~~Describe \code{n.trees} here~~
The n.trees for the nonparametric GBM approach, if using nonparametric estimation; default is 20000

}
  \item{interaction.depth}{
%%     ~~Describe \code{interaction.depth} here~~
The interaction depth for the nonparametric GBM approach, if using nonparametric estimation; default is 4
}
  \item{shrinkage}{
%%     ~~Describe \code{shrinkage} here~~
The shrinkage for the nonparametric GBM approach, if using nonparametric estimation; default is 0.005
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
TRUE or FALSE  for the nonparametric GBM approach, if using nonparametric estimation; default is FALSE
}
  \item{stop.method}{
%%     ~~Describe \code{stop.method} here~~
The stopping method  for the nonparametric GBM approach, if using nonparametric estimation; default is "es.max"
}
  \item{cv.folds}{
%%     ~~Describe \code{cv.folds} here~~
The cross-validation folds for the nonparametric GBM approach, if using nonparametric estimation; default is 0
}
  \item{standard.error}{
TRUE or FALSE; whether the user wants standard error estimates for the resulting estimates, obtained using boostrapping; default is FALSE; note that if using nonparametric estimation, obtaining the standard error estimates will be very time-intensive
}
  \item{boot.rep}{
The number of bootstrap replications to be used when calculating the standard error; default is 500
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Details can be found in: Parast, L and Griffin, BA (2020). "Quantifying the Bias due to Observed Individual Confounders in Causal Treatment Effect Estimates". Statistics in Medicine, In press.  Please email parast@rand.org for a current version. 
}
\value{
A list is returned:
 \item{delta.naive}{Estimated naive treatment effect}
 \item{p.value.delta.naive}{P-value for estimated naive treatment effect}
 \item{conf.int.delta.naive}{95\% confidence interval for estimated naive treatment effect}
 \item{delta.fully.adjusted}{Estimated treatment effect after adjusting for all variables}
 \item{p.value.delta.fully.adjusted}{P-value for estimated treatment effect after adjusting for all variables}
 \item{conf.int.delta.fully.adjusted}{95\% confidence interval for estimated treatment effect after adjusting for all variables}
  \item{B}{The proportion of the selection bias explained by each variable}
   \item{estimated.selection.bias}{The estimated selection bias, only if Bonly = FALSE}
    \item{lambda}{The defined selection bias movement when including/removing a variable, only if Bonly = FALSE}
    \item{delta.each}{The adjusted delta including only each variable or dropping each variable (depending on type), only if Bonly = FALSE}
     \item{balance.naive.mean}{mean of the standardized effect size differences across all variables between the two groups;  only if balance = TRUE}
      \item{balance.naive.max}{max of the standardized effect size differences across all variables between the two groups;  only if balance = TRUE}
       \item{balance.fully.adjusted.mean}{mean of the standardized effect size differences across all variables between the two groups, after weighting with propensity score weights calculated using all variables;  only if balance = TRUE}
        \item{balance.fully.adjusted.max}{max of the standardized effect size differences across all variables between the two groups, after weighting with propensity score weights calculated using all variables;  only if balance = TRUE}
         \item{balance.mean}{mean of the standardized effect size differences across all variables in the relevant propensity score model between the two groups, after weighting with relevant propensity score weights; if using type = "inclusion" and variable included is a factor, output will be the mean effect size difference across all levels;  only if balance = TRUE}
          \item{balance.max}{max of the standardized effect size differences across all variables in the relevant propensity score model between the two groups, after weighting with relevant propensity score weights;  if using type = "inclusion" and variable included is a factor, output will be the max effect size difference across all levels; only if balance = TRUE}
  \item{B.standard.error}{The estimated standard error for the estimated proportion of the selection bias explained by each variable; only if standard.error = TRUE}
  \item{lambda.standard.error}{The estimated standard error for the estimated selection bias movement when including/removing a variable; only if standard.error = TRUE }
}
\references{
%% ~put references to the literature/web site here ~
Parast, L and Griffin, BA (2020). "Quantifying the Bias due to Observed Individual Confounders in Causal Treatment Effect Estimates." Statistics in Medicine, In press.  
}
\author{
Layla Parast
}
\examples{
data(petsdata)
\donttest{sbdecomp(outcome = petsdata$genhealth, treatment = petsdata$gotdog, confounders = 
as.data.frame(petsdata[,c(2:13)]), type = "inclusion", estimation = "parametric",
balance = FALSE)}

\donttest{sbdecomp(outcome = petsdata$genhealth, treatment = petsdata$gotdog, confounders = 
as.data.frame(petsdata[,c(2:13)]), type = "inclusion", estimation = "parametric", 
Bonly =FALSE, balance = TRUE) }

\donttest{sbdecomp(outcome = "genhealth", treatment = "gotdog", confounders = c("age",
"ismale", "race_coll","hhsize","ownhome", "married",  "ontanf", "hhincome",
"fulltime","spouse_fulltime" ,"liveinhouse", "ruralurban"), data = petsdata, 
type = "inclusion", estimation = "parametric", Bonly =FALSE, balance = TRUE)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}
