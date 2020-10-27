\name{SBdecomp-package}
\alias{SBdecomp-package}
\docType{package}
\title{
\packageTitle{SBdecomp}
}
\description{
\packageDescription{SBdecomp}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{SBdecomp}
\packageIndices{SBdecomp}
This packge provides a function that decomposes the estimated selection bias to quantify what proportion of the estimated
 selection bias is explained by each observed confounder used in the propensity score model; the function is sbdecomp. 
 The function offers two approaches - confounder inclusion or removal, and offers two estimation approaches - parametric 
 or nonparametric. These methods allow one to identify the most important confounder when estimating a propensity score
  weighted treatment effect in the presence of selection bias. 
}
\author{
\packageAuthor{SBdecomp}

Maintainer: \packageMaintainer{SBdecomp}
}
\references{
Parast, L and Griffin, BA (2020). "Quantifying the Bias due to Observed Individual Confounders in Causal Treatment Effect Estimates." Statistics in Medicine, In press. 
}

\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}
\examples{
data(petsdata)
\donttest{sbdecomp(outcome = petsdata$genhealth, treatment = petsdata$gotdog, confounders = 
as.data.frame(petsdata[,c(2:13)]), type = "inclusion", estimation = "parametric")}

\donttest{sbdecomp(outcome = petsdata$genhealth, treatment = petsdata$gotdog, confounders = 
as.data.frame(petsdata[,c(2:13)]), type = "inclusion", estimation = "parametric", 
Bonly =FALSE, balance = TRUE) }

\donttest{sbdecomp(outcome = "genhealth", treatment = "gotdog", confounders = c("age",
"ismale", "race_coll","hhsize","ownhome", "married",  "ontanf", "hhincome",
"fulltime","spouse_fulltime" ,"liveinhouse", "ruralurban"), data = petsdata, 
type = "inclusion", estimation = "parametric", Bonly =FALSE, balance = TRUE)}

}
