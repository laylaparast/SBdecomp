\name{bar.sbdecomp}
\alias{bar.sbdecomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a Bar Plot
}
\description{
This function creates a bar plot to illustrate the estimated proportion of the observed selection bias explained by each observed confounder used in the propensity score model. 
}
\usage{
bar.sbdecomp(output.list, main)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.list}{
List outputed from main sbdecomp function
}
 \item{main}{
Title for plot, optional
}
}
\value{
A bar plot showing the estimated proportions of the estimated selection bias explained by the observed confounders used in the propensity score model is created.
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
output = sbdecomp(outcome = petsdata$genhealth, treatment = petsdata$gotdog, confounders = 
as.data.frame(petsdata[,c(2:13)]), type = "inclusion", estimation = "parametric",
balance = FALSE) 

bar.sbdecomp(output)

}


