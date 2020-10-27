\name{petsdata}
\alias{petsdata}
\docType{data}
\title{
Dog ownership dataset
}
\description{
This dataset is a random subset of the publicly available 2003 California Health Interview Survey data; the data consists responses from 2,102 adults. For the purpose of illustrating the functions in this package, the goal is to investigate the effect of dog ownership on general health.  Dog ownership was assessed with the question ``Do you have any dogs that you allow inside your home?"; 29.0\% of respondents owned a dog. General health status of the individual was measured as the self-reported response to the question ``Would you say that in general your health is excellent, very good, good, fair or poor?" Responses were coded from 1 through 5 with 5 indicating ``Excellent." Available individual characteristics i.e., confounders, in this dataset include age, gender, race/ethnicity, household size, marriage status, whether the individual received TANF (Temporary Assistance for Needy Families), household annual income, whether the individual worked full time, whether the individual had a spouse that worked full time, whether the individual lived in a house, and a rural/urban measure (1= urban; 2= 2nd city; 3 = suburban; 4 = town and rural) for the individual's address. 
}
\usage{data("petsdata")}
\format{
  A data frame with 2102 observations on the following 14 variables.
  \describe{
    \item{\code{gotdog}}{indicates whether individual reported owning a dog where 1 is yes and 0 is no (binary)}
    \item{\code{age}}{age of the individual (numeric)}
    \item{\code{ismale}}{indicates the sex of the individual where 1 is yes and 0 is no (binary)}
    \item{\code{race_coll}}{ace/ethnicity of the individual, collapsed into categories where 1 is Latino, 4 is Asian, 6 is White, and 7 is a combined category of Pacific Islander/Native American/African American/other race/multiple race (factor/categorical)}
\item{\code{hhsize}}{household size (numeric count)}
\item{\code{ownhome}}{indicates whether the individual owns the home they live in where 1 is yes and 0 is no (binary)}
\item{\code{married}}{indicates whether the individual is married where 1 is yes and 0 is no (binary)}
\item{\code{ontanf}}{indicates whether the individual received TANF where 1 is yes and 0 is no (binary)}    
\item{\code{hhincome}}{household income, log-transformed (numeric)}
\item{\code{fulltime}}{indicates whether the individual works full time where 1 is yes and 0 is no (binary)}
\item{\code{spouse_fulltime}}{indicates whether the individual has a spouse that works full time where 1 is yes and 0 is no (binary)}
\item{\code{liveinhouse}}{indicates whether the individual lives in a house (as opposed to an apartment, duplex, or mobile home) where 1 is yes and 0 is no (binary)}
\item{\code{ruralurban}}{the rural/urban measure of the home address of the individual where 1= urban, 2= 2nd city, 3 = suburban, 4 = town and rural (numeric)} 
\item{\code{genhealth}}{self-reported general health of the individual, numeric 1-5 scale where 5 is excellent (numeric)} 
  }
}
\examples{
data(petsdata)
dim(petsdata)
names(petsdata)
}
\references{
	California Health Interview Survey (2005). CHIS 2003 Methodology Series. Los Angeles: UCLA Center for Health Policy Research.
}
\keyword{datasets}
