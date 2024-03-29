% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grouped_perm_glmm.R
\name{grouped_perm_glmm}
\alias{grouped_perm_glmm}
\title{Permutation test for linear mixed-effects models}
\usage{
grouped_perm_glmm(tbl, formla, var_to_perm, permNum = 1000, seed = 42)
}
\arguments{
\item{tbl}{A data frame or a data frame extension (e.g. a tibble) containing
the variables in the model.}

\item{formla}{An object of class "formula" containing a symbolic description
of the linear mixed-effects model to be fitted. The object include both the
fixed-effects and random-effects part of the model.See
\code{\link[lme4]{lmer}} for more details.}

\item{var_to_perm}{The columns of variables to permute.}

\item{permNum}{The number of permutations to generate.}

\item{seed}{A single value interpreted as an integer for initializing a
random permutation generator. See \code{\link[base]{set.seed}} for more
details.}
}
\value{
It returns a tibble containing the following components:
\tabular{ll}{\code{term}\tab  The name of the regression term.\cr
\tab \cr
\code{effect}\tab  The type of the tested effect of the
regression term.\cr
\tab \cr
\code{estimate} \tab The estimated value of the regression term.\cr
\tab \cr
\code{statistic} \tab The value of a T-statistic to use in a hypothesis
that the regression term is non-zero.\cr
\tab \cr
\code{p.perm} \tab The likelihood of observing the test statistic of the
original data among that of permutations.}
}
\description{
\code{grouped_perm_glmm} is a nonparametric test of linear
mixed-effects models. It assesses the significance of coefficients by
permutation tests, which calculate the distribution of the test statistic by
randomly rearranging the observed data.
}
\examples{
library(sdamr)
data("anchoring")
grouped_perm_glmm(anchoring,everest_feet ~ anchor+sex + (1|referrer),  
"everest_feet")
grouped_perm_glmm(anchoring,everest_feet ~ anchor+sex + (1|referrer),  
"everest_feet", permNum = 500, seed = 2)

}
