\name{reverseMAthread}
\alias{reverseMAthread}
\title{reverseMAthread}
\usage{
reverseMAthread(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0.1,
  varM = 1, beta0 = 0, betaX = 1, betaM = c(0, 0.1, 0.2),
  varY = 1, nSim = 100, nSimImai = 1000, SEED = 1, plot.pdf = T,
  plot.name = "reverseMAsim.pdf", alpha_level = 0.05,
  use_multi_processing = F, use_cpp = F, num_jobs = 1)
}
\arguments{
\item{n}{is the sample size.}

\item{pX}{is the minor allele frequency}

\item{gamma0}{is the intercept for M}

\item{gammaX}{is the association of X with M}

\item{varM}{is the variance of M}

\item{beta0}{is the intercept for Y}

\item{betaX}{is the direct effect of X on Y}

\item{betaM}{is a vector of different associations of M with Y}

\item{varY}{is the variance of Y}

\item{nSim}{is the number of simulations to run}

\item{nSimImai}{is the number of simulations to run in mediate from the mediation package}

\item{SEED}{is the seed}

\item{plot.pdf}{is T to output a plot, is F to not output a plot}

\item{plot.name}{is the name of the plot}

\item{alpha_level}{is the significance level}

\item{use_multi_processing}{use multi_processing instead of threading for a speed boost}

\item{use_cpp}{use Rcpp with Eigen}

\item{num_jobs}{the number of cores to use, i.e. the number of parallel procesess or threads to spawn}
}
\value{
a matrix of the power of the mediate method from the mediation package to detect an effect of the mediator M on the outcome Y when M and Y are correctly specified and also when they are incorrectly specified (the true mediator is Y and the true outcome is M)
}
\description{
A function to simulate the performance of the mediate function from the mediation package in scenarios of reverse causality.
}
\author{
Annie Thwing, Sharon Lutz
}
