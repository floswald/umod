\name{ufun_Attanasio}
\alias{ufun_Attanasio}
\title{Utility Function from Attanasio et al}
\usage{
  ufun_Attanasio(ResR, sR, par)
}
\arguments{
  \item{Res}{resources aka consumption}

  \item{s}{vector of house sizes}

  \item{par}{list of parameters}
}
\description{
  computes utility over consumption and housing
}
\examples{
n = 5    # number of states
m = 7    # number of savings choices
cash   <- matrix(1:(n*m),n,m)
hsize  <- sample(0:2,size=n,replace=TRUE)
pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6)
res <- ufun_Attanasio(ResR=cash, sR=hsize, par=pars)
}

