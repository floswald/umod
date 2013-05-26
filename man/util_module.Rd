\name{util_module}
\alias{util_module}
\title{C++ module for computation of V = U + b*EV}
\usage{
  util_module(cashR, saveR, EVR, hsizeR, laborR, par, b)
}
\arguments{
  \item{cashR}{numeric matrix \code{(n,m)} of cash holdings
  conditional on labor supply (that's why \code{m}
  columns)}

  \item{saveR}{numeric matrix \code{(n,k)} of savings
  options. k < n.}

  \item{EVR}{numeric matrix \code{(n,k)} representing
  expected future value at each state,choice combination}

  \item{hsizeR}{vector \code{(n,1)}}

  \item{laborR}{vector \code{(m,1)}, basically
  \code{seq(from=0,to=1,length=m)}}

  \item{b}{boolean indicator equal TRUE if you have NAs in
  savings matrix imposing some borrowing constraint. if
  b==FALSE, no NA checking is done and results will be
  incorrect.}

  \item{theta}{elasticity of substitution between c and h}

  \item{phival}{value of relative utility difference flat
  vs house}

  \item{mu}{weight on additive utility premium}

  \item{gamma}{coefficient of relative risk aversion}

  \item{cutoff}{minimum level of consumption. below cutoff,
  u(c) is quadratically approximated.}

  \item{alpha}{coefficient on labor}
}
\value{
  list with elements \item{values}{\code{(n,m)} matrix of
  conditional value functions. column i is V_i.}
  \item{saving}{\code{(n,m)} matrix of conditional savings
  functions. column i is save_i.} \item{cons}{\code{(n,m)}
  matrix of conditional value functions. column i is
  cons_i.} \item{dchoiceL}{\code{(n,1)} vector indexing
  discrete choice at each state.} \item{maxL}{\code{(n,1)}
  vector of maximal value. maxL = max_d V_d.}
}
\description{
  computes value functions for \code{m} discrete labor
  supply choices. no time separability assumed, i.e. labor
  supply is not implied by current resources. computes V on
  \code{n} states.
}
\examples{
n = 5    # number of states
k = 5    # number of savings choices by state
m = 3    # number of discrete labor choices by state
cash   <- matrix(1:n,n,m)
cash   <- cash + matrix(0:2,n,m,byrow=TRUE)
labo   <- seq(from=0,to=1,length=m)
saving <- matrix(seq(from=0,to=8,length=k),n,k,byrow=TRUE)
EV     <- log(outer(1:n,1:n))
hsize  <- sample(0:2,size=n,replace=TRUE)
pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6)
res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars, b=FALSE)
##
## work with borrowing constraint in savings matrix: all saving less than x inadmissible
##
saving[1,1:3] <-NA	# all savings less than element c(1,4) are illegal in row 1
res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars, b=FALSE) # WRONG
res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars, b=TRUE) # RIGHT
}
\author{
  Florian Oswald <florian.oswald@gmail.com>
}

