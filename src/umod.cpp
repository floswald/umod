

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "ufuns.h"

using namespace Rcpp;

//' C++ module for computation of V = max_s U + b*EV(s)
//'
//' computes value functions for \code{m} discrete labor supply choices. no time separability assumed, i.e. labor supply is not implied by current resources. computes V on \code{n} states. 
//' The current implementation computes the following functional form for the utility function:
//' \deqn{\begin{array}{ll}
//' u(c,l_j,h) &= \frac{\left(c \exp( \alpha l_j ) \right)^{1-\gamma}}{1-\gamma} \exp( \theta \phi(h) ) + \mu \phi(h) \\
//' \phi(h) &= \left\{
//' 	\begin{array}{ll}
//' 	0                       &   h=0 \\
//' 	\mbox{phival} \in (0,1) &   h=1 \\
//' 	1                       &   h=2 
//' 	\end{array}
//' 	\right. \\
//' \mbox{params} &= \left( \gamma, \alpha, \theta, \mbox{phival}, \mu \right)\\
//'  & \hspace{1em}j = 1,2,\dots,m \\
//'  & \hspace{1em}l_1 = 0, l_m=1, l_j < l_{j+1} 
//' \end{array}}{ u(c,l,h) = (c * exp(alpha * l) )^1-gamma / (1-gamma) * exp( theta * phi(h) ) + mu * phi(h), phi(h) = 0 if h=0, phi(h) = phival if h=1, phi(h) = 1 if h=2 }
//'
//' @param cashR numeric matrix \code{(n,m)} of cash holdings conditional on labor supply (that's why \code{m} columns)
//' @param saveR numeric matrix \code{(n,k)} of savings options. k < n.
//' @param EVR numeric matrix \code{(n,k)} representing expected future value at each state,choice combination
//' @param hsizeR vector \code{(n,1)}
//' @param laborR vector \code{(m,1)}, basically \code{seq(from=0,to=1,length=m)}
//' @param b boolean indicator equal TRUE if you have NAs in savings matrix imposing some borrowing constraint. if b==FALSE, no NA checking is done and results will be incorrect.
//' @param theta elasticity of substitution between c and h
//' @param phival value of relative utility difference flat vs house
//' @param mu weight on additive utility premium
//' @param gamma coefficient of relative risk aversion
//' @param cutoff minimum level of consumption. below cutoff, u(c) is quadratically approximated.
//' @param alpha coefficient on labor
//' @return list with elements
//' \item{values}{\code{(n,m)} matrix of conditional value functions. column i is V_i.}
//' \item{saving}{\code{(n,m)} matrix of conditional savings functions. column i is save_i.}
//' \item{cons}{\code{(n,m)} matrix of conditional value functions. column i is cons_i.}
//' \item{dchoiceL}{\code{(n,1)} vector indexing discrete choice at each state.}
//' \item{maxL}{\code{(n,1)} vector of maximal value. maxL = max_d V_d.}
//' @author Florian Oswald <florian.oswald@@gmail.com>
//' @examples
//' n = 5    # number of states
//' k = 5    # number of savings choices by state
//' m = 3    # number of discrete labor choices by state
//' cash   <- matrix(1:n,n,m)
//' cash   <- cash + matrix(0:2,n,m,byrow=TRUE)
//' labo   <- seq(from=0,to=1,length=m)
//' saving <- matrix(seq(from=0,to=8,length=k),n,k,byrow=TRUE)
//' EV     <- log(outer(1:n,1:n))
//' hsize  <- sample(0:2,size=n,replace=TRUE)
//' pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6)
//' res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars, b=FALSE)
//' ##
//' ## work with borrowing constraint in savings matrix: all saving less than x inadmissible
//' ##
//' saving[1,1:3] <-NA	# all savings less than element c(1,4) are illegal in row 1
//' res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars, b=FALSE) # WRONG
//' res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars, b=TRUE) # RIGHT
// [[Rcpp::export]]
List util_module(NumericMatrix cashR, NumericMatrix saveR, NumericMatrix EVR, NumericVector hsizeR, NumericVector laborR, List par, bool b){

    const int n  = cashR.nrow();   // number of states
	const int m  = cashR.ncol();	// number of discrete choices
	const int k  = saveR.ncol();   // number of savings choices

	// some input checks
	if (n != saveR.nrow()){
		throw Rcpp::exception( "util_module: saveR not equal rows as cashR");
		return R_NilValue;
	} 
	if (n != EVR.nrow()){
		throw Rcpp::exception( "util_module: EVR not equal rows as cashR");
		return R_NilValue;
	} 
	if (k != EVR.ncol()){
		throw Rcpp::exception( "util_module: EVR not equal cols as saveR");
		return R_NilValue;
	} 
  if (n != hsizeR.size()){
		throw Rcpp::exception( "util_module: hsizeR not same number of rows as cashR");
		return R_NilValue;
	} 
  // remark: current structure forces us to have n = k. if that changes in the future, this module can handle it.

	// map to R to arma
	arma::mat cash(cashR.begin(), n, m, false);	// advanced constructor. no copying.
	arma::mat EV(EVR.begin(), n, k, false);	
	arma::mat save(saveR.begin(), n, k, false);	
	arma::vec labor(laborR.begin(),m,false);
	arma::vec hsize(hsizeR.begin(),n,false);
	// Rcpp::List par( par_ ) ;

	// allocate tmpcash, utility, and W matrices
	arma::mat tmpcash(n,k);
	arma::mat util(n,k);
	arma::mat W(n,k);

	// allocate indicator objects if need to subset
	arma::uvec u(k);
	arma::rowvec tmpvec2;
	arma::rowvec tmpvec3(m);	// for discrete choice

	// allocate out objects 
	arma::mat tmpy(n,m);	// maximal value for each discrete choice (value)
	arma::mat tmps(n,m);	// optimizing savings choice (value)
	arma::mat tmpc(n,m);	// implied optimal consumption (value)
	arma::uvec retiD(n);		// index of optimal discrete choice (index)
    arma::vec retV(n);		// value of optimal discrete choice 
  //arma::vec retC(n);		// cons at discrete choice 
	//arma::vec retS(n);		// save at discrete choice
  

	// allocate for maximization loop
	arma::uword j;
//	arma::uvec::fixed<n> iy;
	arma::vec y(n);
	arma::rowvec tmpvec(k);

	
	// loop over discrete choices
	for (int i=0;i<m;i++){

		y.zeros();

		// fill tmp with copies of cash
		tmpcash = repmat(cash.col(i),1,k);
		// get consumption at each savings choice
		tmpcash = tmpcash - save;
		// get utility
		util    = ufun_discreteL( tmpcash, par, hsize, labor(i) );
		// get lifetime utility
		W       = util + EV;

		// maximize
		// loop over rows of W and find maximal value and it's index for each row.
		//
		// insert switch about restricted savings space here
		// need to construct an index matrix FALSE where is not finite
		if (b) {
			// there are NAs in save matrix - need to account for that
			for (int irow=0;irow<n;irow++){
			// get current row of savings matrix to see which choices are illegal
			   tmpvec = save.row(irow);
			   // for this row, find all infinite indexes
			   for (int r=0;r<k;r++){
			       u(r) = arma::is_finite(tmpvec(r));
			   }
			   // overwrite with values of W
			   tmpvec = W.row(irow);
				// max over row
			   // while is.na(u[j]), increase j to get to first non-na u[j]
			   int r = 0;
			   while ( !u(r) ){
			   	r += 1;
			   }
			   tmpvec2 = tmpvec.subvec(r,k-1);
			   y(irow) = tmpvec2.max(j);	// j is the uword recording the position of the maximal element
			   tmpc(irow,i) = tmpcash(irow,j + r);
			   tmps(irow,i) = save(irow,j + r);	// savings matrix is the same for all discrete labor choices
			}
		} else {
			// there are no NAs in savings matrix, unrestricted choice
			for (int irow=0; irow<n; irow++) {
				tmpvec       = W.row(irow);
				y(irow)      = tmpvec.max( j );
				tmpc(irow,i) = tmpcash(irow,j);
				tmps(irow,i) = save(irow,j);	// savings matrix is the same for all discrete labor choices
		//		iy(irow)    = j + 1; // go to 1 based indices we don't even need those indices anymore.
			}
		}
		// put into return object at discrete choice column i
		tmpy.col(i)  = y;
	}

	// find optimal discrete choice
	for (int irow=0; irow<n; irow++){
		tmpvec3     = tmpy.row(irow);
		retV(irow)  = tmpvec3.max( j );
    //retC(irow)  = tmpc(irow,j);
    //retS(irow)  = tmps(irow,j);
		retiD(irow) = j + 1; // record discrete choice as 1-based index.
	}

  // return vectors at optimal discrete choice: don't return conditional policy functions
  
	Rcpp::List list = Rcpp::List::create( _["values"] = tmpy, _["saving"] = tmps, _["cons"] = tmpc, _["dchoiceL"] = retiD, _["maxL"] = retV);
	return list;
}



//' C++ module for computation of V = U + b*EV(0)
//'
//' computes value functions for \code{m} discrete labor supply choices when there is no savings choice.
//' The current implementation computes the following functional form for the utility function:
//' \deqn{\begin{array}{ll}
//' u(c,l_j,h) &= \frac{\left(c \exp( \alpha l_j ) \right)^{1-\gamma}}{1-\gamma} \exp( \theta \phi(h) ) + \mu \phi(h) \\
//' \phi(h) &= \left\{
//' 	\begin{array}{ll}
//' 	0                       &   h=0 \\
//' 	\mbox{phival} \in (0,1) &   h=1 \\
//' 	1                       &   h=2 
//' 	\end{array}
//' 	\right. \\
//' \mbox{params} &= \left( \gamma, \alpha, \theta, \mbox{phival}, \mu \right)\\
//'  & \hspace{1em}j = 1,2,\dots,m \\
//'  & \hspace{1em}l_1 = 0, l_m=1, l_j < l_{j+1} 
//' \end{array}}{ u(c,l,h) = (c * exp(alpha * l) )^1-gamma / (1-gamma) * exp( theta * phi(h) ) + mu * phi(h), phi(h) = 0 if h=0, phi(h) = phival if h=1, phi(h) = 1 if h=2 }
//'
//' @param cashR numeric matrix \code{(n,m)} of cash holdings conditional on labor supply (that's why \code{m} columns)
//' @param EVR numeric matrix \code{(n,1)} representing expected future value at tomorrow's assets = 0
//' @param hsizeR vector \code{(n,1)}
//' @param laborR vector \code{(m,1)}, basically \code{seq(from=0,to=1,length=m)}
//' @param theta elasticity of substitution between c and h
//' @param phival value of relative utility difference flat vs house
//' @param mu weight on additive utility premium
//' @param gamma coefficient of relative risk aversion
//' @param cutoff minimum level of consumption. below cutoff, u(c) is quadratically approximated.
//' @param alpha coefficient on labor
//' @return list with elements
//' \item{values}{\code{(n,m)} matrix of conditional value functions. column i is V_i.}
//' \item{dchoiceL}{\code{(n,1)} vector indexing discrete choice at each state.}
//' \item{maxL}{\code{(n,1)} vector of maximal value. maxL = max_d V_d.}
//' @author Florian Oswald <florian.oswald@@gmail.com>
//' @examples
//' n = 5    # number of states
//' m = 3    # number of discrete labor choices by state
//' cash   <- matrix(1:n,n,m)
//' cash   <- cash + matrix(0:2,n,m,byrow=TRUE)
//' labo   <- seq(from=0,to=1,length=m)
//' EV     <- log(1:n)
//' hsize  <- sample(0:2,size=n,replace=TRUE)
//' pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6)
//' res <- util_module_file(cashR=cash, EVR=EV, hsizeR=hsize, laborR=labo, par=pars)
// [[Rcpp::export]]
List util_module_file(NumericMatrix cashR, NumericMatrix EVR, NumericVector hsizeR, NumericVector laborR, List par){

    const int n  = cashR.nrow();   // number of states
	const int m  = cashR.ncol();	// number of discrete choices

	// some input checks
	if (n != EVR.nrow()){
		throw Rcpp::exception( "util_module: EVR not equal rows as cashR");
		return R_NilValue;
	} 
  if (n != hsizeR.size()){
		throw Rcpp::exception( "util_module: hsizeR not same number of rows as cashR");
		return R_NilValue;
	} 
  // remark: current structure forces us to have n = k. if that changes in the future, this module can handle it.

	// map to R to arma
	arma::mat cash(cashR.begin(), n, m, false);	// advanced constructor. no copying.
	arma::mat EV(EVR.begin(), n, 1, false);	
	arma::vec labor(laborR.begin(),m,false);
	arma::vec hsize(hsizeR.begin(),n,false);
	// Rcpp::List par( par_ ) ;

	// allocate tmpcash, utility, and W matrices
	arma::mat tmpcash(n,1);
	arma::mat util(n,1);
	arma::mat W(n,1);

	// allocate indicator objects if need to subset
	arma::rowvec tmpvec3(m);	// for discrete choice
	arma::uword j;

	// allocate out objects 
	arma::mat tmpy(n,m);	// maximal value for each discrete choice (value)
	arma::uvec retiD(n);		// index of optimal discrete choice (index)
    arma::vec retV(n);		// value of optimal discrete choice 
  

	// loop over discrete choices
	for (int i=0;i<m;i++){

		// get consumption at each current labor supply choice
		tmpcash = cash.col(i);
		// get utility
		util    = ufun_discreteL( tmpcash, par, hsize, labor(i) );
		// get lifetime utility
		W       = util + EV;

		// put into return object at discrete choice column i
		tmpy.col(i)  = W;
	}

	// find optimal discrete choice
	for (int irow=0; irow<n; irow++){
		tmpvec3     = tmpy.row(irow);
		retV(irow)  = tmpvec3.max( j );
		retiD(irow) = j + 1; // record discrete choice as 1-based index.
	}

  // return vectors at optimal discrete choice: don't return conditional policy functions
  
	Rcpp::List list = Rcpp::List::create( _["values"] = tmpy, _["dchoiceL"] = retiD, _["maxL"] = retV);
	return list;
}

