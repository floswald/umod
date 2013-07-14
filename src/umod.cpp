

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

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
//' @param theta elasticity of substitution between c and h
//' @param phival value of relative utility difference flat vs house
//' @param mu weight on additive utility premium
//' @param gamma coefficient of relative risk aversion
//' @param cutoff minimum level of consumption. below cutoff, u(c) is quadratically approximated.
//' @param alpha coefficient on labor
//' @param quad boolean of whether neg cons quadratically approximated or not
//' @param borrconst boolean of whether there are borrowing constraints built into the savings matrix
//' @param myNA numerical value to be assigned to values with negative consumption (if quad == FALSE)
//' @param tau numerical value of proportional consumption scaling. a value in [0,infty). used in welfare experiments
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
//' pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6,quad=FALSE,borrconst=FALSE,myNA=-1e9,tau=1)
//' res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars)
//' ##
//' ## work with borrowing constraint in savings matrix: all saving less than x inadmissible
//' ##
//' saving[1,1:3] <-NA	# all savings less than element c(1,4) are illegal in row 1
//' pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6,quad=TRUE,borrconst=TRUE,myNA=-1e9,tau=1)
//' res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars) 
// [[Rcpp::export]]
List util_module(NumericMatrix cashR, NumericMatrix saveR, NumericMatrix EVR, NumericVector hsizeR, NumericVector laborR, List par){

    signal(SIGSEGV, handler);   // install our handler
   // start timer
   Timer timer; 

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


  timer.step("initial checks");

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

	// get consumption tax
	double tau = Rcpp::as<double>(par["tau"]);
	// get borrowing constraint switch
	bool b   = Rcpp::as<bool>(par["borrconst"]);


	timer.step("allocate memory");
	
	// loop over discrete choices
	for (int i=0;i<m;i++){

   // in-loop timer
		y.zeros();

		// fill tmp with copies of cash
		tmpcash = repmat(cash.col(i),1,k);
		// get consumption at each savings choice
		tmpcash = tmpcash - save;

		// apply the consumption tax
		tmpcash = tmpcash * tau;
		if (i==0) timer.step("allocate inloop-memory");
		// get utility
		util    = ufun_discreteL( tmpcash, par, hsize, labor(i) );
		// get lifetime utility
		W       = util + EV;
		if (i==0) timer.step("compute utility inloop");


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
		if (i==0) timer.step("maximize by row");
		}
		// put into return object at discrete choice column i
		tmpy.col(i)  = y;
	}

	timer.step("loop dchoice and maximize");

	// find optimal discrete choice
	for (int irow=0; irow<n; irow++){
		tmpvec3     = tmpy.row(irow);
		retV(irow)  = tmpvec3.max( j );
    //retC(irow)  = tmpc(irow,j);
    //retS(irow)  = tmps(irow,j);
		retiD(irow) = j + 1; // record discrete choice as 1-based index.
	}

	timer.step("find optimal discrete choice");

  // return vectors at optimal discrete choice: don't return conditional policy functions
  
	Rcpp::List list = Rcpp::List::create( _["values"] = tmpy, _["saving"] = tmps, _["cons"] = tmpc, _["dchoiceL"] = retiD, _["maxL"] = retV, _["timer"]=timer);
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
//' @param quad boolean of whether neg cons quadratically approximated or not
//' @param borrconst boolean of whether there are borrowing constraints built into the savings matrix
//' @param myNA numerical value to be assigned to values with negative consumption (if quad == FALSE)
//' @param tau numerical value of proportional consumption scaling. a value in [0,infty). used in welfare experiments
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
//' pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6,quad=FALSE,borrconst=FALSE,myNA=-1e9,tau=1)
//' res <- util_module_file(cashR=cash, EVR=EV, hsizeR=hsize, laborR=labo, par=pars)
// [[Rcpp::export]]
List util_module_file(NumericMatrix cashR, NumericMatrix EVR, NumericVector hsizeR, NumericVector laborR, List par){

    signal(SIGSEGV, handler);   // install our handler
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
	
	// get consumption tax
	double tau = Rcpp::as<double>(par["tau"]);
	cash = cash * tau;

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

		// get consumption at each current labor supply choice and apply cons tax
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



//' Utility Function from Attanasio et al
//'
//' computes utility over consumption and housing
//' @param Res resources aka consumption
//' @param s vector of house sizes
//' @param par list of parameters
//' @examples
//' n = 5    # number of states
//' m = 7    # number of savings choices
//' cash   <- matrix(1:(n*m),n,m)
//' hsize  <- sample(0:2,size=n,replace=TRUE)
//' pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6)
//' res <- ufun_Attanasio(ResR=cash, sR=hsize, par=pars)
// [[Rcpp::export]]
NumericMatrix ufun_Attanasio( NumericMatrix ResR, NumericVector sR, List par){

    signal(SIGSEGV, handler);   // install our handler
	const int n = ResR.nrow();
	const int m = ResR.ncol();

	if (sR.size() != n){
		throw Rcpp::exception( "umod::ufun_Attanasio: hsize and Res are not equal rows!" );
		return R_NilValue;
	}

	// map R to arma
	mat Res(ResR.begin(),n,m,false);
	vec hsize(sR.begin(),n,false);


	//par_ = list(gamma,theta,phival,mu,cutoff)

	// now get parameters out of lists par 
	double gamma = Rcpp::as<double>(par["gamma"]);
	double theta = Rcpp::as<double>(par["theta"]);
	double cutoff = Rcpp::as<double>(par["cutoff"]);
	double phival = Rcpp::as<double>(par["phival"]);
	double mu = Rcpp::as<double>(par["mu"]);
	double mgamma = 1-gamma;
	double imgamma = 1/mgamma;

	vec phivals;
	phivals << 0 << phival << 1 << endr;
	vec phivec(hsize.size());
	for (int i=0; i<hsize.size(); i++) {
		phivec(i) = phivals( hsize( i ) );
	}

	//setup output matrix
	mat ret(Res);
	ret.zeros();
	double diff, g, u, grad, hess;
	rowvec tmpvec(m);

	// compute utility from non-durable consumption
	for (int i=0; i<n; i++){
		if( Res(i,m-1) < cutoff ){	// check if the last entry in row i is below cutoff. if last is not, none is. Res = c - savings, savings increase to the right.
			for (int j=0; j<m; j++){
				if( Res(i,j) < cutoff ){ // compute approximated utility at all entries of row i below cutoff
					g        = pow(cutoff,mgamma);	 
					grad     = g / cutoff;   
					hess     = -gamma * grad / cutoff;	
					diff     = Res(i,j) - cutoff;
					ret(i,j) = imgamma*g + (grad * diff) + 0.5 * hess * pow(diff,2);
				} else {
					g = pow(Res(i,j),mgamma);	  
					ret(i,j) = imgamma * g;	    
				}
			}
		} else {
			tmpvec = pow(Res.row(i),mgamma);   // equation (6)
			ret.row(i) = imgamma * tmpvec;    // equation (11)
		}
	}
	// add additive premium if houseing is of right size
	mat phimat = repmat(phivec,1,m);
	mat hfac = exp( theta * phimat);
	ret = ret % hfac;
	ret = ret + mu * phimat;
	return wrap(ret);
}


//' Segfault Test Function boom
//'
//' test function produces a C++ segfault. Calls function baz which
//' allocates a wrong pointer. If you compiled this code with 
//' \code{CXXFLAGS=-g3 -rdynamic} the installed function \code{handler}
//' will print a traceback of the stack that contains the name of the offending function.
//' Without this compiler flag, you miss the function name.
//' you should place a call to \code{signal(SIGSEGV, handler);} at the beginning of each
//' function you want to check for segfaults.
//' @return R will crash with a segfault but you will see a traceback. ONLY run in console.
// [[Rcpp::export]]
void boom(){
    signal(SIGSEGV, handler);   // install our handler
	foo();
}
