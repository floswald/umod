#ifndef _RcppUtils_UFUNS_H
#define _RcppUtils_UFUNS_H

// [[Rcpp::depends(Rcpp,RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

using namespace arma;
using namespace Rcpp;


//////////////////////////////////////
// start C++ error handling
// tracing out a segmentation fault 
// to stdout

// Define Error Handler
inline void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  //backtrace_symbols(array, size);
  exit(1);
}

// Define Segfault Test function boom
void baz() {
 int *foo = (int*)-1; // make a bad pointer
  printf("%d\n", *foo);       // causes segfault
}

// Define foo and bar
void bar() { baz(); }
void foo() { bar(); }


// end C++ error handling
/////////////////////////
 
  
  
  

// utility with positive and negative consumption
vec u_pos(vec c,double lab,vec h,double alph,double gamm){
    signal(SIGSEGV, handler);   // install our handler
  double mgamm = 1-gamm;
  double imgamm = 1/mgamm;
	double g,z,y;
	vec ret(c);
	for (int i=0;i<c.n_elem;i++){
		y      = pow( exp( alph * lab ), mgamm);
		z      = y * h[i];
		g      = pow( c[i], mgamm);
		ret[i] = imgamm * g * z;
	}
	return ret;
}

// matrix utility
mat u_mat(mat c,double lab,vec h,double alph,double gamm){
    signal(SIGSEGV, handler);   // install our handler
    double mgamm = 1-gamm;
    double imgamm = 1/mgamm;
	double g,z,y;
	y      = pow( exp( alph * lab ), mgamm);
	mat ret(c.n_rows,c.n_cols);
	for (int i=0;i<c.n_rows;i++){

		z      = y * h[i];

		for (int j=0; j<c.n_cols;j++){

			g        = pow( c(i,j) , mgamm);
			ret(i,j) = imgamm * g * z;

		}
	}
	return ret;
}

// quadratic approximation function for case of no work for vector
vec u_neg(vec c,double cuto,double lab,vec h,double alph,double gamm){
    signal(SIGSEGV, handler);   // install our handler
	vec diff;
	vec ret(c);
	double grad,hess,y,z;
	double mgamm = 1-gamm;
	double imgamm = 1/mgamm;
	diff = c - cuto;
	for (int i=0;i<c.n_elem;i++){
		y      = pow( exp( alph * lab ), mgamm);
		z      = y * h[i];
		grad   = pow( cuto, -gamm) * z;
		hess   = -(gamm * grad) / cuto;
		ret[i] = imgamm * pow( cuto, 1-gamm) * z  + (grad * diff[i]) + 0.5 * hess * pow(diff[i],2);
	}
	return ret;
}


//' CRRA utility augmented with discrete labor supply and utility from housing
//'
//' function not exposed to R. used to compute CRRA utility at discrete labor supply choices.
//' @param Resources numeric matrix of consumption levels
//' @param hsize numeric vector of house sizes
//' @param labor numeric vector of discrete labor supply choices
//' @param params list of parameters 
//' \itemize{
//' \item{theta}{elasticity of substitution between c and h}
//' \item{phival}{value of relative utility difference flat vs house}
//' \item{mu}{weight on additive utility premium}
//' \item{gamma}{coefficient of relative risk aversion}
//' \item{cutoff}{minimum level of consumption. below cutoff, u(c) is quadratically approximated.}
//' \item{alpha}{coefficient on labor}
//' \item{quad}{boolean TRUE if quadratic approx for neg cons required.}
//' \item{myNA}{numerical value for infeasible choice}
//' }
//' @return numeric matrix of utility values 
mat ufun_discreteL(mat cons, Rcpp::List par, vec hsize, double labor){

    signal(SIGSEGV, handler);   // install our handler
	uword n = cons.n_rows;
	uword m = cons.n_cols;

	if ( hsize.n_elem != cons.n_rows){
		throw std::runtime_error("ufun_discreteL:::error. hsize and cons are not conformable");
	}

	// extract elements from list par
	double theta     = Rcpp::as<double>(par["theta"]);
	double phival    = Rcpp::as<double>(par["phival"]);
	double mu        = Rcpp::as<double>(par["mu"]);
	double gamma     = Rcpp::as<double>(par["gamma"]);
	double cutoff    = Rcpp::as<double>(par["cutoff"]);
	double alpha     = Rcpp::as<double>(par["alpha"]);
	double myNA      = Rcpp::as<double>(par["myNA"]);

	// construct the multiplicative housing factor from house size
	vec phivals;
	phivals << 0 << phival << 1 << endr;
	vec phivec(hsize.size());
	for (int i=0; i<hsize.n_elem; i++) {
		phivec(i) = phivals( hsize( i ) );
	}
	vec hfac  = exp( theta * phivec );

	// initiate return objects as zero
	mat util(n,m);
	util.zeros();

	// compute multiplicative part of utility at each consumption level
	util      = u_mat( cons, labor, hfac, alpha, gamma );

	// find indices of negative consumptions
	uvec ineg = find( cons < cutoff );
	
	// add additive premium 
	util = util + mu * repmat(phivec,1,m);

	// if there is negative consumption 
	if ( !(ineg.is_empty()) ){
		vec x(ineg.n_elem);
		x.fill( myNA );
 		// set utility to myNA value
		util.elem(ineg) = x;
	}
	
	return util;
}


#endif
