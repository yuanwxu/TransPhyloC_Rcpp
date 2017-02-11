/* Rewrite probTTree with allowTransPostSamp = TRUE
   Calculates the log-probability of a transmission tree.
   @param ttree Transmission tree
   @param off.r First parameter of the negative binomial distribution for offspring number
   @param off.p Second parameter of the negative binomial distribution for offspring number
   @param pi probability of sampling an infected individual
   @param w.shape Shape parameter of the Gamma probability density function representing the generation time
   @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
   @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
   @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
   @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
   @param allowTransPostSamp Whether or not to allow transmission after sampling of a host
   @return Probability of the transmission tree */


// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <limits>
using namespace Rcpp;

struct wstar_functor
{ // Functor also returning 1st derivative.
  wstar_functor(double const& pi, double const& p, double const& r) : pi(pi), p(p), r(r){};
  
  std::pair<double, double> operator()(double const& x)
  {
    // Return both f(x) and f'(x).
    double temp = pow((1-p)/(1-p*x),r);
    double fx = x - (1-pi)*temp; 
    double dx = 1 - (1-pi)*p*r/(1-p*x)*temp;
    return std::make_pair(fx, dx);   
  }
private:
  double pi;
  double p;
  double r;
};


double wstar_rootFinder(double pi, double p, double r)
{
  using namespace boost::math::tools;
  
  const int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type T.
  int get_digits = static_cast<int>(digits * 0.6);    // Accuracy doubles with each step, so stop when we have
                                                      // just over half the digits correct.
  const boost::uintmax_t maxit = 20;
  boost::uintmax_t it = maxit;
  double result = newton_raphson_iterate(wstar_functor(pi, p, r), 0.5, 0.0, 1.0, get_digits, it);
  return result;
}


double alphastar(int d, double p, double r, double wstar)
{
  if(abs(r-1.0)<1e-6) // Exact solution available
    return (1-p)/(1-p*wstar)*pow(p/(1-p*wstar), d);

  int k = d;
  std::vector<double> toSum;

  /* For scalar calculation, compute pmf of NegBinom 
     using c++ Boost "binomial_coefficient".  
     Rcpp sugar version of dnbinom is vectorised, 
     so not convinient to work with here. */
  
  // Function dnbinom("dnbinom");
  
  while(true){
    // double term = as<double>(dnbinom(k, r, p))*pow(wstar,k);
    double dnb = boost::math::binomial_coefficient<double>(k+r-1,k) * pow(p,k) * pow(1-p, r);
    double term = dnb * pow(wstar,k);

    toSum.push_back(term);
    if(term < 1e-8) break; // Achieved desired accuracy
    
    k++;
  }

  NumericVector toSumR = wrap(toSum); // Convert toSum to Rcpp NumericVector
  NumericVector v(k-d+1);
  for(int i=0; i<v.size(); i++) v[i] = i+d;

  return sum(choose(v, d)*toSumR)/pow(wstar,d);
}
    
       
 
// [[Rcpp::export]]
double probTTreeC(NumericMatrix ttree, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double dateT){

  int numCases = ttree.nrow();
  boost::math::gamma_distribution<double> genGamma(shGen, scGen);
  
  if(dateT == INFINITY){ // finished outbreak
    double wstar = wstar_rootFinder(pi, pOff, rOff);

    NumericVector sstatus = ifelse(is_na(ttree(_,1)), 1-pi, pi*dgamma(ttree(_,1)-ttree(_,0),shSam,scSam));
    NumericVector lsstatus = log(sstatus);

    std::map<int, std::vector<int> > infMap; // Map from infector to infected
    std::vector<std::vector<int> > progeny(numCases);
    for(int i=0; i<numCases; ++i){
      if(ttree(i,2) == 0) continue; // Found root node i 

      progeny[ttree(i,2)-1].push_back(i); // C++ index starts from 0
      infMap[ttree(i,2)-1] = progeny[ttree(i,2)-1]; 
    }

    double accum = 0.0;
    for(int i=0; i<numCases; ++i){
      accum += log(alphastar(progeny[i].size(), pOff, rOff, wstar));
     
      for(int j=0; j<progeny[i].size(); ++j)
	accum += log(pdf(genGamma, ttree(progeny[i][j],0) - ttree(i,0)));
    }
      
    return sum(lsstatus) + accum;
  }

  // Ongoing case not implemented.
  return 0;
}
