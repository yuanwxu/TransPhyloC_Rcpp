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
#include <boost/math/distributions/negative_binomial.hpp>
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
  if(std::abs(r-1.0)<1e-6) // Exact solution available
    return (1-p)/(1-p*wstar)*pow(p/(1-p*wstar), d);

  int k = d;
  std::vector<double> toSum;
  
  // Function dnbinom("dnbinom");
  boost::math::negative_binomial_distribution<double> nbinom(r,p);
  while(true){
    // double term = as<double>(dnbinom(k, r, p))*pow(wstar,k);
    
    double dnb = pdf(nbinom,k);
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

// [[Rcpp::export]]
List inferTTreeC(List ptree, double dateT, double shGen, double scGen, double shSam, double scSam, int mcmcIterations=1000, int thinning=1, double startNeg=100.0/365.0, double start_rOff=1.0, double start_pOff=0.5, double startPi=0.5, bool updateNeg=true, bool update_rOff=true, bool update_pOff=false, bool updatePi=true, bool updateTTree=true, bool optiStart=true){

  /* To adopt the interface used in Xavier's inferTTree, 
     need to extract from the list "ptree" its matrix 
     component, do some operation, and then reconstruct
     back to a list. */
  
  NumericMatrix ptreeM = as<NumericMatrix>(ptree["ptree"]);
  ptreeM(_,0) = ptreeM(_,0) + runif(ptreeM.nrow())*1e-10; // Ensure that all leaves have unique times
  List ptreeNew = List::create(Named("ptree")=ptreeM, Named("nam")=ptree["nam"]);
  
  double neg = startNeg, rOff = start_rOff, pOff = start_pOff, pi = startPi;
  
  Environment transPhylo("package:TransPhylo"); // Obtain R environment for TransPhylo
  Function makeCtreeFromPTree = transPhylo["makeCtreeFromPTree"]; // Use existing functions in XD's TransPhylo package
  Function extractTTree = transPhylo["extractTTree"];
  Function probPTreeGivenTTree = transPhylo["probPTreeGivenTTree"];

  /* If we want to use the ".proposal" function defined in 
     proposal.R in TransPhylo, we need a trick. We have to
     first source("proposal.R") into R global environment, and 
     then make it available to call here. This is because 
     .proposal() is not availabe as package function. 
     Since .proposal() in turn depends on .computeHost(),
     need also source("computeHost.R").
     
     XD, consider if ok to make proposal function available as package function?
  */

  Environment env = Environment::global_env();
  Function proposal = env[".proposal"]; 
  
  List ctree;
  if(optiStart)
    ctree = makeCtreeFromPTree(ptreeNew, rOff, pOff, neg, pi, shGen, scGen, shSam, scSam, dateT, _["allowTransPostSamp"]=true);
  else
    ctree = makeCtreeFromPTree(ptreeNew, NA_REAL, pOff, neg, pi, shGen, scGen, shSam, scSam, dateT, _["allowTransPostSamp"]=true);

  List ttree = extractTTree(ctree);
  NumericMatrix ttreeM = ttree["ttree"];
  double pTTree = probTTreeC(ttreeM, rOff, pOff, pi, shGen, scGen, shSam, scSam, dateT); // YX's C++ implementation of probTTree
  double pPTree = as<double>(probPTreeGivenTTree(ctree, neg));
  
  List out(mcmcIterations/thinning);
  for(int i=1; i<=mcmcIterations; ++i){ // Now the MCMC loop
    if(i % thinning == 0)
      out[i/thinning-1] = List::create(Named("ctree")=ctree, Named("pTTree")=pTTree, Named("pPTree")=pPTree, Named("neg")=neg, Named("off.r")=rOff, Named("off.p")=pOff, Named("pi")=pi, Named("w.shape")=shGen, Named("w.scale")=scGen, Named("ws.shape")=shSam, Named("ws.scale")=scSam);
    
    if(updateTTree){
      List prop = proposal(ctree["ctree"]);
      List ctree2 = List::create(Named("ctree")=prop["tree"], Named("nam")=ctree["nam"]);
      List ttree2 = extractTTree(ctree2);
      NumericMatrix ttree2M = ttree2["ttree"];
      double pTTree2 = probTTreeC(ttree2M, rOff, pOff, pi, shGen, scGen, shSam, scSam, dateT); // YX's C++ implementation of probTTree
      double pPTree2 = as<double>(probPTreeGivenTTree(ctree2, neg));
      
      if(log(R::runif(0,1)) < log(as<double>(prop["qr"]))+pTTree2+pPTree2-pTTree-pPTree){ // Use scalar mode runif
	ctree = ctree2;
	ttree = ttree2;
	pTTree = pTTree2;
	pPTree = pPTree2;
      }
    }
    if(updateNeg){
      double neg2 = std::abs(neg + (R::runif(0,1)-0.5)*0.5);
      double pPTree2 = as<double>(probPTreeGivenTTree(ctree,neg2));
      if(log(R::runif(0,1)) < pPTree2-pPTree-neg2+neg){
	neg = neg2;
	pPTree = pPTree2;
      }
    }
    if(update_rOff){
      double r2Off = std::abs(rOff + (R::runif(0,1)-0.5)*0.5);
      double pTTree2 = probTTreeC(ttree["ttree"], r2Off, pOff, pi, shGen, scGen, shSam, scSam, dateT); // YX's C++ implementation of probTTree
      if(log(R::runif(0,1)) < pTTree2-pTTree-r2Off+rOff){
	rOff = r2Off;
	pTTree = pTTree2;
      }
    }
    if(update_pOff){
      double p2Off = std::abs(pOff + (R::runif(0,1)-0.5)*0.1);
      if(p2Off > 1) p2Off = 2 - p2Off;
      double pTTree2 = probTTreeC(ttree["ttree"], rOff, p2Off, pi, shGen, scGen, shSam, scSam, dateT); // YX's C++ implementation of probTTree
      if(log(R::runif(0,1)) < pTTree2-pTTree){
	pOff = p2Off;
	pTTree = pTTree2;
      }
    }
    if(updatePi){
      double pi2 = pi + (R::runif(0,1)-0.5)*0.1;
      if(pi2 < 0.01) pi2 = 0.02-pi2;
      if(pi2 > 1) pi2 = 2-pi2;
      double pTTree2 = probTTreeC(ttree["ttree"], rOff, pOff, pi2, shGen, scGen, shSam, scSam, dateT); // YX's C++ implementation of probTTree
      if(log(R::runif(0,1)) < pTTree2-pTTree){
	pi = pi2;
	pTTree = pTTree2;
      }
    }
  } // End MCMC loop
    
  return out;
}

