/*

   Rewrite inferTTree using Rcpp.

   - Omitting startCTree and allowTransPostSamp

   - For finished outbreak, must specify dateT = Inf since not all c++
   default arguments can be parsed into R equivalents. It is
   repositioned in the secound argument to inferTTreeC.

   - Must specify shSam and scSam, shape and scale parameter of
   sampling time density since C++ language rules on default
   arguments.

   - Ignoring "source" component of MCMC output, since need
     understanding of original R code

#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param startCTree Optional combined tree to start from
#' @param updateTTree Whether or not to update the transmission tree
#' @param optiStart Whether or not to optimise the MCMC start point
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @param allowTransPostSamp Whether or not to allow transmission after sampling of a host
#' @return posterior sample set of transmission trees
#' @export */

#include "probTTree.h"
#include <Rcpp.h>
using namespace Rcpp;

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
