// Spatiotemporal Bayesian ETAS
// Author: Max Schneider, based on Gordon Ross's bayesianETAS code

#include <Rcpp.h>
using namespace Rcpp;
// use 'g++ -std=c++11 -o cxx12_random cxx12_random.cpp'

#include <R.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <list>
#include <algorithm> 
#include <string.h> 

#include <stdio.h>
#include <stdlib.h>

using namespace std;



extern "C" {
  double distanceCalculate(double x1, double y1, double x2, double y2) {
    double x = x1 - x2; //calculating number to square in next step
    double y = y1 - y2;
    double dist, dist_sq;
    
    dist_sq = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
    dist = sqrt(dist_sq);                  
    
    return dist;
  }
}

extern "C" {
  double sqDistanceCalculate(double x1, double y1, double x2, double y2) {
    double x = x1 - x2; //calculating number to square in next step
    double y = y1 - y2;
    double dist_sq;
    
    dist_sq = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
    // dist = sqrt(dist_sq);                  
    
    return dist_sq;
  }
}


// function found online at acmath
// TODO: track down URL
inline double 
  BilinearInterpolation(double q11, double q12, double q21, double q22, double x1, double x2, double y1, double y2, double x, double y) 
  {
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    // std::cout << "Inputs " << x1 << " "  << x2 << " " << y1 << " " << y2 << "\n";
    // std::cout << "Components " << x2x1 << " "  << y2y1 << " " << x2x << " " << y2y << " "  << yy1 << " " << xx1 << " " << "\n";
    return 1.0 / (x2x1 * y2y1) * (
        q11 * x2x * y2y +
          q21 * xx1 * y2y +
          q12 * x2x * yy1 +
          q22 * xx1 * yy1
    );
  }

// [[Rcpp::export]]
// with or without & in argument?
double interpolateGevalsRCpp(double this_d, double this_q, NumericMatrix Gevals_grid, int index){
  // Perform 2D linear interpolation between nearest d/q elements to approximate the Geval for this_d/this_q
  // sizeof array / sizeof array[0] is the best way to get nrow of an array
  // https://stackoverflow.com/questions/10274162/how-to-find-2d-array-size-in-c/24984668
  // int num_rows = sizeof(Gevals_grid)/sizeof(Gevals_grid[0]);
  // int num_cols = sizeof(Gevals_grid[0]) / sizeof(int); 
  // int num_rows = Gevals_grid.size();
  // int num_cols = Gevals_grid[0].size();
  int num_rows = Gevals_grid.nrow();
  int col_eqk = index;
  // std::cout << "Number rows/cols: " << num_rows << " / " << col_eqk << "\n";
  
  // Pull out the d and q values in this grid 
  std::vector<double> grid_ds (num_rows);
  std::vector<double> grid_qs (num_rows);
  for (int i = 0; i < num_rows; i++) {
    grid_ds[i] = Gevals_grid(i, 0);
    grid_qs[i] = Gevals_grid(i, 1);
  }
  
  // find the unique d and q values
  // vector<int> v{1,2,3,1,2,3,3,4,5,4,5,6,7};
  sort(grid_ds.begin(), grid_ds.end());
  auto lastd = unique(grid_ds.begin(), grid_ds.end());
  grid_ds.erase(lastd, grid_ds.end());
  sort(grid_qs.begin(), grid_qs.end());
  auto lastq = unique(grid_qs.begin(), grid_qs.end());
  grid_qs.erase(lastq, grid_qs.end());
  
  // for (int i = 0; i < 4; i++) {
  //   std::cout << grid_ds[i] << " " << grid_qs[i] << "\n";
  // }
  
  // Find the bounding d and q values from the grid to this_d/this_q
  // rbegin/rend will go in backwards order
  // datatype auto will let the compiler decide what the datatype should be. avoids errors of the type
  // conversion from ‘X’ to non-scalar type ‘Y’ requested (with X and Y related to iterators)
  auto d_below = std::find_if(grid_ds.rbegin(), grid_ds.rend(), [&](double val){ return val < this_d;});
  auto d_above = std::find_if(grid_ds.begin(), grid_ds.end(), [&](double val){ return val > this_d;});
  auto q_below = std::find_if(grid_qs.rbegin(), grid_qs.rend(), [&](double val){ return val < this_q;});
  auto q_above = std::find_if(grid_qs.begin(), grid_qs.end(),  [&](double val){ return val > this_q;});
  
  // std::cout << "d_below is " << *d_below << "\n";
  // std::cout << "d_above is " << *d_above << "\n";
  // std::cout << "q_below is " << *q_below << "\n";
  // std::cout << "q_above is " << *q_above << "\n";
  
  
  // // Interpolate Geval for this_d
  // frac_d = this_d-d_below/(d_above-d_below);
  // // Interpolate Geval for this_d
  // frac_q = this_q-q_below/(q_above-q_below);
  //
  // Find the rows corresponding to this_d/this_q
  int r_d_below_q_below, r_d_above_q_below,  r_d_below_q_above, r_d_above_q_above;
  for (int r = 0; r < num_rows; r++) {
    if ((Gevals_grid(r, 0) == *d_below) & (Gevals_grid(r, 1) == *q_below)) {
      r_d_below_q_below = r;
    } else if ((Gevals_grid(r, 0) == *d_above) & (Gevals_grid(r, 1) == *q_below)) {
      r_d_above_q_below = r;
    } else if ((Gevals_grid(r, 0) == *d_below) & (Gevals_grid(r, 1) == *q_above)) {
      r_d_below_q_above = r;
    } else if ((Gevals_grid(r, 0) == *d_above) & (Gevals_grid(r, 1) == *q_above)) {
      r_d_above_q_above = r;
    }
  }
  
  // Setup the BilinearInterpolation function
  //
  // std::cout << "Q11 is " << Gevals_grid(r_d_below_q_below,col_eqk) << "\n";
  // std::cout << "Q21 is " << Gevals_grid(r_d_above_q_below,col_eqk) << "\n";
  // std::cout << "Q12 is " << Gevals_grid(r_d_below_q_above,col_eqk) << "\n";
  // std::cout << "Q22 is " << Gevals_grid(r_d_above_q_above,col_eqk) << "\n";
  // 
  // std::cout << "x1 is " << Gevals_grid(r_d_below_q_below,0) << "\n";
  // std::cout << "x2 is " << Gevals_grid(r_d_above_q_above,0) << "\n";
  // std::cout << "y1 is " << Gevals_grid(r_d_below_q_below,1) << "\n";
  // std::cout << "y2 is " << Gevals_grid(r_d_above_q_above,1) << "\n";
  
  double q11, q12, q21, q22, x1, x2, y1, y2, x, y;
  // x variable for d, y variable for q
  // double this_interp_Geval =
  //   BilinearInterpolation(q11=Gevals_grid(r_d_below_q_below,last_col), q12=Gevals_grid(r_d_above_q_below,last_col),
  //                         q21=Gevals_grid(r_d_below_q_above,last_col), q22=Gevals_grid(r_d_above_q_above,last_col),
  //                         x1=Gevals_grid(r_d_below_q_below,0), x2=Gevals_grid(r_d_above_q_above,0),
  //                         y1=Gevals_grid(r_d_below_q_below,1), y2=Gevals_grid(r_d_above_q_above,1), x=this_d, y=this_q);
  // double this_interp_Geval = 5;
  double this_interp_Geval = BilinearInterpolation(q11=Gevals_grid(r_d_below_q_below,col_eqk), q12=Gevals_grid(r_d_above_q_below,col_eqk),
                                                   q21=Gevals_grid(r_d_below_q_above,col_eqk), q22=Gevals_grid(r_d_above_q_above,col_eqk),
                                                   x1=*d_below, x2=*d_above,
                                                   y1=*q_below, y2=*q_above, x=this_d, y=this_q);
  // std::cout << "interpolated Gval " << this_interp_Geval << "\n";
  
  return(this_interp_Geval);
}

// This function just returns the branching structure only, for testing in R
// [[Rcpp::export]]
std::vector<int> STsampleBranching2(std::vector<double> &ts, std::vector<double> &lons, std::vector<double> &lats, 
                                    std::vector<double> &marks, double M0, double mu, double logK, double alpha, 
                                    double c, double p, double d, double q,  double maxT, 
                                    NumericMatrix &cat_Gevals_grid,
                                    std::vector<int> &branching) {
  // Initialize some objects.
  int n = ts.size();
  int parent;
  double temp;
  double K = exp(logK);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  // Initialize branching vector. Clear what was in there before, make it size n and add an element "0"
  branching.clear();
  branching.reserve(n);
  branching.push_back(0);
  
  // Create probs and kappaevals vectors of size n
  std::vector<double> probs;
  probs.reserve(n);
  
  std::vector<double> kappaevals;
  kappaevals.reserve(n);
  
  // Rcout << "sampleBranching2" << endl;
  // Rcout << B_prior_days << endl;
  
  
  //(GR: precompute)
  // For each eqk...
  for (int i = 0; i < n; i++) {
    // Compute the kappa fn for each eqk
    kappaevals.push_back(K*exp(alpha*(marks[i]-M0)));
  }
  
  // For each eqk...
  for (int i = 1; i < n; i++) {
    // Clear the probs vector and insert in the mu
    probs.clear();
    probs.push_back(mu);
    // Rcout <<  i << endl;
    
    // For each previous eqk before this one...
    for (int j = 0; j < i; j++) {//(GR: check this iteratures through...)
      // Calculate the aftershock triggering of that prev eqk (j) from the current eqk (i)
      // This is kappa * the temporal triggering fn. Add this to the probs vector
      // Rcout << "dist: " << i << " and " << j << ": "<< sqDistanceCalculate(lons[i], lats[i], lons[j], lats[j]) << endl;
      temp = kappaevals[j]*pow(ts[i]-ts[j]+c,-p)*pow(sqDistanceCalculate(lons[i], lats[i], lons[j], lats[j])+d,-q); // Unnorm Omori
      // temp = kappaevals[j]*(p-1)*pow(c, p-1)*pow(ts[i]-ts[j]+c,-p)*pow( sqDistanceCalculate(lons[i], lats[i], lons[j], lats[j])  +d , -q); // ITO, Ross (21)
      probs.push_back(temp);
      // Rcout <<  temp << endl;
    }
    // double Gintegral = interpolateGevalsRCpp(d, q, cat_Gevals_grid, i);
    // Rcout <<  Gintegral << endl;
    
    // Create a discrete distn, with probability weights given by probs
    // (aftershock triggering values for each preceding eqk, as well as background value)
    std::discrete_distribution<> d(probs.begin(), probs.end()); //correct?
    
    // Draw a single integer from this discrete distn. This is (this iteration's) branching
    // value - the index for the event (stochastically) most likely to have triggered it. 
    // Add this to the branching vec
    parent = d(gen);
    branching.push_back(parent); 
  }
  return(branching);
}

// void sampleBranching(std::vector<double> &ts, std::vector<double> &marks, double M0,   double mu, double logK, double alpha, double c, double p, std::vector<int> &branching);

//GR: maybe change this to avoid using push_back/clear. Instead, use [] to put elements into vector. See http://stackoverflow.com/questions/20168051/why-push-back-is-slower-than-operator-for-an-previous-allocated-vector
// Inputs: eqk catalog,
//    initial values for mu, logK, alpha, c and p
// Outputs: randomly generated branching structure for the catalog, given the current
// param values
// [[Rcpp::export]]
std::vector<int> STsampleBranchingMS(std::vector<double> &ts, std::vector<double> &lons, std::vector<double> &lats, 
                                     std::vector<double> &marks, double M0, 
                                     double mu, double K, double alpha, double c, double p, double d, double q) {
  // Initialize some objects.
  int n = ts.size();
  int parent;
  double temp;
  // double K = exp(logK);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  // Initialize branching vector. Clear what was in there before, make it size n and add an element "0"
  // (the first eqk in a catalog must be bkgd)
  std::vector<int> branching;
  branching.reserve(n);
  branching.push_back(0);
  
  // Create probs and kappaevals vectors of size n
  std::vector<double> probs;
  // probs.reserve(n);
  
  std::vector<double> kappaevals;
  kappaevals.reserve(n);
  
  //(GR: precompute)
  // For each eqk...
  for (int i = 0; i < n; i++) {
    // Compute the kappa fn for each eqk
    kappaevals.push_back(K*exp(alpha*(marks[i]-M0)));
  }
  
  // Rcout <<  "mu: " << mu << endl;
  double pi = 2 * acos(0.0); 
  double utsucoef = (q-1) / (pi*pow(d, 1-q));
  // Rcout << "utsucoef " <<  utsucoef  << endl;
  
  // For each eqk...
  for (int i = 1; i < n; i++) {
    // Clear the probs vector and insert in the mu
    probs.clear();
    probs.push_back(mu);
    // Rcout <<  i << endl;
    
    // For each previous eqk before this one...
    for (int j = 0; j < i; j++) { //(GR: check this iteratures through...)
      // Calculate the aftershock triggering of that prev eqk (j) from the current eqk (i)
      // This is kappa * the temporal triggering fn. Add this to the probs vector
      // temp = kappaevals[j]*(p-1)*pow(c,p-1)*pow(ts[i]-ts[j]+c,-p);
      // temp = kappaevals[j]*pow(ts[i]-ts[j]+c,-p);
      // Unnormed Omori using squared distance (correct)
      // Rcout <<  j <<  endl;
      temp = kappaevals[j]*pow(ts[i]-ts[j]+c,-p)*utsucoef*pow( sqDistanceCalculate(lons[i], lats[i], lons[j], lats[j]) + d,-q); 
      // Rcout <<  "kappa*temp trig part: " << kappaevals[j]*pow(ts[i]-ts[j]+c,-p) << " , kappa*temp*spatial trig part: " << temp << " , bkgd part: " << mu << endl;
      
      // Unnormed Omori using distance (incorrect)
      // temp = kappaevals[j]*pow(ts[i]-ts[j]+c,-p)*pow( distanceCalculate(lons[i], lats[i], lons[j], lats[j]) + d,-q); 
      // temp = kappaevals[j]*(p-1)*pow(c, p-1)*pow(ts[i]-ts[j]+c,-p)*pow( sqDistanceCalculate(lons[i], lats[i], lons[j], lats[j])  +d , -q); // ITO, Ross (21)
      probs.push_back(temp);
    }
    // Create a discrete distn, with probability weights given by probs
    // (aftershock triggering values for each preceding eqk, as well as background value)
    std::discrete_distribution<> d(probs.begin(), probs.end()); //(GR: correct?)
    // 
    //     for (int j = 0; j < probs.size(); j++) {
    //      cout << probs[j] << " ";
    //     }
    //      cout << "\n";
    
    // Draw a single number from this discrete distn. This is (this iteration's) branching
    // value - the event (stochastically) most likely to have triggered it. Add this to the
    // branching vec
    parent = d(gen);
    // Rcout <<  "For eqk i: " << i << " probablistic parent is: " <<  parent << endl;
    
    branching.push_back(parent);
  }
  return(branching);
}

double SThBranchingPosteriorMSPrior(std::vector<double> &ts, std::vector<double> &marks, 
                                    std::vector<double> &z, double maxT, 
                                    std::vector<double> &kappaevals, std::vector<double> &Gevals, 
                                    double c, bool cpgamma, double chyper1, double chyper2, 
                                    double p, bool ppgamma, double phyper1, double phyper2) {
  
  // Uniform priors on c/p
  // if (c <= 0 || p <= 0) {
  //   return(-9999999);
  // }
  
  double hPrior;
  // Rcout <<  "This is cpgamma: " <<  cpgamma <<  "This is chypers: " <<  chyper1 <<  " " <<   chyper2 << " This is phypers: " << phyper1 <<   " " << phyper2 << endl;
  if(!cpgamma){
    // Rcout << "Not gamma priors for c/p! " << endl;
    if (c <= chyper1 || p <= phyper1 || c > chyper2 || p > phyper2) {
      // Rcout << "This is c: " <<  c << " This is p: " << p  << " and its too big" << endl;
      return(-9999999);
    }
    else{
      hPrior = 0;
    }
  } else {
    // Rcout << "Gamma priors for c/p! " << endl;
    // Rcout <<  "This is cpgamma: " <<  cpgamma <<  "This is chypers: " <<  chyper1 <<  " " <<   chyper2 << " This is phypers: " << phyper1 <<   " " << phyper2 << endl;
    hPrior = (chyper1 - 1)*log(c) - chyper2*c - log(phyper2 - phyper1);
    // Rcout <<  "This is cpgamma: " <<  cpgamma <<  "This is chypers: " <<  chyper1 <<  " " <<   chyper2 << " This is phypers: " << phyper1 <<   " " << phyper2 << endl;
    
  } 
  
  
  
  
  // if (c <= 0) {
  //   c = 0.001;
  // }
  // 
  // if (c > 8) {
  //   c = 8;
  // }
  // 
  // if (p <= 1) {
  //   p = 1.001;
  // }
  // 
  // if (p > 8) {
  //   p = 8;
  // }
  
  
  // Prior
  
  int n = ts.size();
  // Initialize values of interest, storage vec
  double hSumTerm, hIntTerm, hPost;
  std::vector<double> hSums;
  hSums.reserve(n);
  std::vector<double> hInts;
  hInts.reserve(n);
  // Rcout << "This is c: " <<  c << " This is p: " << p  << endl;
  
  
  // Rcout << "z1 " <<  z[0] << ", z2 " << z[1] << " z5 " <<  z[4] << " z10 " <<  z[10] << endl;
  
  for (int i = 0; i < z.size(); i++){
    hSums.push_back( -1 * p*log(z[i]+c) );
  }
  
  // Finish compute the sum term
  // For each triggered eqk...
  for (int i = 0; i < ts.size(); i++) {
    
    // Collect sum parts for each eqk
    // 6.11.21. Unnormalized O
    // if(z[i] > 0){
    //   // hSums.push_back( -1 * p*log(z[i]+c) );
    //   hSums[i] = -1 * p*log(z[i]+c);
    // } else{
    //   // hSums.push_back(0);
    //   hSums[i] = 0;
    // }
    
    
    // Collect integral parts for each eqk
    // 6.11.21. Unnormalized O
    double Hevali = (pow(maxT-ts[i] +c, 1-p) - pow(c, 1-p)) / (1-p); // H(T-t_i)
    hInts.push_back(-1 * kappaevals[i]*Gevals[i]*Hevali );
  }
  
  // Calculate the sum terms and integral terms (both are sums)
  hSumTerm = std::accumulate(hSums.begin(), hSums.end(), 0.0);
  hIntTerm = std::accumulate(hInts.begin(), hInts.end(), 0.0);
  // Calculate the posterior (sum of sum and integral terms)
  hPost = hPrior + hSumTerm + hIntTerm;
  // Rcout << "This is c: " <<  c << " This is p: " << p  << " with sum term: " << hSumTerm << " int term: " << hIntTerm << " and post: " << hPost << endl;
  
  // Return the log posterior after each sum has been calculated
  return(hPost);
}



// (GR: this shows up when we cant do the infinite time approximatoin...)
// Inputs: eqk catalog,
//    dist: distance between each triggered event and its triggering event
//    kappaevals: computed in STestimateETASBranchingInteractionMS fn
//    given values for d and q
// Outputs: log marginal posterior value for (d,q), based on the given values
// [[Rcpp::export]]
double STgBranchingPosteriorMSPrior(std::vector<double> &lons, std::vector<double> &lats, 
                                    std::vector<double> &marks, 
                                    std::vector<double> &dists, double maxT, 
                                    std::vector<double> &kappaevals, 
                                    std::vector<double> &Hevals, 
                                    std::vector<double> &Gevals,
                                    double d, bool dpgamma, double dhyper1, double dhyper2, 
                                    double q, bool qpgamma, double qhyper1, double qhyper2) {
  
  // if (d <= 0 || q <= 1) {
  //   return(-9999999);
  // }
  
  // Uniform priors on d/q - general
  double gPrior;
  // Rcout << "This is dhypers: " <<  dhyper1 <<  " " <<  dhyper2 << " This is qhypers: " << qhyper1 <<  " " <<  qhyper2 << endl;
  if(!dpgamma){
    if (d <= dhyper1 || q <= qhyper1 || d > dhyper2 || q > qhyper2) {
      // Rcout << "This is d: " <<  d << " This is q: " << q  << " and its too big" << endl;
      return(-9999999);
    }
    else{
      gPrior = 0;
    }
    // Gamma prior on d, uniform prior on q
  } else {
    gPrior = (dhyper1 - 1)*log(d) - dhyper2*d - log(qhyper2 - qhyper1);
  }
  
  
  // Initialize values of interest, storage vec
  // log post starts off at 0 = log(prior for uniform priors)
  // TODO: adjust this when changing prior
  double pi = 2 * acos(0.0); 
  // int n = dists.size();
  
  int n = lats.size();
  // Initialize values of interest, storage vec
  double gSumTerm, gIntTerm, gPost;
  std::vector<double> gSums;
  gSums.reserve(n);
  std::vector<double> gInts;
  gInts.reserve(n);
  
  for (int i = 0; i < dists.size(); i++){
    gSums.push_back( log(q-1) - log(pi) + ((q-1)*log(d)) - q*(log( dists[i] + d )) );
  }
  
  // Finish compute the sum term
  // For each triggered eqk...
  for (int i = 0; i < lats.size(); i++) {
    
    // // Collect sum parts for each eqk
    // // 6.11.21. Unnormalized O
    // if(dists[i] > 0){
    //   gSums.push_back( log(q-1) - log(pi) + ((q-1)*log(d)) - q*(log( dists[i] + d )) );
    // } else{
    //   gSums.push_back(0);
    // }
    
    // Collect integral parts for each eqk
    // 6.11.21. Unnormalized O
    gInts.push_back(-1 * kappaevals[i]*Gevals[i]*Hevals[i] );
  }
  
  // Calculate the sum terms and integral terms (both are sums)
  gSumTerm = std::accumulate(gSums.begin(), gSums.end(), 0.0);
  gIntTerm = std::accumulate(gInts.begin(), gInts.end(), 0.0);
  // Calculate the posterior (sum of sum and integral terms)
  gPost = gPrior + gSumTerm + gIntTerm;
  // Rcout << "This is d: " <<  d << " This is q: " << q  << " with sum term: " << gSumTerm << " int term: " << gIntTerm << " and post: " << gPost << endl;
  
  // Return the log posterior after each sum has been calculated
  return(gPost);
  
}

// Inputs: eqk catalog,
//    numtriggered: number of eqks triggered by each event in catalog (0 if its a bkgd event)
//    Hevals: computed in etasBranchingInteraction fn
//    Gevals: input by user
//    given values for logK and alpha, M0
// Outputs: log marginal posterior value for (K, alpha), based on the given values
// [[Rcpp::export]]
double STkappaBranchingPosteriorMSPrior(std::vector<double> &ts, std::vector<double> &marks, 
                                        std::vector<int> &numtriggered, double M0,
                                        std::vector<double> &Hevals, std::vector<double> &Gevals, 
                                        double K, bool Kpgamma, double Khyper1, double Khyper2, 
                                        double alpha, bool alphapgamma, double alphahyper1, double alphahyper2) {
  // // Get K off log scale and initialize values, vectors
  // double K = exp(logK);
  
  if (K <= 0 || alpha <= 0|| K > 10 || alpha > 10) {
    return(-9999999);
  }
  
  // Uniform priors on K/alpha - general
  double kPrior;
  // Rcout << "This is Khypers: " <<  Khyper1 << " " <<  Khyper2 << " This is alphahypers: " << alphahyper1 <<  " " << alphahyper2 << endl;
  if(!Kpgamma){
    // Rcout << "This is !Kpgamma: " << !Kpgamma << endl;
    if (K <= Khyper1 || alpha <= alphahyper1 || K > Khyper2 || alpha > alphahyper2) {
      // Rcout << "This is K: " <<  K << " This is alpha: " << alpha  << " and its too big" << endl;
      return(-9999999);
    }
    else{
      kPrior = 0;
    }
    // Gamma prior on K, uniform prior on alpha
  } else {
    // Rcout << "This is Kpgamma: " << Kpgamma << endl;
    kPrior = (Khyper1 - 1)*log(K) - Khyper2*K - log(alphahyper2 - alphahyper1);
  }
  
  int n = ts.size();
  // Initialize the vectors storing the sum and int components (over all eqks)
  std::vector<double> kSums;
  kSums.reserve(n);
  std::vector<double> kInts;
  kInts.reserve(n);
  // Initialize sum and int terms of the post, and post itself
  double kSumTerm, kIntTerm, kPost;
  
  // Finish compute the sum term
  // For each triggered eqk...
  for (int i = 0; i < ts.size(); i++) {
    double  kevali = K*exp(alpha*(marks[i]-M0));
    
    // Collect sum parts for each eqk
    kSums.push_back( numtriggered[i] * log(kevali) );
    
    // Collect integral parts for each eqk
    // 6.11.21. Unnormalized O
    kInts.push_back(-1 * kevali*Gevals[i]*Hevals[i] );
  }
  
  // Calculate the sum terms and integral terms (both are sums)
  kSumTerm = std::accumulate(kSums.begin(), kSums.end(), 0.0);
  kIntTerm = std::accumulate(kInts.begin(), kInts.end(), 0.0);
  // Calculate the posterior (sum of sum and integral terms)
  kPost = kPrior + kSumTerm + kIntTerm;
  // Rcout << "For K and alpha: " <<  K << "  " << alpha << " Sum term: " << kSumTerm  << " Int term: " << kIntTerm << " Post: " << kPost  << endl;
  
  // Return the log posterior after each sum has been calculated
  return(kPost);
}

// Inputs: eqk catalog,
//    numtriggered: number of eqks triggered by each event in catalog
//    Hevals: computed in etasBranchingInteraction fn
//    given values for logK and alpha
// Outputs: log marginal posterior value for (K, alpha), based on the given values
//                                           // TODO: NumericMatrix &cat_Gevals_grid should be an input argument
//[[Rcpp::export]]
List STestimateETASBranchingFreeBMSListPrior(std::vector<double> &ts, std::vector<double> &marks,
                                             std::vector<double> &lons, std::vector<double> &lats,
                                             std::vector<int> &branching, double maxT, double M0,
                                             int sims, int numMCMCSamples,
                                             double mu, double logK,
                                             double spatialArea,
                                             double alpha, double c, double p, double d, double q,
                                             std::vector<double> & initval,
                                             NumericMatrix &cat_Gevals_grid, bool alphaFixed,
                                             bool mupgamma, bool Kpgamma, bool alphapgamma,
                                             bool cpgamma, bool ppgamma, bool dpgamma, bool qpgamma,
                                             double muhyper1, double muhyper2,
                                             double Khyper1, double Khyper2,
                                             double alphahyper1, double alphahyper2,
                                             double chyper1, double chyper2,
                                             double phyper1, double phyper2,
                                             double dhyper1, double dhyper2,
                                             double qhyper1, double qhyper2) {
  
  //cout << "Interaction sampler\n";
  // Initialize objects to store calculated values
  int n = ts.size(), numbackground;
  // double currposterior,newposterior, newc,newp, newd, newq, newlogK,newalpha;
  double posteriorcpNew,posteriorcpCurr;
  double posteriordqNew,posteriordqCurr;
  double posteriorKalphaNew,posteriorKalphaCurr;
  double cNew, pNew, dNew, qNew, logKNew,KNew, alphaNew;
  
  
  // double mualpha = muHypers[0], mubeta = muHypers[1]; //prior params - change for mu_ST?
  // double mu_prior_a = 0.1, mu_prior_b = 0.1; //prior params - change for mu_ST?
  // double mualpha = 0.3, mubeta = 0.03; //prior params - much longer prior, max at 60*4
  // double mualpha = 2, mubeta = 3000000; //prior params - Ends at 1 shock/day*km^2 for a 500^2 area
  // double mualpha = 2, mubeta = 1000000; //prior params - Ends at 3 shock/day*km^2 for a 500^2 area
  // double mualpha = 1, mubeta = 700000; //prior params - Ends at 3.75 shock/day*km^2 for a 500^2 area (= corresponds to Ga(0.1, 0.1) for 1000^2, converted to 500^2 area )
  // double mualpha = 1.1, mubeta = 580000; //prior params - Ends at 4.5 shock/day*km^2 for a 500^2 area. Seemed to work best on 8/1
  // double mualpha = 0.5, mubeta = 35000; //This comes closest to what Ga(0.1, 0.1) spread over 500^2 ends at (60 shocks)
  // double mualpha = 1, mubeta = 500000; //prior params - Ends at 5 shock/day*km^2 for a 500^2 area
  // double mualpha = 0.8, mubeta = 200000; //prior params - changed for mu_ST
  // double mualpha = 0.8, mubeta = 1; //prior params - change for mu_ST? 7/31
  // double mualpha = 1, mubeta = 4; //prior params - change for mu_ST? 7/31
  // double mualpha = 0.01, mubeta = 1; //prior params - changed on 6.14 for mu_ST
  // double mualpha = 0.001, mubeta = 0.001; //prior params for mu_ST
  // double mualpha = 1, mubeta = 5000000; //prior params - Ends at 0.5 shock/day*km^2 for a 500^2 area
  // double mualpha = 1, mubeta = 600000000; //prior params - come out to posterior on the scale we need
  
  // Initialize kappaevals (= mag-based triggering) and Hevals (=)
  // and z (= time between triggered event and triggering event, for each triggered event)
  std::vector<double> kappaevals;
  kappaevals.reserve(n);
  std::vector<double> Hevals;
  Hevals.reserve(n);
  std::vector<double> z;
  z.reserve(n);
  // std::vector<double> sqdists;
  // sqdists.reserve(n);
  std::vector<double> dists;
  dists.reserve(n);
  std::vector<double> Gevals;
  Gevals.reserve(n);
  
  std::vector<int> nbkgds;
  nbkgds.reserve(sims);
  
  
  std::vector<double> mus;
  std::vector<double> Ks;
  std::vector<double> alphas;
  std::vector<double> cs;
  std::vector<double> ps;
  std::vector<double> ds;
  std::vector<double> qs;
  mus.reserve(sims);
  Ks.reserve(sims);
  alphas.reserve(sims);
  cs.reserve(sims);
  ps.reserve(sims);
  ds.reserve(sims);
  qs.reserve(sims);
  
  // Call random num generators we need
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> rnorm(0, 0.1);
  std::normal_distribution<> rnorm_mu(0, 0.05);
  std::normal_distribution<> rnorm_K(0, 0.05);
  // std::normal_distribution<> rnorm_c(0, 0.25);
  // std::normal_distribution<> rnorm_p(0, 0.3);
  // std::normal_distribution<> rnorm_d(0, 0.25); // previously same as c, p
  // std::normal_distribution<> rnorm_q(0, 0.3);
  // std::normal_distribution<> rnorm_K(0, 0.1);
  // std::normal_distribution<> rnorm_c(0, 0.1);
  // std::normal_distribution<> rnorm_p(0, 0.1);
  // std::normal_distribution<> rnorm_d(0, 0.1); // previously same as c, p
  // std::normal_distribution<> rnorm_q(0, 0.1);
  std::uniform_real_distribution<> runif(0, 1);
  // std::normal_distribution<> rnorm_K(0, 0.5);
  // std::normal_distribution<> rnorm_c(0, 0.5);
  // std::normal_distribution<> rnorm_p(0, 0.5);
  // std::normal_distribution<> rnorm_d(0, 0.5); // previously same as c, p
  // std::normal_distribution<> rnorm_q(0, 0.5);
  std::normal_distribution<> rnorm_c(0, 0.05);
  std::normal_distribution<> rnorm_p(0, 0.05);
  std::normal_distribution<> rnorm_d(0, 0.05);
  std::normal_distribution<> rnorm_q(0, 0.05);
  
  
  // Initialize the current parameters at their initial vals
  double muCurr = initval[0];
  double KCurr = initval[1];
  double alphaCurr = initval[2];
  double cCurr = initval[3];
  double pCurr = initval[4];
  double dCurr = initval[5];
  double qCurr = initval[6];
  // Rcout << "^^^^^^^ Starting mu: " <<  muCurr << ", mu scaled: " << muCurr * pow(1000,2) << ", K: " <<  KCurr << ", alpha: " <<  alphaCurr << ", c: " <<  cCurr << ", p: " <<  pCurr << ", d: " <<  dCurr << ", q: " <<  qCurr  << endl;
  
  // double K = exp(logK);
  // For each simulation iteration...
  for (int s=0; s<sims; s++) {
    // Rcout << "DRAWDRAWDRAWDRAWDRAWDRAWDRAWDRAWDRAWDRAWDRAWDRAW" << s << endl;
    Rcout << s << endl;
    // Rcout << "mu: " <<  mu << ", K: " <<  K << ", alpha: " <<  alpha << ", c: " <<  c << ", p: " <<  p << ", d: " <<  d << ", q: " <<  q  << endl;
    // Rcout << "mu: " <<  muCurr << ", mu scaled: " << muCurr * pow(1000,2) << ", K: " <<  KCurr << ", alpha: " <<  alphaCurr << ", c: " <<  cCurr << ", p: " <<  pCurr << ", d: " <<  dCurr << ", q: " <<  qCurr  << endl;
    // Rcout << "mu: " <<  muCurr << ", mu scaled: " << muCurr * pow(900, 2) << ", K: " <<  KCurr << ", alpha: " <<  alphaCurr << ", c: " <<  cCurr << ", p: " <<  pCurr << ", d: " <<  dCurr << ", q: " <<  qCurr  << endl;
    // Rcout << "mu: " <<  muCurr << ", mu scaled: " << muCurr * pow((Starget[2] - Starget[1]), 2) << ", K: " <<  KCurr << ", alpha: " <<  alphaCurr << ", c: " <<  cCurr << ", p: " <<  pCurr << ", d: " <<  dCurr << ", q: " <<  qCurr  << endl;
    Rcout << "mu: " <<  muCurr << ", mu scaled: " << muCurr * spatialArea << ", K: " <<  KCurr << ", alpha: " <<  alphaCurr << ", c: " <<  cCurr << ", p: " <<  pCurr << ", d: " <<  dCurr << ", q: " <<  qCurr  << endl;
    
    branching.clear();
    // Run the STsampleBranchingMS fn with the current param values
    branching = STsampleBranchingMS(ts, lons, lats, marks, M0, muCurr,KCurr, alphaCurr, cCurr,
                                    pCurr, dCurr, qCurr);
    
    // TODO: Can we remove this z.clear?
    z.clear();
    dists.clear();
    // Initialize vector with 0 to put the number of eqks triggered by every event
    // Thus, bkgd events will just have 0 here
    std::vector<int> numtriggered(n,0);
    // Initialize counter of the number of eqks not triggered by any eqk
    numbackground = 0;
    
    // For each eqk i...
    for (int i=0; i<n; i++) {
      
      // If B_i is not 0 (this eqk was triggered by some other eqk)
      if (branching[i] > 0) {
        // Raise the numtriggered for that triggering eqk's index by 1.
        //branching[i]-1 is the index of the triggering event since C++ starts counting at 0
        numtriggered[branching[i]-1]++;
        // Calculate time diff between the event and its trigger, add to z
        z.push_back(ts[i]-ts[branching[i]-1]);
        dists.push_back(sqDistanceCalculate(lons[i], lats[i],
                                            lons[branching[i]-1], lats[branching[i]-1]));
        // numtriggered[branching[i]]++;
        // z.push_back(ts[i]-ts[branching[i]]);
        // dists.push_back(distanceCalculate(lons[i], lats[i], lons[branching[i]], lats[branching[i]]));
        // Otherwise, its a bkgd event so raise that counter
      } else {
        numbackground++;
      }
    }
    Rcout <<  "background is " << numbackground << endl;
    nbkgds.push_back(numbackground);
    
    // Rcout << "c and p " <<  c << "  " << p << " chosen with log posterior " << currposterior << endl;
    // Rcout << "mupgamma is " <<  mupgamma << endl;
    
    // Generate a random mu, add to the vector of mus
    // auto lons_range = std::minmax_element (lons.begin(), lons.end());
    // auto lats_range = std::minmax_element (lats.begin(), lats.end());
    // double lons_dim = *lons_range.second - *lons_range.first;
    // double lats_dim = *lats_range.second - *lats_range.first;
    // std::cout << "range is  " << *lons_range.second - *lons_range.first << '\n';
    
    // std::gamma_distribution<> rgamma(mualpha+numbackground, 1/(mubeta+maxT*(500^2)));
    // mu = rgamma(gen);
    std::gamma_distribution<> rgamma(muhyper1+numbackground, 1/(muhyper2+maxT));
    // if(mupgamma){ // For now, mu always has a gamma prior
    //   std::gamma_distribution<> rgamma(muhyper1+numbackground, 1/(muhyper2+maxT));
    // } else{
    //   std::gamma_distribution<> rgamma(muhyper1+numbackground, 1/(muhyper2+maxT));
    // }
    // int saveit = 500^2;
    // Rcout <<  "500^2 is " << saveit << "pow(500,2) is " << pow(500,2) << " and what I'm calculating is" << (pow(1000,2)/pow(500,2) * pow(500,2)) << endl;
    
    // Rcout <<  "divisor is " << pow(1000,2)/pow(500,2) << endl;
    // mu = rgamma(gen) / (pow(1000,2)/pow(500,2) * pow(500,2)); // from before 7/31
    // mu = rgamma(gen)/pow(1000,2) * (pow(500,2)/pow(1000,2)); // 7/31
    double mutempCurr = rgamma(gen); // 8/3
    // muCurr = mutempCurr/pow(1000,2) * (pow((Starget[2] - Starget[1]), 2)/pow(1000,2));
    // muCurr = mutempCurr/pow(1000,2); // for simulations
    // muCurr = mutempCurr/730767.1; // for PNW
    // muCurr = mutempCurr/pow(900, 2); // for PNW
    // muCurr = mutempCurr/pow((Starget[2] - Starget[1]), 2); // for PNW
    muCurr = mutempCurr/spatialArea; // for PNW
    // mu = rgamma(gen);
    // mu = 0.1/(pow(1000,2) * 4); // fix mu to truth
    mus.push_back(muCurr);
    // Rcout <<  "mutemp is " << mutempCurr <<  "mu is " << muCurr << endl;
    
    Hevals.clear();
    Gevals.clear();
    kappaevals.clear();
    // K = exp(logK);
    // For each eqk...
    for (int i=0; i<n; i++) {
      // Calculate the last term of Gordon's eqn 4 (H)
      // Hevals.push_back(1 - (pow(c,p-1) / pow(maxT-ts[i]+c,p-1)));
      // Changed 5/25 to be for infinite-time normalized Omori
      // Hevals.push_back((pow(maxT-ts[i] +c, 1-p) - pow(c, 1-p)) / (-pow(c, 1-p)));
      // Changed 4/13 to be for finite-time normalized Omori
      // Hevals.push_back((pow(maxT-ts[i] +c, 1-p) - pow(c, 1-p)) / (pow(maxT+c, 1-p) - pow(c, 1-p)));
      // Changed 6.11.21 to be for unnormalized Omori
      // Hevals.push_back( (pow(maxT-ts[i] +c, 1-p) - pow(c, 1-p)) / (1-p) );
      
      // Changed 6.14.21 to be for Ross (21)'s H for his ITO
      // Hevals.push_back( 1 - ( pow(c,p-1) / pow(maxT-ts[i]+c,p-1) ) ); // H(T-t_i)  eqn 4 in Ross (21)  (typo??)
      // Hevals.push_back( 1 - (pow( c / maxT-ts[i]+c, p-1)) ); // H(T-t_i)  eqn 5 in Ross (21)  (typo??)
      
      // Holschneider's Omori
      // Hevals.push_back( (pow(1+ ((maxT-ts[i])/c), 1-p) -1) / (pow(1+ (maxT/c), 1-p) -1) );
      
      // double Gintegral = interpolateGevalsRCpp(d, q, cat_Gevals_grid, i);
      double Gintegral = interpolateGevalsRCpp(dCurr, qCurr, cat_Gevals_grid, i);
      // double Gintegral = 0.99;
      // double Gintegral = 1;
      // Rcout <<  "CHECKMEOUT CHECKMEOUT CHECKMEOUT CHECKMEOUT CHECKMEOUT Gintegral for eqk " << i << " is " << Gintegral << endl;
      // Rcout <<  "CHECKMEOUT CHECKMEOUT CHECKMEOUT CHECKMEOUT CHECKMEOUT kappaevals for eqk " << i << " is " << K*exp(alpha*(marks[i]-M0)) << endl;
      Gevals.push_back(Gintegral);
      
      // kappaevals.push_back(K*exp(alpha*(marks[i]-M0)));
      kappaevals.push_back(KCurr*exp(alphaCurr*(marks[i]-M0)));
      
    }
    // Rcout <<  "20th kappaevals is " << kappaevals[20] <<  " and Gevals is " << Gevals[20] << endl;
    
    
    // Rcout <<  "finished Gevalsfinished Gevalsfinished Gevalsfinished Gevalsfinished Gevalsfinished Gevalsfinished Gevalsfinished Gevals" << endl;
    
    
    // Get the posterior value for p(c, p) using the current c and p
    // currposterior = SThBranchingPosteriorMS(ts,marks,z,maxT,kappaevals, Gevals, c,p);
    posteriorcpCurr = SThBranchingPosteriorMSPrior(ts,marks,z,maxT,kappaevals, Gevals,
                                                   cCurr, cpgamma, chyper1, chyper2,
                                                   pCurr, ppgamma, phyper1, phyper2);
    // Rcout << "Starting c and p " <<  c << "  " << p  << endl;
    // For each MCMC sample...
    for (int i=0; i < numMCMCSamples; i++) {
      // Rcout << "MCMC sample: " << i << endl;
      // Generate proposal c and p value with the 1-d RW
      // newc = c + rnorm_c(gen);
      // newp = p + rnorm_p(gen);
      cNew = cCurr + rnorm_c(gen);
      pNew = pCurr + rnorm_p(gen);
      while(cNew <= 0){
        cNew = cCurr + rnorm_c(gen);
      }
      
      while(pNew <= 0){
        pNew = pCurr + rnorm_p(gen);
      }
      
      // if (newc <= 0) {
      //   newc = 0.001;
      // }
      //
      // if (newc > 8) {
      //   newc = 8;
      // }
      //
      // if (newp <= 1) {
      //   newp = 1.001;
      // }
      //
      // if (newp > 8) {
      //   newp = 8;
      // }
      //
      // if(newc < 0.01){
      //   newc = 0.01;
      //   Rcout <<  "hey! c small" << endl;
      // }
      // Generate a value from the posterior distn using the proposal c and p
      // newposterior = hBranchingPosteriorInteraction(ts,marks,z,maxT,kappaevals,newc,newp);
      // Rcout << "New c and p " <<  newc << "  " << newp  << endl;
      // newposterior = SThBranchingPosteriorMS(ts,marks,z, kappaevals, Hevals, Gevals, newc,newp);
      // newposterior = SThBranchingPosteriorMS(ts,marks,z,maxT,kappaevals, Gevals, newc,newp);
      if(cNew > 0){
        if(pNew > 0){
          posteriorcpNew = SThBranchingPosteriorMSPrior(ts, marks, z, maxT, kappaevals, Gevals,
                                                        cNew, cpgamma, chyper1, chyper2,
                                                        pNew, ppgamma, phyper1, phyper2);
        }
      }
      
      
      // Rcout << "Proposal posterior: " <<  posteriorcpNew <<  "Current posterior: " <<  posteriorcpCurr << ", comparison to old posterior: " << exp(posteriorcpNew-posteriorcpCurr) << endl;
      // newposterior = hBranchingPosteriorMH(ts,marks,z,maxT,kappaevals,newc,newp);
      // Metropolis-Hastings step.
      // Generate a random unif number and accept the proposal c and p and posterior if it
      // is less than the ratio of the new posterior and current posterior
      // ratio of logged posteriors = exponential of differences of logged posteriors
      // if (runif(gen) < exp(newposterior-currposterior)) {
      //   // Rcout << "Proposal posterior: " <<  newposterior << "for new c and p " <<  newc << "  " << newp << endl;
      //   c = newc;
      //   p = newp;
      //   currposterior = newposterior;
      // }
      if (runif(gen) < exp(posteriorcpNew-posteriorcpCurr)) {
        // Rcout << "Proposal posterior: " <<  posteriorcpNew << " for new c and p " <<  cNew << "  " << pNew << " accepted over current posterior " << posteriorcpCurr << " for old c and p " << cCurr  << "  " <<   pCurr << endl;
        cCurr = cNew;
        pCurr = pNew;
        posteriorcpCurr = posteriorcpNew;
      }
    }
    // cCurr = 0.05;
    // pCurr = 1.08;
    // Add the current (thinned out) c, p to their vectors
    // cs.push_back(c);
    // ps.push_back(p);
    cs.push_back(cCurr);
    ps.push_back(pCurr);
    // Rcout << "c and p " <<  c << "  " << p << " chosen with log posterior " << currposterior << endl;
    
    // Set up Hevals with the current c and p
    Hevals.clear();
    // For each eqk...
    for (int i=0; i<n; i++) {
      // Calculate the kappa function (GR aftershock contribution) for each eqk
      // Hevals.push_back( (pow(maxT-ts[i] +c, 1-p) - pow(c, 1-p)) / (1-p) );
      Hevals.push_back( (pow(maxT-ts[i] +cCurr, 1-pCurr) - pow(cCurr, 1-pCurr)) / (1-pCurr) );
      
    }
    
    // Get the posterior value for p(d, q) using the current d and q
    // currposterior = STgBranchingPosteriorMS(lons, lats, marks, dists, maxT, kappaevals, Hevals, cat_Gevals_grid, d, q);
    // currposterior = STgBranchingPosteriorMS(lons, lats, marks, dists, maxT, kappaevals, Hevals, Gevals, d, q);
    posteriordqCurr = STgBranchingPosteriorMSPrior(lons, lats, marks, dists, maxT, kappaevals, Hevals, Gevals,
                                                   dCurr, dpgamma, dhyper1, dhyper2,
                                                   qCurr, qpgamma, qhyper1, qhyper2);
    // // currposterior = STgBranchingPosteriorMS(lons, lats, marks, dists, maxT, kappaevals, Hevals, d, q);
    // // Rcout << "Current posterior: " <<  currposterior << endl;
    //
    // // Rcout <<  "currposterior d/q posterior currposterior d/q posterior" << currposterior << endl;
    // For each MCMC sample...
    for (int i=0; i < numMCMCSamples; i++) {
      // Rcout << "MCMC sample:MCMC sample:MCMC sample:MCMC sample:MCMC sample: " << i << endl;
      // Generate proposal c and p value with the 1-d RW
      // newd = d + rnorm_d(gen);
      // newq = q + rnorm_q(gen);
      dNew = dCurr + rnorm_d(gen);
      qNew = qCurr + rnorm_q(gen);
      
      // if (newd <= 0) {
      //   newd = 0.001;
      // }
      //
      // if (newd > 8) {
      //   newd = 8;
      // }
      //
      // if (newq <= 1) {
      //   newq = 1.001;
      // }
      //
      // if (newq > 8) {
      //   newq = 8;
      // }
      // Rcout << "New d and q " <<  newd << "  " << newq  << endl;
      
      // Generate a value from the posterior distn using the proposal c and p
      // newposterior = hBranchingPosteriorInteraction(ts,marks,z,maxT,kappaevals,newd,newq);
      // newposterior = STgBranchingPosteriorMS(lons, lats, marks, dists, maxT, kappaevals, Hevals, cat_Gevals_grid, newd, newq);
      // newposterior = STgBranchingPosteriorMS(lons, lats, marks, dists, maxT, kappaevals, Hevals, Gevals, newd, newq);
      
      if(dNew > 0){
        if(qNew > 0){
          posteriordqNew = STgBranchingPosteriorMSPrior(lons, lats, marks, dists, maxT, kappaevals, Hevals, Gevals,
                                                        dNew, dpgamma, dhyper1, dhyper2,
                                                        qNew, qpgamma, qhyper1, qhyper2);
        }
      }
      
      // Rcout <<  "new d/q posterior new d/q posterior" << newposterior << endl;
      // newposterior = STgBranchingPosteriorMS(lons, lats, marks, dists, maxT, kappaevals, Hevals, newd, newq);
      // Rcout << "mu: " <<  mu << ", K: " <<  K << ", alpha: " <<  alpha << ", c: " <<  c << ", p: " <<  p << ", d: " <<  d << ", q: " <<  q << "Proposal posterior:Proposal posterior:Proposal posterior:Proposal posterior: " <<  newposterior <<  ", using newd/newq: " << newd << ", " << newq <<  ", comparison to old posterior: " << exp(newposterior-currposterior) << " for MCMC sample " << i << " in post draw " << s << endl;
      // newposterior = SThBranchingPosteriorMS(ts,marks,z, Hevals,kappaevals, Gevals, newc,newp);
      // Rcout << "Proposal posterior: " <<  posteriordqNew << " for new d and q " <<  dNew << "  " << qNew << " accepted over current posterior " << posteriordqCurr << " for old d and q " << dCurr  << "  " <<   qCurr << endl;
      
      // Metropolis-Hastings step.
      // Generate a random unif number and accept the proposal c and p and posterior if it
      // is less than the ratio of the new posterior and current posterior
      // ratio of logged posteriors = exponential of differences of logged posteriors
      // if (runif(gen) < exp(newposterior-currposterior)) {
      //   // Rcout << "Proposal posterior: " <<  newposterior << "for new d and q " <<  newd << "  " << newq << endl;
      //   d = newd;
      //   q = newq;
      //   currposterior = newposterior;
      //   // Rcout << "EUREKA  for MCMC sample " << i << " in post draw " << s << endl;
      //
      // }
      if (runif(gen) < exp(posteriordqNew-posteriordqCurr)) {
        // Rcout << "Proposal posterior: " <<  posteriordqNew << "for new d and q " <<  dNew << "  " << qNew << "accepted over current posterior" << posteriordqCurr << " for old d and q " << dCurr<< "  "  <<  qCurr << endl;
        dCurr = dNew;
        qCurr = qNew;
        posteriordqCurr = posteriordqNew;
        // Rcout << "EUREKA  for MCMC sample " << i << " in post draw " << s << endl;
      }
    }
    // dCurr = 1;
    // qCurr = 2;
    // Add the current (thinned out) c, p to their vectors
    // ds.push_back(d);
    // qs.push_back(q);
    ds.push_back(dCurr);
    qs.push_back(qCurr);
    
    // Reinterpolate Gevals based on the updated d and q
    Gevals.clear();
    // For each eqk...
    for (int i=0; i<n; i++) {
      // Calculate the kappa function (GR aftershock contribution) for each eqk
      // double Gintegral = interpolateGevalsRCpp(d, q, cat_Gevals_grid, i);
      double Gintegral = interpolateGevalsRCpp(dCurr, qCurr, cat_Gevals_grid, i);
      // double Gintegral = 0.99;
      // double Gintegral = 1;
      // Rcout <<  "Gintegral for eqk " << i << " is " << Gintegral << endl;
      Gevals.push_back(Gintegral);
    }
    
    // Generate a draw from the (K, alpha) posterior using these Hevals and values
    // currposterior = STkappaBranchingPosteriorMS(ts,marks,numtriggered,M0, Hevals, Gevals, logK,alpha);
    posteriorKalphaCurr = STkappaBranchingPosteriorMSPrior(ts,marks,numtriggered,M0, Hevals, Gevals,
                                                           KCurr, Kpgamma, Khyper1, Khyper2,
                                                           alphaCurr, alphapgamma, alphahyper1, alphahyper2);
    double logKCurr = log(KCurr);
    // Rcout <<  "got K alpha posteriorgot K alpha posteriorgot K alpha posteriorgot K alpha posteriorgot K alpha posterior" << endl;
    // For each MCMC sample...
    for (int i=0; i < numMCMCSamples; i++) {
      // Generate proposal logK and alpha value with the 1-d RW
      // Rcout << i << endl;
      // newlogK = logK + rnorm_K(gen);
      // newalpha = alpha + rnorm_K(gen);
      logKNew = logKCurr + rnorm_K(gen);
      // logKNew = logKCurr; // Fix K
      KNew = exp(logKNew);
      if(alphaFixed){
        // Rcout << << "!! alphaFixed is = " << alphaFixed << endl;
        alphaNew = alphaCurr; // fix alpha
      } else{
        alphaNew = alphaCurr + rnorm_K(gen); // alpha free
      }
      
      // if (newlogK <= -10000000) {
      //   newlogK = -100000;
      // }
      //
      // if (newalpha <= 0) {
      //   newalpha = 0.001;
      // }
      
      if(KNew > 0){
        if(alphaNew > 0){
          posteriorKalphaNew = STkappaBranchingPosteriorMSPrior(ts,marks,numtriggered,M0, Hevals, Gevals,
                                                                KNew, Kpgamma, Khyper1, Khyper2,
                                                                alphaNew, alphapgamma, alphahyper1, alphahyper2);
          
        }
      }
      // Rcout << "Proposal posterior: " <<  posteriorKalphaNew << " for new K and alpha " <<  KNew << "  " << alphaNew << " accepted over current posterior " << posteriorKalphaCurr << " for old K and alpha " << KCurr  << "  " <<   alphaCurr << endl;
      
      // newposterior = STkappaBranchingPosteriorMS(ts,marks,numtriggered,M0, Hevals, Gevals, newlogK,newalpha);
      // Rcout <<  "new K alpha posteriornew K alpha posteriornew K alpha posteriornew K alpha posteriornew K alpha posteriornew K alpha posteriornew K alpha posterior: " << newposterior  << "for logK" << newlogK << "and alpha: "<<  newalpha <<endl;
      // Metropolis-Hastings step.
      // Generate a random unif number and accept the proposal logK and alpha and posterior if it
      // is less than the ratio of the new posterior and current posterior
      // Need to take exponential because the posterior valeus come back logged?
      // if (runif(gen) < exp(newposterior-currposterior)) {
      //   logK = newlogK;
      //   alpha = newalpha;
      //   currposterior = newposterior;
      //   // Rcout <<  "replaced K alpha posteriorreplaced K alpha posteriorreplaced K alpha posteriorreplaced K alpha posteriorreplaced K alpha posterior" << endl;
      // }
      if (runif(gen) < exp(posteriorKalphaNew-posteriorKalphaCurr)) {
        // Rcout << "Proposal posterior: " <<  posteriorKalphaNew << "for new K and alpha " <<  KNew << "  " << alphaNew << " accepted over current posterior " << posteriorKalphaCurr << "for old K and alpha" << KCurr<< "  "  <<  alphaCurr << endl;
        logKCurr = logKNew;
        KCurr = KNew;
        alphaCurr = alphaNew;
        posteriorKalphaCurr = posteriorKalphaNew;
        // Rcout <<  "replaced K alpha posteriorreplaced K alpha posteriorreplaced K alpha posteriorreplaced K alpha posteriorreplaced K alpha posterior" << endl;
      }
    }
    // Add the current (thinned out) logK, alpha to their vectors
    // KCurr = 0.02;
    // alphaCurr = 1.7;
    Ks.push_back(KCurr);
    alphas.push_back(alphaCurr);
    
    // Rcout <<  "Finished kappa evalsFinished kappa evalsFinished kappa evalsFinished kappa evalsFinished kappa evals" << endl;
    
    // Rcout << "d and q " <<  d << "  " << q << " chosen with log posterior " << currposterior << endl;
    // Rcout << "mu: " <<  mu * pow(500,2)/pow(1000,2) * pow(500,2) << ", K: " <<  K << ", alpha: " <<  alpha << ", c: " <<  c << ", p: " <<  p << ", d: " <<  d << ", q: " <<  q  << endl;
    // Rcout << "mu: " <<  mu << ", mu scaled: " << mu * pow(1000,2)/pow(500,2) * pow(1000,2) << ", K: " <<  K << ", alpha: " <<  alpha << ", c: " <<  c << ", p: " <<  p << ", d: " <<  d << ", q: " <<  q  << endl;
    // Rcout << "mu: " <<  mu << ", mu scaled1: " << mu * pow(1000,2)/pow(500,2) * pow(1000,2) << " mu scaled2:" << mu * 4 <<  ", K: " <<  K << ", alpha: " <<  alpha << ", c: " <<  c << ", p: " <<  p << ", d: " <<  d << ", q: " <<  q  << endl;
    // Rcout << "mu: " <<  mu * pow(1000,2)/pow(500,2) * pow(1000,2) << ", K: " <<  K << ", alpha: " <<  alpha << ", c: " <<  c << ", p: " <<  p << ", d: " <<  d << ", q: " <<  q  << endl;
    
    // 8/2 attempts
    if (s % 100 == 0) {
      Rprintf("Generated %d samples so far...Generated %d samples so far...Generated %d samples so far...\n",s);
    }
  }
  List L = List::create(_["mu"] = mus, _["K"] = Ks, _["alpha"] = alphas ,
                        _["c"] = cs, _["p"] = ps,  _["d"] = ds , _["q"] = qs,
                        _["nbkgd"] = nbkgds,  _["branching"] = branching,
                        _["z"] = z, _["dists"] = dists);
  return L;
  
}
