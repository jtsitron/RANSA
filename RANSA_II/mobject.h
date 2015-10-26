/* mobject.h */

// AM & JT, 10/01/2009

#ifndef _MOBJECT_H_
#define _MOBJECT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <math.h>
#include <float.h>
#include "input.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

//#define UNIFORM ((rand()+0.5)/(RAND_MAX+1.0))  // Uniform inside (0,1)
#define PI 3.14159
#define LOG10 2.302585
#define RT 0.592
//#define RT 0.1
//#define RT 0.6

#define FALSE_RAN -1000000.0 // some fixed # outside of the (0,1) range

// Lighthouse:
#define XLH_LO -2.0
#define XLH_HI 2.0
#define YLH_LO 0.0
#define YLH_HI 2.0
// Gas_sensors:
#define V0_LO -1
#define V0_HI 1

#define D_LO -1
#define D_HI 1

#define p_LO -1
#define p_HI 2

#define C_LO -1
#define C_HI 2

#define a_LO -1
#define a_HI 2

#define b_LO -1
#define b_HI 2

#define aPr_LO -1
#define aPr_HI 2

#define bPr_LO -1
#define bPr_HI 2

#define cPr_LO -1
#define cPr_HI 2


#define TOTAL_LO 0.0001
#define TOTAL_HI 1000.0
// Receptors:
#define dG_LO -20.0
#define dG_HI 5.0


#define A_LO 0.0
#define A_HI 1.0


#define B_LO 0.0
#define B_HI 2.0


#define ALPHA_LO 0.0001
#define ALPHA_HI 100.0
#define SIGMA_LO 0.0001
#define SIGMA_HI 100.0
#define MuTOT_LO -10.0
#define MuTOT_HI -2.0

// This class basically evaluates and stores logL for a given probabilistic model, along with current model prms
// (either fixed or sampled from the prior) and other model-specific
// info (e.g. the log weight of each model).
class Mobject {
 private:
  friend class Nested;
  gsl_matrix *Hessian; // Hessian matrix
  gsl_matrix *MuHessian; // MuHessian matrix
  gsl_matrix *FullHessian; // FullHessian matrix
  gsl_rng *rng; // random number generator (copy of the pointer passed from Nested)
  std::string model_type; // type of the model, determines the formula for computing logL
  double *prms; // set of model prms (whose number and order vary by model type)
  //double *x_min; // minimum value of x coordinate;  needed to compute model y(x')
  double *prms_ran; // set of corresponding (0,1) random numbers
  double *prms_lo; // [lo,hi] prior range for each prm
  double *prms_hi;
  char *prior_type; // 'u' for uniform, 'j' for jeffreys
  double *Ik; // array of Ik values (one for each input datapoint)
  unsigned int Ndat; // Number of datapoints
  unsigned int Nprms; // Number of prms in the current model
  unsigned int Nrec; // Number of receptors
  unsigned int Nlig; // Number of ligands
  double  logL;  // logLikelihood = ln Prob(data | position)
  double  logWt;  // log(Weight), adding to SUM(Wt) = Evidence Z
  void SetUpPrms(DATAtable &dt,DELTAtable &dxt, Mtable &mt, Ptable &pt, int Ngab, int Na, int Ns, int Nmu );
  void EvaluateModel(DATAtable &dt,DELTAtable &dxt, Mtable &mt, Ptable &pt);
  double uniform_prior(double low, double hi, double uf);
  double jeffreys_prior(double low, double hi, double uf);
  double normal_prior(double mu, double gf);
  void logLhood(DATAtable &dt, DELTAtable &dxt, double *prms_loc);
  void writeIk(std::string const &filename, DATAtable &dt, DELTAtable &dxt, double *prms_loc); // prints predicted Ik values into a file in a datafile format
 public:
  // logL & logWt made public for conciseness & efficiency
  //double  logL;  // logLikelihood = ln Prob(data | position)
  //double  logWt; // log(Weight), adding to SUM(Wt) = Evidence Z
  Mobject();
  Mobject(DATAtable &dt, DELTAtable &dxt, Mtable &mt, Ptable &pt, std::string const &mtype, gsl_rng *r);
  ~Mobject();
  void copy(Mobject const &obj); // copy Mobject's datafields into 'this' one
  void Evolve(DATAtable &dt, DELTAtable &dxt, double logLstar, int mcmc); // evolve Mobject by MC with a logL constraint
};

double Determinant(gsl_matrix *m);
void printMatrix(gsl_matrix *m);

//gsl_matrix *Inverse(gsl_matrix *m,int N);

#endif

