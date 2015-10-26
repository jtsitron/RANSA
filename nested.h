/* nested.h */

// AM & JT, 10/01/2009

#ifndef _NESTED_H_
#define _NESTED_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <math.h>
#include "mobject.h"

#define PLUS(x,y) (x>y ? x+log(1+exp(y-x)) : y+log(1+exp(x-y))) // logarithmic addition log(exp(x)+exp(y))

// This class is used to carry out nested sampling. It uses an Mobject class which basically evaluates and stores logL
// for a given model, along with current model prms (either fixed or sampled from the prior) and other model-specific
// info (e.g. the log weight of each model). The Nested class stores posterior Samples, information H and evidence logZ,
// and can evaluate various estimated quantities from the posterior sample such as the mean and stdev of each (non-fixed)
// model prm.
class Nested {
 private:
  //gsl_matrix *inverse; //Inverse matrix
  gsl_rng *r; // random number generator
  Mobject *Obj; // nested "model" object set
  Mobject *Samples; // weighted "model" objects stored for posterior results
  Mobject *Post; // unweighted posterior sample of size D
  unsigned int Obj_size, Samples_size, Post_size; // dims of the arrays above
  unsigned int MCMC; // # MC steps passed into Mobject::Evolve
  double H; // Information, initially 0
  double logZ; // ln(Evidence Z), initially Z=0
  void MCPrint(std::string const &fileout);
  void PrintPrms(std::string const &modeltype, int Nprms, int Ngab, int Na, int Ns, int Nmu, double *means, double *stdev, double *log_means, double *log_stdev);
  void ParamStats(DATAtable &dt, DELTAtable &dxt, Mobject *mobj, unsigned int mobj_size, std::string const &mode);
 public:
  Nested();
  Nested(DATAtable &dt,DELTAtable &dxt, Mtable &mt, Ptable &pt, std::string const &mtype, int nset, int maxnest, int mcmc, unsigned long int rs);
  ~Nested();
  void Results(DATAtable &dt,DELTAtable &dxt);
  void MCPosterior(DATAtable &dt, DELTAtable &dxt, unsigned int D, std::string const &fileout, std::string const &mode);
};

gsl_matrix *Inverse(gsl_matrix *m,int N);
#endif
