/* sample.c */

// AM & JT, 10/01/2009
// This program is designed to evaluate parameters of the model, evidence logZ and information H by nested sampling.
// The models are built for evaluating dG,a,b for receptor-ligand sensors and also for estimating mixing ratios alpha,beta,..
// in multiple-ligand mixtures. However, the nesting sampling approach is general and can be easily adapted to other
// applications.

#include "options.h"
#include "nested.h"

int main(int argc, char **argv) {
  // (A) ==== Read in command-line options: ==== //
  if (argc < 2) {
    fprintf(stderr,"usage: %s [-f <Input data file>] [-model <model type>] [-gapfile <file with concentration dilution steps.] [-mfile <file with model parameters>] [-pfile <file with priors>] [-n <number of nested objects>] [-max <max number of nesting iterations>] [-mc <number of MCMC trials for finding a new object>] [-D <size of unweighted posterior sample>] [-out <file with unweighted posterior sample>] [-rs <random seed>] \n",argv[0]);
    exit(1);
  }
  int argbase = 1;
  Options *opt = ReadArguments(argbase, argc, argv);
  
  // (B) ==== Read in data and model parameters from files ==== //
  DATAtable dt(opt->file);
  DELTAtable dxt(opt->gapfile);
  Mtable mt(opt->mfile);
  Ptable pt(opt->pfile);
 
  // (C) ==== Carry out nested sampling ==== //
  Nested nst(dt, dxt, mt, pt, opt->model, opt->n, opt->max, opt->mc, opt->rs);
  nst.Results(dt, dxt);
  if (opt->outfile != "none") {
    nst.MCPosterior(dt, dxt, opt->D, opt->outfile, "adaptive"); // adaptive/fixed
  }
  
}

