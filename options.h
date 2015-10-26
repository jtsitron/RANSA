/* options.h */

// AM & JT, 10/01/2009

#ifndef _UTIL_H_
#define _UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

// This header file declares an Options object which parses command-line args.

// WARNING: option name should NOT be a substring of any other option name (e.g. -PDB is a substring of -PDB_id so -PDB
// is not allowed). Otherwise the ReadArguments function will mix them up!!

class Options {
public:
  std::string file; // input data filename
  std::string gapfile; // input filename with 'delta' values
  std::string mfile; // filename with model parameters
  std::string pfile; // filename with parameters of the priors (ranges, etc.)
  std::string outfile; // output file for the unweighted posterior sample
  std::string model; // model type (could be inferred from model prms but we choose to provide it explicitly as a check)
  int max; // max number of nesting iterations
  int mc;  // number of MCMC steps for finding a new object with L > L^{\star}
  int n; // number of objects in the nested set
  int D; // number of unweighted posterior samples
  unsigned long int rs; //random seed form the command line
  Options() { file = "none"; gapfile = "none"; mfile = "none"; pfile = "none"; outfile = "none"; model = "rec1_lig1"; mc = 20; max = 1000; n = 100; D = 1000; rs = 0; }
};

Options *ReadArguments(int &argbase, int argc, char **argv);

#endif

