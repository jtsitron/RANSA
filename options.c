/* options.c */

// AM & JT, 10/01/2009

#include "options.h"

/////////////////////////////////////////////////////////////////////
///// SYNOPSIS //////////////////////////////////////////////////////
///// Command-line arguments are read here //////////////////////////
/////////////////////////////////////////////////////////////////////

Options *ReadArguments(int &argbase, int argc, char **argv)
{
  Options *opt = new Options();
  while (argbase < argc && argv[argbase][0]=='-') {
    if (strstr(argv[argbase],"-f")) {
      argbase++;
      if (argbase >= argc) { // input data file
	printf("Error: -f must be followed by a data file!\n");
	exit(1);
      }
      opt->file = argv[argbase];
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-gapfile")) { // filename with 'delta' values
      //argbase++;
      /*if (argbase >= argc) { // input data file
	printf("Error: -f must be followed by a data file!\n");
	exit(1);
      }*/
      opt->gapfile = argv[++argbase];
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-mfile")) { // filename with model parameters
      opt->mfile = argv[++argbase];
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-pfile")) { // filename with parameters of the priors (ranges, etc.)
      opt->pfile = argv[++argbase];
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-out")) { // output file with unweighted posterior samples
      opt->outfile = argv[++argbase];
      argbase++;
      continue;
    }
    
    if (strstr(argv[argbase],"-model")) { // model type
      opt->model = argv[++argbase];
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-max")) { // max number of nesting iterations
      opt->max = atoi(argv[++argbase]);
      if (opt->max <= 0) {
	printf("Error: -max must be followed by a positive integer!\n");
	exit(1);
      }
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-n")) { // number of objects in the nested set
      opt->n = atoi(argv[++argbase]);
      if (opt->n <= 0) {
	printf("Error: -n must be followed by a positive integer!\n");
	exit(1);
      }
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-mc")) { // number of MCMC steps for finding a new object with L > L^{\star}
      opt->mc = atoi(argv[++argbase]);
      if (opt->mc <= 0) {
	printf("Error: -mc must be followed by a positive integer!\n");
	exit(1);
      }
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-D")) { // number of objects in the unweighted posterior sample
      opt->D = atoi(argv[++argbase]);
      if (opt->D <= 0) {
	printf("Error: -D must be followed by a positive integer!\n");
	exit(1);
      }
      argbase++;
      continue;
    }
    if (strstr(argv[argbase],"-rs")) { // random seed
      int new_rs = atol(argv[++argbase]);
      if (new_rs < 0) {
	printf("Error: -rs must be followed by a non-negative integer!\n");
	exit(1);
      }
      opt->rs = (unsigned long int) new_rs;
      argbase++;
      continue;
    }
    printf("Error: option %s not recognized\n",argv[argbase]);
    exit(1);
  }
  
  return opt;
}

