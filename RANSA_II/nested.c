/* nested.c */

// AM & JT, 10/01/2009

#include "nested.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

gsl_matrix *Inverse(gsl_matrix *m, int N){
	gsl_matrix *inverse = gsl_matrix_alloc(N, N);
        gsl_permutation * permut = gsl_permutation_alloc(N);
        int signum;
        gsl_matrix * HessianLU = gsl_matrix_alloc(N, N);  
        gsl_matrix_memcpy(HessianLU, m);
        gsl_linalg_LU_decomp(HessianLU, permut, &signum);
        gsl_linalg_LU_invert(HessianLU, permut, inverse);
        gsl_matrix_free(HessianLU);
    	gsl_permutation_free(permut);
    	return inverse;	
}


Nested::Nested(void) {
  r = NULL;
  Obj = NULL;
  Samples = NULL;
  Post = NULL;
  Obj_size = 0;
  Samples_size = 0;
  Post_size = 0;
  MCMC = 0;
  H = 0.0;
  logZ = -DBL_MAX;
}

Nested::~Nested() {
  delete [] Obj;
  delete [] Samples;
  delete [] Post;
  gsl_rng_free (r); // free random number generator
}

Nested::Nested(DATAtable &dt, DELTAtable &dxt,Mtable &mt,Ptable &pt, std::string const &mtype, int nset, int maxnest, int mcmc, unsigned long int rs) {

  // Set up a random number generator:
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  //const gsl_rng_type *T = gsl_rng_randu;  // can use any of the 20+ rng's here, see GSL manual for details
  r = gsl_rng_alloc (T);
  gsl_rng_set (r, rs); // use a different random seed if necessary
  
  MCMC = mcmc;
  
  Obj = new class Mobject[nset]; // objects in the current nested set
  Obj_size = nset;
  for (int j=0; j<nset; ++j) {
    Mobject mobj(dt,dxt,mt,pt,mtype,r);
    Obj[j].copy(mobj);
  }
  Samples = new class Mobject[maxnest]; // all objects kept for posterior analysis
  Samples_size = maxnest;
  Post = NULL; // do not allocate until MCPosterior() is called
  Post_size = 0;
  ////
  H = 0.0;
  logZ = -DBL_MAX;
  ////
  double logwidth; // ln(width in prior mass)
  double logLstar; // ln(Likelihood constraint)
  double logZnew; // updated logZ
  ////
  // Outermost interval of prior mass
  logwidth = log(1.0 - exp(-1.0/nset) );
  // NESTED SAMPLING LOOP:
  for (int nest = 0; nest < maxnest; ++nest) {
    if (nest%1000 == 0) {
    	printf("nested iteration # %d\n",nest+1);
    }
    int worst = 0; // worst object in collection
    for (int i = 1; i < nset; ++i) {
      if (Obj[i].logL < Obj[worst].logL)  worst = i;
    }
    Obj[worst].logWt = logwidth + Obj[worst].logL;
   
    // Update Evidence Z and Information H:
    logZnew = PLUS(logZ, Obj[worst].logWt);
    H = exp(Obj[worst].logWt - logZnew) * Obj[worst].logL
      + exp(logZ - logZnew) * (H + logZ) - logZnew;
    logZ = logZnew;
     
    // Posterior Samples:
    Samples[nest].copy(Obj[worst]);
    // Kill worst object in favor of a copy of different survivor:
    int copy;
    do {
      //copy = (int) (nset * UNIFORM) % nset; // force 0 <= copy < nset
      copy = (int) (nset * gsl_rng_uniform_pos(r) ) % nset; // force 0 <= copy < nset
    } while( copy == worst && nset > 1 );  // don't kill if nset is only 1
    logLstar = Obj[worst].logL; // new likelihood constraint
    Obj[worst].copy(Obj[copy]); // overwrite worst object
    Obj[worst].Evolve(dt,dxt,logLstar,mcmc); // evolve copied object within constraint
    // Shrink interval:
    logwidth -= 1.0 / nset;
    
    } // END OF NESTED SAMPLING LOOP
}

// This function outputs results of the nested sampling run:
// H and logZ estimates plus mean +- stdev for each model prm.

void Nested::Results(DATAtable &dt, DELTAtable &dxt) {

  std::string mtype = Samples[0].model_type;
  printf("======================\n");
  printf("Model type: %s\n",mtype.c_str() );
  printf("======================\n");
  printf("max # nested iterations = %d, nH = %g\n", Samples_size, Obj_size*H);
  printf("# nested objects = %d\n", Obj_size);
  printf("# mcmc iterates for object evolution = %d\n", MCMC);
  printf("Evidence: ln(Z) = %g +- %g\n", logZ, sqrt(H/Obj_size) );
  printf("Information: H = %g [or %g bits]\n", H, H/log(2.) );
  
  ParamStats(dt, dxt, Samples, Samples_size, "weighted");
}

//This function outputs dG, A, B prm values onto screen
void Nested::PrintPrms(std::string const &modeltype, int Nprms, int Ngab, int Na, int Ns, int Nmu, double *means, double *stdev, double *log_means, double *log_stdev){
    if (Nprms != Ngab + Na + Ns + Nmu) {
      fprintf(stderr,"!!!Number of model prms: %d does not match the model type: %s!\nExiting ..\n",Nprms, modeltype.c_str());
      //cerr << "Number of model prms: " << Nprms << " does not match the model type: " << mode << "!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    if (modeltype.substr(0,6) != "NONLIN" && modeltype.substr(0,6) != "LINEAR" ) {
      for (unsigned int j=0; j< Ngab; j+=3) {
        printf("dG%1d = %g +- %g\n", j/3+1, means[j], stdev[j]);
        printf("A%1d = %g +- %g\n", j/3+1, means[j+1], stdev[j+1]);
        printf("B%1d = %g +- %g\n", j/3+1, means[j+2], stdev[j+2]);
      }
    }else if (modeltype.substr(0,6) == "LINEAR") {
      for (unsigned int j=0; j< Ngab; j+=3) {
        printf("V_0%1d = %g +- %g\n", j/3+1, means[j], stdev[j]);
        printf("A%1d = %g +- %g\n", j/3+1, means[j+1], stdev[j+1]);
        printf("B%1d = %g +- %g\n", j/3+1, means[j+2], stdev[j+2]);
      }
    }
    else if (modeltype.substr(0,12) == "NONLIN_CALIB") {
      for (unsigned int j=0; j< Ngab; j+=6) {
        printf("V_0%1d = %g +- %g\n", j/6+1, means[j], stdev[j]);
        printf("D%1d = %g +- %g\n", j/6+1, means[j+1], stdev[j+1]);
        printf("p%1d = %g +- %g\n", j/6+1, means[j+2], stdev[j+2]);
        printf("C%1d = %g +- %g\n", j/6+1, means[j+3], stdev[j+3]);
        printf("B%1d = %g +- %g\n", j/6+1, means[j+4], stdev[j+4]);
        printf("q%1d = %g +- %g\n", j/6+1, means[j+5], stdev[j+5]);
       // printf("V_0 = %g +- %g\n", means[j], stdev[j]);
      }
      //for (unsigned int j=0; j< Ngab; j+=6) {
        //printf("B%1d = %g +- %g\n", j/6+1, means[j+2], stdev[j+2]);
      //}
    }
    else if (modeltype.substr(0,11) == "NONLIN_PRED") {
      for (unsigned int j=0; j< Ngab; j+=6) {
        printf("a%1d = %g +- %g\n", j/6+1, means[j], stdev[j]);
        printf("b%1d = %g +- %g\n", j/6+1, means[j+1], stdev[j+1]);
        printf("p%1d = %g +- %g\n", j/6+1, means[j+2], stdev[j+2]);
        printf("a'%1d = %g +- %g\n", j/6+1, means[j+3], stdev[j+3]);
        printf("b'%1d = %g +- %g\n", j/6+1, means[j+4], stdev[j+4]);
        printf("c'%1d = %g +- %g\n", j/6+1, means[j+5], stdev[j+5]);
       // printf("V_0 = %g +- %g\n", means[j], stdev[j]);
      }
      //for (unsigned int j=0; j< Ngab; j+=6) {
        //printf("B%1d = %g +- %g\n", j/6+1, means[j+2], stdev[j+2]);
      //}
    }
    /*for (unsigned int j=0; j< Ngab; j+=3) {
      //printf("dG%1d = %g +- %g\n", j/3+1, means[j], stdev[j]);
      printf("A%1d = %g +- %g\n", j/3+1, means[j+1], stdev[j+1]);
      printf("B%1d = %g +- %g\n", j/3+1, means[j+2], stdev[j+2]);
    }*/
    
    for (unsigned int j=Ngab; j< Ngab+Na; ++j) {
      //printf("test!!\n");
      printf("alpha%1d = %g +- %g\n", j-Ngab+1, means[j], stdev[j]);
      printf("log(alpha%1d) = %g +- %g\n", j-Ngab+1, log_means[j], log_stdev[j]);
    }
    if (modeltype.substr(0,9) != "FULL_HESS") {
     for (unsigned int j=Ngab+Na; j< Ngab+Na+Ns; ++j) {
       printf("sigma%1d = %g +- %g\n",j-(Ngab+Na)+1, means[j], stdev[j]);
       printf("log(sigma%1d) = %g +- %g\n",j-(Ngab+Na)+1, log_means[j], log_stdev[j]);
     }
    }
    if (modeltype.substr(0,5) == "MuTOT" || modeltype.substr(0,9) == "FULL_HESS" ) {
   //if (modeltype.substr(0,5) == "MuTOT" || modeltype.substr(0,9) == "FULL_HESS" || modeltype.substr(0,10) == "gas_sensor") {
      for (unsigned int j=Ngab+Na+Ns; j<Ngab+Na+Ns+Nmu; j++) {
        printf("mutot = %g +- %g\n", means[j], stdev[j]);  
      }
    }
    else if (modeltype.substr(0,6) == "NONLIN" || modeltype.substr(0,6) == "LINEAR" ) {
      for (unsigned int j=Ngab+Na+Ns; j<Ngab+Na+Ns+Nmu; j++) {
        printf("[TOTAL]_max = %g +- %g\n", means[j], stdev[j]);  
      }
    }
}




void Nested::ParamStats(DATAtable &dt, DELTAtable &dxt, Mobject *mobj, unsigned int mobj_size, std::string const &mode) {

  gsl_matrix *inverse=NULL; 
  
  unsigned int mode_flag;
  if (mode == "weighted") {
    mode_flag = 1;
  }else if (mode == "unweighted") {
    mode_flag = 2;
  }else {
    fprintf(stderr,"Unknown mode: %s in Nested::ParamStats!\nExiting ..\n",mode.c_str() );
    exit(EXIT_FAILURE);
  }
  double det = mobj[mobj_size-1].logL;
  unsigned int Nprms = mobj[0].Nprms; // Nprms is the same in all input Mobjects
  unsigned int Nrec = mobj[0].Nrec; // Nrec is the same in all input Mobjects
  unsigned int Nlig = mobj[0].Nlig; // Nlig is the same in all input Mobjects
  std::string mtype = mobj[0].model_type;
  
  // Compute mu & sigma for all prms:
  double *means = new double[Nprms];
  double *stdev = new double[Nprms];
  // Compute mu & sigma for log(alpha) and log(sigma):
  double *log_means = new double[Nprms]; // Nrec alphas + Nrec sigma
  double *log_stdev = new double[Nprms];
  
  for (unsigned int i = 0; i < Nprms; ++i) {
    means[i] = 0.0;
    stdev[i] =  0.0;
    log_means[i] = 0.0;
    log_stdev[i] =  0.0;
  }

  
  //debug:  (make sure weights sum to 1.0)
  double ww;
    ww = 0;
  for(int i = 0; i < mobj_size; ++i) {
  
    double w;
    
    if (mode_flag == 1) {
      //w = exp(mobj[i].logWt);
      w = exp(mobj[i].logWt - logZ);
      //printf("%d, w = %g\n",i,w);
     // ww += w;
    }else if (mode_flag == 2) {
      w = 1.0/mobj_size;
    }
    ww += w;
  }  
  printf("\n\nww = %e\n\n",ww);
  
    
  for(int i = 0; i < mobj_size; ++i) {
    double w;
    
    if (mode_flag == 1) {
     //w = exp(mobj[i].logWt);
     w = exp(mobj[i].logWt - logZ);
      //printf("%d, w = %g\n",i,w);
      //ww += w;
    }else if (mode_flag == 2) {
      w = 1.0/mobj_size;
    }
    
    for (unsigned int j = 0; j < Nprms; ++j) {
      means[j] += w/ww * mobj[i].prms[j];
      stdev[j] += w/ww * mobj[i].prms[j] * mobj[i].prms[j];
      double logval = log(mobj[i].prms[j]);
      log_means[j] += w/ww * logval;
      log_stdev[j] += w/ww * logval * logval;
     
    }
  }
  //debug:
  /*for (unsigned int j = 0; j < Nprms; ++j) {
   printf("means[%d] = %g\n",j,means[j]);
   }
   printf("[NH3]_max = %g +- %g\n", mobj[mobj_size-1].prms[8]);
  */
  for (unsigned int j = 0; j < Nprms; ++j) {
    double var = stdev[j] - means[j] * means[j];
    if (var < 1e-08) {
      printf("WARNING: negative or small variance for prm # %d: %g!\nResetting to 0.0 ..\n",j+1,var);
      var = 0.0;
    }
    stdev[j] = sqrt(var);
  
    double logvar = log_stdev[j] - log_means[j] * log_means[j];
    if (logvar < 1e-08) {
      printf("WARNING: negative or small log-variance for prm # %d: %g!\nResetting to 0.0 ..\n",j+1,logvar);
      logvar = 0.0;
    }
    log_stdev[j] = sqrt(logvar);
  }

  
  // Model_type dependent output:
  if (mtype == "lighthouse") {
    if (Nprms != 2) {
      fprintf(stderr,"Number of model prms: %d does not match the model type: %s!\nExiting ..\n",Nprms,mtype.c_str() );
      exit(EXIT_FAILURE);
    }
    printf("x = %g +- %g\n", means[0], stdev[0]);
    printf("y = %g +- %g\n", means[1], stdev[1]);
    
  }else if (mtype.substr(0,6) == "LINEAR" && Nlig == 1){
    printf("logL = %g\n\n",det);
    PrintPrms(mtype, Nprms, 3*Nrec*Nlig, Nlig, Nrec, 1, means, stdev, log_means, log_stdev);
  
  }else if (mtype.substr(0,6) == "LINEAR" && Nlig == 2){
    printf("logL = %g\n\n",det);
    PrintPrms(mtype, Nprms, 3*Nrec*Nlig, Nlig-1, Nrec, 1, means, stdev, log_means, log_stdev);
  
  }else if (mtype.substr(0,6) == "NONLIN"){
    printf("logL = %g\n\n",det);
    PrintPrms(mtype, Nprms, 3*Nrec*Nlig, Nlig-1, Nrec, 1, means, stdev, log_means, log_stdev);
  
  }else if (mtype.substr(0,3) == "rec" && mtype.substr(5,3) == "lig") {
    printf("logL = %g\n\n",det);
    PrintPrms(mtype, Nprms, 3*Nrec*Nlig, Nlig-1, 1, 0, means, stdev, log_means, log_stdev);

  }else if (mtype.substr(0,8) == "MULT_SIG") { 
    printf("logL = %g\n\n",det);
    PrintPrms(mtype, Nprms, 3*Nrec*Nlig, Nlig-1, Nrec, 0, means, stdev, log_means, log_stdev);
 
  }else if (mtype.substr(0,5) == "MuTOT") { 
    printf("logL = %g\n\n",det);
    PrintPrms(mtype, Nprms, 3*Nrec*Nlig, Nlig-1, Nrec, 1, means, stdev, log_means, log_stdev);
   
  }else if (mtype.substr(0,9) == "FULL_HESS" && mtype.substr(10,3) == "rec" && mtype.substr(16,3) == "lig") {
    PrintPrms(mtype, Nprms, 3*Nrec*Nlig, Nlig-1, Nrec, 1, means, stdev, log_means, log_stdev);
  
    inverse = Inverse(mobj[mobj_size-1].FullHessian,Nlig);
    
    printf("FullHessian:\n");
    printMatrix(mobj[mobj_size-1].FullHessian);
    printf("determinant: \n%g\n", det);
    printf("Inverse:\n");
    printMatrix(inverse);
    
  }
  
  delete [] means;
  delete [] stdev;
  delete [] log_means;
  delete [] log_stdev;
  gsl_matrix_free(inverse);
}


// This function builds a posterior sample of size D by MC, and outputs it into a file.
void Nested::MCPosterior(DATAtable &dt, DELTAtable &dxt, unsigned int D, std::string const &fileout, std::string const &mode) {
  double step = 0.1;   // Initial guess suitable step-size, given as a fraction of the prior's range
  unsigned int accept = 0;   // # MCMC acceptances
  unsigned int reject = 0;   // # MCMC rejections
  
  unsigned int mode_flag;
  if (mode == "adaptive") {
    mode_flag = 1;
  }else if (mode == "fixed") {
    mode_flag = 2;
  }else{
    fprintf(stderr,"Unknown mode: %s in Nested::MCPosterior!\nExiting ..\n",mode.c_str() );
    exit(EXIT_FAILURE);
  }
  
  // Reset a random number generator if necessary:
  //gsl_rng_set(r, 314159265); // use a different random seed from the NESTED constructor
  
  printf("Building an unweighted posterior sample of size %d ..\n",D);
  if (D == 0) return;
  
  Post = new class Mobject[D]; // allocate memory for the unweighted posterior sample
  Post_size = D;
  Post[0].copy(Obj[gsl_rng_uniform_int(r, Obj_size)]); // init. MC chain with one of the randomly chosen objects from the latest nested set
  accept++;
  
  unsigned int Nprms = Post[0].Nprms;
  
  while (accept < D) {
    Mobject TrialPost;
    TrialPost.copy(Post[accept-1]); // copy the last accepted one
    for (unsigned int i=0; i<Nprms; ++i) {
      if (TrialPost.prms_ran[i] != FALSE_RAN) {
	if (TrialPost.prior_type[i] == 'u' || TrialPost.prior_type[i] == 'j') {
	  TrialPost.prms_ran[i] += step * TrialPost.uniform_prior(-1.0, 1.0, gsl_rng_uniform_pos(r) );  // |move| < step
	  TrialPost.prms_ran[i] -= floor(TrialPost.prms_ran[i]); // wraparound to stay within (0,1)
	}else if (TrialPost.prior_type[i] == 'g') {
	  do {
	    TrialPost.prms_ran[i] += step * gsl_ran_gaussian(r, TrialPost.prms_hi[i]);
	  } while (TrialPost.prms_ran[i] == FALSE_RAN);
	}else if (TrialPost.prior_type[i] == 'l') {
	  do {
	    TrialPost.prms_ran[i] += step * gsl_ran_lognormal(r, 0.0, TrialPost.prms_hi[i]);
	  } while (TrialPost.prms_ran[i] == FALSE_RAN);
	}else{
	  std::cerr << "Unknown prior type: "<<TrialPost.prior_type[i]<<" in Nested::MCPosterior!\nExiting ..\n";
	  exit(EXIT_FAILURE);
	}
      }
    }
 
    for (unsigned int j=0; j<Nprms; ++j) {
      if (TrialPost.prms_ran[j] != FALSE_RAN) {
	if (TrialPost.prior_type[j] == 'u') {
	  TrialPost.prms[j] = TrialPost.uniform_prior(TrialPost.prms_lo[j], TrialPost.prms_hi[j], TrialPost.prms_ran[j]);
	}else if (TrialPost.prior_type[j] == 'j') {
	 // printf("%d\t%g\t%g\n",j,TrialPost.prms_lo[j],TrialPost.prms_hi[j]);
	  TrialPost.prms[j] = TrialPost.jeffreys_prior(TrialPost.prms_lo[j], TrialPost.prms_hi[j], TrialPost.prms_ran[j]);
	}else if (TrialPost.prior_type[j] == 'g') {
	  TrialPost.prms[j] = TrialPost.normal_prior(TrialPost.prms_lo[j], TrialPost.prms_ran[j]); // prms_lo contains Gaussian mean (mu) for x
	}else if (TrialPost.prior_type[j] == 'l') {
	  TrialPost.prms[j] = TrialPost.normal_prior(TrialPost.prms_lo[j], TrialPost.prms_ran[j]); // prms_lo contains Gaussian mean (mu) for log(x)
	}else{
	  std::cerr << "Unknown prior type: "<<TrialPost.prior_type[j]<<" in Nested::MCPosterior!\nExiting ..\n";
	  exit(EXIT_FAILURE);
	}
      }
    }
    //exit(0);
    TrialPost.logLhood(dt, dxt, TrialPost.prms); // this resets logL in TrialPost
    if (TrialPost.logL > Post[accept-1].logL + log(gsl_rng_uniform_pos(r) ) ) { // Accept this object
      Post[accept].copy(TrialPost);
      //printf("accept = %d\n",accept);
      accept++;
    }else{
      reject++;
    }
    // In adaptive mode, refine step-size to let acceptance ratio converge around 50%:
    if (mode_flag == 1) {
      if (accept > reject)  step *= exp(1.0 / accept);
      if (accept < reject)  step /= exp(1.0 / reject);
    }
  }
  
  // Report accept/reject ratio:
  printf("# accepted = %d, # rejected = %d\n",accept,reject);
  
  // Reset the log-weights:
  for (unsigned int i=0; i<Post_size; ++i) {
    Post[i].logWt = 0.0;
  }
  
  // Print means & stdevs estimated from the unweighted sample:
  printf("Parameter estimates from the unweighted posterior sample:\n");
  ParamStats(dt, dxt, Post, Post_size, "unweighted");
  
  // Output unweighted posterior sample into a file:
  MCPrint(fileout);
}

// This function outputs unweighted posterior sample into a file:
void Nested::MCPrint(std::string const &fileout) {
  printf("Outputting posterior sample into file: %s ..\n",fileout.c_str() );
  
  FILE *fp = fopen(fileout.c_str(),"w");
  if (fp == NULL) {
    fprintf(stderr,"Error writing to file %s\n",fileout.c_str() );
    exit(1);
  }
  
  if (Post[0].model_type == "lighthouse") { // lighthouse model
    fprintf(fp,"   %8s   %6s   %6s\n","logL","x","y"); // header
  }else if (Post[0].model_type.substr(0,3) == "rec" && Post[0].model_type.substr(5,3) == "lig") {
    fprintf(fp,"   %8s","logL"); // header
    for (unsigned int j=0; j<Post[0].Nrec*Post[0].Nlig; ++j) {
      fprintf(fp,"   %5s%1d   %5s%1d   %5s%1d","dG",j+1,"A",j+1,"B",j+1);
    }
    for (unsigned int j=0; j<Post[0].Nlig - 1; ++j) {
      fprintf(fp,"   %5s%1d","alpha",j+1);
    }
    fprintf(fp,"   %6s\n","sigma"); // end of header
  }else if (Post[0].model_type.substr(0,9) == "FULL_HESS" && Post[0].model_type.substr(10,3) == "rec" && Post[0].model_type.substr(16,3) == "lig") {
    fprintf(fp,"   %8s","determinant"); // header
    for (unsigned int j=0; j<Post[0].Nrec*Post[0].Nlig; ++j) {
      fprintf(fp,"   %5s%1d   %5s%1d   %5s%1d","dG",j+1,"A",j+1,"B",j+1);
    }
    for (unsigned int j=0; j<Post[0].Nlig - 1; ++j) {
      fprintf(fp,"   %5s%1d","alpha",j+1);
    }
    fprintf(fp,"   %6s\n","sigma"); // end of header
  }else if (Post[0].model_type.substr(0,5) == "MuTOT" || Post[0].model_type.substr(0,4) == "MULT" ) {
    fprintf(fp,"   %8s","logL"); // header
    for (unsigned int j=0; j<Post[0].Nrec*Post[0].Nlig; ++j) {
      fprintf(fp,"   %5s%1d   %5s%1d   %5s%1d","dG",j+1,"A",j+1,"B",j+1);
    }
    for (unsigned int j=0; j<Post[0].Nlig - 1; ++j) {
      fprintf(fp,"   %5s%1d","alpha",j+1);
    }
    fprintf(fp,"   %6s\n","sigma"); // end of header
  }
  /*
  else if (Post[0].model_type == "rec1_lig1") { // one receptor one ligand, a model with common (dG,A,B,sigma) assumed for all datapoints
    fprintf(fp,"logL     dG    A    B    sigma\n"); // header
  }
  */
  for (unsigned int i=0; i<Post_size; ++i) {
    fprintf(fp,"   %8.3f",Post[i].logL);
    for (unsigned int j=0; j<Post[0].Nprms; ++j) {
      fprintf(fp,"   %6.3f",Post[i].prms[j]);
    }
    fprintf(fp,"\n");
  }
  
  fclose(fp);
}

