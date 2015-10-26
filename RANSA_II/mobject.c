/* mobject.c */

// AM & JT, 10/01/2009

#include "mobject.h"

double Determinant(gsl_matrix *m){
	int n = m->size1; //assuming square
	gsl_permutation * p = gsl_permutation_alloc(n);
	gsl_matrix *tmp = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(tmp, m);
	int signum;
	gsl_linalg_LU_decomp(tmp, p, &signum);
    	double determ = gsl_linalg_LU_det(tmp, signum);
    	gsl_matrix_free(tmp);
    	gsl_permutation_free(p);
    	return determ;	
}


void printMatrix(gsl_matrix *m){
	int i,j;
	for(i = 0; i < m->size1; i++){
		for(j = 0; j < m->size2; j++){
			printf("%g\t", gsl_matrix_get(m, i, j));
		}
		printf("\n");
	}

}
 

Mobject::Mobject(void) {
  model_type = "none";
  prms = NULL;
  //x_min = NULL;
  prms_ran = NULL;
  prms_lo = NULL;
  prms_hi = NULL;
  prior_type = NULL;
  Ik = NULL;
  Nprms = 0;
  Nrec = 0;
  Nlig = 0;
  Ndat = 0;
  logL = -DBL_MAX;
  logWt = -DBL_MAX;
  rng = NULL; // do NOT delete; allocated outside of Mobjects
  Hessian = NULL;
  MuHessian = NULL;
  FullHessian = NULL;
}

Mobject::~Mobject(void) {
  delete [] prms;
  //delete [] x_min;
  delete [] prms_ran;
  delete [] prms_lo;
  delete [] prms_hi;
  delete [] prior_type;
  delete [] Ik;
  if(Hessian != NULL){
	  gsl_matrix_free (Hessian);
  }
  if(MuHessian != NULL){
  	gsl_matrix_free (MuHessian);
  }
  if(FullHessian != NULL){
  	gsl_matrix_free (FullHessian);
  }
}

// Copies data fields from obj into 'this' instantiation of Mobject.
void Mobject::copy(Mobject const &obj) {
  model_type = obj.model_type;
  Nprms = obj.Nprms;
  Nrec = obj.Nrec;
  Nlig = obj.Nlig;
  Ndat = obj.Ndat;
  
  if(Hessian == NULL){
    if (Nlig > 1){
      Hessian = gsl_matrix_alloc (Nlig-1, Nlig-1);
    }
  }
  
  if (Nlig > 1){
    gsl_matrix_memcpy(Hessian, obj.Hessian);
  }
  
  if(MuHessian == NULL){
    MuHessian = gsl_matrix_alloc(1,1);
  }
  
  gsl_matrix_memcpy(MuHessian, obj.MuHessian);


  if(FullHessian == NULL){
    FullHessian = gsl_matrix_alloc(Nlig,Nlig);
  }
  
  gsl_matrix_memcpy(FullHessian, obj.FullHessian);
  
  if (prms == NULL) {
    prms = new double[Nprms];
  }
  //if (x_min == NULL) {
    //x_min = new double[Nrec];
  //}
  if (prms_ran == NULL) {
    prms_ran = new double[Nprms];
  }
  if (prms_lo == NULL) {
    prms_lo = new double[Nprms];
  }
  if (prms_hi == NULL) {
    prms_hi = new double[Nprms];
  }
  if (prior_type == NULL) {
    prior_type = new char[Nprms];
  }
  if (Ik == NULL) {
    Ik = new double[Ndat];
  }
  for (unsigned int i=0; i<Nprms; ++i) {
    prms[i] = obj.prms[i];
    prms_ran[i] = obj.prms_ran[i];
    prms_lo[i] = obj.prms_lo[i];
    prms_hi[i] = obj.prms_hi[i];
    prior_type[i] = obj.prior_type[i];
  }
  for (unsigned int i=0; i<Ndat; ++i) {
    Ik[i] = obj.Ik[i];
  }
  //for (unsigned int i=0; i<Nrec; ++i) {
    //x_min[i] = obj.x_min[i];
  //}  
  logL = obj.logL;
  logWt = obj.logWt;
  rng = obj.rng;
}




Mobject::Mobject(DATAtable &dt, DELTAtable &dxt, Mtable &mt, Ptable &pt, std::string const &mtype, gsl_rng *r) {
  rng = r; // copy state of the GSL random number generator, to be used instead of UNIFORM
  model_type = mtype;
  prms = NULL;
  //x_min = NULL;
  prms_ran = NULL;
  prms_lo = NULL;
  prms_hi = NULL;
  prior_type = NULL;
  Ik = NULL;
  Nprms = 0;
  Nrec = 0;
  Nlig = 0;
  Ndat = 0;
  logL = -DBL_MAX;
  logWt = -DBL_MAX;
  Hessian = NULL;
  MuHessian = NULL;
  FullHessian = NULL;
  
  EvaluateModel(dt,dxt,mt,pt); // will set up all prms and compute logL and logWt
}


void Mobject::SetUpPrms(DATAtable &dt, DELTAtable &dxt, Mtable &mt, Ptable &pt, int Ngab, int Na, int Ns, int Nmu ){

  using std::cerr;
  
  if (pt.get_name() == "none") { // use default prior prms
    if  (model_type.substr(0,6) != "LINEAR" && model_type.substr(0,6) != "NONLIN"){
      for (unsigned int j=0; j<Ngab; j+=3) { // this is for dG,A,B
      //if  (model_type.substr(0,6) != "LINEAR" && model_type.substr(0,6) != "NONLIN"){ 
        prms_lo[j] = dG_LO;
        prms_hi[j] = dG_HI;
        prms_lo[j+1] = A_LO;
        prms_hi[j+1] = A_HI;
        prms_lo[j+2] = B_LO;
        prms_hi[j+2] = B_HI;
      }  
    }else if ( model_type.substr(0,6) == "LINEAR"){
      for (unsigned int j=0; j<Ngab; j+=3) { // this is for V_0,A (slope and y-intercept)
        prms_lo[j] = V0_LO;
        prms_hi[j] = V0_HI;
        prms_lo[j+1] = A_LO;
        prms_hi[j+1] = A_HI; 
      }
    }else if ( model_type.substr(0,12) == "NONLIN_CALIB"){
      for (unsigned int j=0; j<Ngab; j+=6) { // this is for V_0,D,p,C
        prms_lo[j] = V0_LO;
        prms_hi[j] = V0_HI;
        prms_lo[j+1] = D_LO;
        prms_hi[j+1] = D_HI;
        prms_lo[j+2] = p_LO;
        prms_hi[j+2] = p_HI;
        prms_lo[j+3] = C_LO;
        prms_hi[j+3] = C_HI; 
      }
    }else if ( model_type.substr(0,12) == "NONLIN_PRED"){
      for (unsigned int j=0; j<Ngab; j+=6) { // this is for a,b,a',b',c'
        prms_lo[j] = a_LO;
        prms_hi[j] = a_HI;
        prms_lo[j+1] = b_LO;
        prms_hi[j+1] = b_HI;
        prms_lo[j+3] = aPr_LO;
        prms_hi[j+3] = aPr_HI;
        prms_lo[j+4] = bPr_LO;
        prms_hi[j+4] = bPr_HI;
        prms_lo[j+5] = cPr_LO;
        prms_hi[j+5] = cPr_HI; 
      }
    }
    for (unsigned int j=0; j<Ngab; ++j) { // this is for dG,A,B
      prior_type[j] = 'u';
    }
    for (unsigned int j=Ngab; j<Ngab+Na; ++j) { // this is for alpha
      prms_lo[j] = ALPHA_LO;
      prms_hi[j] = ALPHA_HI;
    }
    for (unsigned int j=Ngab+Na; j<Ngab+Na+Ns; j++){ // this is for sigma
      prms_lo[j] = SIGMA_LO;
      prms_hi[j] = SIGMA_HI;
    }
    for (unsigned int j=Ngab+Na+Ns; j<Ngab+Na+Ns+Nmu; ++j) { // this is for mutot/TOTAL
      if (model_type.substr(0,6) == "LINEAR" || model_type.substr(0,12) == "NONLIN_CALIB"){ 
        prms_lo[j] = TOTAL_LO;
        prms_hi[j] = TOTAL_HI;
      }else if (model_type.substr(0,11) == "NONLIN_PRED"){
        int high_delt = 0;
        //int lo_x;
        for (unsigned int k = 0; k < dt.length(); ++k){ // loop over datapoints
          if (high_delt < dxt.get_dx(k)) {
            high_delt = dxt.get_dx(k);
          }
        } 
        //Ik[datapoint_cnt] -= prms_loc[3]*pow((prms_loc[Nprms-1]-dxt.get_dx(k)),prms_loc[2]);
        prms_lo[j] = high_delt;
        //printf("high_delt = %d\n",high_delt);
        //prms_lo[j] = MuTOT_LO;
        prms_hi[j] = TOTAL_HI;
      }else {
        prms_lo[j] = MuTOT_LO;
        prms_hi[j] = MuTOT_HI;
      }
    }
    for (unsigned int j=Ngab; j<Ngab+Na+Ns; ++j) { // this is for alpha,sigma
      //printf("%d\n",Ngab);
      //printf("%d\n",Na);
      //printf("%d\n",Ns);
      //printf("%d\n",j);
      prior_type[j] = 'j';
    }
    for (unsigned int j=Ngab+Na+Ns; j<Ngab+Na+Ns+Nmu; j++){ //this is for mutot
        prior_type[j] = 'u';
    }
    
    //for (unsigned int j=0; j<Nprms; j++){
    //printf("%1d\n",j);
    //fprintf("prior_type[%d] = %s\n",j, prior_type[j]);
    //}
  }else{ // use prior prms from the pfile
    // Range check:
    if (pt.getNtot() != Nprms || pt.getN() != Nrec*Nlig) {
    //if (pt.getNtot() != Nprms) {
      cerr << "ERROR: file with priors must have prior info for exactly "<<Nprms<<" parameters which describe "<<Nrec<<" receptor(s) and "<<Nlig<<" ligands!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    unsigned int cnt_loc = 0;
    for (unsigned int j=0; j<Ngab; j+=3) {
      prms_lo[j] = pt.get_prior_prm("ddG",cnt_loc,0); // 0 for lo, 1 for hi
      prms_hi[j] = pt.get_prior_prm("ddG",cnt_loc,1);
      prms_lo[j+1] = pt.get_prior_prm("A",cnt_loc,0);
      prms_hi[j+1] = pt.get_prior_prm("A",cnt_loc,1);
      prms_lo[j+2] = pt.get_prior_prm("B",cnt_loc,0);
      prms_hi[j+2] = pt.get_prior_prm("B",cnt_loc,1);
      ////
      prior_type[j] = pt.get_prior_type("ddG",cnt_loc);
      prior_type[j+1] = pt.get_prior_type("A",cnt_loc);
      prior_type[j+2] = pt.get_prior_type("B",cnt_loc);
      cnt_loc++;
    }
    cnt_loc = 0;
    for (unsigned int j=Ngab; j<Ngab+Na; ++j) {
      prms_lo[j] = pt.get_prior_prm("alpha",cnt_loc,0);
      prms_hi[j] = pt.get_prior_prm("alpha",cnt_loc,1);
      prior_type[j] = pt.get_prior_type("alpha",cnt_loc);
      cnt_loc++;
    }
    cnt_loc = 0;
    for (unsigned int j=Ngab+Na; j<Ngab+Na+Ns; ++j) {
      prms_lo[j] = pt.get_prior_prm("sigma",cnt_loc,0);
      prms_hi[j] = pt.get_prior_prm("sigma",cnt_loc,1);
      prior_type[j] = pt.get_prior_type("sigma",cnt_loc);
      cnt_loc++;
    }
    if (Nmu == 1){
      prms_lo[Ngab+Na+Ns+Nmu-1] = pt.get_mutot(0);
      prms_hi[Ngab+Na+Ns+Nmu-1] = pt.get_mutot(1);
      prior_type[Ngab+Na+Ns+Nmu-1] = pt.get_mutot_type();
    }
  }
  if (mt.get_name() != "none") { // some prms are explicitly listed in the model file, otherwise all prms are to be fitted
    if (mt.getNmax() > Nrec*Nlig) { // N can be 0 (only sigma listed in the file) or 1 (some of the model prms are listed in the file)
      cerr << "Too many receptor/ligand combinations found in Mtable::Mtable "<<mt.get_name()<<" for "<<model_type<<" model type!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    unsigned int cnt = 0;
    for (unsigned int j=0; j<3*Ngab; j+=3) {
      if (mt.get_ddG_onoff(cnt) == true) { // fixed dG
	prms_ran[j] = FALSE_RAN;
	prms[j] = mt.get_ddG(cnt);
	//printf("V_o = %d\n",prms[j]);
      }
      if (mt.get_A_onoff(cnt) == true) { // fixed A
	prms_ran[j+1] = FALSE_RAN;
	prms[j+1] = mt.get_A(cnt);
	//printf("D = %d\n",prms[j+1]);
      }
      if (mt.get_B_onoff(cnt) == true) { // fixed B
	prms_ran[j+2] = FALSE_RAN;
	prms[j+2] = mt.get_B(cnt);
	//printf("p = %d\n",prms[j+2]);
      }
      cnt++;
      if (cnt >= mt.getNmax() ) break;
    }
    cnt = 0;
    for (unsigned int j=Ngab; j<Ngab+Na; ++j) {
      if (mt.get_alpha_onoff(cnt) == true) { // fixed alpha
	prms_ran[j] = FALSE_RAN;
	prms[j] = mt.get_alpha(cnt);
      }
      cnt++;
      if (cnt >= mt.getNmax() ) break;
    }
    cnt = 0;
    for (unsigned int j=Ngab+Na; j<Ngab+Na+Ns; ++j) {
      if (mt.get_sigma_onoff(cnt) == true) { // fixed sigma
	prms_ran[j] = FALSE_RAN;
	prms[j] = mt.get_sigma(cnt);
      }
      cnt++;
      if (cnt >= mt.getNmax() ) break;
    }
    if (Nmu == 1){
      if (mt.get_mutot_onoff() == true) { // fixed mutot
        prms_ran[Ngab+Na+Ns+Nmu-1] = FALSE_RAN;
        prms[Ngab+Na+Ns+Nmu-1] = mt.get_mutot();
      }
    }
  }
}

void Mobject::EvaluateModel(DATAtable &dt, DELTAtable &dxt, Mtable &mt, Ptable &pt) {
  using std::cerr;
  if (model_type.substr(0,6) != "LINEAR" && model_type.substr(0,12) != "NONLIN_CALIB" && model_type.substr(0,11) != "NONLIN_PRED" && model_type.substr(0,3) != "rec" && model_type.substr(0,5) != "MuTOT" && model_type.substr(0,8) != "MULT_SIG" && model_type.substr(0,9) != "FULL_HESS"){
   cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    
  if (model_type.substr(0,9) != "FULL_HESS") {
   if (dt.get_name() == "none") { // no data available, quit
     cerr << "No data found to evaluate the model with, please check the DATAtable::DATAtable object passed into Mobject::EvaluateModel!\nExiting ..\n";
     exit(EXIT_FAILURE);
   }
  }
  
  if ( model_type.substr(0,5) == "MuTOT" && model_type.find("lig",11) != std::string::npos) { 
    //&& model_type.find("lig",11) != std::string::npos) {
    //model_type.substr(0,6) == "LINEAR" || model_type.substr(0,12) == "NONLIN_CALIB" || model_type.substr(0,11) == "NONLIN_PRED" || model_type.substr(0,3) == "rec" ||
    if (dxt.get_name() == "none") { // no delta values available, quit
      cerr << "No delta values found to evaluate the model with, please check the DELTAtable::DELTAtable object passed into Mobject::EvaluateModel!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }  
  }
  
  
  Ndat = dt.length();
  
  
  // Choose # prms acc. to model_type:
  if (model_type == "lighthouse") {
    Nprms = 2; // order: x y
    Nrec = 0;
    Nlig = 0;
    
  }else if (model_type.substr(0,6) == "LINEAR") {
    if (model_type.length() != 16) {
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one-digit numbers to enumerate receptors!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
   
    Nrec = atoi(model_type.substr(7,1).c_str() ); // one-digit number!!
    Nlig = atoi(model_type.substr(12,1).c_str() ); // one-digit number!!
    //Nlig = 2;
    
    if (Nrec == 0 || Nlig == 0) { // atoi returns 0 if confused
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    if (Nlig == 1) { //calibration
      Nprms = 3*Nrec*Nlig + (Nlig) + Nrec + 1; // order: V_0, A, B, alpha, sigma, [TOTAL]_max  
    }else if (Nlig == 2) {  
      Nprms = 3*Nrec*Nlig + (Nlig-1) + Nrec + 1; // order: (V_0,A,B) for each receptor/ligand pair + Nlig-1 alphas + Nrec sigmas + [TOTAL]_max  
    }
  }else if (model_type.substr(0,12) == "NONLIN_CALIB") {
    if (model_type.length() != 17) {
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one-digit numbers to enumerate receptors!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
   
    Nrec = atoi(model_type.substr(13,1).c_str() ); // one-digit number!!
    Nlig = 2;
    
    if (Nrec == 0 || Nlig == 0) { // atoi returns 0 if confused
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    
    Nprms = 3*Nrec*Nlig + (Nlig-1) + Nrec + 1; // order: (V_0,D,p,C,B,q) for each receptor + Nlig-1 alphas + Nrec sigmas + [TOTAL]_max  
     
  }else if (model_type.substr(0,11) == "NONLIN_PRED") {
    if (model_type.length() != 16) {
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one-digit numbers to enumerate receptors!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
   
    Nrec = atoi(model_type.substr(12,1).c_str() ); // one-digit number!!
    Nlig = 2;
    
    if (Nrec == 0 || Nlig == 0) { // atoi returns 0 if confused
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    
    Nprms = 3*Nrec*Nlig + Nlig-1 + Nrec + 1; // order: (a,b,p,a',b',c') for each receptor + Nlig-1 alphas + Nrec sigmas + [TOTAL]_max  
     
  }else if (model_type.substr(0,3) == "rec" && model_type.find("lig",5) != std::string::npos) {
    if (model_type.length() != 9) {
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one-digit numbers to enumerate receptors and ligands!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    Nrec = atoi(model_type.substr(3,1).c_str() ); // one-digit number!!
    Nlig = atoi(model_type.substr(8,1).c_str() ); // one-digit number!!
    if (Nrec == 0 || Nlig == 0) { // atoi returns 0 if confused
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    
    Nprms = 3 * Nrec * Nlig + (Nlig - 1) + 1;  // order: (dG,A,B) for each receptor-ligand combination + (Nlig - 1) alphas + sigma
    
  }else if (model_type.substr(0,8) == "MULT_SIG" && model_type.find("lig",14) != std::string::npos) {
     if (model_type.length() != 18) {
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one-digit numbers to enumerate receptors and ligands!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    Nrec = atoi(model_type.substr(12,1).c_str() ); // one-digit number!!
    Nlig = atoi(model_type.substr(17,1).c_str() ); // one-digit number!!
    if (Nrec == 0 || Nlig == 0) { // atoi returns 0 if confused
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    
    Nprms = 3 * Nrec * Nlig + (Nlig - 1) + Nrec;  // order: (dG,A,B) for each receptor-ligand combination + (Nlig - 1) alphas + (Nrec) sigmas
   
  }else if (model_type.substr(0,5) == "MuTOT" && model_type.find("lig",11) != std::string::npos) {
     if (model_type.length() != 15) {
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one-digit numbers to enumerate receptors and ligands!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    Nrec = atoi(model_type.substr(9,1).c_str() ); // one-digit number!!
    Nlig = atoi(model_type.substr(14,1).c_str() ); // one-digit number!!
    if (Nrec == 0 || Nlig == 0) { // atoi returns 0 if confused
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    
    Nprms = 3 * Nrec * Nlig + (Nlig - 1) + Nrec + 1;  // order: (dG,A,B) for each receptor-ligand combination + (Nlig - 1) alphas + (Nrec) sigmas + 1 muTOT
    
  }else if (model_type.substr(0,9) == "FULL_HESS" && model_type.find("rec",10) != std::string::npos) {
    if (model_type.length() != 21) {
      cerr << "!!Unknown model type: "<<model_type<<"!\nPlease use two-digit numbers to enumerate receptors and ligands!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    Nrec = atoi(model_type.substr(13,2).c_str() ); // two-digit number!!
    
    Nlig = atoi(model_type.substr(19,2).c_str() ); // two-digit number!!
    if (Nrec == 0 || Nlig == 0) { // atoi returns 0 if confused
      cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nLINEAR_Xrec_Ylig\nNONLIN_CALIB_Xrec\nNONLIN_PRED_Xrec\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    
    Nprms = 3 * Nrec * Nlig + (Nlig - 1) + Nrec + 1;  // order: (dG,A,B) for each receptor-ligand combination + (Nlig - 1) alphas + (Nrec) sigmas + 1 mutot
    
  }else{
    cerr << "Unknown model type: "<<model_type<<"!\nPlease use one of the following modes:\nrecX_ligY\nMULT_SIG_recX_ligY\nMuTOT_recX_ligY\nFULL_HESS_recXX_ligYY\n\nExiting ..\n";
    exit(EXIT_FAILURE);
  }
  
  // Allocate arrays:
  prms = new double[Nprms];
  //x_min = new double[Nrec];
  prms_ran = new double[Nprms];
  prms_lo = new double[Nprms];
  prms_hi = new double[Nprms];
  prior_type = new char[Nprms];
  for (unsigned int j=0; j<Nprms; ++j) {
    prms[j] = 0.0;
    prms_ran[j] = 0.0;
    prms_lo[j] = 0.0;
    prms_hi[j] = 0.0;
    prior_type[j] = '-';
  }
  
  if (Nlig>1){
    Hessian = gsl_matrix_alloc (Nlig-1,Nlig-1);
  }
  
  MuHessian = gsl_matrix_alloc (1, 1);
  
  FullHessian = gsl_matrix_alloc(Nlig,Nlig);
  
  // Set up model parameters:
  if (model_type == "lighthouse") {
    // Problem:
    //  Lighthouse at (x,y) emitted n flashes observed at D[.] on coast.
    // Inputs:
    //  2D position  -2 < x < 2, 0 < y < 2 with flat prior
    //  Likelihood  is L(x,y) = PRODUCT[k] (y/pi) / ((D[k] - x)^2 + y^2)
    // Outputs:
    //  Evidence is Z = INTEGRAL L(x,y) Prior(x,y) dxdy
    //  Posterior is P(x,y) = L(x,y) / Z estimating lighthouse position
    //  Information is H = INTEGRAL P(x,y) log(P(x,y)/Prior(x,y)) dxdy  
    // ************************** //
    // Toy model from Sivia & Skilling: 2 coordinates (x,y) for the lighthouse position
    // NOTE: unlike "real" receptor-ligand models, prms of this problem are not given by an Mtable object
    if (mt.get_name() != "none") { // some file with model prms has been passed in, although "lighthouse" mode does not accept any
      cerr << "Mtable::Mtable "<<mt.get_name()<<" was passed into Mobject::EvaluateModel in the "<<model_type<<" mode!\n";
      cerr << "Please remove the -mfile option from the command line!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    if (pt.get_name() != "none") { // some file with prior prms has been passed in, although "lighthouse" mode does not accept any
      cerr << "Ptable::Ptable "<<pt.get_name()<<" was passed into Mobject::EvaluateModel in the "<<model_type<<" mode!\n";
      cerr << "Please remove the -pfile option from the command line!\nExiting ..\n";
      exit(EXIT_FAILURE);
    }
    //// assign parameters
    prms_lo[0] = XLH_LO; // (-2,2) uniform random number (x)
    prms_hi[0] = XLH_HI;
    prms_lo[1] = YLH_LO; // (0,2) uniform random number (y)
    prms_hi[1] = YLH_HI;
    for (unsigned int j=0; j<Nprms; ++j) {
      prior_type[j] = 'u';
    }
    
  }else if (model_type.substr(0,6) == "LINEAR" && Nlig==1){
    int Ngab = 3*Nrec*Nlig;
    int Na = Nlig;
    int Ns = Nrec;
    int Nmu = 1;
    SetUpPrms(dt,dxt,mt,pt,Ngab,Na,Ns,Nmu);
    
  }else if (model_type.substr(0,6) == "LINEAR" && Nlig==2){
    int Ngab = 3*Nrec*Nlig;
    int Na = Nlig-1;
    int Ns = Nrec;
    int Nmu = 1;
    SetUpPrms(dt,dxt,mt,pt,Ngab,Na,Ns,Nmu);
    
  }else if (model_type.substr(0,6) == "NONLIN"){
    int Ngab = 3*Nrec*Nlig;
    int Na = Nlig-1;
    int Ns = Nrec;
    int Nmu = 1;
    SetUpPrms(dt,dxt,mt,pt,Ngab,Na,Ns,Nmu);
    
  }else if (model_type.substr(0,3) == "rec" && model_type.substr(5,3) == "lig") {
    int Ngab = 3*Nrec*Nlig; // Number of Delta G's, A's, and B's
    int Na = Nlig -1; // Number of \alpha's
    int Ns = 1; // Number of \sigma's
    int Nmu = 0; // Number of \mu's
    SetUpPrms(dt,dxt,mt,pt,Ngab,Na,Ns,Nmu);
    
    
  }else if (model_type.substr(0,8) == "MULT_SIG"){
    int Ngab = 3*Nrec*Nlig; // Number of Delta G's, A's, and B's
    int Na = Nlig -1; // Number of \alpha's
    int Ns = Nrec; // Number of \sigma's
    int Nmu = 0; // Number of \mu's
    SetUpPrms(dt,dxt,mt,pt,Ngab,Na,Ns,Nmu);
    
    
  }else if (model_type.substr(0,5) == "MuTOT" && model_type.substr(6,3) == "rec" && model_type.substr(11,3) == "lig"){
    int Ngab = 3*Nrec*Nlig; // Number of Delta G's, A's, and B's
    int Na = Nlig -1; // Number of \alpha's
    int Ns = Nrec; // Number of \sigma's
    int Nmu = 1; // Number of \mu's
    SetUpPrms(dt,dxt,mt,pt,Ngab,Na,Ns,Nmu);
    
    
  }else if (model_type.substr(0,9) == "FULL_HESS" && model_type.substr(10,3) == "rec" && model_type.substr(16,3) == "lig") {
    int Ngab = 3*Nrec*Nlig; // Number of Delta G's, A's, and B's
    int Na = Nlig -1; // Number of \alpha's
    int Ns = Nrec; // Number of \sigma's
    int Nmu = 1; // Number of \mu's
    SetUpPrms(dt,dxt,mt,pt,Ngab,Na,Ns,Nmu);
    
    
  }
  
  for (unsigned int j=0; j<Nprms; ++j) {
    if (prms_ran[j] != FALSE_RAN) {
      if (prior_type[j] == 'u') { // uniform
	prms_ran[j] = gsl_rng_uniform_pos(rng);
        //printf("TEST%f\n%f\t%f\n",j,prms_lo[j],prms_hi[j]);
	prms[j] = uniform_prior(prms_lo[j], prms_hi[j], prms_ran[j]);
      }else if (prior_type[j] == 'j') { // jeffreys
	prms_ran[j] = gsl_rng_uniform_pos(rng);
	prms[j] = jeffreys_prior(prms_lo[j], prms_hi[j], prms_ran[j]);
      }else if (prior_type[j] == 'g') { // normal
	do {
	  prms_ran[j] = gsl_ran_gaussian(rng, prms_hi[j]); // prms_hi contains Gaussian sigma for this prior_type
	} while (prms_ran[j] == FALSE_RAN);
	prms[j] = normal_prior(prms_lo[j], prms_ran[j]); // prms_lo contains Gaussian mean (mu) for this prior_type
      }else if (prior_type[j] == 'l') { // log-normal, needs mu and sigma in log space
	do {
	  prms_ran[j] = gsl_ran_lognormal(rng, 0.0, prms_hi[j]); // prms_hi contains Gaussian sigma for log(x)
	} while (prms_ran[j] == FALSE_RAN);
	prms[j] = normal_prior(prms_lo[j], prms_ran[j]); // prms_lo contains Gaussian mean (mu) for log(x)
      }else{
	cerr << "Unknown prior type: "<<prior_type[j]<<" for param: "<<j<<" in Mobject::EvaluateModel!\nExiting ..\n";
	exit(EXIT_FAILURE);
      }
    }
  }

  
  // Evaluate logL:
  Ik = new double[Ndat];
  logLhood(dt,dxt,prms);
}



void Mobject::logLhood(DATAtable &dt, DELTAtable &dxt, double *prms_loc) {
  logL = 0.0;
  if (model_type == "lighthouse") {
    for (unsigned int k = 0; k < dt.length();  ++k) {
      logL += log((prms_loc[1]/PI) / ((dt.get_y(k)-prms_loc[0])*(dt.get_y(k)-prms_loc[0]) + prms_loc[1]*prms_loc[1]) );
    }
    
  }else if (model_type.substr(0,6) == "LINEAR"){ 
    unsigned int datapoint_cnt = 0;
    double lin_coeff[Nlig];
    for (unsigned int j = 0;  j < Nrec;  ++j) { // loop over receptor types
      for (unsigned int k = 0; k < dt.length(); ++k){ // loop over datapoints
	
        if (dxt.length() != dt.length()){
          printf("dt.length = %d  dxt.length = %d\n\n",dt.length(), dxt.length());
          std::cerr << "Length of gapfile '"<<dxt.get_name()<<"' does not equal length of datafile '"<<dt.get_name()<<"'!! PLEASE CHECK INPUT FILES!\n\nExiting..\n"; 
          exit(EXIT_FAILURE);
        }
        if (dt.get_recN(k) == j+1) { // only take datapoints for the current receptor type
          double V_0 = 0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    V_0 += prms_loc[3*Nlig*j + 3*l];
	  }
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    if (l == 0){
	      lin_coeff[l] =  prms_loc[3*Nlig*j + 3*l + 1];
	    }
	    else {
              //printf("TEST\n");
	      lin_coeff[l] =  prms_loc[3*Nlig*j + 3*l + 1]*prms_loc[3*Nlig*Nrec];
	    }
	  }
	  Ik[datapoint_cnt] = V_0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    //printf("%g\n",prms_loc[3*Nlig*Nrec]);
	    Ik[datapoint_cnt] += (lin_coeff[l]/(1+prms_loc[3*Nlig*Nrec]))*(prms_loc[Nprms-1]-dxt.get_dx(k));
	  }  
	     //          Ik[datapoint_cnt] = prms_loc[3*j] + prms_loc[3*j+1]*(prms_loc[Nprms-1]-dxt.get_dx(k)) + prms_loc[3*j+2]*prms_loc[3*Nrec*Nlig]*(prms_loc[Nprms-1]-dxt.get_dx(k));
          logL -= (dt.get_y(k) - Ik[datapoint_cnt])*(dt.get_y(k) - Ik[datapoint_cnt]) / (2 * prms_loc[3*Nrec*Nlig+1+j] * prms_loc[3*Nrec*Nlig+1+j]) + log(prms_loc[3*Nrec*Nlig+1+j]);
         // Ik[datapoint_cnt] = prms_loc[2*j]* (prms_loc[Nprms-1]-dxt.get_dx(k)) + prms_loc[2*j +1]*prms_loc[2*Nrec]*(prms_loc[Nprms-1]-dxt.get_dx(k));
         // logL -= ( (dt.get_y(k) - Ik[datapoint_cnt])*(dt.get_y(k) - Ik[datapoint_cnt]) / (2 * prms_loc[2*Nrec+j+1] * prms_loc[2*Nrec+j+1]) ) + log(prms_loc[2*Nrec+j+1]);
          datapoint_cnt++;
        }
      }
    }
        
    
    if (datapoint_cnt == 0) {
      std::cerr << "No data found in Mobject::logLhood, please check recN assignment in "<<dt.get_name()<<"!\n";
      exit(EXIT_FAILURE);       }
    if (datapoint_cnt != Ndat) {
      std::cerr << "Not all data found in "<<dt.get_name()<<" was used, please check recN assignment in the input file!\n";
      exit(EXIT_FAILURE);
    }
        
    logL -= 0.5 * datapoint_cnt * log(2 * PI ); 
    
  }else if (model_type.substr(0,12) == "NONLIN_CALIB"){
    unsigned int datapoint_cnt = 0;
    double lin_coeff[Nlig];
    for (unsigned int j = 0;  j < Nrec;  ++j) { // loop over receptors
      for (unsigned int k = 0; k < dt.length(); ++k){ // loop over datapoints
	
        if (dxt.length() != dt.length()){
          printf("dt.length = %d  dxt.length = %d\n\n",dt.length(), dxt.length());
          std::cerr << "Length of gapfile '"<<dxt.get_name()<<"' does not equal length of datafile '"<<dt.get_name()<<"'!! PLEASE CHECK INPUT FILES!\n\nExiting..\n"; 
          exit(EXIT_FAILURE);
        }
        if (dt.get_recN(k) == j+1) { // only take datapoints for the current receptor
	  Ik[datapoint_cnt] = prms_loc[6*j]; // V_0
	  
	  //                      D * [NH3]_max + C * [NH3]_max^p
	  Ik[datapoint_cnt] += (prms_loc[3*Nlig*j + 1])*(prms_loc[Nprms-1]-dxt.get_dx(k)) + (prms_loc[6*j+3]/pow((prms_loc[3*Nlig*Nrec+j] +  1),prms_loc[2]))*pow((prms_loc[Nprms-1]-dxt.get_dx(k)),prms_loc[2]);
	  
	  //                       C * [NH3]_max^p
	  //Ik[datapoint_cnt] += prms_loc[6*j+3]*pow((prms_loc[Nprms-1]-dxt.get_dx(k)),prms_loc[2]);
	  	
	  logL -= (dt.get_y(k) - Ik[datapoint_cnt])*(dt.get_y(k) - Ik[datapoint_cnt]) / (2 * prms_loc[3*Nrec*Nlig+Nrec+j] * prms_loc[3*Nrec*Nlig+Nrec+j]) + log(prms_loc[3*Nrec*Nlig+Nrec+j]);
	  
	  
         
          datapoint_cnt++;
        }
      }
    }
        
    
    if (datapoint_cnt == 0) {
      std::cerr << "No data found in Mobject::logLhood, please check recN assignment in "<<dt.get_name()<<"!\n";
      exit(EXIT_FAILURE);       }
    if (datapoint_cnt != Ndat) {
      std::cerr << "Not all data found in "<<dt.get_name()<<" was used, please check recN assignment in the input file!\n";
      exit(EXIT_FAILURE);
    }
        
    logL -= 0.5 * datapoint_cnt * log(2 * PI ); 
    
  }else if (model_type.substr(0,11) == "NONLIN_PRED"){ 
    unsigned int datapoint_cnt = 0;
    //double coeff[Nlig];
    
    for (unsigned int j = 0;  j < Nrec;  ++j) { // loop over values of ALPHA/or over receptors
      for (unsigned int k = 0; k < dt.length(); ++k){ // loop over datapoints
	
        if (dxt.length() != dt.length()){
          printf("dt.length = %d  dxt.length = %d\n\n",dt.length(), dxt.length());
          std::cerr << "Length of gapfile '"<<dxt.get_name()<<"' does not equal length of datafile '"<<dt.get_name()<<"'!! PLEASE CHECK INPUT FILES!\n\nExiting..\n"; 
          exit(EXIT_FAILURE);
        }
        if (dt.get_recN(k) == j+1) { // only take datapoints for the current receptor
          Ik[datapoint_cnt] = 0.0;
          
	  Ik[datapoint_cnt] += ((prms_loc[3*Nlig*j + 3] + prms_loc[3*Nlig*j + 3 + 1]*prms_loc[3*Nlig*Nrec] + prms_loc[3*Nlig*j + 3 + 2]*pow(prms_loc[3*Nlig*Nrec],2))/(pow((1+prms_loc[3*Nlig*Nrec]),(prms_loc[3*Nlig*j] + prms_loc[3*Nlig*j + 1]*prms_loc[3*Nlig*Nrec]))))* pow((prms_loc[Nprms-1]-dxt.get_dx(k)),prms_loc[3*Nlig*j] + prms_loc[3*Nlig*j + 1]*prms_loc[3*Nlig*Nrec] + prms_loc[3*Nlig*j + 2]*pow(prms_loc[3*Nlig*Nrec],2));
	  
	  	
	  logL -= (dt.get_y(k) - Ik[datapoint_cnt])*(dt.get_y(k) - Ik[datapoint_cnt]) / (2 * prms_loc[3*Nrec*Nlig+Nlig-1+j] * prms_loc[3*Nrec*Nlig+Nlig-1+j]) + log(prms_loc[3*Nrec*Nlig+Nlig-1+j]);
         
          datapoint_cnt++;
        }
      }
    }
        
    
    if (datapoint_cnt == 0) {
      std::cerr << "No data found in Mobject::logLhood, please check recN assignment in "<<dt.get_name()<<"!\n";
      exit(EXIT_FAILURE);       }
    if (datapoint_cnt != Ndat) {
      std::cerr << "Not all data found in "<<dt.get_name()<<" was used, please check recN assignment in the input file!\n";
      exit(EXIT_FAILURE);
    }
        
    logL -= 0.5 * datapoint_cnt * log(2 * PI ); 
    
  }else if (model_type.substr(0,3) == "rec" && model_type.substr(5,3) == "lig") { // recX_ligY class of models
    unsigned int datapoint_cnt = 0;
    for (unsigned int j = 0;  j < Nrec;  ++j) { // loop over receptor types
      for (unsigned int k = 0; k < dt.length();  ++k) { // loop over datapoints
	if (dt.get_recN(k) == j+1) { // only take datapoints for the current receptor type
	  double mu_k = dt.get_x(k) * LOG10; // this is log(n)
	  
	  double mu[Nlig], pk[Nlig];
	  ////
	  for  (unsigned int l1 = 0; l1 < Nlig;  ++l1) {
	    double alpha_tot = 0.0;
	    for  (unsigned int l2 = 0; l2 < Nlig-1;  ++l2) {
	      //alpha_tot += prms_loc[Nprms-Nlig+l2];
	      alpha_tot += prms_loc[3*Nrec*Nlig+l2];
	    }
	    double log_alpha_tot = log(1 + alpha_tot);
	    if (l1 == 0) {
	      mu[l1] = mu_k - log_alpha_tot;
	    }else{
	      //mu[l1] = mu_k - log_alpha_tot + log(prms_loc[Nprms-Nlig+l1-1]);
	      mu[l1] = mu_k - log_alpha_tot + log(prms_loc[3*Nlig*Nrec+l1-1]);
	    }
	  }
	  double Z = 1.0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    Z += exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT);
	  }
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    pk[l] = exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT)/Z;
	  }
	  Ik[datapoint_cnt] = 0.0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    Ik[datapoint_cnt] += prms_loc[3*Nlig*j+ 3*l + 1] * pk[l] + (1.0/Nlig) * prms_loc[3*Nlig*j + 3*l + 2];  //this is Equation 2
	  }
          logL -= (dt.get_y(k) - Ik[datapoint_cnt]) * (dt.get_y(k) - Ik[datapoint_cnt]) / (2 * prms_loc[3*Nrec*Nlig + (Nlig-1) + 1 -1] * prms_loc[3*Nrec*Nlig + (Nlig-1) + 1 -1]); //first term of Equation SI_1
	  datapoint_cnt++;
	}
      }
    }
    if (datapoint_cnt == 0) {
      std::cerr << "No data found in Mobject::logLhood, please check recN assignment in "<<dt.get_name()<<"!\n";
      exit(EXIT_FAILURE);
    }
    if (datapoint_cnt != Ndat) {
      std::cerr << "Not all data found in "<<dt.get_name()<<" was used, please check recN assignment in the input file!\n";
      exit(EXIT_FAILURE);
    }
    logL -= 0.5 * datapoint_cnt * log(2 * PI * prms_loc[3*Nrec*Nlig + (Nlig-1) + 1 -1] * prms_loc[3*Nrec*Nlig + (Nlig-1) + 1 -1]); //Equation SI_1
  
    
  }else if (model_type.substr(0,8) == "MULT_SIG" && model_type.substr(9,3) == "rec" && model_type.substr(14,3) == "lig"){// MULT_SIG_recX_ligY class of models
   //int n = 1;
   unsigned int datapoint_cnt = 0;
   //double log_sigma_k = 0;
    for (unsigned int j = 0;  j < Nrec;  ++j) { // loop over receptor types
      for (unsigned int k = 0; k < dt.length();  ++k) { // loop over datapoints
	if (dt.get_recN(k) == j+1) { // only take datapoints for the current receptor type
	
	  double mu_k = dt.get_x(k) * LOG10; // this is log(n)
	  
	  double mu[Nlig], pk[Nlig];
	  ////
	  for  (unsigned int l1 = 0; l1 < Nlig;  ++l1) {
	    double alpha_tot = 0.0;
	    for  (unsigned int l2 = 0; l2 < Nlig-1;  ++l2) {
	      alpha_tot += prms_loc[3*Nlig*Nrec+l2];
	    }
	    double log_alpha_tot = log(1 + alpha_tot);
	    if (l1 == 0) {
	      mu[l1] = mu_k - log_alpha_tot;
	    }else{
	      mu[l1] = mu_k - log_alpha_tot + log(prms_loc[3*Nlig*Nrec+l1-1]);
	    }
	  }
	  double Z = 1.0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    Z += exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT);
	  }
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    pk[l] = exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT)/Z;
	  }
	  Ik[datapoint_cnt] = 0.0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    Ik[datapoint_cnt] += prms_loc[3*Nlig*j+ 3*l + 1] * pk[l] + (1.0/Nlig) * prms_loc[3*Nlig*j + 3*l + 2];
	  }
	  logL -= ( (dt.get_y(k) - Ik[datapoint_cnt]) * (dt.get_y(k) - Ik[datapoint_cnt]) / (2 * prms_loc[3*Nlig*Nrec + (Nlig-1) +j] * prms_loc[3*Nlig*Nrec + (Nlig-1) +j]) ) + log(prms_loc[3*Nlig*Nrec + (Nlig-1) +j]); //first term of Equation SI_1
	  datapoint_cnt++;
	  
	}
      }
    }
    if (datapoint_cnt == 0) {
      std::cerr << "No data found in Mobject::logLhood, please check recN assignment in "<<dt.get_name()<<"!\n";
      exit(EXIT_FAILURE);
    }
    if (datapoint_cnt != Ndat) {
      std::cerr << "Not all data found in "<<dt.get_name()<<" was used, please check recN assignment in the input file!\n";
      exit(EXIT_FAILURE);
    }
    
    logL -= 0.5 * datapoint_cnt * log(2 * PI);
	 

  }else if (model_type.substr(0,5) == "MuTOT" && model_type.substr(6,3) == "rec" && model_type.substr(11,3) == "lig"){// MuTOT_recX_ligY class of models
   //int n = 1;
   unsigned int datapoint_cnt = 0;
   //double log_sigma_k = 0;
   //int mod = 0;
    for (unsigned int j = 0;  j < Nrec;  ++j) { // loop over receptor types
      for (unsigned int k = 0; k < dt.length();  ++k) { // loop over datapoints
        if (dxt.length() != dt.length()){
          printf("dt.length = %d  dxt.length = %d\n\n",dt.length(), dxt.length());
          std::cerr << "Length of gapfile '"<<dxt.get_name()<<"' does not equal length of datafile '"<<dt.get_name()<<"'!! PLEASE CHECK INPUT FILES!\n\nExiting..\n"; 
          exit(EXIT_FAILURE);
        }
	if (dt.get_recN(k) == j+1) { // only take datapoints for the current receptor type
	  double mu[Nlig], pk[Nlig];
	  ////
	  for  (unsigned int l1 = 0; l1 < Nlig;  ++l1) {
	    double alpha_tot = 0.0;
	    for  (unsigned int l2 = 0; l2 < Nlig-1;  ++l2) {
	      alpha_tot += prms_loc[3*Nlig*Nrec+l2];
	    }
	    double log_alpha_tot = log(1 + alpha_tot);
	      if (l1 == 0) {
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - dxt.get_dx(k)*LOG10 - log_alpha_tot;
	      }else{
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - dxt.get_dx(k)*LOG10- log_alpha_tot + log(prms_loc[3*Nlig*Nrec+l1-1]);
	     }
	  }
	    
	   
	    
	  double Z = 1.0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    Z += exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT);
	  }
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    pk[l] = exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT)/Z;
	  }
	  Ik[datapoint_cnt] = 0.0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    Ik[datapoint_cnt] += prms_loc[3*Nlig*j+ 3*l + 1] * pk[l] + (1.0/Nlig) * prms_loc[3*Nlig*j + 3*l + 2];
	  }
	  
	  logL -= ( (dt.get_y(k) - Ik[datapoint_cnt]) * (dt.get_y(k) - Ik[datapoint_cnt]) / (2 * prms_loc[3*Nlig*Nrec + (Nlig-1) +j] * prms_loc[3*Nlig*Nrec + (Nlig-1) +j]) ) + log(prms_loc[3*Nlig*Nrec + (Nlig-1) +j]);
	  datapoint_cnt++;
	  
	}
      }
    }
    if (datapoint_cnt == 0) {
      std::cerr << "No data found in Mobject::logLhood, please check recN assignment in "<<dt.get_name()<<"!\n";
      exit(EXIT_FAILURE);
    }
    if (datapoint_cnt != Ndat) {
      std::cerr << "Not all data found in "<<dt.get_name()<<" was used, please check recN assignment in the input file!\n";
      exit(EXIT_FAILURE);
    }
    
   logL -= (0.5 * datapoint_cnt * log(2 * PI) );



  }else if (model_type.substr(0,9) == "FULL_HESS" && model_type.substr(10,3) == "rec" && model_type.substr(16,3) == "lig") { // FULL_HESS_recXX_ligYY class of models
    for (unsigned int x = 0; x < Nlig; x++){
      for (unsigned int y= 0; y < Nlig; y++){
        gsl_matrix_set (FullHessian, x, y, 0);
      }
    }
    double alpha_tot = 0.0;
    unsigned int datapoint_cnt = 0;
    for (unsigned int j = 0;  j < Nrec;  ++j) { // loop over receptor types     
    //REAL DATA:    
    for (unsigned int r = 0; r < 4; r++){
      for (unsigned int k = 0; k < 9; k++){
      
      //for (unsigned int r = 0; r < 7; r++){  //7 replicates   
        //for (unsigned int k = 0; k < 13; k++){  //12 data points each
     
	  double mu[Nlig], pk[Nlig];
	  ////
	  for  (unsigned int l1 = 0; l1 < Nlig;  ++l1) {
	    alpha_tot = 0.0;
	    //for  (unsigned int l2 = 0; l2 < Nlig-1;  ++l2) {
            for  (unsigned int l2 = 3*Nrec*Nlig; l2 < (3*Nrec*Nlig)+(Nlig-1);  ++l2) {
	      alpha_tot += prms_loc[l2];
	    }
	    double log_alpha_tot = log(1 + alpha_tot);
	    //for Receptor 2211:
	  ///*  
	        if (r == 3){
	    	if (l1 == 0) {
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - (2+k*(.5))*LOG10 - log_alpha_tot;
	      }else{
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - (2+k*(.5))*LOG10 - log_alpha_tot + log(prms_loc[3*Nlig*Nrec+l1-1]);
	     }
	    }else{ 
            //*/
	    
	    
	    //
	      if (k != 8){
	      if (l1 == 0) {
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - k*(.5)*LOG10 - log_alpha_tot;
	      }else{
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - k*(.5)*LOG10 - log_alpha_tot + log(prms_loc[3*Nlig*Nrec+l1-1]);
	     }
	    
	    //
	      }
	    ///*
	    else if (k == 8){
	      if (l1 == 0) {
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - 12*(.5)*LOG10 - log_alpha_tot;
	      }else{
	        mu[l1] = prms_loc[3*Nlig*Nrec + (Nlig-1) + Nrec]*LOG10 - 12*(.5)*LOG10 - log_alpha_tot + log(prms_loc[3*Nlig*Nrec+l1-1]);
	      }
	    }
	    //*/
	    //
	   }
	    
	  }
	  ////
	  double Z = 1.0;
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    Z += exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT);
	  }
	  ////
	  for  (unsigned int l = 0; l < Nlig;  ++l) {
	    pk[l] = exp(mu[l]-prms_loc[3*Nlig*j + 3*l]/RT)/Z;
	  }
	  
	    
	  double sum_Al_pl =0.0;
	  for (unsigned int l= 0; l < Nlig; ++l) {
	     sum_Al_pl += prms_loc[3*Nlig*j + 3*l +1] * pk[l];
	  }
	  
	  for (unsigned int x = 0; x < 1; x++){
            for (unsigned int y= 0; y < 1; y++){
              gsl_matrix_set (FullHessian, x, y, gsl_matrix_get(FullHessian, x, y) + (1/RT)*(1/RT)*( sum_Al_pl / Z )*( sum_Al_pl / Z ) );
            }
          }
          
          for (unsigned int x = 1; x < Nlig; x++){
            for (unsigned int y= 0; y < 1; y++){
              gsl_matrix_set (FullHessian, x, y, gsl_matrix_get(FullHessian, x, y) + (1/RT)*( sum_Al_pl / Z )*( ((pk[x]/prms_loc[3*Nrec*Nlig+x-1]) * (prms_loc[3*Nlig*j + 3*(x) +1]-sum_Al_pl)) - (sum_Al_pl)/(Z * (1+alpha_tot))) );           
            }
          }
          
          for (unsigned int x = 0; x < 1; x++){
            for (unsigned int y= 1; y < Nlig; y++){
              gsl_matrix_set (FullHessian, x, y, gsl_matrix_get(FullHessian, x, y) + (1/RT)*( sum_Al_pl / Z )*( ((pk[y]/prms_loc[3*Nrec*Nlig+y-1]) * (prms_loc[3*Nlig*j + 3*(y) +1]-sum_Al_pl)) - (sum_Al_pl)/(Z * (1+alpha_tot))) );           
            }
          }
          
          
	  
          for (unsigned int x = 1; x < Nlig; x++){
            for (unsigned int y= 1; y < Nlig; y++){
              gsl_matrix_set (FullHessian, x, y, gsl_matrix_get(FullHessian,x,y)+( ((pk[x]/prms_loc[3*Nrec*Nlig+x-1]) * (prms_loc[3*Nlig*j + 3*(x) +1]-sum_Al_pl)) - (sum_Al_pl)/(Z * (1+alpha_tot))) * ( ((pk[y]/prms_loc[3*Nrec*Nlig+y-1]) * (prms_loc[3*Nlig*j + 3*(y) +1]-sum_Al_pl)) - (sum_Al_pl)/(Z * (1+alpha_tot)))) ;
            }
          }
          
	  datapoint_cnt++;
	  
	  
	}
	
      }
      
    }
 
    if (datapoint_cnt == 0) {
      std::cerr << "No data found in Mobject::logLhood, please check recN assignment in "<<dt.get_name()<<"!\n";
      exit(EXIT_FAILURE);
    }
 
   
   logL = Determinant(FullHessian)*(pow(10,16));//scaling factor of 10^16
    //logL = Determinant(FullHessian);
 
  }
}



// Returns a random number sampled uniformly in the (low,hi) interval, GIVEN a random number in the (0,1) interval.
// So essentially this function rescales a (0,1) random number into an appropriate (low,hi) range.
double Mobject::uniform_prior(double low, double hi, double uf) {
  if (low >= hi) {
    std::cerr << "Unable to generate a uniform random number in the ("<<low<<","<<hi<<") range!\nExiting ..\n";
    exit(EXIT_FAILURE);
  }
  if (uf <= 0.0 || uf >= 1.0) {
    std::cerr << "Uniform random number "<<uf<<"is not in the (0,1) range!\nExiting ..\n";
    exit(EXIT_FAILURE);
  }
  //double uf = UNIFORM;
  double scale = hi - low;
  double shift = low;
  double ran = scale * uf + shift;
  return ran;
}

// Returns a random number sampled as Prob(x|I) ~ 1/x in the (low,hi) interval, GIVEN a uniform random number in the (0,1) interval.
// Essentially, this function returns a  random number x such that Prob(log(x)|I) = const.
double Mobject::jeffreys_prior(double low, double hi, double uf) {
  double ran_log = uniform_prior(log(low), log(hi), uf); // this will check ranges
  double ran = exp(ran_log);
  return ran;
}

// Returns a Gaussian-distributed (with mu & sigma) random number in the (-\infty,+\infty) interval given
// a Gaussian-distributed random number with mu=0.0 & sigma. So, basically just shifts x (or log(x)) by mu.
double Mobject::normal_prior(double mu, double gf) {
  return mu + gf;
}

// Evolves Mobject by MC with a logL constraint (logL > logLstar).
void Mobject::Evolve(DATAtable &dt, DELTAtable &dxt, double logLstar, int mcmc) {
  double step = 0.1;   // Initial guess suitable step-size, given as a fraction of the prior's range
  int accept = 0;   // # MCMC acceptances
  int reject = 0;   // # MCMC rejections
  
  double *trial_prms = new double[Nprms];
  double *trial_prms_ran = new double[Nprms];
  double *Ik_old = new double[Ndat];
  
  for(int m=0; m < mcmc; ++m) {
    for (unsigned int i=0; i<Nprms; ++i) {
      if (prms_ran[i] != FALSE_RAN) {
	if (prior_type[i] == 'u' || prior_type[i] == 'j') {
	  //trial_prms_ran[i] = prms_ran[i] + step * uniform_prior(-1.0, 1.0, UNIFORM);  // |move| < step
	  trial_prms_ran[i] = prms_ran[i] + step * uniform_prior(-1.0, 1.0, gsl_rng_uniform_pos(rng) );  // |move| < step
	  trial_prms_ran[i] -= floor(trial_prms_ran[i]); // wraparound to stay within (0,1)
	}else if (prior_type[i] == 'g') {
	  do {
	    trial_prms_ran[i] = prms_ran[i] + step * gsl_ran_gaussian(rng, prms_hi[i]);
	  } while (trial_prms_ran[i] == FALSE_RAN);
	}else if (prior_type[i] == 'l') {
	  do {
	    trial_prms_ran[i] = prms_ran[i] + step * gsl_ran_lognormal(rng, 0.0, prms_hi[i]);
	  } while (trial_prms_ran[i] == FALSE_RAN);
	}else{
	  std::cerr << "Unknown prior type: "<<prior_type[i]<<" in Mobject::Evolve!\nExiting ..\n";
	  exit(EXIT_FAILURE);
	}
      }else{
	trial_prms_ran[i] = prms_ran[i];
      }
    }
    
    for (unsigned int j=0; j<Nprms; ++j) {
      if (prms_ran[j] != FALSE_RAN) {
	if (prior_type[j] == 'u') {
	  //printf("%d\n%d\t%d\n",j,prms_lo[j],prms_hi[j]);
	  trial_prms[j] = uniform_prior(prms_lo[j], prms_hi[j], trial_prms_ran[j]);
	}else if (prior_type[j] == 'j') {
	  trial_prms[j] = jeffreys_prior(prms_lo[j], prms_hi[j], trial_prms_ran[j]);
	}else if (prior_type[j] == 'g') {
	  trial_prms[j] = normal_prior(prms_lo[j], trial_prms_ran[j]); // prms_lo contains Gaussian mean (mu) for this prior_type
	}else if (prior_type[j] == 'l') {
	  trial_prms[j] = normal_prior(prms_lo[j], trial_prms_ran[j]); // prms_lo contains Gaussian mean (mu) for this prior_type
	}else{
	  std::cerr << "Unknown prior type: "<<prior_type[j]<<" in Mobject::Evolve!\nExiting ..\n";
	  exit(EXIT_FAILURE);
	}
      }else{
	trial_prms[j] = prms[j];
      }
    }
    
    double logL_old = logL;
    for (unsigned int i=0; i<Ndat; ++i) {
      Ik_old[i] = Ik[i];
    }
    if (model_type.substr(0,6) == "LINEAR" ||model_type.substr(0,6) == "NONLIN" || model_type.substr(0,3) == "rec" && model_type.substr(5,3) == "lig" || model_type.substr(0,8) == "MULT_SIG" || model_type.substr(0,5) == "MuTOT" && model_type.substr(6,3) == "rec" && model_type.substr(11,3) == "lig" ){
      
      logLhood(dt, dxt, trial_prms); // this resets logL & Ik
    
      if (logL > logLstar) { // Accept if and only if within hard likelihood constraint
        for (unsigned int i=0; i<Nprms; ++i) { // accept prms
	  if (prms_ran[i] != FALSE_RAN) {
	    prms_ran[i] = trial_prms_ran[i];
	    prms[i] = trial_prms[i];
	  }
        }
        accept++;
      }else{
           
        logL = logL_old;
        for (unsigned int i=0; i<Ndat; ++i) {
	  Ik[i] = Ik_old[i];
        }
        reject++;
      }
      // Refine step-size to let acceptance ratio converge around 50%:
      if (accept > reject)  step *= exp(1.0 / accept);
      if (accept < reject)  step /= exp(1.0 / reject);
    
      
    }else if(model_type.substr(0,9) == "FULL_HESS" && model_type.substr(10,3) == "rec" && model_type.substr(16,3) == "lig") { // FullHessian-recX_ligY class of models

      gsl_matrix *tmpFullHessian = gsl_matrix_alloc(Nlig, Nlig);
      
      gsl_matrix_memcpy(tmpFullHessian, FullHessian);

      logLhood(dt, dxt, trial_prms); // this resets logL & Ik
    
      if (logL > logLstar) { // Accept if and only if within hard likelihood constraint
        for (unsigned int i=0; i<Nprms; ++i) { // accept prms
	  if (prms_ran[i] != FALSE_RAN) {
	    prms_ran[i] = trial_prms_ran[i];
	    prms[i] = trial_prms[i];
	  }
        }
        accept++;
      }else{
      
        gsl_matrix_memcpy(FullHessian, tmpFullHessian);        
        logL = logL_old;
        for (unsigned int i=0; i<Ndat; ++i) {
	  Ik[i] = Ik_old[i];
        }
        reject++;
      }
      // Refine step-size to let acceptance ratio converge around 50%:
      if (accept > reject)  step *= exp(1.0 / accept);
      if (accept < reject)  step /= exp(1.0 / reject);
      
      gsl_matrix_free(tmpFullHessian);  
    } 
  }
  
  // Report accept/reject ratio:
  //printf("# accepted = %d, # rejected = %d\n",accept,reject);
  
  //if (accept == 0) {
  //  std::cerr << "Could not find a single acceptable object - please increase number of MCMC trials or decrease step_size!\nExiting ..\n";
  //  exit(EXIT_FAILURE);
  //}
  
  delete [] trial_prms;
  delete [] trial_prms_ran;
  delete [] Ik_old;
}

void Mobject::writeIk(std::string const &filename,DATAtable &dt, DELTAtable &dxt, double *prms_loc) {
  // Recompute logL & Ik:
  logLhood(dt, dxt, prms_loc); // resets logL & Ik
  // Output results into a file:
  FILE *fp = fopen(filename.c_str(),"w");
  if (fp == NULL) {
    fprintf(stderr,"Error writing to file %s\n",filename.c_str() );
    exit(1);
  }
  fprintf(fp,"mu        I         Ik       recN\n"); // header
  for (unsigned int i=0; i<Ndat; ++i) {
    fprintf(fp,"%6.3f   %6.3f   %6.3f   %2d\n",dt.get_x(i),dt.get_y(i),Ik[i],dt.get_recN(i) );
  }
  fclose(fp);
}

