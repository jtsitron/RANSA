/* input.c */

// AM & JT, 10/01/2009

#include "input.h"

DATAtable::DATAtable(void) {
  name = "none";
  x = NULL;
  y = NULL;
  rec_num = NULL;
}

DATAtable::~DATAtable() {
  delete [] x;
  delete [] y;
  delete [] rec_num;
}

DATAtable::DATAtable(std::string const & filename) {
  name = filename;
  if (name == "none") {
    //std::cerr << "Attempted to read a file with reserved name="<<name<<" in DATAtable::DATAtable!\nExiting ..\n";
    //exit(EXIT_FAILURE);
    x = NULL;
    y = NULL;
    rec_num = NULL;
    N = 0;
  }else{
  ReadDATAtable(filename);
  }
}


void DATAtable::ReadDATAtable(std::string const & filename) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::ifstream;
  
  ifstream file(filename.c_str() );
  
  char buffer[1024];
  
  if (file.fail() ) {
    cerr << "*****Could not open file: "<<filename<<" for reading!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  
  // Count data points in file:
  int cnt = 0; // data point counter
  double tmp[2];
  unsigned int tmp_int;
  while (!file.eof() ) {
    file.getline(buffer,1024,'\n');
    if (file.fail() && !file.eof() ) file.clear();
    if (sscanf(buffer,"%lf %lf %d",
	       &tmp[0],&tmp[1],&tmp_int) == 3) {
      ++cnt;
    }
  }
  file.clear(); // forget we hit the end of file
  file.seekg(0, std::ios::beg); // move to the start of the file
  
  // Allocate arrays:
  N = cnt;
  x = new double[N];
  y = new double[N];
  rec_num = new unsigned int[N];
  
  cnt = -1; // reset basestep counter
  while (!file.eof() ) {
    file.getline(buffer,512,'\n');
    if (file.fail() && !file.eof() ) file.clear();
    if (cnt==-1) {
      ++cnt;
      continue; // skip header
    }
    if (sscanf(buffer,"%lf %lf %d",&x[cnt],&y[cnt],&rec_num[cnt]) == 3) {
      ++cnt;
    }
  }
  file.close();
  
  cout << "Read "<<cnt<<" data points from file: "<<filename<<" ..\n";

}

const unsigned int DATAtable::length(void) {
  return N;
}

// No range check for speed!
const double DATAtable::get_x(int i) {
  return x[i];
}

// No range check for speed!
const double DATAtable::get_y(int i) {
  return y[i];
}

// No range check for speed!
const unsigned int DATAtable::get_recN(int i) {
  return rec_num[i];
}

const std::string DATAtable::get_name(void) {
  return name;
}

void DATAtable::WriteDATAtable(std::string const & filename) {
  FILE *fp = fopen(filename.c_str(),"w");
  if (fp == NULL) {
    fprintf(stderr,"Error writing to file %s\n",filename.c_str() );
    exit(1);
  }
  fprintf(fp,"mu        I        recN\n"); // header
  for (unsigned int i=0; i<N; ++i) {
    fprintf(fp,"%6.3f   %6.3f   %2d\n",x[i],y[i],rec_num[i]);
  }
  fclose(fp);
}

/* ==================== */
/* =======DELTAtable=====*/
/* ====================*/

DELTAtable::DELTAtable(void) {
  name = "none";
  dx = NULL;
}

DELTAtable::~DELTAtable() {
  delete [] dx;
  //delete [] y;
  //delete [] rec_num;
}

DELTAtable::DELTAtable(std::string const & filename) {
  name = filename;
  if (name == "none") {
    //std::cerr << "Attempted to read a file with reserved name="<<name<<" in DELTAtable::DELTAtable!\nExiting ..\n";
    //exit(EXIT_FAILURE);
    dx = NULL;
  }else{
    ReadDELTAtable(filename);
  }
  
}

void DELTAtable::ReadDELTAtable(std::string const & filename) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::ifstream;
  
  ifstream file(filename.c_str() );
  
  char buffer[1024];
  
  if (file.fail() ) {
    cerr << "!!Could not open file: "<<filename<<" for reading!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  
  // Count data points in file:
  int cnt = 0; // data point counter
  double tmp;
  unsigned int tmp_int;
  while (!file.eof() ) {
    file.getline(buffer,1024,'\n');
    if (file.fail() && !file.eof() ) file.clear();
    if (sscanf(buffer,"%lf ",
	       &tmp) == 1) {
      ++cnt;
    }
  }
  file.clear(); // forget we hit the end of file
  file.seekg(0, std::ios::beg); // move to the start of the file
  
  // Allocate arrays:
  N = cnt;
  dx = new double[N];

  cnt = -1; // reset basestep counter
  while (!file.eof() ) {
    file.getline(buffer,512,'\n');
    if (file.fail() && !file.eof() ) file.clear();
    if (cnt==-1) {
      ++cnt;
      continue; // skip header
    }
    if (sscanf(buffer,"%lf ",&dx[cnt]) == 1) {
      ++cnt;
    }
  }
  file.close();
  
  cout << "Read "<<cnt<<" 'delta' values from file: "<<filename<<" ..\n";
  
}

const unsigned int DELTAtable::length(void) {
  return N;
}

// No range check for speed!
const double DELTAtable::get_dx(int i) {
  return dx[i];
}


const std::string DELTAtable::get_name(void) {
  return name;
}

/*void DATAtable::WriteDATAtable(std::string const & filename) {
  FILE *fp = fopen(filename.c_str(),"w");
  if (fp == NULL) {
    fprintf(stderr,"Error writing to file %s\n",filename.c_str() );
    exit(1);
  }
  fprintf(fp,"mu        I        recN\n"); // header
  for (unsigned int i=0; i<N; ++i) {
    fprintf(fp,"%6.3f   %6.3f   %2d\n",x[i],y[i],rec_num[i]);
  }
  fclose(fp);
}*/


/* ===================== */
/* ========Mtable========= */
/* ===================== */

Mtable::Mtable(void) {
  name = "none";
  ddG = NULL;
  A = NULL;
  B = NULL;
  alpha = NULL;
  Sigma = NULL; //line added 3/8/10
  //mu_tot = NULL; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ddG_onoff = NULL;
  A_onoff = NULL;
  B_onoff = NULL;
  alpha_onoff = NULL;
  Sigma_onoff = NULL; //line added 3/8/10
  //mutot_onoff = NULL; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  N = 0;
  Nmax = 0;
}

Mtable::~Mtable() {
  delete [] ddG;
  delete [] A;
  delete [] B;
  delete [] alpha;
  delete [] Sigma; //line added 3/8/10
  //delete [] mu_tot; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  delete [] ddG_onoff;
  delete [] A_onoff;
  delete [] B_onoff;
  delete [] alpha_onoff;
  delete [] Sigma_onoff; //line added 3/8/10
  //delete [] mutot_onoff; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

Mtable::Mtable(std::string const & filename) {
  name = filename;
  if (name == "none") {
    //std::cerr << "Attempted to read a file with reserved name="<<name<<" in Mtable::Mtable!\nExiting ..\n";
    //exit(EXIT_FAILURE);
    ddG = NULL;
    A = NULL;
    B = NULL;
    alpha = NULL;
    Sigma = NULL; //line added 3/8/10
    ddG_onoff = NULL;
    A_onoff = NULL;
    B_onoff = NULL;
    alpha_onoff = NULL;
    Sigma_onoff = NULL; //line added 3/8/10
    N = 0;
    Nmax = 0;
  }else{
    ReadMtable(filename);
  }
}

void Mtable::ReadMtable(std::string const & filename) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::ifstream;
  
  ifstream file(filename.c_str() );
  
  char buffer[1024];
  char tmp_string[1024];
  unsigned int tmp_ind,tmp_onoff;
  double tmp_val;
  
  if (file.fail() ) {
    cerr << "Could not open file: "<<filename<<" for reading!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  
  // Allocate arrays (N_RECLIG is static):
  ddG = new double[N_RECLIG];
  A = new double[N_RECLIG];
  B = new double[N_RECLIG];
  alpha = new double[N_LIG-1+1];
  Sigma = new double[N_REC]; 
  int *rec_types = new int[N_RECLIG];
  for (unsigned int i=0; i<N_RECLIG; ++i) { // Init. prms
    ddG[i] = FALSE_VALUE;
    A[i] = FALSE_VALUE;
    B[i] = FALSE_VALUE;
    rec_types[i] = (int) FALSE_VALUE;
  }
  for (unsigned int i=0; i<N_LIG-1; ++i) {
    alpha[i] = FALSE_VALUE;
  }
  for (unsigned int i=0; i<N_REC; ++i){
    Sigma[i] = FALSE_VALUE;
  }						
  mu_tot = FALSE_VALUE;
  ////
  ddG_onoff = new bool[N_RECLIG];
  A_onoff = new bool[N_RECLIG];
  B_onoff = new bool[N_RECLIG];
  alpha_onoff = new bool[N_LIG-1];
  Sigma_onoff = new bool[N_REC];

  for (unsigned int i=0; i<N_RECLIG; ++i) { // Init. prms
    ddG_onoff[i] = false;
    A_onoff[i] = false;
    B_onoff[i] = false;
  }
  for (unsigned int i=0; i<N_LIG-1; ++i) {
    alpha_onoff[i] = false;
  }
  for (unsigned int i=0; i<N_REC; ++i) {
    Sigma_onoff[i] = false;
  }						
  mutot_onoff = false;
  
  int cnt = -1; // model prm counter
  while (!file.eof() ) {
    file.getline(buffer,1024,'\n');
    if (file.fail() && !file.eof() ) file.clear();
    if (cnt==-1) {
      ++cnt;
      continue; // skip header
    }
    if (sscanf(buffer,"%s %lf %d %d",tmp_string,&tmp_val,&tmp_onoff,&tmp_ind) == 4) { // read off dG,A,B,alpha_i,Sigma_i
      if (tmp_ind <= 0) {
	cerr << "Receptor index (given as "<<tmp_ind<<") must be a positive integer!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      if (tmp_onoff != 0 && tmp_onoff != 1) {
	cerr << "Parameter on/off state (given as "<<tmp_onoff<<") must be a 0/1 integer!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      // Make sure there are no clashes with FALSE_VALUE:
      if (tmp_val == FALSE_VALUE) {
	cerr << "Parameter value (given as "<<tmp_val<<") cannot be equal to FALSE_VALUE!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      if (strcmp(tmp_string,"ddG") == 0) {
	ddG[tmp_ind-1] = tmp_val;
	ddG_onoff[tmp_ind-1] = (bool) tmp_onoff;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"A") == 0) {
	A[tmp_ind-1] = tmp_val;
	A_onoff[tmp_ind-1] = (bool) tmp_onoff;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"B") == 0) {
	B[tmp_ind-1] = tmp_val;
	B_onoff[tmp_ind-1] = (bool) tmp_onoff;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"alpha") == 0) {
	alpha[tmp_ind-1] = tmp_val;
	alpha_onoff[tmp_ind-1] = (bool) tmp_onoff;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"sigma") == 0) { 
	Sigma[tmp_ind-1] = tmp_val;
	Sigma_onoff[tmp_ind-1] = (bool) tmp_onoff;
	rec_types[tmp_ind-1] = tmp_ind;             
      }
      else{
	cerr << "Unknown parameter type: "<<tmp_string<<"!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      ++cnt;
    }else if  (sscanf(buffer,"%s %lf %d",tmp_string,&tmp_val,&tmp_onoff) == 3) { // read off mutot
      if (tmp_onoff != 0 && tmp_onoff != 1) {
	cerr << "Parameter on/off state (given as "<<tmp_onoff<<") must be a 0/1 integer!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      // Make sure there are no clashes with FALSE_VALUE:
      if (tmp_val == FALSE_VALUE) {
	cerr << "Parameter value (given as "<<tmp_val<<") cannot be equal to FALSE_VALUE!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      if (strcmp(tmp_string,"mutot") == 0) { 
	mu_tot = tmp_val;
	mutot_onoff = (bool) tmp_onoff;
      }else if (strcmp(tmp_string,"ddG") == 0 || strcmp(tmp_string,"A") == 0 || strcmp(tmp_string,"B") == 0 || strcmp(tmp_string,"alpha") == 0 || strcmp(tmp_string,"sigma") == 0) {
        cerr << "The parameter "<<tmp_string<<" must have an index value!\nExiting..\n";
	exit(EXIT_FAILURE);
      }else{
	cerr << "!!!Unknown parameter type: "<<tmp_string<<"!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      ++cnt;
    }
  }
  file.close();
  
  // Compute # different receptor/ligand combinations:
  N = 0;
  Nmax = 0;
  for (unsigned int i=0; i<N_RECLIG; ++i) {
    if (rec_types[i] != FALSE_VALUE) {
      N++;
      if (rec_types[i] > Nmax) Nmax = rec_types[i];
    }
  }
  
  cout << "Read "<<cnt<<" model parameter(s) from file: "<<filename<<" ..\n";
  cout << "Found one or more model parameter for "<<N<<" receptor/ligand combination(s) ..\n";
  
  delete [] rec_types;
}

void Mtable::WriteMtable(std::string const & filename) {
  FILE *fp = fopen(filename.c_str(),"w");
  if (fp == NULL) {
    fprintf(stderr,"Error writing to file %s\n",filename.c_str() );
    exit(1);
  }
  fprintf(fp,"Prm      Val    OnOff  Ind\n"); // header
  if (name != "none") {
    for (unsigned int i=0; i<N_RECLIG; ++i) {
      if (ddG[i] != FALSE_VALUE) {
	fprintf(fp,"%5s   %6.3f   %2d   %2d\n","ddG",ddG[i],ddG_onoff[i],i+1);
      }
      if (A[i] != FALSE_VALUE) {
	fprintf(fp,"%5s   %6.3f   %2d   %2d\n","A",A[i],A_onoff[i],i+1);
      }
      if (B[i] != FALSE_VALUE) {
	fprintf(fp,"%5s   %6.3f   %2d   %2d\n","B",B[i],B_onoff[i],i+1);
      }
    }
    for (unsigned int i=0; i<N_LIG-1; ++i) {
      if (alpha[i] != FALSE_VALUE) {
	fprintf(fp,"%5s   %6.3f   %2d   %2d\n","alpha",alpha[i],alpha_onoff[i],i+1);
      }
    }
    for (unsigned int i=0; i<N_REC; ++i) {                                           
      if (Sigma[i] != FALSE_VALUE) {
        fprintf(fp,"%5s   %6.3f   %2d	%2d\n","sigma",Sigma[i],Sigma_onoff[i],i+1);  
      }
    }                                                                                 
  }
  fclose(fp);
}

// No range check for speed!
const double Mtable::get_ddG(int i) {
  return ddG[i];
}

// No range check for speed!
const double Mtable::get_A(int i) {
  return A[i];
}

// No range check for speed!
const double Mtable::get_B(int i) {
  return B[i];
}

// No range check for speed!
const double Mtable::get_alpha(int i) {
  return alpha[i];
}

/*const double Mtable::get_sigma(void) {
  return sigma;
}*/

const double Mtable::get_sigma(int i) {		
  return Sigma[i];
}                                               

const double Mtable::get_mutot(void) {	
  return mu_tot;
}                                               
// No range check for speed!
const bool Mtable::get_ddG_onoff(int i) {
  return ddG_onoff[i];
}

// No range check for speed!
const bool Mtable::get_A_onoff(int i) {
  return A_onoff[i];
}

// No range check for speed!
const bool Mtable::get_B_onoff(int i) {
  return B_onoff[i];
}

// No range check for speed!
const bool Mtable::get_alpha_onoff(int i) {
  return alpha_onoff[i];
}

/*const bool Mtable::get_sigma_onoff(void) {
  return sigma_onoff;
}*/

const bool Mtable::get_sigma_onoff(int i) {		
  return Sigma_onoff[i];
}                                               

const bool Mtable::get_mutot_onoff(void) {	
  return mutot_onoff;
}                                              

const std::string Mtable::get_name(void) {
  return name;
}

const unsigned int Mtable::getN(void) {
  return N;
}

const unsigned int Mtable::getNmax(void) { 
  return Nmax;
}

const unsigned int Mtable::getNreclig(void) {
  return N_RECLIG;
}

/* ===================== */
/* ========Ptable========= */
/* ===================== */

Ptable::Ptable(void) {
  name = "none";
  ddG = NULL;
  A = NULL;
  B = NULL;
  alpha = NULL;
  Sigma = NULL;  
  mu_tot = NULL;  
  ddG_prior_type = NULL;
  A_prior_type = NULL;
  B_prior_type = NULL;
  alpha_prior_type = NULL;
  Sigma_prior_type = NULL; 
  mutot_prior_type = '\0';  
  //mutot_prior_type = 'NULL';  
  N = 0;
  Ntot = 0;
}

Ptable::~Ptable() {
  ////
  if (ddG != NULL) {
    for (unsigned int i=0; i<N_RECLIG; ++i) {
      delete [] ddG[i];
    }
    delete [] ddG;
  }
  ////
  if (A != NULL) {
    for (unsigned int i=0; i<N_RECLIG; ++i) {
      delete [] A[i];
    }
    delete [] A;
  }
  ////
  if (B != NULL) {
    for (unsigned int i=0; i<N_RECLIG; ++i) {
      delete [] B[i];
    }
    delete [] B;
  }
  ////
  if (alpha != NULL) {
    for (unsigned int i=0; i<N_LIG-1; ++i) {
      delete [] alpha[i];
    }
    delete [] alpha;
  }
  ////
  if (Sigma != NULL) {                             
    for (unsigned int i=0; i<N_REC; ++i) {
      delete [] Sigma[i];
    }
    delete [] Sigma;
  }                                               
  
  delete [] mu_tot;  
  
  delete [] ddG_prior_type;
  delete [] A_prior_type;
  delete [] B_prior_type;
  delete [] alpha_prior_type;
  delete [] Sigma_prior_type;                 
}

Ptable::Ptable(std::string const & filename) {
  name = filename;
  if (name == "none") {
    //std::cerr << "Attempted to read a file with reserved name="<<name<<" in Ptable::Ptable!\nExiting ..\n";
    //exit(EXIT_FAILURE);
    ddG = NULL;
    A = NULL;
    B = NULL;
    alpha = NULL;
    Sigma = NULL;                              
    mu_tot = NULL; 
    ddG_prior_type = NULL;
    A_prior_type = NULL;
    B_prior_type = NULL;
    alpha_prior_type = NULL;           
    Sigma_prior_type = NULL;                    
    mutot_prior_type = '\0';  
    N = 0;
    Ntot = 0;
  }else{
    ReadPtable(filename);
  }
}

void Ptable::ReadPtable(std::string const & filename) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::ifstream;
  
  ifstream file(filename.c_str() );
  
  
  
  char buffer[1024];
  char tmp_string[1024];
  unsigned int tmp_ind;
  double tmp_val1,tmp_val2;
  char tmp_char;
  
  if (file.fail() ) {
    cerr << "Could not open file: "<<filename<<" for reading!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  
  // Allocate arrays (N_RECLIG is static):
  int *rec_types = new int[N_RECLIG];
  ////
  ddG = new double*[N_RECLIG];
  A = new double*[N_RECLIG];
  B = new double*[N_RECLIG];
  for (unsigned int i=0; i<N_RECLIG; ++i) {
    ddG[i] = new double[N_PRI];
    A[i] = new double[N_PRI];
    B[i] = new double[N_PRI];
  }
  ddG_prior_type = new char[N_RECLIG];
  A_prior_type = new char[N_RECLIG];
  B_prior_type = new char[N_RECLIG];
  ////
  alpha = new double*[N_LIG-1];
  for (unsigned int i=0; i<N_LIG-1; ++i) {
    alpha[i] = new double[N_PRI];
  }
  alpha_prior_type = new char[N_LIG-1];
  ////
  Sigma = new double*[N_REC];                
  for (unsigned int i=0; i<N_REC; i++){
    Sigma[i] = new double[N_PRI];            
  }
  Sigma_prior_type = new char[N_REC];        
  ////
  mu_tot = new double[N_PRI]; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ////
  for (unsigned int i=0; i<N_RECLIG; ++i) { // Init. prms
    rec_types[i] = (int) FALSE_VALUE;
    ddG_prior_type[i] = '-';
    A_prior_type[i] = '-';
    B_prior_type[i] = '-';
    for (unsigned int j=0; j<N_PRI; ++j) {
      ddG[i][j] = FALSE_VALUE;
      A[i][j] = FALSE_VALUE;
      B[i][j] = FALSE_VALUE;
    }
  }
  ////
  for (unsigned int i=0; i<N_LIG-1; ++i) {
    alpha_prior_type[i] = '-';
    for (unsigned int j=0; j<N_PRI; ++j) {
      alpha[i][j] = FALSE_VALUE;
    }
  }
  ////
  for (unsigned int i=0; i<N_REC; ++i) {          
    Sigma_prior_type[i] = '-';
    for (unsigned int j=0; j<N_PRI; ++j) {
      Sigma[i][j] = FALSE_VALUE;
    }
  }                                              

  ////
  for (unsigned int j=0; j<N_PRI; ++j) { 
    mu_tot[j] = FALSE_VALUE;
  }
  mutot_prior_type = '-';    
  //// 
  int cnt = -1; // model prm counter
  while (!file.eof() ) {
    file.getline(buffer,1024,'\n');
    if (file.fail() && !file.eof() ) file.clear();
    if (cnt==-1) {
      ++cnt;
      continue; // skip header
    }
    if (sscanf(buffer,"%s %lf %lf %c %d",tmp_string,&tmp_val1,&tmp_val2,&tmp_char,&tmp_ind) == 5) { // read off ddG,A,B,alpha,Sigma prior prms
      if (tmp_ind <= 0) {
	cerr << "Receptor index (given as "<<tmp_ind<<") must be a positive integer!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      if (tmp_char == '-') {
	cerr << "PriorType "<<tmp_char<<"is invalid!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      if (strcmp(tmp_string,"ddG") == 0) {
	ddG[tmp_ind-1][0] = tmp_val1;
	ddG[tmp_ind-1][1] = tmp_val2;
	ddG_prior_type[tmp_ind-1] = tmp_char;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"A") == 0) {
	A[tmp_ind-1][0] = tmp_val1;
	A[tmp_ind-1][1] = tmp_val2;
	A_prior_type[tmp_ind-1] = tmp_char;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"B") == 0) {
	B[tmp_ind-1][0] = tmp_val1;
	B[tmp_ind-1][1] = tmp_val2;
	B_prior_type[tmp_ind-1] = tmp_char;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"alpha") == 0) {
	alpha[tmp_ind-1][0] = tmp_val1;
	alpha[tmp_ind-1][1] = tmp_val2;
	alpha_prior_type[tmp_ind-1] = tmp_char;
	rec_types[tmp_ind-1] = tmp_ind;
      }else if (strcmp(tmp_string,"sigma") == 0) {  
	Sigma[tmp_ind-1][0] = tmp_val1;
	Sigma[tmp_ind-1][1] = tmp_val2;
	Sigma_prior_type[tmp_ind-1] = tmp_char;
	rec_types[tmp_ind-1] = tmp_ind;             
      }
      else{
	cerr << "Unknown parameter type: "<<tmp_string<<"!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      ++cnt;
    }else if  (sscanf(buffer,"%s %lf %lf %c",tmp_string,&tmp_val1,&tmp_val2,&tmp_char) == 4) { // read off mutot prior prms
      if (tmp_char == '-') {
	cerr << "PriorType "<<tmp_char<<"is invalid!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      if (strcmp(tmp_string,"mutot") == 0) {  
	mu_tot[0] = tmp_val1;
	mu_tot[1] = tmp_val2;
	mutot_prior_type = tmp_char; 
      }else{
	cerr << "Unknown parameter type: "<<tmp_string<<"!\nExiting..\n";
	exit(EXIT_FAILURE);
      }
      ++cnt;
    }
  }
  file.close();
  
  // Compute # different receptor/ligand combinations:
  Ntot = cnt;
  N = 0;
  for (unsigned int i=0; i<N_RECLIG; ++i) {
    if (rec_types[i] != FALSE_VALUE) {
      N++;
    }
  }
  
  cout << "Read "<<Ntot<<" model parameter(s) from file: "<<filename<<" ..\n";
  cout << "Found one or more model parameter for "<<N<<" receptor/ligand combination(s) ..\n";
  
  delete [] rec_types;
}

void Ptable::WritePtable(std::string const & filename) {
  FILE *fp = fopen(filename.c_str(),"w");
  if (fp == NULL) {
    fprintf(stderr,"Error writing to file %s\n",filename.c_str() );
    exit(1);
  }
  fprintf(fp,"Prm        Val1       Val2     PriorType    Ind\n"); // header
  if (name != "none") {
    for (unsigned int i=0; i<N_RECLIG; ++i) {
      if (ddG_prior_type[i] != '-') {
	fprintf(fp,"%5s   %9.3f   %9.3f   %1c      %2d\n","ddG",ddG[i][0],ddG[i][1],ddG_prior_type[i],i+1);
      }
      if (A_prior_type[i] != '-') {
	fprintf(fp,"%5s   %9.3f   %9.3f   %1c      %2d\n","A",A[i][0],A[i][1],A_prior_type[i],i+1);
      }
      if (B_prior_type[i] != '-') {
	fprintf(fp,"%5s   %9.3f   %9.3f   %1c      %2d\n","B",B[i][0],B[i][1],B_prior_type[i],i+1);
      }
    }
    for (unsigned int i=0; i<N_LIG-1; ++i) {
      if (alpha_prior_type[i] != '-') {
	fprintf(fp,"%5s   %6.3f   %6.3f   %1c      %2d\n","alpha",alpha[i][0],alpha[i][1],alpha_prior_type[i],i+1);
      }
    }
  }
  fclose(fp);
}

const std::string Ptable::get_name(void) {
  return name;
}

const unsigned int Ptable::getN(void) { // returns total # of receptor-ligand combinations read from the file
  return N;
}

const unsigned int Ptable::getNtot(void) { // returns total # of prior prm sets read from the file
  return Ntot;
}

const double Ptable::get_prior_prm(std::string const & prm_type, unsigned int rl_ind, unsigned int pr_ind) {
  // Range check:
  if (rl_ind >= N_RECLIG) {
    std::cerr << "receptor_ligand index "<<rl_ind<<"is out of range in Ptable::get_prior_prm!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  if (pr_ind >= N_PRI) {
    std::cerr << "prior prm index "<<pr_ind<<"is out of range in Ptable::get_prior_prm!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  if (prm_type == "ddG") {
    return ddG[rl_ind][pr_ind];
  }else if (prm_type == "A") {
    return A[rl_ind][pr_ind];
  }else if (prm_type == "B") {
    return B[rl_ind][pr_ind];
  }else if (prm_type == "alpha") {
    // Additional range check:
    if (rl_ind >= N_LIG-1) {
      std::cerr << "receptor_ligand index "<<rl_ind<<"is out of range in Ptable::get_prior_prm!\nExiting..\n";
      exit(EXIT_FAILURE);
    }
    return alpha[rl_ind][pr_ind];
  }else if (prm_type == "sigma") {                                                     
    // Additional range check:
    if (rl_ind >= N_REC) {
      std::cerr << "receptor_ligand index "<<rl_ind<<"is out of range in Ptable::get_prior_prm!\nExiting..\n";
      exit(EXIT_FAILURE);
    }
    return Sigma[rl_ind][pr_ind];                                                      
  }
  else{
    std::cerr << "Unknown parameter type: "<<prm_type<<" in Ptable::get_prior_prm!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
}


const double Ptable::get_mutot(unsigned int pr_ind) {  
  // Range check:
  if (pr_ind >= N_PRI) {
    std::cerr << "prior prm index "<<pr_ind<<"is out of range in Ptable::get_mutot!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  return mu_tot[pr_ind];
}                                                     
const char Ptable::get_prior_type(std::string const & prm_type, unsigned int rl_ind) {
  // Range check:
  if (rl_ind >= N_RECLIG) {
    std::cerr << "receptor_ligand index "<<rl_ind<<"is out of range in Ptable::get_prior_type!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
  if (prm_type == "ddG") {
    return ddG_prior_type[rl_ind];
  }else if (prm_type == "A") {
    return A_prior_type[rl_ind];
  }else if (prm_type == "B") {
    return B_prior_type[rl_ind];
  }else if (prm_type == "alpha") {
    // Additional range check:
    if (rl_ind >= N_LIG-1) {
      std::cerr << "receptor_ligand index "<<rl_ind<<"is out of range in Ptable::get_prior_type!\nExiting..\n";
      exit(EXIT_FAILURE);
    }
    return alpha_prior_type[rl_ind];
  }else if (prm_type == "sigma") {                                                      
    // Additional range check:
    if (rl_ind >= N_REC) {
      std::cerr << "receptor_ligand index "<<rl_ind<<"is out of range in Ptable::get_prior_type!\nExiting..\n";
      exit(EXIT_FAILURE);
    }
    return Sigma_prior_type[rl_ind];                                                   
  }
  else{
    std::cerr << "Unknown parameter type: "<<prm_type<<" in Ptable::get_prior_type!\nExiting..\n";
    exit(EXIT_FAILURE);
  }
}


const char Ptable::get_mutot_type(void) {  
  return mutot_prior_type;
}                                          
