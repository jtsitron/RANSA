/* input.h */

// AM & JT, 10/01/2009

#ifndef _INPUT_H_
#define _INPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>

# define FALSE_VALUE -1000000.0

// This header file declares objects used for reading data from input files.

static const unsigned int N_RECLIG = 9801; // max number of receptor/ligand combinations
static const unsigned int N_REC = 99; // max number of receptor types ***********(line added 3/8/10)*************
static const unsigned int N_LIG = 99; // max number of ligand types
static const unsigned int N_PRI = 2; // max number of prior prms in each receptor/ligand combination

// This class reads in a three-column table with input data:   mu    I    receptor_number.
class DATAtable {
 private:
  void ReadDATAtable(std::string const & filename);
  double *x; // values of mu
  double *y; // values of I
  unsigned int *rec_num; // positive integers indicating receptor (or experiment) type
  unsigned int N; // number of data points in the file
  std::string name; // name of the dataset (=filename)
 public:
  DATAtable();
  DATAtable(std::string const & filename);
  ~DATAtable();
  void WriteDATAtable(std::string const & filename);
  const unsigned int length(void);
  const double get_x(int i);
  const double get_y(int i);
  const unsigned int get_recN(int i);
  const std::string get_name(void);
  //gsl_matrix* Export(void); // export data into a gsl_matrix (allocated inside!)
};

// This class reads in a one-column table with 'delta' values.  It is used for the muTOT mode where we are fitting for the log[total ligand concentration]
class DELTAtable {
 private:
  void ReadDELTAtable(std::string const & filename);
  double *dx; // values of the diffs of concentrations between data points
  unsigned int N; // number of 'delta' values in the file
  std::string name; // name of the dataset (=filename)
 public:
  DELTAtable();
  DELTAtable(std::string const & filename);
  ~DELTAtable();
  //void WriteDELTAtable(std::string const & filename);
  const unsigned int length(void);
  const double get_dx(int i);
  const std::string get_name(void);
  //gsl_matrix* Export(void); // export data into a gsl_matrix (allocated inside!)
};


// This class reads in a file with model parameters.
class Mtable {
 private:
  void ReadMtable(std::string const & filename);
  double *ddG; // ddG
  double *A; // A
  double *B; // B
  double *alpha; // mixture coeffs
  double *Sigma; // Gaussian noise prm   *********(line added on 3/8/10)**********
  //double sigma; // Gaussian noise prm
  double mu_tot; // log[total concentration] (a.k.a. chem potential)  
  bool *ddG_onoff; // 0/1 flag for fitting ddG
  bool *A_onoff; // 0/1 flag for fitting A
  bool *B_onoff; // 0/1 flag for fitting B
  bool *alpha_onoff; // 0/1 flag for fitting mixture coeffs
  bool *Sigma_onoff; // 0/1 flag for Gaussian noise       ********(line added on 3/8/10)**********
  //bool sigma_onoff; // 0/1 flag for Gaussian noise
  bool mutot_onoff; // 0/1 flag for log[total ligand concentration]
  std::string name; // name of the file with prms
  unsigned int N; // actual number of receptor/ligand types read from the file
  unsigned int Nmax; // one-based max rec_lig index actually found in the file
 public:
  Mtable();
  Mtable(std::string const & filename);
  ~Mtable();
  void WriteMtable(std::string const & filename);
  //const unsigned int length(void);
  ////
  const double get_ddG(int i);
  const double get_A(int i);
  const double get_B(int i);
  const double get_alpha(int i);
  //const double get_sigma(void);
  const double get_sigma(int i);
  const double get_mutot(void); 
  ////
  const bool get_ddG_onoff(int i);
  const bool get_A_onoff(int i);
  const bool get_B_onoff(int i);
  const bool get_alpha_onoff(int i);
  //const bool get_sigma_onoff(void);
  const bool get_sigma_onoff(int i);
  const bool get_mutot_onoff(void); 
  ////
  const std::string get_name(void);
  const unsigned int getN(void);
  const unsigned int getNmax(void);
  const unsigned int getNreclig(void);
};

// This class reads in a file with prior parameters (ranges etc.)
class Ptable {
 private:
  void ReadPtable(std::string const & filename);
  double **ddG;        // ddG
  double **A;          // A
  double **B;          // B
  double **alpha;      // mixture coeffs
  double **Sigma;      // Gaussian noise prms for MULT_SIG case  
  double *mu_tot;      //MuTOT  
  char *ddG_prior_type;
  char *A_prior_type;
  char *B_prior_type;
  char *alpha_prior_type;
  char *Sigma_prior_type; 
  char mutot_prior_type; 
  std::string name;    // name of the file with prior prms
  unsigned int N;      // actual number of receptor/ligand prior types read from the file
  unsigned int Ntot;   // actual number of prior prm sets read from the file
 public:
  Ptable();
  Ptable(std::string const & filename);
  ~Ptable();
  void WritePtable(std::string const & filename);
  const std::string get_name(void);
  const unsigned int getN(void);
  const unsigned int getNtot(void);
  const double get_prior_prm(std::string const & prm_type, unsigned int rl_ind, unsigned int pr_ind);
  const char get_prior_type(std::string const & prm_type, unsigned int rl_ind);
  //const double get_sigma(unsigned int pr_ind);
  //const char get_sigma_type(void);
  const double get_mutot(unsigned int pr_ind);  
  const char get_mutot_type(void); 
};

#endif
