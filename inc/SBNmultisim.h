#ifndef SBNMULTISIM_H_
#define SBNMULTISIM_H_

#include <cmath>
#include <vector>
#include <iostream>

#include "SBNspec.h"
#include "SBNconfig.h"
#include "SBNchi.h"

#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TLine.h"

#include "TROOT.h"
#include "TRint.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THnSparse.h"

#include <map>
#include <ctime>
#include "params.h"

namespace sbn{


class SBNmultisim : public SBNconfig{
	bool is_small_negative_eigenvalue;
	

	public:
		
	SBNspec spec_CV;	

	std::vector<TH2D> h_spec2d_CV;
	std::vector<TH2D> h_spec2d_sig;
	
	
	SBNmultisim(std::string xmlname);

	int formCovarianceMatrix(std::string tag);
	int writeOut();
	int printMatricies(std::string tag);

	int qualityTesting();
	virtual bool eventSelection(int file);
	virtual int fillHistograms(int file, int uni, double wei);
	

	std::vector<SBNspec> multi_sbnspec;
	std::vector<std::vector<double>> multi_vecspec;

	TMatrixD full_covariance;
	TMatrixD frac_covariance;
	TMatrixD full_correlation;

	std::vector<TMatrixD> vec_full_covariance;
	std::vector<TMatrixD> vec_frac_covariance;
	std::vector<TMatrixD> vec_full_correlation;

	//for plotting, hsa been superseeded by printMatricies a bit
	TH2D * hist_frac_cov;
	TH2D * hist_full_cor;
	TH2D * hist_full_cov;

	//Some checks on multisims
	double tolerence_positivesemi;
	int universes_used;
	double abnormally_large_weight;


	//Multisim input variables
	std::vector<std::string> variations;
	std::vector<std::vector<double> > vars;

	std::vector<int> num_universes_per_variation;
	std::map<int, std::string> map_universe_to_var;



	int Nfiles;
	std::vector<TFile *> files;
	std::vector<TTree *> trees;

	std::vector<int> nentries;
	std::vector< TBranch *> * bWeight;
	std::vector< std::map<std::string, std::vector<double> > * > * fWeights;

	std::vector<std::vector<int> > vars_i;
	std::vector<std::vector<double> > vars_d;




	// In testing 2D fits with 4D covariance..
	SBNspec template_spec;	
	SBNspec spec_CV2;	
	THnSparseD *frac_covariance_4d;
	TMatrixD full_covariance2D;
	TMatrixD frac_covariance2D;
	TMatrixD full_correlation2D;
	std::vector<std::vector<double>> multi_vecspec2D;
	std::vector<double> vecspec2DCV;



};


};
#endif
