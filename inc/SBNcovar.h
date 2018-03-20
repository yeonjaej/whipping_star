#ifndef SBNCOVAR_H_
#define SBNCOVAR_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNconfig.h"

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


class SBNcovar : public SBNconfig{
	bool is_small_negative_eigenvalue;
	

	public:
	double tolerence_positivesemi;
	int universes_used;
	int Nfiles;

	SBNspec spec_CV;	
	SBNspec spec_sig;

	//these are vectors so 1 for each file
	std::vector<TH2D> h_spec2d_CV;
	std::vector<TH2D> h_spec2d_sig;
	

	SBNspec template_spec;	
	SBNspec spec_CV2;	
	SBNcovar(std::string xmlname);

	// a vector of num_multisim vectors, with a vector of subchannel*bin histograms in each	
	std::vector<SBNspec> multi_sbnspec;
	std::vector<std::vector<double>> multi_vecspec;
	std::vector<std::vector<double>> multi_vecspec2D;
	std::vector<double> vecspec2DCV;
	int formCovarianceMatrix();
	int writeOut();

	TMatrixD full_covariance;
	TMatrixD frac_covariance;
	TMatrixD full_correlation;

	TMatrixD full_covariance2D;
	TMatrixD frac_covariance2D;
	TMatrixD full_correlation2D;



	TH2D * hist_frac_cov;
	TH2D * hist_full_cor;// = new TH2D("Corr","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
	TH2D * hist_full_cov;// = new TH2D("Full Cov","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);




	std::vector<std::vector<int> > vars_i;
	std::vector<std::vector<double> > vars_d;

	std::vector<TFile *> files;
	std::vector<TTree *> trees;


	std::vector<int> nentries;

	std::vector< TBranch *> * bWeight;
	std::vector< TBranch *> * bLepMom;

	std::vector< std::map<std::string, std::vector<double> > * > * fWeights;
	std::vector<TLorentzVector * > *fLepMom;



	virtual bool eventSelection(int file);
	virtual int fillHistograms(int file, int uni, double wei);

	//one for each file



	THnSparseD *frac_covariance_4d;

};


};
#endif
