#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main(int argc, char* argv[])
{
	std::string xml = "default.xml";
	std::string bkg = "../../data/precomp/SBN_bkg_all";
	std::string sig = "../../data/precomp/SBN_bkg_all";
	int iarg = 0;
	opterr=1;
	int index;
	
	// uboone pot scaling factor (explained below)
	double uboonepot = 0.5;	
	
	const struct option longopts[] = 
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"bkg",			required_argument,	0, 'b'},
		{"sig",			required_argument,	0, 's'},
		{0,			no_argument, 		0,  0},
	};

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:b:s:", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'b':
				bkg = optarg;
				uboonepot = 1.0;
				break;
			case 's':
				sig = optarg;
				break;	
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}

	
	TFile *fin = new TFile("covariance_matrix_ALL.root","read");
	TMatrixD * m = (TMatrixD*)fin->Get("TMatrixT<double>;2");

	SBNspec cv_spec("SBN_CV.root",xml);
	SBNspec bkg_spec = cv_spec;
	bkg_spec.Scale("nu_uBooNE_elike_signal",0.0);
	
	cv_spec.compressVector();
	bkg_spec.compressVector();

	SBNchi mychi(cv_spec, *m);


	std::cout<<"CHI^2: "<<mychi.CalcChi(&bkg_spec)<<std::endl;


	mychi.printMatricies("test_coll.root");


}
