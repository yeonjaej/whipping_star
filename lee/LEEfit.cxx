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

	// Load up background. Keep in mind that UBooNE was considered for double POT for SBN, so we'll amend that here
	SBNspec bkg_spec(bkg.c_str(),xml);
	bkg_spec.Scale("uboone",uboonepot);
	bkg_spec.compressVector();

	SBNchi test_chi(bkg_spec);

	SBNspec sig_spec(sig.c_str(),xml);
	sig_spec.compressVector();

 	double tchi = test_chi.CalcChi(sig_spec);
	std::cout << "The two spectra agree with a chi2 of " << tchi << std::endl;
}
