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
	std::string xml = "build_uboone_covar.xml";
	int iarg = 0;
	opterr=1;
	int index;
	
	const struct option longopts[] = 
	{
		{"xml", 		required_argument, 	0, 'x'},
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
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}

	std::string tag = "LEEtest";

	//Load up the calculated (or another) spectrum.
	SBNspec central_value_spec("LEEtest.SBNspec.root",xml);

	//Create a second SBNspec for a BKG-only cv spectrum
	SBNspec bkg_only_spec = central_value_spec;
	//But remove the signal componants
	bkg_only_spec.Scale("elike_signal",0.0);

	/******************LEE signal done seperately example ***********************/
	double bins[4] = {0,500,1000,3000};
	double made_up_number = 20;

	TH1D  * lee_signal = new TH1D("lee","lee",3, bins);
	lee_signal->SetBinContent(1,made_up_number); //Fill this anyway you want. Should be scaled to the correct POT.

	//The following line would then add this ad-hoc TH1 to the "nu_uBooNE_elike_intrinsic" spectra
	//It will break here as the example XML has more than 3 bins! But is here for syntax reasons.
	//central_value_spec.Add("nu_uBooNE_elike_intrinsic",lee_signal);	

	//WARNING! you must check yourself to make sure POT, scaling, etc..  is correct with this TH1, its outside the SBNfit .xml system so bit more dangerous.
	/*****************************************/

	//Load up the calculated covariance matrix	
	TFile *covar_file = new TFile("LEEtest.SBNcovar.root","read");
	TMatrixD *full_covariance = (TMatrixD*)covar_file->Get("full_covariance_LEEtest");
	
	SBNchi uboone_chi(central_value_spec, *full_covariance);

	std::cout<<"CHI^2: "<<uboone_chi.calcChi(&bkg_only_spec)<<std::endl;

	uboone_chi.printMatricies(tag);

}
