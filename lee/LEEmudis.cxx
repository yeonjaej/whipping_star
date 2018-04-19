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
#include "SBNgenerate.h"
#include "prob.h"

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

	std::string tag = "LEEgen";

	if(true){

	
	//model mass, ue4 um4 
	neutrinoModel testModel(1.4, 0.1, 0.1);
	// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin62 amp, 1 CV oscilatted 
	SBNgenerate gen(xml,testModel);
	// Write them
	gen.writePrecomputedOscSpecs(tag);
	
	//can get unoscillate spectra from sbn_CV

	return 0;
	}	

	SBNosc bkg("LEEgen.SBNspec.root",xml, testModel);

	//Load up the calculated covariance matrix	
	TMatrixD statonly(bkg.num_bins_total, bkg.num_bins_total);
	statonly.Zero();
	
	
	SBNchi uboone_chi(bkg, statonly);

	if(false){

		for(dm){
		for(Um4){
			SBNosc osc("LEEgen.SBNspec.root",xml, testModel);
			osc.OscillateThis(tag);
			//osc.writeOut("osctest");
	
			double ans = uboone_chi.calcChi(&osc);

			std::cout<<dm<<" "<<um4<<" "<<ans<<std::endl;
			// for plotting/looking 
			//	bkg.compareSBNspecs(&osc,tag);
		}
		return 0;
	}
}
