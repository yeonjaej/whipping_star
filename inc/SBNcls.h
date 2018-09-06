#ifndef SBNCLS_H_
#define SBNCLS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNchi.h"
#include "SBNconfig.h"

#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"

#include <ctime>
#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

namespace sbn{


class SBNcls{

	public:

	SBNspec * H0;
	SBNspec * H1;
	
	SBNchi chi;
	TMatrixD covariance_matrix;

	TRandom3 * rangen;


	SBNcls(SBNspec *inH0, SBNspec * inH1, TMatrixD matin) : H0(inH0), H1(inH1), covariance_matrix(matin), chi(*inH0, matin){

		rangen= new TRandom3(0);
	}

	int calcCLS(int);



};


};
#endif
