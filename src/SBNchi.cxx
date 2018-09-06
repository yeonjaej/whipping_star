#include "SBNchi.h"
using namespace sbn;


/***********************************************
 *		Constructors
 * ********************************************/

SBNchi::SBNchi(std::string xml) : SBNconfig(xml,false){};

SBNchi::SBNchi(SBNspec in, TMatrixT<double> Msysin) : SBNconfig(in.xmlname), bkgSpec(in){
	lastChi = -9999999;
	stat_only= false;

	matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
	Msys.ResizeTo(Msysin.GetNrows(), Msysin.GetNcols());
	MfracCov.ResizeTo(Msysin.GetNrows(), Msysin.GetNcols());

	MfracCov = Msysin;
	Msys.Zero();

	this->reload_core_spec(&bkgSpec);
}

//Alternative constrctors
SBNchi::SBNchi(SBNspec in, std::string newxmlname) : SBNconfig(newxmlname), bkgSpec(in){
	stat_only = false;

	matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
	MfracCov.ResizeTo(num_bins_total, num_bins_total);

	if(fullnames.size() !=in.fullnames.size()){
		std::cerr<<"ERROR: SBNchi::SBNchi | Selected covariance matrix and background spectrum are different sizes!"<<std::endl;
		exit(EXIT_FAILURE);
	}else{
		for(int i=0; i< fullnames.size(); i++){
			if(fullnames[i]!=in.fullnames[i]){
				std::cerr<<"ERROR: SBNchi::SBNchi | Spectrum and Covariance matrix have different (or different order) subchannels!"<<std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	MfracCov = sys_fill_direct();
	lastChi = -9999999;
	bkgSpec.collapseVector();
	Msys.ResizeTo(num_bins_total,num_bins_total);
	Msys.Zero();
	Msys=MfracCov;

	this->reload_core_spec(&in);
}

SBNchi::SBNchi(SBNspec in) : SBNchi(in,false){}


SBNchi::SBNchi(SBNspec in, bool is_stat_only): SBNconfig(in.xmlname), bkgSpec(in), stat_only(is_stat_only){
	lastChi = -9999999;

	matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
	Msys.ResizeTo(num_bins_total, num_bins_total);
	MfracCov.ResizeTo(num_bins_total, num_bins_total);


	if(is_stat_only){
		MfracCov.Zero();
		Msys.Zero();

	}else{
		MfracCov = sys_fill_direct();
		Msys.Zero();
	}


	this->reload_core_spec(&bkgSpec);

}


/***********************************************
 *		Rest for now
 * ********************************************/

int SBNchi::reload_core_spec(SBNspec *bkgin){
	otag = "SBNchi::reload_core_spec\t|| ";

	bool is_fractional = true;

	if(isVerbose)std::cout<<otag<<"Begininning to reload core spec! First set new core spec"<<std::endl;
	bkgSpec = *bkgin;
	bkgSpec.collapseVector();

	if(isVerbose)std::cout<<otag<<" || Clear all previous chi^2 data"<<std::endl;
	lastChi_vec.clear();
	lastChi_vec.resize(num_bins_total_compressed, std::vector<double>( num_bins_total_compressed,0) );

	//Reset Msys to fractional
	if(isVerbose) std::cout<<otag<<" Reseting Msys to MfracCov"<<std::endl;
	Msys = MfracCov;

	if(Msys.GetNcols()!=num_bins_total ){
		std::cout<<otag<<"ERROR: trying to pass a matrix to SBNchi that isnt the right size"<<std::endl;
		std::cout<<otag<<"ERROR: num_bins_total: "<<num_bins_total<<" and matrix is: "<<Msys.GetNcols()<<std::endl;
		exit(EXIT_FAILURE);
	}

	if(isVerbose)std::cout<<otag<<"Go from fracCovariance to fullCovariance. Msys.GetNcols(): "<<Msys.GetNcols()<<" Msys.GetNrows(): "<<Msys.GetNrows()<<" core->fullvec.size(): "<<bkgSpec.fullVec.size()<<std::endl;
	// systematics per scaled event
	for(int i =0; i<Msys.GetNcols(); i++)
	{
		for(int j =0; j<Msys.GetNrows(); j++)
		{
			if(is_fractional){
				if(isnan(Msys(i,j)))
					Msys(i,j) = 0;
				else
					Msys(i,j)=Msys(i,j)*bkgSpec.fullVec.at(i)*bkgSpec.fullVec.at(j);
			}
		}
	}

	if(isVerbose)std::cout<<otag<<"Filling stats into cov matrix"<<std::endl;
	// Fill stats from the back ground vector
	TMatrixT <double> Mstat(num_bins_total, num_bins_total);
	stats_fill(Mstat, bkgSpec.fullVec);

	if(Mstat.IsSymmetric()){
		if(isVerbose)std::cout<<otag<<"Stat matrix is symmetric (it is just diagonal core)"<<std::endl;
	}else{
		std::cout<<otag<<"ERROR: SBNchi::formCovarianceMatrix, stats  is not symmetric!"<<std::endl;
		exit(EXIT_FAILURE);
	}

	//And then define the total covariance matrix in all its glory
	TMatrixT <double> Mtotal(num_bins_total,num_bins_total);
	Mtotal.Zero();

	if(stat_only){
		if(isVerbose)std::cout<<otag<<"Using stats only in covariance matrix"<<std::endl;
		Mtotal = Mstat;
	}else{
		if(isVerbose)std::cout<<otag<<" Using stats+sys in covariance matrix"<<std::endl;
		Mtotal = Mstat + Msys;
	}

	if(isVerbose)std::cout<<otag<<"Mstat: "<<Mstat.GetNrows()<<" x "<<Mstat.GetNcols()<<std::endl;
	if(isVerbose)std::cout<<otag<<"Msys: "<<Msys.GetNrows()<<" x "<<Msys.GetNcols()<<std::endl;
	if(isVerbose)std::cout<<otag<<"Mtotal: "<<Mtotal.GetNrows()<<" x "<<Mtotal.GetNcols()<<std::endl;

	if(Mtotal.IsSymmetric() ){
		if(isVerbose)	std::cout<<otag<<"Total Mstat +Msys is symmetric"<<std::endl;
	}else{

		double tol = 1e-13;
		double biggest_deviation = 0;
		int bi =0;
		int bj=0;

		if(isVerbose)std::cout<<otag<<"WARNING: Stats + sys result appears to be not symmetric!"<<std::endl;
		for(int i=0; i<Mtotal.GetNrows(); i++){
			for(int j=0; j<Mtotal.GetNcols(); j++){
				double dev = fabs(Mtotal(i,j)-Mtotal(j,i));
				if(dev>biggest_deviation){
					biggest_deviation = 2*dev/(fabs(Mtotal(i,j))+fabs(Mtotal(j,i)));
					bi=i;
					bj=j;
				}
			}
		}

		if(isVerbose)	std::cout<<otag<<"WARNING: Biggest Relative Deviation from symmetry is i:"<<bi<<" j: "<<bj<<" of order "<<biggest_deviation<<" M(j,i)"<<Mtotal(bj,bi)<<" M(i,j)"<<Mtotal(bi,bj)<<std::endl;

		if(biggest_deviation >tol){

			std::cout<<"ERROR: Thats too unsymettric, killing process. Better check your inputs."<<std::endl;

			exit(EXIT_FAILURE);
		}else{

			if(isVerbose)	std::cout<<otag<<"WARNING: Thats within tolderence. Continuing."<<std::endl;
		}
	}

	TMatrixT<double > Mctotal(num_bins_total_compressed,num_bins_total_compressed);

	collapse_layer3(Mtotal, Mctotal);

	matrix_collapsed = Mctotal;

	vMc = to_vector(Mctotal);
	double invdet=0;

	TMatrixD McI(num_bins_total_compressed,num_bins_total_compressed);
	McI.Zero();

	std::cout<<otag<<" About to do a SVD decomposition"<<std::endl;
	TDecompSVD svd(Mctotal);
	if (!svd.Decompose()) {
		std::cout <<otag<<"Decomposition failed, matrix not symettric?, has nans?" << std::endl;
		std::cout<<otag<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;
		exit(EXIT_FAILURE);
	} else {
		McI = svd.Invert();
	}
	if( !McI.IsValid()){
		std::cout<<otag<<"ERROR: The inverted matrix isnt valid! Something went wrong.."<<std::endl;
		exit(EXIT_FAILURE);

	}

	std::cout<<otag<<"SUCCESS! Inverted."<<std::endl;
	//McI.Print();
	vMcI = to_vector(McI);


	// test for validity
	bool is_small_negative_eigenvalue = false;
	double tolerence_positivesemi = 1e-5;


	//if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
	TMatrixDEigen eigen (Mtotal);
	TVectorD eigen_values = eigen.GetEigenValuesRe();


	for(int i=0; i< eigen_values.GetNoElements(); i++){
		if(eigen_values(i)<0){
			is_small_negative_eigenvalue = true;
			if(fabs(eigen_values(i))> tolerence_positivesemi ){
				std::cout<<otag<<"covariance matrix contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}


	if(is_small_negative_eigenvalue){
		if(isVerbose)	std::cout<<otag<<"Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
	}else{
		if(isVerbose)	std::cout<<otag<<"Generated covariance matrix is also positive semi-definite."<<std::endl;
	}

	bkgSpec.collapseVector();

	return 0;
}



/*********************************************
 *		Different Chi^2 calculations
 * ******************************************/

//Standard chi^2 calculation
double SBNchi::calcChi(SBNspec *sigSpec){
	double tchi = 0;

	if(sigSpec->compVec.size()==0){
		if(isVerbose)	std::cout<<"WARNING: SBNchi::calcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		sigSpec->collapseVector();
	}

	int k=0;

	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			k++;

			if(i==j && vMcI.at(i).at(j)<0){
				std::cout<<"ERROR: SBNchi::calcChi || diagonal of inverse covariance is negative!"<<std::endl;
			}
			lastChi_vec.at(i).at(j) =(bkgSpec.compVec.at(i)-sigSpec->compVec.at(i))*vMcI.at(i).at(j)*(bkgSpec.compVec.at(j)-sigSpec->compVec.at(j) );
			tchi += lastChi_vec.at(i).at(j);
		}
	}

	lastChi = tchi;
	return tchi;
}


//same as above but passing in a vector instead of whole SBNspec
double SBNchi::calcChi(std::vector<double> sigVec){
	double tchi = 0;

	if(sigVec.size() != num_bins_total_compressed ){
		std::cerr<<"ERROR: SBNchi::calcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
		std::cerr<<"sigVec.size(): "<<sigVec.size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
		exit(EXIT_FAILURE);
	}

	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			tchi += (bkgSpec.compVec[i]-sigVec[i])*vMcI[i][j]*(bkgSpec.compVec[j]-sigVec[j] );
		}
	}

	lastChi = tchi;
	return tchi;
}

//A log-lilihood based one used @ MiniBooNE
double SBNchi::calcChiLog(SBNspec *sigSpec){
	double tchi = 0;

	if(sigSpec->compVec.size()==0){
		std::cout<<"WARNING: SBNchi::calcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		sigSpec->collapseVector();
	}

	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			lastChi_vec.at(i).at(j) =(bkgSpec.compVec[i]-sigSpec->compVec[i])*vMcI[i][j]*(bkgSpec.compVec[j]-sigSpec->compVec[j] );
			tchi += lastChi_vec.at(i).at(j);
		}
	}

	double absDetM = log(fabs(matrix_collapsed.Determinant()));

	lastChi = tchi+absDetM;
	return tchi+absDetM;
}


double SBNchi::calcChi(SBNspec *sigSpec, SBNspec *obsSpec){
	double tchi=0;
	if(sigSpec->compVec.size()==0){
		std::cout<<"WARNING: SBNchi::calcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		sigSpec->collapseVector();
	}

	if(obsSpec->compVec.size()==0){
		std::cout<<"WARNING: SBNchi::calcChi, inputted obsSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		obsSpec->collapseVector();
	}


	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			tchi += (obsSpec->compVec[i]-sigSpec->compVec[i])*vMcI[i][j]*(obsSpec->compVec[j]-sigSpec->compVec[j] );
		}
	}

	lastChi = tchi;
	return tchi;

}



/**************************************************************************
 *			Collapsing code
 * ************************************************************************/


//This is the powerhouse, takes each detector matrix filled with num_channels channels of num_subchannels[i] subchannels, and collapses it.
void SBNchi::collapse_layer1(TMatrixT <double> & M, TMatrixT <double> & Mc){
	bool debug = false;
	if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<115<<std::endl;
	if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<30<<std::endl;

	std::vector<std::vector<TMatrixT<double>>> Summed(num_channels, std::vector<TMatrixT<double>>(num_channels) );	//Initialise a matrix of matricies, to ZERO.
	for(int ic = 0; ic < num_channels; ic++){
		for(int jc =0; jc < num_channels; jc++){
			Summed[ic][jc].ResizeTo(num_bins[jc],num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
			Summed[ic][jc] = 0.0;
		}
	}

	int mrow = 0.0;
	int mcol = 0.0;

	for(int ic = 0; ic < num_channels; ic++){ 	 //Loop over all rows
		for(int jc =0; jc < num_channels; jc++){ //Loop over all columns

			if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;

			for(int m=0; m < num_subchannels[ic]; m++){
				for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing together
					Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
				}
			}


			mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
		}//end of column loop

		mrow = 0; // as we end this row, reset row count, but jump down 1 column
		mcol += num_subchannels[ic]*num_bins[ic];
	}//end of row loop



	///********************************* And put them back together! ************************//
	Mc.Zero();
	mrow = 0;
	mcol = 0;

	//Repeat again for Contracted matrix
	for(int ic = 0; ic < num_channels; ic++){
		for(int jc =0; jc < num_channels; jc++){

			Mc.SetSub(mrow,mcol,Summed[ic][jc]);
			mrow += num_bins[jc];
		}

		mrow = 0;
		mcol +=num_bins[ic];
	}

	return;
}



//This is the detector layer, Take a given mode and run over each detector V detector sub matrix
void SBNchi::collapse_layer2(TMatrixT <double> & M, TMatrixT <double> & Mc){

	Mc.Zero();
	int nrow = num_bins_detector_block;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
	int crow = num_bins_detector_block_compressed; //N_e_bins+N_m_bins;

	for(int m =0; m< num_detectors; m++){
		for(int n =0; n< num_detectors; n++){
			TMatrixT<double> imat(nrow,nrow);
			TMatrixT<double> imatc(crow,crow);

			imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
			collapse_layer1(imat,imatc);
			Mc.SetSub(n*crow,m*crow,imatc);
		}
	}

	return;
}

//This is the Mode layer, Take a given full matrix and runs over each Mode V Mode sub matrix
void SBNchi::collapse_layer3(TMatrixT <double> & M, TMatrixT <double> & Mc){

	Mc.Zero();
	int nrow = num_bins_mode_block;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
	int crow=  num_bins_mode_block_compressed;// (N_e_bins+N_m_bins)*N_dets;

	for(int m =0; m< num_modes ; m++){
		for(int n =0; n< num_modes; n++){

			TMatrixT<double> imat(nrow,nrow);
			TMatrixT<double> imatc(crow,crow);

			imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);

			collapse_layer2(imat,imatc);
			Mc.SetSub(n*crow,m*crow,imatc);

		}
	}

	return;
}





/**************************************************************************
 *			Misc
 * ************************************************************************/


int SBNchi::fillCollapsedCovarianceMatrix(TMatrixT<double>*in){
	in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
	for(int i=0; i<num_bins_total_compressed;i++){
		for(int j=0; j<num_bins_total_compressed;j++){
			(*in)(i,j) = vMc.at(i).at(j);
		}
	}

	return 0;
}


int SBNchi::fillCollapsedCorrelationMatrix(TMatrixT<double>*in){
	in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
	for(int i=0; i<num_bins_total_compressed;i++){
		for(int j=0; j<num_bins_total_compressed;j++){
			(*in)(i,j) = vMc.at(i).at(j)/(sqrt(vMc.at(j).at(j))*sqrt(vMc.at(i).at(i)));
		}
	}

	return 0;
}

int SBNchi::fillCollapsedFractionalMatrix(TMatrixT<double>*in){
	in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
	for(int i=0; i<num_bins_total_compressed;i++){
		for(int j=0; j<num_bins_total_compressed;j++){
			(*in)(i,j) = vMc.at(i).at(j)/(bkgSpec.compVec.at(i)*bkgSpec.compVec.at(j));
		}
	}

	return 0;
}


TMatrixT<double> * SBNchi::getCompressedMatrix(){
	TMatrixT<double> * tmp = new TMatrixT<double>(num_bins_total_compressed,num_bins_total_compressed);
	for(int i=0; i<num_bins_total_compressed;i++){
		for(int j=0; j<num_bins_total_compressed;j++){
			(*tmp)(i,j) = vMc.at(i).at(j);
		}
	}

	return tmp;
}




void SBNchi::fake_fill(TMatrixT <double> &M){
	//Fills a square matrix of dim matrix_size with random numbers for now.
	TRandom3 *rangen = new TRandom3(0);
	int matrix_size=M.GetNrows();
	if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}
	for(int i=0; i<matrix_size; i++){
		for (int j = i;j<matrix_size;j++){
			M(i,j)=rangen->Uniform(0,1);
			M(j,i)=M(i,j);
		}
	}
	return ;
}


std::vector<std::vector<double >> SBNchi::to_vector(TMatrixT <double > Min)
{
	int dimension =  Min.GetNrows();

	std::vector<std::vector<double >>  ans(dimension, std::vector<double>(dimension));

	for(int i = 0; i< dimension; i++){
		for(int k = 0; k< dimension; k++){
			ans[i][k]=Min(i,k);
			if(ans[i][k]==-0){
				ans[i][k]=0;
			}
		}
	}
	return ans;
}

void SBNchi::stats_fill(TMatrixT <double> &M, std::vector<double> diag){
	int matrix_size = M.GetNrows();

	if(matrix_size != diag.size()){std::cout<<"#ERROR: stats_fill, matrix not equal to diagonal"<<std::endl;}
	if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}

	M.Zero();

	for(int i=0; i<matrix_size; i++)
	{

		//This NEEDS to be removed soon
		//This was just for wierd MiniBooNE run
		//if(i>=11 && i< 30) continue;
		//if(i>=41) continue;
		M(i,i) = diag.at(i);

	}

	return ;
}



TMatrixT<double> SBNchi::sys_fill_direct(){
	return sys_fill_direct(correlation_matrix_rootfile, correlation_matrix_name);
}


TMatrixT<double > SBNchi::sys_fill_direct(std::string rootname, std::string matname){
	//Pretty much obsolete now, should fill directly really.
	std::cout<<"SBNchi::sys_fill_direct || filling from "<<rootname<<std::endl;

	TMatrixT<double> temp2(num_bins_total,num_bins_total);
	TFile *fm= new TFile(rootname.c_str());

	TMatrixT<float> * temp = (TMatrixT <float>* )fm->Get(matname.c_str());
	//TMatrixT<double> * temp = (TMatrixT <double>* )fm->Get(matname.c_str());

	std::vector<std::vector<double>> mcont;

	for(int p:used_bins){
		std::vector<double> tvec;
		for(int u:used_bins){
			tvec.push_back( (*temp)(p,u) );
		}
		mcont.push_back(tvec);
	}

	for(int i =0; i<num_bins_total; i++)
	{
		for(int j =0; j<num_bins_total; j++)
		{
			temp2(i,j)=mcont[i][j];
		}
	}
	delete temp;

	std::cout<<"SBNchi::sys_fill_direct || loaded with dim : "<<temp2.GetNcols()<<" "<<temp2.GetNrows()<<std::endl;

	fm->Close();
	delete fm;

	if(temp2.IsSymmetric()){
		if(isVerbose)std::cout<<"Inputted fracCov covariance matrix is symmetric"<<std::endl;
	}else{
		std::cerr<<"ERROR: SBNchi::sys_fill_direct, Msys input is not symmetric!"<<std::endl;
		//exit(EXIT_FAILURE);
	}

	return temp2;

}





TH2D* SBNchi::getChiogram(){
	TH2D *tmp = new TH2D("chi-o-gram","chi-o-gram",num_bins_total_compressed,0, num_bins_total_compressed ,num_bins_total_compressed,0, num_bins_total_compressed);

	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){

			tmp->SetBinContent(i+1, j+1, lastChi_vec.at(i).at(j));
		}
	}

	return tmp;
}

int SBNchi::printMatricies(std::string tag){
	TFile* fout = new TFile(("SBNfit_collapsed_matrix_plots_"+tag+".root").c_str(),"recreate");
	fout->cd();


	gStyle->SetOptStat(0);

	TMatrixD full, frac, corr;
	this->fillCollapsedCovarianceMatrix(&full);
	this->fillCollapsedFractionalMatrix(&frac);
	this->fillCollapsedCorrelationMatrix(&corr);

	corr.Write();
	TH2D h2_corr(corr);
	h2_corr.SetName("corr");
	//h2_corr.Write();
	TCanvas *c_corr = new TCanvas("collapsed correlation matrix");
	c_corr->cd();
	c_corr->SetFixedAspectRatio();
	h2_corr.Draw("colz");
	h2_corr.SetTitle("Collapsed correlation matrix");
	h2_corr.GetXaxis()->SetTitle("Reco Bin i");
	h2_corr.GetYaxis()->SetTitle("Reco Bin j");

	c_corr->SetRightMargin(0.150);

	int use_corr =0;
	for(int im =0; im<num_modes; im++){
		for(int id =0; id<num_detectors; id++){
			for(int ic = 0; ic < num_channels; ic++){
				TLine *lv = new TLine(0, num_bins.at(ic)+use_corr, num_bins_total_compressed, num_bins.at(ic)+use_corr);
				TLine *lh = new TLine(num_bins.at(ic)+use_corr,0, num_bins.at(ic)+use_corr, num_bins_total_compressed);
				lv->SetLineWidth(1.5);
				lh->SetLineWidth(1.5);
				use_corr+=num_bins.at(ic);
				lv->Draw();
				lh->Draw();
			}
		}
	}
	c_corr->Write();


	frac.Write();
	TH2D h2_frac(frac);
	//h2_frac.Write();
	h2_frac.SetName("frac");
	TCanvas *c_frac = new TCanvas("collapsed fractional covariance matrix");
	c_frac->cd();
	c_frac->SetFixedAspectRatio();
	h2_frac.Draw("colz");
	h2_frac.SetTitle("Collapsed fractional covariance matrix");
	h2_frac.GetXaxis()->SetTitle("Reco Bin i");
	h2_frac.GetYaxis()->SetTitle("Reco Bin j");

	c_frac->SetRightMargin(0.150);

	int use_frac =0;
	for(int im =0; im<num_modes; im++){
		for(int id =0; id<num_detectors; id++){
			for(int ic = 0; ic < num_channels; ic++){
				TLine *lv = new TLine(0, num_bins.at(ic)+use_frac, num_bins_total_compressed, num_bins.at(ic)+use_frac);
				TLine *lh = new TLine(num_bins.at(ic)+use_frac,0, num_bins.at(ic)+use_frac, num_bins_total_compressed);
				lv->SetLineWidth(1.5);
				lh->SetLineWidth(1.5);
				use_frac+=num_bins.at(ic);
				lv->Draw();
				lh->Draw();
			}
		}
	}
	c_frac->Write();


	full.Write();
	TH2D h2_full(full);
	//h2_full.Write();
	h2_corr.SetName("full");
	TCanvas *c_full = new TCanvas("collapsed covariance matrix");
	c_full->cd();
	c_full->SetFixedAspectRatio();
	h2_full.Draw("colz");
	h2_full.SetTitle("Collapsed covariance matrix");
	h2_full.GetXaxis()->SetTitle("Reco Bin i");
	h2_full.GetYaxis()->SetTitle("Reco Bin j");

	c_full->SetRightMargin(0.150);

	int use_full =0;
	for(int im =0; im<num_modes; im++){
		for(int id =0; id<num_detectors; id++){
			for(int ic = 0; ic < num_channels; ic++){
				TLine *lv = new TLine(0, num_bins.at(ic)+use_full, num_bins_total_compressed, num_bins.at(ic)+use_full);
				TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total_compressed);
				lv->SetLineWidth(1.5);
				lh->SetLineWidth(1.5);
				use_full+=num_bins.at(ic);
				lv->Draw();
				lh->Draw();
			}
		}
	}
	c_full->Write();


	TCanvas *chiogram = new TCanvas("Chi-o-gram","Chi-o-gram");
	chiogram->cd();
	chiogram->SetFixedAspectRatio();

	TH2D * h_chiogram = (TH2D*)this->getChiogram();
	h_chiogram->Draw("colz");
	h_chiogram->SetTitle("Reco Bin i");
	h_chiogram->SetTitle("Reco Bin j");

	int use_chio =0;
	for(int im =0; im<num_modes; im++){
		for(int id =0; id<num_detectors; id++){
			for(int ic = 0; ic < num_channels; ic++){
				TLine *lv = new TLine(0, num_bins.at(ic)+use_chio, num_bins_total_compressed, num_bins.at(ic)+use_chio);
				TLine *lh = new TLine(num_bins.at(ic)+use_chio,0, num_bins.at(ic)+use_chio, num_bins_total_compressed);
				lv->SetLineWidth(1.5);
				lh->SetLineWidth(1.5);
				use_chio+=num_bins.at(ic);
				lv->Draw();
				lh->Draw();
			}
		}
	}
	chiogram->Write();

	fout->Close();
	return 0;
}


//This one varies the input comparative spectrum, and as sucn has  only to calculate the Msys once
TH1D SBNchi::toyMC_varyInput(SBNspec *specin, int num_MC){
	double center = this->calcChi(specin);
	int nlower=0;

	TRandom3 *rangen = new TRandom3(0);

	TH1D ans("","",100,0,100);
	//So save the core one that we will sample for
	ans.GetXaxis()->SetCanExtend(kTRUE);
	isVerbose = false;
	for(int i=0; i < num_MC;i++){

		SBNspec tmp = *specin;
		tmp.poissonScale(rangen);
		tmp.collapseVector(); //this line important isnt it!
		//tmp.printfull_vectorFullVec();

		double thischi = this->calcChi(&tmp);
		ans.Fill(thischi);
		if(thischi<=center) nlower++;

		if(i%100==0) std::cout<<"SBNchi::toyMC_varyInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
	}
	std::cout<<"pval: "<<nlower/(double)num_MC<<std::endl;

	isVerbose = true;

	return ans;


}


//This one varies the input comparative spectrum, and as sucn has  only to calculate the Msys once
std::vector<double> SBNchi::toyMC_varyInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
	std::vector<int> nlower(chival.size(),0);

	TRandom3 *rangen = new TRandom3(0);

	TH1D ans("","",100,0,100);
	//So save the core one that we will sample for
	ans.GetXaxis()->SetCanExtend(kTRUE);
	isVerbose = false;
	for(int i=0; i < num_MC;i++){

		SBNspec tmp = *specin;
		tmp.poissonScale(rangen);
		tmp.collapseVector(); //this line important isnt it!
		//tmp.printFullVec();

		double thischi = this->calcChi(&tmp);
		ans.Fill(thischi);

		for(int j=0; j< chival.size(); j++){
			if(thischi>=chival.at(j)) nlower.at(j)++;
		}

		if(i%1000==0) std::cout<<"SBNchi::toyMC_varyInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
	}
	std::vector<double> pval;
	for(auto n: nlower){
		pval.push_back(n/(double)num_MC);

	}

	isVerbose = true;
	return pval;


}


//This one varies the core spectrum, and as sucn has to recalculate the Msys each stem
TH1D SBNchi::toyMC_varyCore(SBNspec *specin, int num_MC){
	double center = this->calcChi(specin);
	int nlower=0;

	TRandom3 *rangen = new TRandom3(0);

	TH1D ans("MCans","MCans",100,center-100,center+200);
	//So save the core one that we will sample for
	SBNspec core  = bkgSpec;

	isVerbose = false;
	for(int i=0; i<num_MC;i++){

		SBNspec tmp = core;
		tmp.poissonScale(rangen);
		this->reload_core_spec(&tmp);
		double thischi = this->calcChi(specin);
		ans.Fill(thischi);
		if(thischi<=center)nlower++;

		if(i%1000==0) std::cout<<"SBNchi::toyMC_varyCore(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
	}
	std::cout<<"pval: "<<nlower/(double)num_MC<<std::endl;

	isVerbose = true;
	return ans;


}
