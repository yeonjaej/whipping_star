#include "SBNchi.h"
#include "openacc.h"
#include "openacc_curand.h"

using namespace sbn;


/***********************************************
 *		Constructors
 * ********************************************/

SBNchi::SBNchi(std::string xml) : SBNconfig(xml,false){};

SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin) : SBNconfig(in.xmlname), core_spectrum(in){
  last_calculated_chi = -9999999;
  is_stat_only= false;

  matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
  matrix_systematics.ResizeTo(num_bins_total, num_bins_total);
  matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);

  TMatrixD m = matrix_systematicsin;
  for (int i = 0; i < used_bins.size(); i++){
    TMatrixDColumn(m,i) = TMatrixDColumn(m,used_bins.at(i));
    m.ResizeTo(used_bins.size(),matrix_systematicsin.GetNcols());
  }
		

  matrix_fractional_covariance = m;
  matrix_systematics.Zero();


  this->ReloadCoreSpectrum(&core_spectrum);
}

//Alternative constrctors
SBNchi::SBNchi(SBNspec in, std::string newxmlname) : SBNconfig(newxmlname), core_spectrum(in){
  is_stat_only = false;

  matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
  matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);

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

  matrix_fractional_covariance = FillSystematicsFromXML();
  last_calculated_chi = -9999999;
  core_spectrum.CollapseVector();
  matrix_systematics.ResizeTo(num_bins_total,num_bins_total);
  matrix_systematics.Zero();
  matrix_systematics=matrix_fractional_covariance;


  this->ReloadCoreSpectrum(&in);
}

SBNchi::SBNchi(SBNspec in) : SBNchi(in,true){}


SBNchi::SBNchi(SBNspec in, bool is_is_stat_only): SBNconfig(in.xmlname), core_spectrum(in), is_stat_only(is_is_stat_only){
  last_calculated_chi = -9999999;

  matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
  matrix_systematics.ResizeTo(num_bins_total, num_bins_total);
  matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);




  if(is_is_stat_only){
    matrix_fractional_covariance.Zero();
    matrix_systematics.Zero();

  }else{
    matrix_fractional_covariance = FillSystematicsFromXML();
    matrix_systematics.Zero();
  }

  this->ReloadCoreSpectrum(&core_spectrum);

}


/***********************************************
 *		Rest for now
 * ********************************************/

int SBNchi::ReloadCoreSpectrum(SBNspec *bkgin){
  otag = "SBNchi::ReloadCoreSpectrum\t|| ";

  bool is_fractional = true;
  cholosky_performed = false;

  if(is_verbose)std::cout<<otag<<"Begininning to reload core spec! First Set new core spec"<<std::endl;
  core_spectrum = *bkgin;
  core_spectrum.CollapseVector();

  if(is_verbose)std::cout<<otag<<" || Clear all previous chi^2 data"<<std::endl;
  vec_last_calculated_chi.clear();
  vec_last_calculated_chi.resize(num_bins_total_compressed, std::vector<double>( num_bins_total_compressed,0) );

  //Reset matrix_systematics to fractional
  if(is_verbose) std::cout<<otag<<" Reseting matrix_systematics to matrix_fractional_covariance"<<std::endl;
  matrix_systematics = matrix_fractional_covariance;

  if(matrix_systematics.GetNcols()!=num_bins_total ){
    std::cout<<otag<<"ERROR: trying to pass a matrix to SBNchi that isnt the right size"<<std::endl;
    std::cout<<otag<<"ERROR: num_bins_total: "<<num_bins_total<<" and matrix is: "<<matrix_systematics.GetNcols()<<std::endl;
    exit(EXIT_FAILURE);
  }

  if(is_verbose)std::cout<<otag<<"Go from fracCovariance to fullCovariance. matrix_systematics.GetNcols(): "<<matrix_systematics.GetNcols()<<" matrix_systematics.GetNrows(): "<<matrix_systematics.GetNrows()<<" core->fullvec.size(): "<<core_spectrum.full_vector.size()<<std::endl;
  // systematics per scaled event
  for(int i =0; i<matrix_systematics.GetNcols(); i++)
    {
      for(int j =0; j<matrix_systematics.GetNrows(); j++)
	{
	  if(is_fractional){
	    if(isnan(matrix_systematics(i,j)))
	      matrix_systematics(i,j) = 0;
	    else
	      matrix_systematics(i,j)=matrix_systematics(i,j)*core_spectrum.full_vector.at(i)*core_spectrum.full_vector.at(j);
	  }
	}
    }

  if(is_verbose)std::cout<<otag<<"Filling stats into cov matrix"<<std::endl;
  // Fill stats from the back ground vector
  TMatrixT <double> Mstat(num_bins_total, num_bins_total);
  FillStatsMatrix(Mstat, core_spectrum.full_vector);

  if(Mstat.IsSymmetric()){
    if(is_verbose)std::cout<<otag<<"Stat matrix is symmetric (it is just diagonal core)"<<std::endl;
  }else{
    std::cout<<otag<<"ERROR: SBNchi::FormCovarianceMatrix, stats  is not symmetric!"<<std::endl;
    exit(EXIT_FAILURE);
  }

  //And then define the total covariance matrix in all its glory
  TMatrixT <double> Mtotal(num_bins_total,num_bins_total);
  Mtotal.Zero();

  if(is_stat_only){
    if(is_verbose)std::cout<<otag<<"Using stats only in covariance matrix"<<std::endl;
    Mtotal = Mstat;
  }else{
    if(is_verbose)std::cout<<otag<<" Using stats+sys in covariance matrix"<<std::endl;
    Mtotal = Mstat + matrix_systematics;
  }

  if(is_verbose)std::cout<<otag<<"Mstat: "<<Mstat.GetNrows()<<" x "<<Mstat.GetNcols()<<std::endl;
  if(is_verbose)std::cout<<otag<<"matrix_systematics: "<<matrix_systematics.GetNrows()<<" x "<<matrix_systematics.GetNcols()<<std::endl;
  if(is_verbose)std::cout<<otag<<"Mtotal: "<<Mtotal.GetNrows()<<" x "<<Mtotal.GetNcols()<<std::endl;

  if(Mtotal.IsSymmetric() ){
    if(is_verbose)	std::cout<<otag<<"Total Mstat +matrix_systematics is symmetric"<<std::endl;
  }else{

    double tol = 1e-13;
    double biggest_deviation = 0;
    int bi =0;
    int bj=0;

    if(is_verbose)std::cout<<otag<<"WARNING: Stats + sys result appears to be not symmetric!"<<std::endl;
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

    if(is_verbose)	std::cout<<otag<<"WARNING: Biggest Relative Deviation from symmetry is i:"<<bi<<" j: "<<bj<<" of order "<<biggest_deviation<<" M(j,i)"<<Mtotal(bj,bi)<<" M(i,j)"<<Mtotal(bi,bj)<<std::endl;

    if(biggest_deviation >tol){

      std::cout<<"ERROR: Thats too unsymettric, killing process. Better check your inputs."<<std::endl;

      exit(EXIT_FAILURE);
    }else{

      if(is_verbose)	std::cout<<otag<<"WARNING: Thats within tolderence. Continuing."<<std::endl;
    }
  }

  TMatrixT<double > Mctotal(num_bins_total_compressed,num_bins_total_compressed);

  CollapseModes(Mtotal, Mctotal);

  matrix_collapsed = Mctotal;

  vec_matrix_collapsed = TMatrixDToVector(Mctotal);
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
  vec_matrix_inverted = TMatrixDToVector(McI);


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
    if(is_verbose)	std::cout<<otag<<"Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
  }else{
    if(is_verbose)	std::cout<<otag<<"Generated covariance matrix is also positive semi-definite."<<std::endl;
  }

  core_spectrum.CollapseVector();

  return 0;
}



/*********************************************
 *		Different Chi^2 calculations
 * ******************************************/

//Standard chi^2 calculation
double SBNchi::CalcChi(SBNspec *sigSpec){
  double tchi = 0;

  if(sigSpec->collapsed_vector.size()==0){
    if(is_verbose)	std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
    sigSpec->CollapseVector();
  }

  int k=0;

  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){
      k++;

      if(i==j && vec_matrix_inverted.at(i).at(j)<0){
	std::cout<<"ERROR: SBNchi::CalcChi || diagonal of inverse covariance is negative!"<<std::endl;
      }
      vec_last_calculated_chi.at(i).at(j) =(core_spectrum.collapsed_vector.at(i)-sigSpec->collapsed_vector.at(i))*vec_matrix_inverted.at(i).at(j)*(core_spectrum.collapsed_vector.at(j)-sigSpec->collapsed_vector.at(j) );
      tchi += vec_last_calculated_chi.at(i).at(j);
    }
  }

  last_calculated_chi = tchi;
  return tchi;
}

double SBNchi::CalcChi(std::vector<double> * sigVec){
  double tchi = 0;

#ifndef _OPENACC
  if(sigVec->size() != num_bins_total_compressed ){
    std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
    std::cerr<<"sigVec.size(): "<<sigVec->size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
    exit(EXIT_FAILURE);
  }
#endif
  
  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){
      tchi += (core_spectrum.collapsed_vector[i]-(*sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(*sigVec)[j] );
    }
  }

  //last_calculated_chi = tchi;
  return tchi;
}

double SBNchi::CalcChi(double* sigVec){
  double tchi = 0;

  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){
      tchi += (core_spectrum.collapsed_vector[i]-(sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(sigVec)[j] );
    }
  }

  //last_calculated_chi = tchi;
  return tchi;
}

double SBNchi::CalcChi(double **invert_matrix, double* core, double *sig){
  double tchi = 0;

  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){
      tchi += (core[i]-sig[i])*invert_matrix[i][j]*(core[j]-sig[j] );
    }
  }

  return tchi;
}




//same as above but passing in a vector instead of whole SBNspec
double SBNchi::CalcChi(std::vector<double> sigVec){
  double tchi = 0;

  if(sigVec.size() != num_bins_total_compressed ){
    std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
    std::cerr<<"sigVec.size(): "<<sigVec.size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
    exit(EXIT_FAILURE);
  }

  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){
      tchi += (core_spectrum.collapsed_vector[i]-sigVec[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-sigVec[j] );
    }
  }

  last_calculated_chi = tchi;
  return tchi;
}

//A log-lilihood based one used @ MiniBooNE
double SBNchi::CalcChiLog(SBNspec *sigSpec){
  double tchi = 0;

  if(sigSpec->collapsed_vector.size()==0){
    std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
    sigSpec->CollapseVector();
  }

  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){
      vec_last_calculated_chi.at(i).at(j) =(core_spectrum.collapsed_vector[i]-sigSpec->collapsed_vector[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-sigSpec->collapsed_vector[j] );
      tchi += vec_last_calculated_chi.at(i).at(j);
    }
  }

  double absDetM = log(fabs(matrix_collapsed.Determinant()));

  last_calculated_chi = tchi+absDetM;
  return tchi+absDetM;
}


double SBNchi::CalcChi(SBNspec *sigSpec, SBNspec *obsSpec){
  double tchi=0;
  if(sigSpec->collapsed_vector.size()==0){
    std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
    sigSpec->CollapseVector();
  }

  if(obsSpec->collapsed_vector.size()==0){
    std::cout<<"WARNING: SBNchi::CalcChi, inputted obsSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
    obsSpec->CollapseVector();
  }


  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){
      tchi += (obsSpec->collapsed_vector[i]-sigSpec->collapsed_vector[i])*vec_matrix_inverted[i][j]*(obsSpec->collapsed_vector[j]-sigSpec->collapsed_vector[j] );
    }
  }

  last_calculated_chi = tchi;
  return tchi;

}



/**************************************************************************
 *			Collapsing code
 * ************************************************************************/


//This is the powerhouse, takes each detector matrix filled with num_channels channels of num_subchannels[i] subchannels, and collapses it.
void SBNchi::CollapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc){
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
	for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing toGether
	  Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
	}
      }


      mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
    }//end of column loop

    mrow = 0; // as we end this row, reSet row count, but jump down 1 column
    mcol += num_subchannels[ic]*num_bins[ic];
  }//end of row loop



  ///********************************* And put them back toGether! ************************//
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
void SBNchi::CollapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc){

  Mc.Zero();
  int nrow = num_bins_detector_block;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
  int crow = num_bins_detector_block_compressed; //N_e_bins+N_m_bins;

  for(int m =0; m< num_detectors; m++){
    for(int n =0; n< num_detectors; n++){
      TMatrixT<double> imat(nrow,nrow);
      TMatrixT<double> imatc(crow,crow);

      imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
      CollapseSubchannels(imat,imatc);
      Mc.SetSub(n*crow,m*crow,imatc);
    }
  }

  return;
}

//This is the Mode layer, Take a given full matrix and runs over each Mode V Mode sub matrix
void SBNchi::CollapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc){

  Mc.Zero();
  int nrow = num_bins_mode_block;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
  int crow=  num_bins_mode_block_compressed;// (N_e_bins+N_m_bins)*N_dets;

  for(int m =0; m< num_modes ; m++){
    for(int n =0; n< num_modes; n++){

      TMatrixT<double> imat(nrow,nrow);
      TMatrixT<double> imatc(crow,crow);

      imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);

      CollapseDetectors(imat,imatc);
      Mc.SetSub(n*crow,m*crow,imatc);

    }
  }

  return;
}





/**************************************************************************
 *			Misc
 * ************************************************************************/


int SBNchi::FillCollapsedCovarianceMatrix(TMatrixT<double>*in){
  in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
  for(int i=0; i<num_bins_total_compressed;i++){
    for(int j=0; j<num_bins_total_compressed;j++){
      (*in)(i,j) = vec_matrix_collapsed.at(i).at(j);
    }
  }

  return 0;
}


int SBNchi::FillCollapsedCorrelationMatrix(TMatrixT<double>*in){
  in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
  for(int i=0; i<num_bins_total_compressed;i++){
    for(int j=0; j<num_bins_total_compressed;j++){
      (*in)(i,j) = vec_matrix_collapsed.at(i).at(j)/(sqrt(vec_matrix_collapsed.at(j).at(j))*sqrt(vec_matrix_collapsed.at(i).at(i)));
    }
  }

  return 0;
}

int SBNchi::FillCollapsedFractionalMatrix(TMatrixT<double>*in){
  in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
  for(int i=0; i<num_bins_total_compressed;i++){
    for(int j=0; j<num_bins_total_compressed;j++){
      (*in)(i,j) = vec_matrix_collapsed.at(i).at(j)/(core_spectrum.collapsed_vector.at(i)*core_spectrum.collapsed_vector.at(j));
    }
  }

  return 0;
}


TMatrixT<double> * SBNchi::GetCollapsedMatrix(){
  TMatrixT<double> * tmp = new TMatrixT<double>(num_bins_total_compressed,num_bins_total_compressed);
  for(int i=0; i<num_bins_total_compressed;i++){
    for(int j=0; j<num_bins_total_compressed;j++){
      (*tmp)(i,j) = vec_matrix_collapsed.at(i).at(j);
    }
  }

  return tmp;
}




void SBNchi::FakeFillMatrix(TMatrixT <double> &M){
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


std::vector<std::vector<double >> SBNchi::TMatrixDToVector(TMatrixT <double > Min)
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

void SBNchi::FillStatsMatrix(TMatrixT <double> &M, std::vector<double> diag){
  int matrix_size = M.GetNrows();

  if(matrix_size != diag.size()){std::cout<<"#ERROR: FillStatsMatrix, matrix not equal to diagonal"<<std::endl;}
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



TMatrixT<double> SBNchi::FillSystematicsFromXML(){
  return FillSystematicsFromXML(correlation_matrix_rootfile, correlation_matrix_name);
}


TMatrixT<double > SBNchi::FillSystematicsFromXML(std::string rootname, std::string matname){
  //Pretty much obsolete now, should fill directly really.
  std::cout<<"SBNchi::FillSystematicsFromXML || filling from "<<rootname<<std::endl;

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

  std::cout<<"SBNchi::FillSystematicsFromXML || loaded with dim : "<<temp2.GetNcols()<<" "<<temp2.GetNrows()<<std::endl;

  fm->Close();
  delete fm;

  if(temp2.IsSymmetric()){
    if(is_verbose)std::cout<<"Inputted fracCov covariance matrix is symmetric"<<std::endl;
  }else{
    std::cerr<<"ERROR: SBNchi::FillSystematicsFromXML, matrix_systematics input is not symmetric!"<<std::endl;
    //exit(EXIT_FAILURE);
  }

  return temp2;

}





TH2D* SBNchi::GetChiogram(){
  TH2D *tmp = new TH2D("chi-o-gram","chi-o-gram",num_bins_total_compressed,0, num_bins_total_compressed ,num_bins_total_compressed,0, num_bins_total_compressed);

  for(int i =0; i<num_bins_total_compressed; i++){
    for(int j =0; j<num_bins_total_compressed; j++){

      tmp->SetBinContent(i+1, j+1, vec_last_calculated_chi.at(i).at(j));
    }
  }

  return tmp;
}

int SBNchi::PrintMatricies(std::string tag){
  TFile* fout = new TFile(("SBNfit_collapsed_matrix_plots_"+tag+".root").c_str(),"recreate");
  fout->cd();


  gStyle->SetOptStat(0);

  TMatrixD full, frac, corr;
  this->FillCollapsedCovarianceMatrix(&full);
  this->FillCollapsedFractionalMatrix(&frac);
  this->FillCollapsedCorrelationMatrix(&corr);

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

  TH2D * h_chiogram = (TH2D*)this->GetChiogram();
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


int SBNchi::PerformCholoskyDecomposition(SBNspec *specin){
  TRandom3 * rangen = new TRandom3(0);
  specin->CalcFullVector();

  double tol = 1e-7;

  TMatrixD U  = matrix_fractional_covariance;

  for(int i =0; i<U.GetNcols(); i++)
    {
      for(int j =0; j<U.GetNrows(); j++)
	{
	  if(isnan(U(i,j)))
	    U(i,j) = 0;
	  else
	    U(i,j)=U(i,j)*specin->full_vector.at(i)*specin->full_vector.at(j);
	}
    }

  for(int i =0; i<U.GetNcols(); i++)
    {
      U(i,i) += specin->full_vector.at(i);
    }

  //First up, we have some problems with positive semi-definite and not positive definite
  TMatrixDEigen eigen (U);
  TVectorD eigen_values = eigen.GetEigenValuesRe();

  int n_t = U.GetNcols();

  for(int i=0; i< eigen_values.GetNoElements(); i++){
    if(eigen_values(i)<=0){
      if(fabs(eigen_values(i))< tol){
	std::cout<<"SBNchi::SampleCovariance\t|| cov has a very small, < "<<tol<<" , negative eigenvalue. Adding it back to diagonal of : "<<eigen_values(i)<<std::endl;

	for(int a =0; a<U.GetNcols(); a++){
	  U(a,a) += eigen_values(i);
	}

      }else{
	std::cout<<"SBNchi::SampleCovariance\t|| 0 or negative eigenvalues! error."<<std::endl;
	exit(EXIT_FAILURE);
      }
    }

    if(fabs(eigen_values(i))< tol){
      //SP_WARNING()<<"U has a very small, < "<<tol<<", eigenvalue which for some reason fails to decompose. Adding 1e9 to diagonal of U"<<std::endl;

      for(int a =0; a<U.GetNcols(); a++){
	U(a,a) += tol;
      }

    }	
  }


  //Seconndly attempt a Cholosky Decomposition
  TDecompChol * chol = new TDecompChol(U,0.1);
  bool worked = chol->Decompose();

  if(!worked){
    std::cout<<"SBNchi::SampleCovariance\t|| Cholosky Decomposition Failed."<<std::endl;
    exit(EXIT_FAILURE);

  }

  TMatrixT<double> upper_trian(n_t,n_t);
  matrix_lower_triangular.ResizeTo(n_t,n_t);
  upper_trian = chol->GetU();
  matrix_lower_triangular = upper_trian;
  matrix_lower_triangular.T();


  vec_matrix_lower_triangular.resize(n_t, std::vector<double>(n_t));
  for(int i=0; i< num_bins_total; i++){
    for(int j=0; j< num_bins_total; j++){
      vec_matrix_lower_triangular[i][j] = matrix_lower_triangular[i][j];
    }
  }

  //cholosky_performed = true;	
  return 0;
}


TH1D SBNchi::SampleCovarianceVaryInput(SBNspec *specin, int num_MC){ 
  std::vector<double>  tmp;
  return SampleCovarianceVaryInput(specin,num_MC,&tmp);
}

TH1D SBNchi::SampleCovarianceVaryInput(SBNspec *specin, int num_MC, std::vector<double> * chival){
  if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

  double** a_vec_matrix_lower_triangular = new double*[num_bins_total];
  double** a_vec_matrix_inverted = new double*[num_bins_total_compressed];
  
  for(int i=0; i < num_bins_total; i++){
    a_vec_matrix_lower_triangular[i] = new double[num_bins_total];
  }

  for(int i=0; i < num_bins_total_compressed; i++){
    a_vec_matrix_inverted[i] = new double[num_bins_total_compressed];
  }

  for(int i=0; i < num_bins_total; i++){
    for(int j=0; j < num_bins_total; j++){
      a_vec_matrix_lower_triangular[i][j] = vec_matrix_lower_triangular[i][j]; 
    }
  }
  
  for(int i=0; i< num_bins_total_compressed; i++){
    for(int j=0; j< num_bins_total_compressed; j++){
      a_vec_matrix_inverted[i][j] = vec_matrix_inverted[i][j]; 
    }
  }

   double *a_specin = new double[num_bins_total];
  double *a_corein = new double[num_bins_total_compressed];

  for(int i=0; i< num_bins_total; i++){
    a_specin[i] = specin->full_vector[i];
  }

  for(int i=0; i< num_bins_total_compressed; i++) {
    a_corein[i] = core_spectrum.collapsed_vector[i];
  }

  TRandom3 * rangen = new TRandom3(0);

  TH1D ans("","",150,0,150);
  ans.GetXaxis()->SetCanExtend(kTRUE);
  is_verbose = false;

  std::vector<double> vec_chis (num_MC, 0.0);

  double* a_vec_chis  = vec_chis.data();
  double* a_chival = chival->data();
  int num_chival = chival->size();
  
  int *nlower = new int[num_chival];
  for(int i=0; i< num_chival; i++){
    nlower[i]=0; 
  }

  std::vector < double > gaus_sample_v(num_bins_total), sampled_fullvector_v(num_bins_total);
  std::vector<double> collapsed_v(num_bins_total_compressed, 0.0);
	
  double gaus_sample[54];
  double sampled_fullvector[54];
  double collapsed[38];

  
#ifdef USE_GPU
  unsigned long long seed[num_MC];
  unsigned long long seq = 0ULL;
  unsigned long long offset = 0ULL;
  curandState_t state;
  
  for(int i=0; i<num_MC; ++i) {
    seed[i] = (int)rangen->Uniform(1e8);
  }
#endif

#ifdef USE_GPU
#pragma acc parallel loop num_gangs(2048) private(gaus_sample[:54],sampled_fullvector[:54],collapsed[:38],state) \
  copyin(this[0:1],							\
	 a_specin[:num_bins_total],					\
	 a_vec_matrix_lower_triangular[:num_bins_total][:num_bins_total],\ 
	 a_corein[:num_bins_total_compressed],				\
	 a_vec_matrix_inverted[:num_bins_total_compressed][:num_bins_total_compressed],	\
	 seed[0:num_MC],						\
	 a_chival[:num_chival],						\
	 this->a_num_bins[:num_channels],				\
	 this->a_num_subchannels[:num_channels])			\    
    copyout(a_vec_chis[:num_MC]) \
    copy(nlower[:num_chival])
#endif
    
    for(int i=0; i < num_MC;i++){

#ifdef USE_GPU
      unsigned long long seed_sd = seed[i];
      curand_init(seed_sd, seq, offset, &state);
      
      for(int a=0; a<num_bins_total; a++) {
	gaus_sample[a]= curand_normal(&state);
      }
#else
      for(int a=0; a<num_bins_total; a++) {
	gaus_sample[a]= (double)rangen->Gaus(0,1);
      }      
#endif

      for(int j = 0; j < num_bins_total; j++){
	sampled_fullvector[j] = a_specin[j];
	for(int k = 0; k < num_bins_total; k++){
	  sampled_fullvector[j] += a_vec_matrix_lower_triangular[j][k] * gaus_sample[k];
	}
	if(sampled_fullvector[j]<0) sampled_fullvector[j]=0.0;
       }
	
	
      this->CollapseVectorStandAlone(sampled_fullvector, collapsed);
	
      a_vec_chis[i] = this->CalcChi(a_vec_matrix_inverted, a_corein, collapsed);

     //Just to get some pvalues that were asked for.
     for(int j=0; j< num_chival; j++){
    #pragma acc atomic update
        if(a_vec_chis[i]>=a_chival[j]) nlower[j]++;
     }


    }


  is_verbose = true;
  
  
  for(int i=0; i<num_MC; i++){
    if (i<(int)1e3) 
      std::cout << "@i=" << a_vec_chis[i] << std::endl;
    ans.Fill(a_vec_chis[i]);
  }
  
  for(int n =0; n< num_chival; n++){
    chival->at(n) = nlower[n]/(double)num_MC;
  }
  
  
  
  delete[] a_corein;
  delete[] a_specin;
  delete[] nlower;
  
  for(int i=0; i < num_bins_total; i++){
    delete[] a_vec_matrix_lower_triangular[i];
  }

  for(int i=0; i < num_bins_total_compressed; i++){
    delete[] a_vec_matrix_inverted[i];  
  }
  
  delete[] a_vec_matrix_lower_triangular;
  delete[] a_vec_matrix_inverted;
    
  return ans;
}

int SBNchi::CollapseVectorStandAlone(std::vector<double> * full_vector, std::vector<double> *collapsed_vector){
  
  for(int im = 0; im < num_modes; im++){
    for(int id =0; id < num_detectors; id++){
      int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
      int out_edge = edge;
      int tmp_chan = 0;
      for(int ic = 0; ic < num_channels; ic++){
	int corner=edge;
						
	for(int j=0; j< num_bins[ic]; j++){

	  double tempval=0;

	  for(int sc = 0; sc < num_subchannels[ic]; sc++){
	    tempval += (*full_vector)[j+sc*num_bins[ic]+corner];
	    edge +=1;	//when your done with a channel, add on every bin you just summed
	  }
	  //we can size this vector beforehand and get rid of all push_back()

	  int collapsed_index = tmp_chan+out_edge;
	  (*collapsed_vector)[collapsed_index] = tempval;
	  tmp_chan++;
	}
      }
    }
  }


  return 0;
}

int SBNchi::CollapseVectorStandAlone(double* full_vector, double *collapsed_vector){

  //int tmp_num_bins[3] = {25,25,6};
  //int tmp_num_subchannels[3] = {2,1,1};
  
  for(int im = 0; im < num_modes; im++){
    for(int id =0; id < num_detectors; id++){
      int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
      int out_edge = edge;
      int chan = 0;
      for(int ic = 0; ic < num_channels; ic++){
	int corner=edge;

	for(int j=0; j< this->a_num_bins[ic]; j++){
	  //for(int j=0; j< tmp_num_bins[ic]; j++){

	  double tempval=0;

	  for(int sc = 0; sc < this->a_num_subchannels[ic]; sc++){
	    //for(int sc = 0; sc < tmp_num_subchannels[ic]; sc++){
	    tempval += full_vector[j+sc*this->a_num_bins[ic]+corner];
	    //tempval += full_vector[j+sc*tmp_num_bins[ic]+corner];
	    edge +=1;	//when your done with a channel, add on every bin you just summed
	  }
	  //we can size this vector beforehand and get rid of all push_back()

	  int collapsed_index = chan+out_edge;
	  collapsed_vector[collapsed_index] = tempval;
	  chan++;
	}
      }
    }
  }


  return 0;
}


SBNspec SBNchi::SampleCovariance(SBNspec *specin){
  if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

  int n_t = specin->full_vector.size();


  TVectorT<double> u(n_t);
  for(int i=0; i<n_t; i++){
    u(i) = specin->full_vector.at(i);
  }

  TRandom3 * rangen = new TRandom3(0);

  is_verbose = false;

  TVectorT<double> gaus_sample(n_t);
  TVectorT<double> multi_sample(n_t);
  for(int a=0; a<n_t; a++){
    gaus_sample(a) = rangen->Gaus(0,1);	
  }

  multi_sample = u + matrix_lower_triangular*gaus_sample;

  std::vector<double> sampled_fullvector(n_t,0.0);
  for(int j=0; j<n_t; j++){
    sampled_fullvector.at(j) = multi_sample(j);
  }
  SBNspec sampled_spectra(sampled_fullvector, specin->xmlname ,false);

  sampled_spectra.CollapseVector(); //this line important isnt it!


  is_verbose = true;



  return sampled_spectra;
}




TH1D SBNchi::SamplePoissonVaryInput(SBNspec *specin, int num_MC){ 
  std::vector<double>  tmp = {10};
  return SamplePoissonVaryInput(specin,num_MC,&tmp);
}
//This one varies the input comparative spectrum, and as sucn has  only to calculate the matrix_systematics once
TH1D SBNchi::SamplePoissonVaryInput(SBNspec *specin, int num_MC, std::vector<double> *chival){
  std::vector<int> nlower(chival->size(),0);

  TRandom3 *rangen = new TRandom3(0);

  TH1D ans("","",150,0,150);
  //So save the core one that we will sample for
  ans.GetXaxis()->SetCanExtend(kTRUE);
  is_verbose = false;
  for(int i=0; i < num_MC;i++){
    SBNspec tmp = *specin;
    tmp.ScalePoisson(rangen);
    tmp.CollapseVector(); //this line important isnt it!
    //tmp.PrintFullVector();

    double thischi = this->CalcChi(&tmp);
    ans.Fill(thischi);

    for(int j=0; j< chival->size(); j++){
      if(thischi>=chival->at(j)) nlower.at(j)++;
    }

    if(i%1000==0) std::cout<<"SBNchi::SamplePoissonVaryInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
  }
  for(int n =0; n< nlower.size(); n++){
    chival->at(n) = nlower.at(n)/(double)num_MC;
  }

  is_verbose = true;
  return ans;


}
/*
  std::vector<double> SBNchi::SampleCovarianceVaryInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
  if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

  int n_t = specin->full_vector.size();
  std::vector<int> nlower(chival.size(),0);

  TVectorT<double> u(n_t);
  for(int i=0; i<n_t; i++){
  u(i) = specin->full_vector.at(i);
  }

  TRandom3 * rangen = new TRandom3(0);


  TH1D ans("","",100,0,100);
  ans.GetXaxis()->SetCanExtend(kTRUE);
  is_verbose = false;
  for(int i=0; i < num_MC;i++){

  TVectorT<double> gaus_sample(n_t);
  TVectorT<double> multi_sample(n_t);
  for(int a=0; a<n_t; a++){
  gaus_sample(a) = rangen->Gaus(0,1);	
  }

  multi_sample = u + matrix_lower_triangular*gaus_sample;

  std::vector<double> sampled_fullvector(n_t,0.0);
  for(int i=0; i<n_t; i++){
  sampled_fullvector.at(i) = multi_sample(i);
  }
  SBNspec sampled_spectra(sampled_fullvector, specin->xmlname ,false);

  sampled_spectra.CollapseVector(); //this line important isnt it!

  double thischi = this->CalcChi(&sampled_spectra);
  ans.Fill(thischi);

  for(int j=0; j< chival.size(); j++){
  if(thischi>=chival.at(j)) nlower.at(j)++;
  }


  if(i%1000==0) std::cout<<"SBNchi::SampleCovarianceVaryInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
  }
  is_verbose = true;

  std::vector<double> pval;
  for(auto n: nlower){
  pval.push_back(n/(double)num_MC);

  }

  return pval;

  }


  //This one varies the input comparative spectrum, and as sucn has  only to calculate the matrix_systematics once
  std::vector<double> SBNchi::SamplePoissonVaryInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
  std::vector<int> nlower(chival.size(),0);

  TRandom3 *rangen = new TRandom3(0);

  TH1D ans("","",100,0,100);
  //So save the core one that we will sample for
  ans.GetXaxis()->SetCanExtend(kTRUE);
  is_verbose = false;
  for(int i=0; i < num_MC;i++){

  SBNspec tmp = *specin;
  tmp.ScalePoisson(rangen);
  tmp.CollapseVector(); //this line important isnt it!
  //tmp.PrintFullVector();

  double thischi = this->CalcChi(&tmp);
  ans.Fill(thischi);

  for(int j=0; j< chival.size(); j++){
  if(thischi>=chival.at(j)) nlower.at(j)++;
  }

  if(i%1000==0) std::cout<<"SBNchi::SamplePoissonVaryInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
  }
  std::vector<double> pval;
  for(auto n: nlower){
  pval.push_back(n/(double)num_MC);

  }

  is_verbose = true;
  return pval;


  }
*/

//This one varies the core spectrum, and as sucn has to recalculate the matrix_systematics each stem
TH1D SBNchi::SamplePoissonVaryCore(SBNspec *specin, int num_MC){
  double center = this->CalcChi(specin);
  int nlower=0;

  TRandom3 *rangen = new TRandom3(0);

  TH1D ans("MCans","MCans",100,center-100,center+200);
  //So save the core one that we will sample for
  SBNspec core  = core_spectrum;

  is_verbose = false;
  for(int i=0; i<num_MC;i++){

    SBNspec tmp = core;
    tmp.ScalePoisson(rangen);
    this->ReloadCoreSpectrum(&tmp);
    double thischi = this->CalcChi(specin);
    ans.Fill(thischi);
    if(thischi<=center)nlower++;

    if(i%1000==0) std::cout<<"SBNchi::SamplePoissonVaryCore(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
  }
  std::cout<<"pval: "<<nlower/(double)num_MC<<std::endl;

  is_verbose = true;
  return ans;


}
