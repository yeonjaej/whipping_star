#include "SBNcovariance.h"
#include <chrono>

using namespace sbn;

SBNcovariance::SBNcovariance(std::string xmlname) : SBNconfig(xmlname) {

  universes_used = 0;
  tolerence_positivesemi = 1e-5;
  is_small_negative_eigenvalue = false;
  abnormally_large_weight = 1e3;


  std::map<std::string, int> parameter_sims;

  //Initialise the central value SBNspec.
  SBNspec tm(xmlname,-1,false);
  spec_central_value = tm;

  int num_files = multisim_file.size();

  for(auto &fn: multisim_file){
    files.push_back(new TFile(fn.c_str()));
  }

  for(int i=0; i<multisim_name.size(); i++){
    trees.push_back((TTree*)files.at(i)->Get(multisim_name.at(i).c_str()) );
  }

  for(int i=0; i<multisim_file.size(); i++){
    if( multisim_file_friend_treename_map.count(multisim_file.at(i))>0){
      for(int k=0; k< multisim_file_friend_treename_map.at(multisim_file.at(i)).size(); k++){

	std::string treefriendname = (multisim_file_friend_treename_map.at(multisim_file.at(i))).at(k);
	std::string treefriendfile = (multisim_file_friend_map.at(multisim_file.at(i))).at(k);

	std::cout<<"SBNcovariance::SBNcovariance\t|| Adding a friend tree  "<< treefriendfile<<" to file "<<multisim_file.at(i)<<std::endl;

	trees.at(i)->AddFriend( treefriendname.c_str()   ,  treefriendfile.c_str()   );
      }
    }
  }

  std::vector<int> nentries(multisim_name.size(),0);
  for(int i=0; i< multisim_name.size(); i++){
    nentries[i]= trees[i]->GetEntries();
  }

  f_weights = new std::vector<std::map<std::string, std::vector<double>>* >(num_files,0);


  for(int i=0; i< num_files; i++){
    delete f_weights->at(i);	f_weights->at(i) = 0;
    trees.at(i)->SetBranchAddress("weights", &f_weights->at(i) );
    delete f_weights->at(i);	f_weights->at(i) = 0;
    for(int k=0; k<branch_variables.at(i).size(); k++){
      trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->name.c_str(), branch_variables.at(i).at(k)->GetValue() );
    }
  }

  variations.clear();

  //This bit will calculate how many "multisims" the file has. if ALL default is the inputted xml value
  int good_event = 0;

  std::vector<int> used_multisims(num_files,0);
  for(int j = 0;j<num_files;j++){
    delete f_weights->at(j);
    f_weights->at(j)=0;

    trees.at(j)->GetEntry(good_event);
    for(std::map<std::string, std::vector<double> >::iterator  it = f_weights->at(j)->begin(); it != f_weights->at(j)->end(); ++it)
      {
	if(it->first == "bnbcorrection_FluxHist") continue;
	used_multisims.at(j) += it->second.size();

	std::cout<<"SBNcovariance::SBNcovariance\t|| "<<it->first<<" has "<<it->second.size()<<" multisims in file "<<j<<std::endl;
	variations.push_back(it->first);
      }
    delete f_weights->at(j);
    f_weights->at(j)=0;
  }
  //Now remove all, duplicates!
  sort( variations.begin(), variations.end() );
  variations.erase( unique( variations.begin(), variations.end() ), variations.end() );

  //make a map and start filling, before filling find if already in map, if it is check size.
  std::cout<<"SBNcovariance::SBNcovariance\t|| We have "<<variations.size()<<" unique variations: "<<std::endl;

  num_universes_per_variation.clear();

  for(auto &v: variations){
    std::cout<<"SBNcovariance::SBNcovariance\t|| "<<v<<std::endl;
    trees.at(0)->GetEntry(good_event);
    int thissize = f_weights->at(0)->at(v).size();


    for(int p=0;p<thissize; p++){
      map_universe_to_var[num_universes_per_variation.size()] = v;
      num_universes_per_variation.push_back(thissize);
    }
  }

  for(int i=1; i<num_files; i++){
    //std::cout<<"File: "<<i-1<<" has "<<used_multisims.at(i-1)<<" multisims"<<std::endl;
    std::cout<<"SBNcovariance::SBNcovariance\t|| File: "<<i<<" has "<<used_multisims.at(i)<<" multisims"<<std::endl;
    if( used_multisims.at(i)!= used_multisims.at(i-1)){
      std::cerr<<"SBNcovariance::SBNcovariance\t|| Warning, number of universes for are different between files"<<std::endl;
      std::cerr<<"SBNcovariance::SBNcovariance\t|| The missing universes are Set to weights of 1. Make sure this is what you want!"<<std::endl;
      for(int j=0; j< num_files; j++){
	if(universes_used < used_multisims.at(j)) universes_used = used_multisims.at(j);
	std::cerr<<"SBNcovariance::SBNmultisom\t|| File "<<j<<" multisims: "<<used_multisims.at(j)<<std::endl;
      }
      //exit(EXIT_FAILURE);
    }
    universes_used = used_multisims.at(0);
  }
  if(num_files ==1){
    universes_used = used_multisims.front();
  }








  std::cout<<"SBNcovariance::SBNcovariance\t|| -------------------------------------------------------------\n";
  std::cout<<"SBNcovariance::SBNcovariance\t|| Initilizing "<<universes_used<<" universes."<<std::endl;
  std::cout<<"SBNcovariance::SBNcovariance\t|| -------------------------------------------------------------\n";

  std::vector<double> base_vec (spec_central_value.num_bins_total,0.0);


  std::cout<<"SBNcovariance::SBNcovariance\t|| Full concatanated vector has : "<<spec_central_value.num_bins_total<<std::endl;
  for(int m=0; m<universes_used; m++){
    multi_vecspec.push_back(base_vec);
  }

  std::cout<<"SBNcovariance::SBNcovariance\t|| multi_vecspec now initilized of size :"<<multi_vecspec.size()<<std::endl;

  for(int j=0;j<num_files;j++){

    delete f_weights->at(j);
    f_weights->at(j)=0;



    for(int i=0; i< std::min(  multisim_maxevents[j]  ,nentries[j]); i++){
      trees.at(j)->GetEntry(i);
      std::map<std::string, std::vector<double>> * thisfWeight = f_weights->at(j);

      if(i%2500==0) std::cout<<"SBNcovariance::SBNcovariance\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<multisim_file[j]<<std::endl;

      std::vector<double> weights;
      std::vector<int> vec_universes;
      double global_weight = 1;

      global_weight = global_weight*multisim_scale.at(j);


      if(thisfWeight->count("bnbcorrection_FluxHist")>0){
	global_weight = global_weight*thisfWeight->at("bnbcorrection_FluxHist").front();
      }

      if(std::isinf(global_weight) || global_weight != global_weight){
	std::cout<<"SBNcovariance::SBNcovariance\t|| ERROR  error @ "<<i<<" in File "<<multisim_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
	exit(EXIT_FAILURE);
      }


      if( this->EventSelection(j) ){

	//else if(parameter_names.at(j)[0]=="ALL") //Currently just do all variations
	{
	  //Loop over all variations!
	  for(std::string &var : variations){

	    //check if variation is in this file, if it isn't: then just push back 1's of appropiate number to keep universes consistent
	    if(thisfWeight->count(var) <=0 ){
	      for(int g=0; g<f_weights->at(0)->at(var).size(); g++){
		weights.push_back(global_weight);
	      }
	      continue;
	    }



	    for(double &wei: thisfWeight->at(var)){
	      bool is_inf = std::isinf(wei);
	      bool is_nan = wei!=wei;

	      if(is_inf || is_nan){
		std::cout<<"SBNcovariance::SBNcovariance\t|| ERROR! Killing :: event # "<<i<<" in File "<<multisim_file.at(j)<<" weight: "<<wei<<" global bnb: "<<global_weight<<" in "<<var<<std::endl;
		exit(EXIT_FAILURE);
	      }

	      if(wei > abnormally_large_weight){
		std::cout<<"ATTENTION! SBNcovariance::SBNcovariance\t|| HUGE weight: "<<wei<<" at "<<var<<" event "<<i<<" file "<<j<<std::endl;
	      }
	      weights.push_back(wei*global_weight);
	    }


	  }//end of variations

	  //So the size of weights must equal global universes ya?
	  if(weights.size() != universes_used || universes_used != multi_vecspec.size() || universes_used != num_universes_per_variation.size()){
	    std::cout<<"SBNcovariance::SBNcovariance\t|| ERROR "<<std::endl;
	    std::cout<<"weights.size() "<<weights.size()<<std::endl;
	    std::cout<<"universes_used "<<universes_used<<std::endl;
	    std::cout<<"multi_vecspec.size() "<<multi_vecspec.size()<<std::endl;
	    std::cout<<"num_universes_per_variation.size() "<<num_universes_per_variation.size()<<std::endl;
	    exit(EXIT_FAILURE);
	  }



	}//end of "ALL" option
	//Begininning of single parameter optifons
	/*		else {

			std::vector<double> this_param_weights = thisfWeight->at(parameter_names.at(j)[0]);



			for(double &wei : this_param_weights){

			if(std::isnan(wei) || wei != wei){
			std::cout<<"Killing. event # "<<i<<" in File "<<multisim_file.at(j)<<" weight: "<<wei<<" global bnb: "<<global_weight<<" in "<<parameter_names.at(j)[0]<<std::endl;
			exit(EXIT_FAILURE);

			weights.push_back(wei*global_weight);

			}//end of this weight loop
			}//end of single parameter scan



			delete f_weights->at(j);f_weights->at(j) = 0;


			//So the size of weights must equal global universes ya?
			if(weights.size() != universes_used || universes_used != multi_vecspec.size()){
			std::cerr<<"num_sim (weights.size()) != universes_used "<<weights.size()<<" "<<universes_used<<std::endl;
			std::cout<<"num_sim (weights.size()) != universes_used "<<weights.size()<<" "<<universes_used<<std::endl;
			//continue;
			exit(EXIT_FAILURE);
			}
			}
	*/

	for(int t=0; t<branch_variables[j].size();t++){
	  //Need the histogram index, the value, the global bin...
	  int ih = spec_central_value.map_hist.at(branch_variables[j][t]->associated_hist);
	  double reco_var = *(static_cast<double*>(branch_variables[j][t]->GetValue()));
	  int reco_bin = spec_central_value.GetGlobalBinNumber(reco_var,ih);
	  spec_central_value.hist[ih].Fill(reco_var,global_weight);

	  for(int m=0; m< weights.size(); m++){
	    if(reco_bin>=0)  multi_vecspec[m][reco_bin]   +=  weights[m];

	    //important check. failure mode
	    if(weights[m]!=weights[m] || std::isinf(weights[m]) ){
	      std::cout<<"SBNcovariance::SBNcovariance\t|| ERROR weight has a value of: "<<weights[m]<<". So I am killing all. on Dim: "<<m<<" global_eright is "<<global_weight<<std::endl;
	      exit(EXIT_FAILURE);
	    }

	  }
	}
      }
    } //end of entry loop
  }//end of file loop

  delete f_weights;

  /***************************************************************
   *		Now some clean-up and Writing
   * ************************************************************/

  for(auto &f: files){
    f->Close();
  }

}
/***************************************************************
 *		Some virtual functions for selection and histogram filling
 * ************************************************************/

bool SBNcovariance::EventSelection(int which_file){
  //from here have access to vars_i  and vars_d  to make a selection
  return true;
}

int SBNcovariance::FillHistograms(int file, int uni, double wei){
  //double en = vars_d.at(file)[0];
  //multi_sbnspec.at(uni).hist.at(file).Fill(en, wei);
  return 0;
}

/***************************************************************
 *		And then form a covariance matrix (well 3)
 * ************************************************************/


int SBNcovariance::FormCovarianceMatrix(std::string tag){

  std::cout<<"SBNcovariance::FormCovariancematrix\t|| Starting.\n";
  full_covariance.ResizeTo(num_bins_total, num_bins_total);
  frac_covariance.ResizeTo(num_bins_total, num_bins_total);
  full_correlation.ResizeTo(num_bins_total, num_bins_total);


  for(auto &v: variations){

    map_var_to_matrix[v] = vec_full_covariance.size();

    vec_full_covariance.push_back(full_covariance);
    vec_frac_covariance.push_back(frac_covariance);
    vec_full_correlation.push_back(full_correlation);

    vec_full_covariance.back().Zero();
    vec_frac_covariance.back().Zero();
    vec_full_correlation.back().Zero();


  }


  for(auto &h: multi_sbnspec){
    h.CalcFullVector();
  }


  spec_central_value.CalcFullVector();
  std::vector<double> CV = spec_central_value.full_vector;

  std::cout<<"SBNcovariance::FormCovariancematrix\t|| Begining to form the "<<num_bins_total<<"X"<<num_bins_total<<" covariance matrix.\n";
  auto a00 = std::chrono::high_resolution_clock::now();

  for(int i=0; i<num_bins_total; i++){
    for(int j=0; j<num_bins_total; j++){

      full_covariance(i,j)=0;
      for(int v=0; v< variations.size();v++){
	vec_full_covariance.at(v)(i,j)=0.0;
      }	


      for(int m=0; m < universes_used; m++){

	std::string var = map_universe_to_var.at(m);
	int which_matrix = map_var_to_matrix.at(var);

	full_covariance(i,j) += 1.0/((double)num_universes_per_variation[m])*(CV[i]-multi_vecspec[m][i])*(CV[j]-multi_vecspec[m][j]);
	vec_full_covariance.at(which_matrix)(i,j) += 1.0/((double)num_universes_per_variation[m])*(CV[i]-multi_vecspec[m][i])*(CV[j]-multi_vecspec[m][j]);

	//if(false) std::cout<<"BinTest: "<<i<<" "<<j<<" universe: "<<m<<" @ "<<var<<" total: "<<num_universes_per_variation.at(m)<<" which_matrix "<<which_matrix<<std::endl;

	if(full_covariance(i,j)!=full_covariance(i,j)){
	  std::cout<<"SBNcovariance::FormCovariancematrix\t|| ERROR: nan : at (i,j):  "<<i<<" "<<j<<" fullcov: "<<full_covariance(i,j)<<" multi hist sise "<<multi_vecspec.size()<<" CV: "<<CV[i]<<" "<<CV[j]<<" multihisg "<<multi_vecspec[m][i]<<" "<<multi_vecspec[m][j]<<" on dim : "<<m<<std::endl;
	  exit(EXIT_FAILURE);
	}

      }

    }
  }

  auto b00 = std::chrono::high_resolution_clock::now();
  auto a00b00 = std::chrono::duration_cast<std::chrono::microseconds>(b00-a00);
  std::cout << "Time: " << a00b00.count() << "µs" << std::endl;


  std::cout<<"SBNcovariance::FormCovariancematrix\t|| Now calculating fractional covariance and correlation matrix from full covariance."<<std::endl;
  for(int i=0; i<num_bins_total; i++){
    for(int j=0; j<num_bins_total; j++){
      frac_covariance(i,j) = full_covariance(i,j)/(spec_central_value.full_vector[i]*spec_central_value.full_vector[j]) ;
      full_correlation(i,j)= full_covariance(i,j)/(sqrt(full_covariance(i,i))*sqrt(full_covariance(j,j)));

      for(int m=0; m< vec_full_covariance.size(); m++){
	vec_frac_covariance.at(m)(i,j) = vec_full_covariance.at(m)(i,j)/(spec_central_value.full_vector[i]*spec_central_value.full_vector[j]) ;
	vec_full_correlation.at(m)(i,j)= vec_full_covariance.at(m)(i,j)/(sqrt(vec_full_covariance.at(m)(i,i))*sqrt(vec_full_covariance.at(m)(j,j)));
      }

    }
  }


  /************************************************************
   *			Saving to file				    *
   * *********************************************************/
  TFile *fout=new TFile((tag+".SBNcovar.root").c_str(),"RECREATE");
  fout->cd();
  full_covariance.Write(("full_covariance_"+tag).c_str(),TObject::kWriteDelete);
  frac_covariance.Write(("frac_covariance_"+tag).c_str(),TObject::kWriteDelete);
  full_correlation.Write(("full_correlation_"+tag).c_str(),TObject::kWriteDelete);


  TDirectory *individualDir = fout->GetDirectory("individualDir"); 
  if (!individualDir) { 
    individualDir = fout->mkdir("individualDir");       
  }
  fout->cd(); 
  individualDir->cd();


  for(int m=0; m< variations.size();m++){
    vec_full_correlation.at(m).Write( (variations.at(m)+"_full_correlation_"+tag).c_str(), TObject::kWriteDelete);
    vec_frac_covariance.at(m).Write( (variations.at(m)+"_frac_covariance_"+tag).c_str(), TObject::kWriteDelete);
    vec_full_covariance.at(m).Write( (variations.at(m)+"_full_covariance_"+tag).c_str(), TObject::kWriteDelete);
  }



  std::vector<TH2D> h2_corr;
  std::vector<TH2D> h2_cov;
  std::vector<TH2D> h2_fcov;

  /*
    for(int m=0; m< variations.size();m++){
    //		vec_frac_covariance.at(m).Write((variations.at(m)+"_frac_covariance_"+tag).c_str() ,TObject::kWriteDelete);
    //		vec_full_covariance.at(m).Write((variations.at(m)+"_full_covariance_"+tag).c_str() ,TObject::kWriteDelete);
    //		vec_full_correlation.at(m).Write((variations.at(m)+"_full_correlation_"+tag).c_str() ,TObject::kWriteDelete);

    h2_corr.push_back(TH2D(vec_full_correlation.at(m)));
    h2_cov.push_back(TH2D(vec_full_covariance.at(m)));
    h2_fcov.push_back(TH2D(vec_frac_covariance.at(m)));

    h2_fcov.back().SetName((variations.at(m)+"_frac_covariance_"+tag).c_str());
    h2_corr.back().SetName((variations.at(m)+"_full_correlation_"+tag).c_str());
    h2_cov.back().SetName((variations.at(m)+"_full_covariance_"+tag).c_str());

    h2_fcov.back().Write();
    h2_cov.back().Write();
    h2_corr.back().Write();

    }
  */

  fout->Close();
  //and also writeout
  spec_central_value.WriteOut(tag);



  this->qualityTesting();
  return 0;
}


int SBNcovariance::qualityTesting(){



  /************************************************************
   *		Quality Testing Suite			    *
   * *********************************************************/
  std::cout<<"SBNcovariance::qualityTesting\t||-----------------------------------------------------\n";
  std::cout<<"SBNcovariance::qualityTesting\t||----------------Quality Testing Suite"<<std::endl;
  std::cout<<"SBNcovariance::qualityTesting\t||-----------------------------------------------------\n";


  std::cout<<"SBNcovariance::qualityTesting\t|| Checking if generated matrix is indeed a valid covariance matrix.\n";
  std::cout<<"SBNcovariance::qualityTesting\t|| First checking if matrix is symmetric.\n";

  double max_sym_violation = 0;
  for(int i=0; i<num_bins_total; i++){
    for(int j=0; j<num_bins_total; j++){
      double tnp = fabs((full_covariance(j,i)-full_covariance(i,j))/(full_covariance(j,i)+full_covariance(i,j)));
      if(tnp > max_sym_violation) max_sym_violation = tnp;
    }
  }



  if(max_sym_violation < 1e-13){
    std::cout<<"SBNcovariance::qualityTesting\t|| PASS: Generated covariance matrix is symmetric"<<std::endl;
  }else{
    std::cout<<"SBNcovariance::qualityTesting\t||  ERROR result is not symmetric! "<<max_sym_violation<<std::endl;
    for(int i=0; i<num_bins_total; i++){
      for(int j=0; j<num_bins_total; j++){
	double tnp = fabs((full_covariance(j,i)-full_covariance(i,j))/(full_covariance(j,i)+full_covariance(i,j)));
	if(full_covariance(i,j) != full_covariance(j,i)) std::cout<<i<<" "<<j<<" "<<full_covariance(i,j)<<" "<<full_covariance(j,i)<<" "<<tnp<<std::endl;
      }
    }
    std::cout<<"SBNcovariance::qualityTesting\t||  ERROR result is not symmetric!"<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout<<"SBNcovariance::qualityTesting\t|| Checking if generated matrix is positive semi-definite by looking at eigenvalues.\n";
  //if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
  TMatrixDEigen eigen (full_covariance);
  TVectorD eigen_values = eigen.GetEigenValuesRe();


  for(int i=0; i< eigen_values.GetNoElements(); i++){
    if(eigen_values(i)<0){
      is_small_negative_eigenvalue = true;
      if(fabs(eigen_values(i))> tolerence_positivesemi ){
	std::cout<<"SBNcovariance::qualityTesting\t|| ERROR contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
	exit(EXIT_FAILURE);
      }
    }
  }


  if(is_small_negative_eigenvalue){
    std::cout<<"SBNcovariance::qualityTesting\t|| PASS: Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
  }else{
    std::cout<<"SBNcovariance::qualityTesting\t|| PASS: Generated covariance matrix is positive semi-definite."<<std::endl;
  }
  std::cout<<"SBNcovariance::qualityTesting\t|| Congratulations, matrix is indeed a valid covariance matrix.\n";

  return 0;
}

int SBNcovariance::PrintVariations(std::string tag){
  TFile *fout = new TFile(("SBNfit_variation_plots_"+tag+".root").c_str(),"recreate");
  fout->cd();

  std::cout<<"SBNcovariance::PrintVariations\t|| Starting to Print all variations, this can take a little. "<<std::endl;

  std::vector<TDirectory*> vec_dir;

  std::vector<std::vector<TCanvas*>> vec_canvas;

  for(auto &v: variations){

    //std::cout<<"SBNcovariance::PrintVariations\t|| Preparing directory and canvases for variation: "<<v<<std::endl;
    fout->cd();
    vec_dir.push_back( fout->GetDirectory(v.c_str()));
    if (!vec_dir.back()) { 
      vec_dir.back() = fout->mkdir(v.c_str());       
    }
    vec_dir.back()->cd();

    std::vector<TCanvas *> tmpc;

    for(int i=0; i< spec_central_value.hist.size(); i++){
      tmpc.push_back(new TCanvas((fullnames.at(i)+"||"+v).c_str()));
      tmpc.back()->cd();
      TH1D * temp_cv_spec = (TH1D*)spec_central_value.hist.at(i).Clone((std::to_string(i)+v).c_str());
      temp_cv_spec->Scale(1,"width");

      tmpc.back()->cd();
      double maxval = temp_cv_spec->GetMaximum();
      if(maxval > 0) 	temp_cv_spec->SetMaximum(maxval*1.45);
      temp_cv_spec->SetStats(false);
      temp_cv_spec->SetLineColor(kBlack);
      temp_cv_spec->SetLineWidth(2);
      temp_cv_spec->GetXaxis()->SetTitle(fullnames.at(i).c_str());
      temp_cv_spec->GetYaxis()->SetTitle("Events/unit");
      temp_cv_spec->SetTitle((v + " || " +fullnames.at(i)).c_str());
      temp_cv_spec->DrawCopy("hist");

      delete temp_cv_spec;
    }
    vec_canvas.push_back(tmpc);	
  }

  std::cout<<"SBNcovariance::PrintVariations\t|| Starting universe loop. "<<std::endl;
  TRandom3 *rangen = new TRandom3(0);
  for(int m=0; m < universes_used; m++){
    std::string var = map_universe_to_var.at(m);
    int which_matrix = map_var_to_matrix.at(var);

    vec_dir.at(which_matrix)->cd();

    SBNspec temp_spec(multi_vecspec.at(m), xmlname,false);

    
    for(int i=0; i< temp_spec.hist.size(); i++){
      vec_canvas.at(which_matrix).at(i)->cd();
      temp_spec.hist.at(i).Scale(1,"width");
      temp_spec.hist.at(i).SetLineColor((int)rangen->Uniform(300,1000));	
      temp_spec.hist.at(i).DrawCopy("same hist");
    }	



  }//end universe loop

  std::cout<<"SBNcovariance::PrintVariations\t|| Finished. Just tidying up and writing TCanvas. "<<std::endl;


  for(int v =0; v< variations.size(); v++){
    fout->cd();
    vec_dir.at(v)->cd();

    for(int i=0; i< spec_central_value.hist.size(); i++){
      vec_canvas.at(v).at(i)->cd();
      TH1D * temp_cv_spec = (TH1D*)spec_central_value.hist.at(i).Clone((std::to_string(i)+variations.at(v)+"tmp2").c_str());
      temp_cv_spec->Scale(1,"width");
      temp_cv_spec->SetLineColor(kBlack);
      temp_cv_spec->SetMarkerStyle(34);
      temp_cv_spec->SetLineWidth(2);
      temp_cv_spec->DrawCopy("same hist p");
		
      vec_canvas.at(v).at(i)->Write();
      delete temp_cv_spec;
    }
  }


  fout->Close();
  return 0;
}


int SBNcovariance::PrintMatricies(std::string tag){
  TFile* fout = new TFile(("SBNfit_covariance_plots_"+tag+".root").c_str(),"recreate");
  fout->cd();


  gStyle->SetOptStat(0);

  //correlation
  TH2D h2_corr(full_correlation);
  h2_corr.SetName("corr");
  //h2_corr.Write();
  TCanvas *c_corr = new TCanvas("full correlation matrix");
  TPad *p_corr = (TPad*)c_corr->cd();
  c_corr->SetFixedAspectRatio();
  h2_corr.Draw("colz");
  h2_corr.SetTitle("Full correlation matrix");
  h2_corr.GetXaxis()->SetTitle("Reco Bin i");
  h2_corr.GetYaxis()->SetTitle("Reco Bin j");

  c_corr->SetFrameFillColor(kWhite);
  c_corr->SetFillColor(kWhite);
  p_corr->SetFillColor(kWhite);

  c_corr->SetRightMargin(0.150);
  int use_corr =0;
  for(int im =0; im<num_modes; im++){
    for(int id =0; id<num_detectors; id++){
      for(int ic = 0; ic < num_channels; ic++){
	for(int isc = 0; isc < num_subchannels.at(ic)-1; isc++){
	  TLine *lscv = new TLine(0, num_bins.at(ic)+use_corr, num_bins_total, num_bins.at(ic)+use_corr);
	  TLine *lsch = new TLine(num_bins.at(ic)+use_corr,0, num_bins.at(ic)+use_corr, num_bins_total);
	  lscv->SetLineWidth(2);
	  lsch->SetLineWidth(2);
	  lscv->SetLineColor(kRed);
	  lsch->SetLineColor(kRed);
	  use_corr+=num_bins.at(ic);
	  lscv->Draw();
	  lsch->Draw();

	}
	TLine *lv = new TLine(0, num_bins.at(ic)+use_corr, num_bins_total, num_bins.at(ic)+use_corr);
	TLine *lh = new TLine(num_bins.at(ic)+use_corr,0, num_bins.at(ic)+use_corr, num_bins_total);
	lv->SetLineWidth(2);
	lh->SetLineWidth(2);
	use_corr+=num_bins.at(ic);
	lv->Draw();
	lh->Draw();

      }
    }
  }
  c_corr->Write();

  //full covariance
  TH2D h2_full(full_covariance);
  h2_full.SetName("full");
  //h2_full.Write();
  TCanvas *c_full = new TCanvas("full covariance matrix");
  c_full->cd();
  c_full->SetFixedAspectRatio();
  h2_full.Draw("colz");
  h2_full.SetTitle("Full covariance matrix");
  h2_full.GetXaxis()->SetTitle("Reco Bin i");
  h2_full.GetYaxis()->SetTitle("Reco Bin j");

  c_full->SetRightMargin(0.150);
  int use_full =0;
  for(int im =0; im<num_modes; im++){
    for(int id =0; id<num_detectors; id++){
      for(int ic = 0; ic < num_channels; ic++){
	for(int isc = 0; isc < num_subchannels.at(ic)-1; isc++){
	  TLine *lscv = new TLine(0, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
	  TLine *lsch = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total);
	  lscv->SetLineWidth(2);
	  lsch->SetLineWidth(2);
	  lscv->SetLineColor(kRed);
	  lsch->SetLineColor(kRed);
	  use_full+=num_bins.at(ic);
	  lscv->Draw();
	  lsch->Draw();

	}
	TLine *lv = new TLine(0, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
	TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total);
	lv->SetLineWidth(2);
	lh->SetLineWidth(2);
	use_full+=num_bins.at(ic);
	lv->Draw();
	lh->Draw();

      }
    }
  }
  c_full->Write();

  //fracelation
  TH2D h2_frac(frac_covariance);
  h2_frac.SetName("frac");
  //h2_frac.Write();
  TCanvas *c_frac = new TCanvas("full fractional covariance matrix");
  c_frac->cd();
  c_frac->SetFixedAspectRatio();
  h2_frac.Draw("colz");
  h2_frac.SetTitle("Full fractional covariance matrix");
  h2_frac.GetXaxis()->SetTitle("Reco Bin i");
  h2_frac.GetYaxis()->SetTitle("Reco Bin j");
  c_frac->SetRightMargin(0.150);

  int use_frac =0;
  for(int im =0; im<num_modes; im++){
    for(int id =0; id<num_detectors; id++){
      for(int ic = 0; ic < num_channels; ic++){
	for(int isc = 0; isc < num_subchannels.at(ic)-1; isc++){
	  TLine *lscv = new TLine(0, num_bins.at(ic)+use_frac, num_bins_total, num_bins.at(ic)+use_frac);
	  TLine *lsch = new TLine(num_bins.at(ic)+use_frac,0, num_bins.at(ic)+use_frac, num_bins_total);
	  lscv->SetLineWidth(2);
	  lsch->SetLineWidth(2);
	  lscv->SetLineColor(kRed);
	  lsch->SetLineColor(kRed);
	  use_frac+=num_bins.at(ic);
	  lscv->Draw();
	  lsch->Draw();

	}
	TLine *lv = new TLine(0, num_bins.at(ic)+use_frac, num_bins_total, num_bins.at(ic)+use_frac);
	TLine *lh = new TLine(num_bins.at(ic)+use_frac,0, num_bins.at(ic)+use_frac, num_bins_total);
	lv->SetLineWidth(2);
	lh->SetLineWidth(2);
	use_frac+=num_bins.at(ic);
	lv->Draw();
	lh->Draw();

      }
    }
  }
  c_frac->Write();



  //Print the collapsed matricies too: Need to fudge this a bit
  SBNchi collapse_chi(xmlname);

  TMatrixT<double > coll_correlation(num_bins_total_compressed,num_bins_total_compressed);
  TMatrixT<double > coll_frac_covariance(num_bins_total_compressed,num_bins_total_compressed);
  TMatrixT<double > coll_covariance(num_bins_total_compressed,num_bins_total_compressed);

  collapse_chi.CollapseModes(full_covariance, coll_covariance);

  for(int i=0; i<num_bins_total_compressed; i++){
    for(int j=0; j<num_bins_total_compressed; j++){
      coll_frac_covariance(i,j) = coll_covariance(i,j)/(spec_central_value.full_vector.at(i)*spec_central_value.full_vector.at(j)) ;
      coll_correlation(i,j)= coll_covariance(i,j)/(sqrt(coll_covariance(i,i))*sqrt(coll_covariance(j,j)));
    }
  }

  TH2D h2_coll_corr(coll_correlation);
  h2_coll_corr.SetName("coll_corr");
  TCanvas *c_coll_corr = new TCanvas("collapsed correlation matrix");
  c_coll_corr->cd();
  c_coll_corr->SetFixedAspectRatio();
  h2_coll_corr.Draw("colz");
  h2_coll_corr.SetTitle("Collapsed correlation matrix");
  h2_coll_corr.GetXaxis()->SetTitle("Reco Bin i");
  h2_coll_corr.GetYaxis()->SetTitle("Reco Bin j");

  c_coll_corr->SetRightMargin(0.150);

  int use_coll_corr =0;
  for(int im =0; im<num_modes; im++){
    for(int id =0; id<num_detectors; id++){
      for(int ic = 0; ic < num_channels; ic++){
	TLine *lv = new TLine(0, num_bins.at(ic)+use_coll_corr, num_bins_total_compressed, num_bins.at(ic)+use_coll_corr);
	TLine *lh = new TLine(num_bins.at(ic)+use_coll_corr,0, num_bins.at(ic)+use_coll_corr, num_bins_total_compressed);
	lv->SetLineWidth(1.5);
	lh->SetLineWidth(1.5);
	use_coll_corr+=num_bins.at(ic);
	lv->Draw();
	lh->Draw();
      }
    }
  }
  c_coll_corr->Write();


  TH2D h2_coll_frac(coll_frac_covariance);
  h2_coll_frac.SetName("coll_frac");
  TCanvas *c_coll_frac = new TCanvas("collapsed fractional covariance matrix");
  c_coll_frac->cd();
  c_coll_frac->SetFixedAspectRatio();
  h2_coll_frac.Draw("colz");
  h2_coll_frac.SetTitle("Collapsed fractional covariance matrix");
  h2_coll_frac.GetXaxis()->SetTitle("Reco Bin i");
  h2_coll_frac.GetYaxis()->SetTitle("Reco Bin j");

  c_coll_frac->SetRightMargin(0.150);

  int use_coll_frac =0;
  for(int im =0; im<num_modes; im++){
    for(int id =0; id<num_detectors; id++){
      for(int ic = 0; ic < num_channels; ic++){
	TLine *lv = new TLine(0, num_bins.at(ic)+use_coll_frac, num_bins_total_compressed, num_bins.at(ic)+use_coll_frac);
	TLine *lh = new TLine(num_bins.at(ic)+use_coll_frac,0, num_bins.at(ic)+use_coll_frac, num_bins_total_compressed);
	lv->SetLineWidth(1.5);
	lh->SetLineWidth(1.5);
	use_coll_frac+=num_bins.at(ic);
	lv->Draw();
	lh->Draw();
      }
    }
  }
  c_coll_frac->Write();


  TH2D h2_coll_full(coll_covariance);
  h2_coll_full.SetName("coll_full");
  TCanvas *c_coll_full = new TCanvas("collapsed covariance matrix");
  c_coll_full->cd();
  c_coll_full->SetFixedAspectRatio();
  h2_coll_full.Draw("colz");
  h2_coll_full.SetTitle("Collapsed covariance matrix");
  h2_coll_full.GetXaxis()->SetTitle("Reco Bin i");
  h2_coll_full.GetYaxis()->SetTitle("Reco Bin j");

  c_coll_full->SetRightMargin(0.150);

  int use_coll_full =0;
  for(int im =0; im<num_modes; im++){
    for(int id =0; id<num_detectors; id++){
      for(int ic = 0; ic < num_channels; ic++){
	TLine *lv = new TLine(0, num_bins.at(ic)+use_coll_full, num_bins_total_compressed, num_bins.at(ic)+use_coll_full);
	TLine *lh = new TLine(num_bins.at(ic)+use_coll_full,0, num_bins.at(ic)+use_coll_full, num_bins_total_compressed);
	lv->SetLineWidth(1.5);
	lh->SetLineWidth(1.5);
	use_coll_full+=num_bins.at(ic);
	lv->Draw();
	lh->Draw();
      }
    }
  }
  c_coll_full->Write();

  for(int m=0; m< variations.size();m++){
    this->plot_one(vec_full_correlation.at(m), variations.at(m)+" Correlation", fout);
    this->plot_one(vec_frac_covariance.at(m), variations.at(m)+" Fractional Covariance", fout);
    this->plot_one(vec_full_covariance.at(m), variations.at(m)+" Full Covariance", fout);

  }



  fout->cd();


  fout->Close();
  return 0;
}


int SBNcovariance::plot_one(TMatrixD matrix, std::string tag, TFile *fin){
  fin->cd();
  TDirectory *individualDir = fin->GetDirectory("individualDir"); 
  if (!individualDir) { 
    individualDir = fin->mkdir("individualDir");       
  }
  fin->cd(); 
  individualDir->cd();



  TH2D h2_full(matrix);
  h2_full.SetName((tag+"_th2d").c_str());
  TCanvas *c_full = new TCanvas((tag+"_canvas").c_str());
  c_full->cd();
  c_full->SetFixedAspectRatio();
  h2_full.Draw("colz");
  h2_full.SetTitle(tag.c_str());
  h2_full.GetXaxis()->SetTitle("Reco Bin i");
  h2_full.GetYaxis()->SetTitle("Reco Bin j");

  c_full->SetRightMargin(0.150);
  int use_full =0;
  for(int im =0; im<num_modes; im++){
    for(int id =0; id<num_detectors; id++){
      for(int ic = 0; ic < num_channels; ic++){
	for(int isc = 0; isc < num_subchannels.at(ic)-1; isc++){
	  TLine *lscv = new TLine(0, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
	  TLine *lsch = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total);
	  lscv->SetLineWidth(2);
	  lsch->SetLineWidth(2);
	  lscv->SetLineColor(kRed);
	  lsch->SetLineColor(kRed);
	  use_full+=num_bins.at(ic);
	  lscv->Draw();
	  lsch->Draw();

	}
	TLine *lv = new TLine(0, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
	TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total);
	lv->SetLineWidth(2);
	lh->SetLineWidth(2);
	use_full+=num_bins.at(ic);
	lv->Draw();
	lh->Draw();

      }
    }
  }
  c_full->Write();

  return 0;
}
