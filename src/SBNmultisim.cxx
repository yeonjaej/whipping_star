#include "SBNmultisim.h"
//#include "MCEventWeight.h"

using namespace sbn;


SBNmultisim::SBNmultisim(std::string xmlname) : SBNconfig(xmlname) {

	universes_used = 0;
	tolerence_positivesemi = 1e-5;
	is_small_negative_eigenvalue = false;
	abnormally_large_weight = 1e3;



	std::map<std::string, int> parameter_sims;

	//Initialise the central value SBNspec.
	SBNspec tm(xmlname,-1,false);
	spec_CV = tm;

	int Nfiles = multisim_file.size();

	for(auto &fn: multisim_file){
		files.push_back(new TFile(fn.c_str()));
	}

	for(int i=0; i<multisim_name.size(); i++){
		trees.push_back((TTree*)files.at(i)->Get(multisim_name.at(i).c_str()) );
	}


	std::vector<int> nentries;
	for(auto &t: trees){
		nentries.push_back(t->GetEntries());
	}

	fWeights = new std::vector<std::map<std::string, std::vector<double>>* >(Nfiles,0);


	for(int i=0; i< Nfiles; i++){
		delete fWeights->at(i);	fWeights->at(i) = 0;
		trees.at(i)->SetBranchAddress("weights", &fWeights->at(i) );
		delete fWeights->at(i);	fWeights->at(i) = 0;
		for(int k=0; k<branch_variables.at(i).size(); k++){
			trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->name.c_str(), branch_variables.at(i).at(k)->getValue() );
		}
	}


	variations.clear();


	//This bit will calculate how many "multisims" the file has. if ALL default is the inputted xml value 
	int good_event = 1;
	if(parameter_names.at(0)[0]!="ALL"){
		/*
		   if(parameter_names.at(0).at(0) == "bnbcorrection_FluxHist"){
		   std::cout<<"ERROR!\tSBNmultisim::SBNmultisim ||  bnbcorrection_FluxHist is not a valid varying parameter!\n";
		   exit(EXIT_FAILURE);	
		   }

		   std::vector<int> used_multisims;
		   for(int j=0; j< Nfiles; j++){
		   delete fWeights->at(j);
		   fWeights->at(j)=0;
		   trees.at(j)->GetEntry(good_event);
		   std::vector<double> num_sim_here = fWeights->at(j)->at(parameter_names.at(j)[0]);
		   std::cout<<"SBNmultisim::SBNmultisim || File: "<<j<<" has: "<<num_sim_here.size()<<" universes for parameter: "<<parameter_names.at(j)[0]<<std::endl; 
		   used_multisims.push_back(num_sim_here.size());
		   delete fWeights->at(j);
		   fWeights->at(j)=0;
		   }

		   for(int i=1; i<Nfiles; i++){
		   std::cout<<"File: "<<i-1<<" has "<<used_multisims.at(i-1)<<" multisims"<<std::endl;
		   std::cout<<"File: "<<i<<" has "<<used_multisims.at(i)<<" multisims"<<std::endl;

		   if( used_multisims.at(i)!= used_multisims.at(i-1)){
		   std::cout<<"ERROR: ERROR! SBNmultisim::SBNmultisim || number of Multisims for "<<parameter_names.at(0)[0]<<" are different between files in "<<"  "<<parameter_names.at(i)[0]<<std::endl;
		   exit(EXIT_FAILURE);
		   }
		   universes_used = used_multisims.at(0);
		   }*/	
	}else {
		//ALL catagory
		std::vector<int> used_multisims(Nfiles,0);
		for(int j = 0;j<Nfiles;j++){
			delete fWeights->at(j);
			fWeights->at(j)=0;

			trees.at(j)->GetEntry(good_event);
			for(std::map<std::string, std::vector<double> >::iterator  it = fWeights->at(j)->begin(); it != fWeights->at(j)->end(); ++it) 
			{	
				if(it->first == "bnbcorrection_FluxHist") continue;
				used_multisims.at(j) += it->second.size();
				std::cout<<"SBNmultisim::SBNmultisim\t|| "<<it->first<<" has "<<it->second.size()<<" multisims in file "<<j<<std::endl;
				variations.push_back(it->first);

			}
			delete fWeights->at(j);
			fWeights->at(j)=0;
		}

		//Now remove all, duplicates!
		sort( variations.begin(), variations.end() );
		variations.erase( unique( variations.begin(), variations.end() ), variations.end() );

		//make a map and start filling, before filling find if already in map, if it is check size.

		std::cout<<"SBNmultisim::SBNmultisim\t|| We have "<<variations.size()<<" unique variations: "<<std::endl;
		for(auto &v: variations){
			std::cout<<"SBNmultisim::SBNmultisim\t|| "<<v<<std::endl;
		}


		for(int i=1; i<Nfiles; i++){
			//std::cout<<"File: "<<i-1<<" has "<<used_multisims.at(i-1)<<" multisims"<<std::endl;
			std::cout<<"SBNmultisim::SBNmultisim\t|| File: "<<i<<" has "<<used_multisims.at(i)<<" multisims"<<std::endl;
			if( used_multisims.at(i)!= used_multisims.at(i-1)){
				std::cerr<<"SBNmultisim::SBNmultisim\t|| Warning, number of Multisims for "<<parameter_names.at(0)[0]<<" are different between files"<<std::endl;
				//exit(EXIT_FAILURE);
			}
			universes_used = used_multisims.at(0);
		}	


	}


	std::cout<<"SBNmultisim::SBNmultisim\t|| -------------------------------------------------------------\n";
	std::cout<<"SBNmultisim::SBNmultisim\t|| Initilizing "<<universes_used<<" universes for "<<parameter_names[0][0]<<std::endl;
	std::cout<<"SBNmultisim::SBNmultisim\t|| -------------------------------------------------------------\n";

	std::vector<double> base_vec (spec_CV.num_bins_total,0.0);


	std::cout<<"SBNmultisim::SBNmultisim\t|| Full concatanated vector has : "<<spec_CV.num_bins_total<<std::endl;
	for(int m=0; m<universes_used; m++){
		multi_vecspec.push_back(base_vec);
	}

	std::cout<<"SBNmultisim::SBNmultisim\t|| multi_vecspec now initilized of size :"<<multi_vecspec.size()<<std::endl;

	for(int j=0;j<Nfiles;j++){

		delete fWeights->at(j);
		fWeights->at(j)=0;



		for(int i=0; i< std::min(  multisim_maxevents.at(j)  ,nentries.at(j)); i++){
			trees.at(j)->GetEntry(i);
			std::map<std::string, std::vector<double>> * thisfWeight = fWeights->at(j);

			if(i%2500==0) std::cout<<"SBNmultisim::SBNmultisim\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<multisim_file[j]<<std::endl;

			std::vector<double> weights;
			double global_weight = 1;


			if(thisfWeight->count("bnbcorrection_FluxHist")>0){
				global_weight = global_weight*thisfWeight->at("bnbcorrection_FluxHist").front();
			}

			if(std::isinf(global_weight) || global_weight != global_weight){
				std::cout<<"SBNmultisim::SBNmultisim\t|| ERROR  error @ "<<i<<" in File "<<multisim_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
				exit(EXIT_FAILURE);
			}


			if( this->eventSelection(j) ){

				if(parameter_names.at(j).size()!=1){
					std::cout<<"SBNmultisim::SBNmultisim\t|| ERROR: Currently can only do either 1 multi_sim parameter or 'ALL'"<<std::endl;
					exit(EXIT_FAILURE);

				}
				else if(parameter_names.at(j)[0]=="ALL")
				{
					//Loop over all variations!
					for(std::string &var : variations){

						//check if variation is in this file, if it isn't: then just push back 1's of appropiate number to keep universes consistent
						if(thisfWeight->count(var) <=0 ){
							for(int g=0; g<fWeights->at(0)->at(var).size(); g++){
								weights.push_back(global_weight);
							}
							continue;					
						}

						for(double &wei: thisfWeight->at(var)){
							bool is_inf = std::isinf(wei);
							bool is_nan = wei!=wei;

							if(is_inf || is_nan){
								std::cout<<"SBNmultisim::SBNmultisim\t|| ERROR! Killing :: event # "<<i<<" in File "<<multisim_file.at(j)<<" weight: "<<wei<<" global bnb: "<<global_weight<<" in "<<var<<std::endl;
								exit(EXIT_FAILURE);
							}

							if(wei > abnormally_large_weight){
								std::cout<<"ATTENTION! SBNmultisim::SBNmultisim\t|| HUGE weight: "<<wei<<" at "<<var<<" event "<<i<<" file "<<j<<std::endl;
							}
							weights.push_back(wei*global_weight);

						}


					}//end of variations

					//So the size of weights must equal global universes ya?
					if(weights.size() != universes_used || universes_used != multi_vecspec.size()){
						std::cout<<"SBNmultisim::SBNmultisim\t|| ERROR "<<std::endl;
						std::cout<<"weights.size() "<<weights.size()<<std::endl;
						std::cout<<"universes_used "<<universes_used<<std::endl;
						std::cout<<"multi_vecspec.size() "<<multi_vecspec.size()<<std::endl;
						exit(EXIT_FAILURE);
					}



				}//end of "ALL" option
				//Begininning of single parameter options
				/*		else {

						std::vector<double> this_param_weights = thisfWeight->at(parameter_names.at(j)[0]);



						for(double &wei : this_param_weights){

						if(std::isnan(wei) || wei != wei){
						std::cout<<"Killing. event # "<<i<<" in File "<<multisim_file.at(j)<<" weight: "<<wei<<" global bnb: "<<global_weight<<" in "<<parameter_names.at(j)[0]<<std::endl;
						exit(EXIT_FAILURE);

						weights.push_back(wei*global_weight);

						}//end of this weight loop
						}//end of single parameter scan



						delete fWeights->at(j);fWeights->at(j) = 0;


				//So the size of weights must equal global universes ya?
				if(weights.size() != universes_used || universes_used != multi_vecspec.size()){
				std::cerr<<"num_sim (weights.size()) != universes_used "<<weights.size()<<" "<<universes_used<<std::endl;
				std::cout<<"num_sim (weights.size()) != universes_used "<<weights.size()<<" "<<universes_used<<std::endl;
				//continue;
				exit(EXIT_FAILURE);
				}
				}
				*/

				for(int t=0; t<branch_variables.at(j).size();t++){
					//Need the histogram index, the value, the global bin...
					int ih = spec_CV.map_hist.at(branch_variables.at(j).at(t)->associated_hist);
					double reco_var = *(static_cast<double*>(branch_variables.at(j).at(t)->getValue()));
					int reco_bin = spec_CV.getGlobalBinNumber(reco_var,ih);
					spec_CV.hist.at(ih).Fill(reco_var,global_weight);

					for(int m=0; m< weights.size(); m++){
						if(reco_bin>=0)  multi_vecspec.at(m).at(reco_bin)   +=  weights.at(m);

						//important check. failure mode	
						if(weights.at(m)!=weights.at(m) || std::isinf(weights.at(m)) ){
							std::cout<<"SBNmultisim::SBNmultisim\t|| ERROR weight has a value of: "<<weights.at(m)<<". So I am killing all. on Dim: "<<m<<" global_eright is "<<global_weight<<std::endl;
							exit(EXIT_FAILURE);
						}

					}
				}



			}	
		} //end of entry loop
	}//end of file loop 

	delete fWeights;

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

bool SBNmultisim::eventSelection(int which_file){
	//from here have access to vars_i  and vars_d  to make a selection
	return true;
}

int SBNmultisim::fillHistograms(int file, int uni, double wei){
	//double en = vars_d.at(file)[0];
	//multi_sbnspec.at(uni).hist.at(file).Fill(en, wei);
	return 0;
}

/***************************************************************
 *		And then form a covariance matrix (well 3)
 * ************************************************************/


int SBNmultisim::formCovarianceMatrix(std::string tag){

	full_covariance.ResizeTo(num_bins_total, num_bins_total);
	frac_covariance.ResizeTo(num_bins_total, num_bins_total);
	full_correlation.ResizeTo(num_bins_total, num_bins_total);

	for(auto &h: multi_sbnspec){
		h.calcFullVector();
	}


	spec_CV.calcFullVector();	
	std::vector<double> CV = spec_CV.fullVec;


	std::cout<<"SBNmultisim::formCovariancematrix\t|| Begining to form the "<<num_bins_total<<"X"<<num_bins_total<<" covariance matrix.\n";
	for(int i=0; i<num_bins_total; i++){
		for(int j=0; j<num_bins_total; j++){

			full_covariance(i,j)=0;

			for(int m=0; m < universes_used; m++){

				//full_covariance(i,j) += (CV[i]-multi_sbnspec.at(m).fullVec.at(i))*(CV[j]-multi_sbnspec.at(m).fullVec.at(j));
				full_covariance(i,j) += (CV[i]-multi_vecspec.at(m).at(i))*(CV[j]-multi_vecspec.at(m).at(j));


				if(full_covariance(i,j)!=full_covariance(i,j)){
					std::cout<<"SBNmultisim::formCovariancematrix\t|| ERROR: nan : at (i,j):  "<<i<<" "<<j<<" fullcov: "<<full_covariance(i,j)<<" multi hist sise "<<multi_vecspec.size()<<" CV: "<<CV[i]<<" "<<CV[j]<<" multihisg "<<multi_vecspec[m][i]<<" "<<multi_vecspec[m][j]<<" on dim : "<<m<<std::endl;
					exit(EXIT_FAILURE);
				}

			}
			full_covariance(i,j) = full_covariance(i,j)/( (double)universes_used-1.0);

		}
	}




	std::cout<<"SBNmultisim::formCovariancematrix\t|| Now calculating fractional covariance and correlation matrix from full covariance."<<std::endl;
	for(int i=0; i<num_bins_total; i++){
		for(int j=0; j<num_bins_total; j++){
			frac_covariance(i,j) = full_covariance(i,j)/(spec_CV.fullVec[i]*spec_CV.fullVec[j]) ;
			full_correlation(i,j)= full_covariance(i,j)/(sqrt(full_covariance(i,i))*sqrt(full_covariance(j,j)));
		}
	}


	this->qualityTesting();



	/************************************************************
	 *			Saving to file				    *
	 * *********************************************************/
	TFile *fout=new TFile((tag+".SBNcovar.root").c_str(),"RECREATE");
	fout->cd();
		full_covariance.Write(("full_covariance_"+tag).c_str(),TObject::kWriteDelete);
		frac_covariance.Write(("frac_covariance_"+tag).c_str(),TObject::kWriteDelete);
		full_correlation.Write(("full_correlation_"+tag).c_str(),TObject::kWriteDelete);
	fout->Close();
	//and also writeout 	
	spec_CV.writeOut(tag);


	return 0;
}


int SBNmultisim::qualityTesting(){



	/************************************************************
	 *		Quality Testing Suite			    *
	 * *********************************************************/
	std::cout<<"SBNmultisim::qualityTesting\t||-----------------------------------------------------\n";
	std::cout<<"SBNmultisim::qualityTesting\t||----------------Quality Testing Suite"<<std::endl;
	std::cout<<"SBNmultisim::qualityTesting\t||-----------------------------------------------------\n";


	std::cout<<"SBNmultisim::qualityTesting\t|| Checking if generated matrix is indeed a valid covariance matrix.\n";
	std::cout<<"SBNmultisim::qualityTesting\t|| First checking if matrix is symmetric.\n";
	if(full_covariance.IsSymmetric()){
		std::cout<<"SBNmultisim::qualityTesting\t|| PASS: Generated covariance matrix is symmetric"<<std::endl;
	}else{
		std::cout<<"SBNmultisim::qualityTesting\t||  ERROR result is not symmetric!"<<std::endl;
		exit(EXIT_FAILURE);
	}

	std::cout<<"SBNmultisim::qualityTesting\t|| Checking if generated matrix is positive semi-definite by looking at eigenvalues.\n";
	//if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
	TMatrixDEigen eigen (full_covariance);
	TVectorD eigen_values = eigen.GetEigenValuesRe();


	for(int i=0; i< eigen_values.GetNoElements(); i++){
		if(eigen_values(i)<0){
			is_small_negative_eigenvalue = true;
			if(fabs(eigen_values(i))> tolerence_positivesemi ){
				std::cout<<"SBNmultisim::qualityTesting\t|| ERROR contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}


	if(is_small_negative_eigenvalue){	
		std::cout<<"SBNmultisim::qualityTesting\t|| PASS: Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
	}else{
		std::cout<<"SBNmultisim::qualityTesting\t|| PASS: Generated covariance matrix is positive semi-definite."<<std::endl;
	}
	std::cout<<"SBNmultisim::qualityTesting\t|| Congratulations, matrix is indeed a valid covariance matrix.\n";

	return 0;
}



int SBNmultisim::printMatricies(std::string tag){
	TFile* fout = new TFile(("SBNfit_covariance_plots_"+tag+".root").c_str(),"recreate");
	fout->cd();


	gStyle->SetOptStat(0);

	//correlation
	TH2D h2_corr(full_correlation);
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



	fout->Close();
	return 0;
}



