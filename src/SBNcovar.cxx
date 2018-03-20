#include "SBNcovar.h"
//#include "MCEventWeight.h"

using namespace sbn;

SBNcovar::SBNcovar(std::string xmlname) : SBNconfig(xmlname) {
	std::string dict_location = "../../src/AutoDict_map_string__vector_double____cxx.so";
	gROOT->ProcessLine("#include <map>");
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <string>");
	//	gSystem->Load("/uboone/app/users/markrl/sbnfit/whipping_star/src/mdict_h.so");

	std::cout<<"Trying to load dictionary: "<<dict_location<<std::endl;
	gSystem->Load(  (dict_location).c_str());
	gStyle->SetOptStat(0);


	universes_used = 0;
	tolerence_positivesemi = 1e-5;
	is_small_negative_eigenvalue = false;
	bool DEBUG_BOOL = false;
	int num_skipped = 0;
	int num_skipped_kaon = 0;

	std::map<std::string, int> parameter_sims;

	SBNspec tm(xmlname,-1);
	spec_CV = tm;
	spec_sig  = tm;


	//Load all files as per xml
	std::vector<TFile *> files;	
	std::vector<TTree *> trees;	

	int Nfiles = multisim_file.size();

	for(auto &fn: multisim_file){
		files.push_back(new TFile(fn.c_str()));
	}

	for(int i=0; i<multisim_name.size(); i++){
		trees.push_back((TTree*)files.at(i)->Get(multisim_name[i].c_str()) );
	}

	std::vector<int> nentries;
	for(auto &t: trees){
		//nentries.push_back(100000);
		nentries.push_back(t->GetEntries());
	}

	//Very brute 2D
	for(int i=0; i<2; i++){
		TH2D tmpcv((channel_names.at(2*i)+"_cv").c_str(), channel_names.at(2*i).c_str(), num_bins.at(0), &bin_edges.at(0)[0], num_bins.at(1), &bin_edges.at(1)[0] );
		h_spec2d_CV.push_back( tmpcv);
		TH2D tmpsig((channel_names.at(2*i)+"_sig").c_str(), channel_names.at(2*i).c_str(), num_bins.at(0), &bin_edges.at(0)[0], num_bins.at(1), &bin_edges.at(1)[0] );
		h_spec2d_sig.push_back( tmpsig);
	}

	//signal crap todo
	//
	//TFile *fS = new TFile("/uboone/app/users/mastbaum/leerw/leerw.root");
	//TTree * tnS =  (TTree*)fS->Get("leerw;1");
	//float lee=0;
	//tnS->SetBranchAddress("leerw",&lee);
	//SBNspec spec_sig(xmlname,-1);


	vars_i= std::vector<std::vector<int>>(Nfiles   , std::vector<int>(branch_names_int.at(0).size(),0));
	vars_d= std::vector<std::vector<double>>(Nfiles   , std::vector<double>(branch_names_double.at(0).size(),0.0));


	fWeights = new std::vector<std::map<std::string, std::vector<double>>* >(Nfiles,0);
	fLepMom = new std::vector<TLorentzVector*>(Nfiles,0);


	for(int i=0; i< Nfiles; i++){
		delete fWeights->at(i);	fWeights->at(i) = 0;
		delete fLepMom->at(i);	fLepMom->at(i) = 0;


		trees.at(i)->SetBranchAddress("mcweight", &fWeights->at(i) );// , &bWeight->at(i) );
		trees.at(i)->SetBranchAddress("leptonMom", &fLepMom->at(i) );//, &bLepMom->at(i)  );

		delete fWeights->at(i);	fWeights->at(i) = 0;
		delete fLepMom->at(i);	fLepMom->at(i) = 0;



		for(auto &bfni: branch_names_int){
			for(int k=0; k< bfni.size();k++){
				trees.at(i)->SetBranchAddress(bfni[k].c_str(), &(vars_i.at(i).at(k)));
			}
		}
		for(auto &bfnd: branch_names_double){
			for(int k=0; k< bfnd.size();k++){
				trees.at(i)->SetBranchAddress(bfnd[k].c_str(), &(vars_d.at(i).at(k)));
			}
		}
	}



	/*************************** test *************************/
	if(false){
		int v0p=0;
		int v1p=0;
		int v0m=0;
		int v1m=0;
		int v0z=0;int v1z=0;	
		for(int i =0; i< 200000; i++){
			trees.at(0)->GetEntry(i);
			trees.at(1)->GetEntry(i);



			if(fWeights->at(0)->at("kminus_PrimaryHadronNormalization").size()!=1000){ v0m++;}
			if(fWeights->at(1)->at("kminus_PrimaryHadronNormalization").size()!=1000){ v1m++;}

			if(fWeights->at(0)->at("kplus_PrimaryHadronFeynmanScaling").size()!=1000){v0p++;}
			if(fWeights->at(1)->at("kplus_PrimaryHadronFeynmanScaling").size()!=1000){v1p++;}

			if(fWeights->at(0)->at("kzero_PrimaryHadronSanfordWang").size()!=1000){v0z++;}
			if(fWeights->at(1)->at("kzero_PrimaryHadronSanfordWang").size()!=1000){v1z++;}

			std::cout<<i<<" "<<"nue kminus_PrimaryHadronNormalization "<<fWeights->at(0)->at("kminus_PrimaryHadronNormalization").size()<<std::endl;
			std::cout<<i<<" "<<"bnb kminus_PrimaryHadronNormalization "<<fWeights->at(1)->at("kminus_PrimaryHadronNormalization").size()<<std::endl;

			std::cout<<i<<" "<<"nue kplus_PrimaryHadronFeynmanScaling "<<fWeights->at(0)->at("kplus_PrimaryHadronFeynmanScaling").size()<<std::endl;
			std::cout<<i<<" "<<"bnb kplus_PrimaryHadronFeynmanScaling "<<fWeights->at(1)->at("kplus_PrimaryHadronFeynmanScaling").size()<<std::endl;

			std::cout<<i<<" "<<"nue kzero_PrimaryHadronSanfordWang "<<fWeights->at(0)->at("kzero_PrimaryHadronSanfordWang").size()<<std::endl;
			std::cout<<i<<" "<<"bnb kzero_PrimaryHadronSanfordWang "<<fWeights->at(1)->at("kzero_PrimaryHadronSanfordWang").size()<<std::endl;
			std::cout<<"kplus: ";
			for(auto wei: fWeights->at(0)->at("kplus_PrimaryHadronFeynmanScaling")){
				std::cout<<wei<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"kmin: ";
			for(auto wei: fWeights->at(0)->at("kminus_PrimaryHadronNormalization")){
				std::cout<<wei<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"kzero: ";
			for(auto wei: fWeights->at(0)->at("kzero_PrimaryHadronSanfordWang")){
				std::cout<<wei<<" ";
			}
			std::cout<<std::endl;


		}

		std::cout<<"1: Min: "<<v1m<<" plus: "<<v1p<<"  zero "<<v1z<<std::endl;
		std::cout<<"0: Min: "<<v0m<<" plus: "<<v0p<<"  zero "<<v0z<<std::endl;
		exit(EXIT_FAILURE);

	}





	//This bit will calculate how many "multisims" the file has. if ALL default is the inputted xml value 
	// use a known good event (2 has been checked)
	int good_event = 2;
	if(parameter_names.at(0)[0]!="ALL"){
		std::vector<int> used_multisims;
		for(int j=0; j< Nfiles; j++){
			delete fWeights->at(j);
			fWeights->at(j)=0;


			trees.at(j)->GetEntry(good_event);
			std::vector<double> num_sim_here = fWeights->at(j)->at(parameter_names.at(j)[0]);
			std::cout<<"File: "<<j<<" has: "<<num_sim_here.size()<<" universes for parameter: "<<parameter_names.at(j)[0]<<std::endl; 
			used_multisims.push_back(num_sim_here.size());
			delete fWeights->at(j);
			fWeights->at(j)=0;
		}

		for(int i=1; i<Nfiles; i++){
			std::cout<<"File: "<<i-1<<" has "<<used_multisims.at(i-1)<<" multisims"<<std::endl;
			std::cout<<"File: "<<i<<" has "<<used_multisims.at(i)<<" multisims"<<std::endl;

			if( used_multisims.at(i)!= used_multisims.at(i-1)){
				std::cerr<<"ERROR: number of Multisims for "<<parameter_names.at(0)[0]<<" are different between files in "<<"  "<<parameter_names.at(i)[0]<<std::endl;
				exit(EXIT_FAILURE);
			}
			universes_used = used_multisims.at(0);
		}	
	}else {






		//warning, currently assumes all the same
		std::vector<int> used_multisims(Nfiles,0);
		for(int j = 0;j<Nfiles;j++){
			delete fWeights->at(j);
			fWeights->at(j)=0;

			trees.at(j)->GetEntry(good_event);
			for(std::map<std::string, std::vector<double> >::iterator  it = fWeights->at(j)->begin(); it != fWeights->at(j)->end(); ++it) 
			{
				if( it->first == "bnbcorrection_FluxHist") continue;


				//	if(it->first == "kplus_PrimaryHadronFeynmanScaling")continue;
				//	|| it->first == "kzero_PrimaryHadronSanfordWang" || it->first== "kminus_PrimaryHadronNormalization")continue;

				used_multisims.at(j) += it->second.size();
				std::cout<<"ALL: "<<it->first<<" has "<<it->second.size()<<" multisims in file "<<j<<std::endl;
			}
			delete fWeights->at(j);
			fWeights->at(j)=0;
		}

		for(int i=1; i<Nfiles; i++){
			std::cout<<"File: "<<i-1<<" has "<<used_multisims.at(i-1)<<" multisims"<<std::endl;
			std::cout<<"File: "<<i<<" has "<<used_multisims.at(i)<<" multisims"<<std::endl;
			if( used_multisims.at(i)!= used_multisims.at(i-1)){
				std::cerr<<"ERROR: number of Multisims for "<<parameter_names.at(0)[0]<<" are different between files"<<std::endl;
				exit(EXIT_FAILURE);
			}
			universes_used = used_multisims.at(0);
		}	


	}


	//Ok now we know now many universes we have, initilize all the sbnspecs
	std::cout<<"Initilizing "<<universes_used<<" universes for "<<parameter_names[0][0]<<std::endl;


		

	SBNspec temp_spec(xmlname,0);
	int num_2d_bins = temp_spec.num_bins[0]*temp_spec.num_bins[1]+temp_spec.num_bins[2]*temp_spec.num_bins[3];
	std::vector<double> base_vec (temp_spec.num_bins_total,0.0);
	std::vector<double> base_vec2D(num_2d_bins, 0);
	std::cout<<"We have : "<<num_2d_bins<<" 2d bins."<<std::endl;

	vecspec2DCV = base_vec2D;

	std::cout<<"Full vector has : "<<temp_spec.num_bins_total<<std::endl;
	for(int m=0; m<universes_used; m++){
		if(m%1000==0)std::cout<<"Initilized : "<<m<<" of "<<universes_used<<std::endl;
		//SBNspec temp_spec(xmlname,m);
		//multi_sbnspec.push_back(temp_spec);

		multi_vecspec.push_back( base_vec);
		multi_vecspec2D.push_back( base_vec2D);
	}

	TFile * flee = new TFile("/home/mark/work/uBooNE/lee_unfolding/git/signal_prediction/sp/Unfold/bin/CCQE_SVD_final_tgraph.root");
	TGraph *lee = (TGraph*)flee->Get("Graph;1");

	TVector3 zaxis(0,0,1);


	for(int j=0;j<Nfiles;j++){
		int nerr=0;

		delete fWeights->at(j);
		fWeights->at(j)=0;

		double pot_factor = pot.at(j)/(pot_scaling.at(j) * (double)nentries.at(j));

		for(int i=0; i< nentries.at(j); i++){
			if(i%1000==0) std::cout<<"Event: "<<i<<" of "<<nentries[j]<<" from File: "<<multisim_file[j]<<" POT factor: "<<pot_factor<<std::endl;

			std::vector<double> weights;
			double global_weight = 1;

			trees.at(j)->GetEntry(i);

			std::map<std::string, std::vector<double>> * thisfWeight = fWeights->at(j);

			//first a check to eliminiate nan/inf in the global bnb correction FluxHist;
			global_weight = thisfWeight->at("bnbcorrection_FluxHist").at(0);


			if(std::isinf(global_weight) || global_weight != global_weight){
				std::cout<<"Skipping event # "<<i<<" in File "<<multisim_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
				num_skipped ++;
				continue;
			}
			global_weight = global_weight*pot_factor;


			//Check to see if event is ok of kaon event-bug
			if(parameter_names.at(j)[0]=="kplus_PrimaryHadronFeynmanScaling" || parameter_names.at(j)[0] == "kminus_PrimaryHadronNormalization" || parameter_names.at(j)[0]== "kzero_PrimaryHadronSanfordWang"){
				if( thisfWeight->at("kplus_PrimaryHadronFeynmanScaling").size()!=1000 || thisfWeight->at("kminus_PrimaryHadronNormalization").size()!=1000 || thisfWeight->at("kzero_PrimaryHadronSanfordWang").size() != 1000){
					std::cout<<"Skipping event # "<<i<<" in File "<<multisim_file.at(j)<<" as one of the kplus/zero/minus is broke"<<std::endl;
					num_skipped_kaon++;
					continue;
				}
			}

			//here we put the selection criteria, for example nuance interaction 1001 == CCQE, virtual bool
			if( this->eventSelection(j) ){

				if(parameter_names.at(j).size()!=1){
					std::cout<<"ERROR: Currently can only do either 1 multi_sim parameter or 'ALL'"<<std::endl;
					exit(EXIT_FAILURE);

				}
				else if(parameter_names.at(j)[0]=="ALL")
				{

					for(std::map<std::string, std::vector<double> >::iterator  it = thisfWeight->begin(); it != thisfWeight->end(); ++it) 
					{

						if(it->first == "bnbcorrection_FluxHist") continue;


						//	std::cout<<"file: "<<j<<" "<<it->first<<" "<<" size "<<it->second.size()<<std::endl;
						for(auto &wei: it->second)
						{
							if(std::isinf(wei) || wei!= wei){
								std::cout<<"Killing. event # "<<i<<" in File "<<multisim_file.at(j)<<" weight: "<<wei<<" global bnb: "<<global_weight<<" in "<<it->first<<std::endl;
								exit(EXIT_FAILURE);
							}

							if(wei> 10){
								std::cout<<"ATTENTION: Large weight: "<<wei<<" at "<<it->first<<" event "<<i<<" file "<<j<<std::endl;
								
							}
							if(wei> 1e5){
								std::cout<<"ATTENTION: HUGE weight: "<<wei<<" at "<<it->first<<" event "<<i<<" file "<<j<<std::endl;
								wei = 1;
							}


							weights.push_back(wei*global_weight);



						}//end of weight vector loop for this parameter


					}//end of all parameter iterator




				}//end of "ALL" option
				//Begininning of single parameter options
				else {

					std::vector<double> this_param_weights = thisfWeight->at(parameter_names.at(j)[0]);

					for(double &wei : this_param_weights){
						//std::cout<<"Weights: "<<wei<<" event: "<<i<<std::endl;

						if(std::isnan(wei) || wei != wei){
							std::cout<<"Killing. event # "<<i<<" in File "<<multisim_file.at(j)<<" weight: "<<wei<<" global bnb: "<<global_weight<<" in "<<parameter_names.at(j)[0]<<std::endl;
							exit(EXIT_FAILURE);
						}

						weights.push_back(wei*global_weight);


					}//end of this weight loop




				}//end of single parameter scan


				delete fWeights->at(j);fWeights->at(j) = 0;

				//So the size of weights must equal global universes ya?
				if(weights.size() > universes_used){
					std::cerr<<"num_sim (weights.size()) != universes_used "<<weights.size()<<" "<<universes_used<<std::endl;
					std::cout<<"num_sim (weights.size()) != universes_used "<<weights.size()<<" "<<universes_used<<std::endl;
					continue;
					nerr++;
					//exit(EXIT_FAILURE);
				}

				double en = vars_d.at(j)[0];
				double lepAngle = zaxis.Angle(fLepMom->at(j)->Vect());		

				int en_bin = temp_spec.getGlobalBinNumber(en, 2*j);
				int cos_bin = temp_spec.getGlobalBinNumber(cos(lepAngle), 2*j+1);
				//std::cout<<"m: "<<m<<"/weights.size() "<<weights.size()<<" EBIN: "<<en_bin<<" CBIN: "<<cos_bin<<" vec: "<<multi_vecspec.at(m).size()<<std::endl;
				//std::cout<<"energy hist: "<<2*j<<" en: "<<en<<" en_bin: "<<en_bin<<std::endl;
				//std::cout<<"cosergy hist: "<<2*j+1<<" cos: "<<cos(lepAngle)<<" cos_bin: "<<cos_bin<<std::endl;
				int en_bin2 = temp_spec.getLocalBinNumber(en,0);
				int cos_bin2 = temp_spec.getLocalBinNumber(cos(lepAngle), 1);



				int bin2D = en_bin2 + num_bins.at(0)*cos_bin2 + j*num_bins.at(0)*num_bins.at(1);
			
				for(int m=0; m<weights.size(); m++){
					//This is the part where we will every histogram in this Universe
					//this->fillHistograms(j, m, weights.at(m) );
					if(en_bin>=0)  multi_vecspec.at(m).at(en_bin)   +=  weights.at(m);
					if(cos_bin>=0) multi_vecspec.at(m).at(cos_bin)  +=  weights.at(m);

					//important check. failure mode	
					if(weights[m]!=weights[m] || std::isinf(weights[m]) ){
						std::cout<<"ERROR: weight has a value of: "<<weights.at(m)<<". So I am killing all. on Dim: "<<m<<" energy"<<vars_d.at(j)[0]<<" global_eright is "<<global_weight<<std::endl;
						exit(EXIT_FAILURE);
					}

					if(en_bin2 >=0&& cos_bin2>=0){
						//std::cout<<"j: "<<j<<" en_bin2: "<<en_bin2<<" cos_bin2: "<<cos_bin2<<" bin2D: "<<bin2D<<" m: "<<m<<std::endl;
						multi_vecspec2D.at(m).at(bin2D) += weights.at(m);
					}



				}

				//blarg, how will I treat this spectrum
				spec_CV.hist.at(2*j).Fill(en,global_weight);
				spec_CV.hist.at(2*j+1).Fill(cos(lepAngle),global_weight);
			
				double lee_weight =1.0;
				if(j==0 && en < 0.8 ){
					lee_weight = 0.5*lee->Eval(1000*en);	
				}


				spec_sig.hist.at(2*j).Fill(en,global_weight*lee_weight);
				spec_sig.hist.at(2*j+1).Fill(cos(lepAngle),global_weight*lee_weight);
				

				h_spec2d_CV.at(j).Fill(en,cos(lepAngle),global_weight);				
				h_spec2d_sig.at(j).Fill(en,cos(lepAngle),global_weight*lee_weight);				
	
				if(en_bin2 >=0&& cos_bin2>=0) vecspec2DCV.at(bin2D) += global_weight;	



				delete fLepMom->at(j);fLepMom->at(j) = 0;


			}//end of CCQE..interaction type check.. OH no this should be in fillHistograms
		} //end of entry loop
		std::cout<<"ERRORed out of : "<<nerr<<std::endl;	
	}//end of file loop 

	delete fWeights;
	delete fLepMom;

	std::cout<<"ATTENTION: Skipped a total of: "<<num_skipped<<" Due to nan/inf bnb weights"<<std::endl;
	std::cout<<"ATTENTION: Skipped a total of: "<<num_skipped_kaon<<" Due to kaon bug"<<std::endl;


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

bool SBNcovar::eventSelection(int which_file){
	//from here have access to vars_i  and vars_d  to make a selection

	bool ans = false;
	if(which_file==0){
		if(vars_i.at(which_file)[3] == 11 || vars_i.at(which_file)[3] == -11  ){
			ans = true;
		}
	} else if (which_file==1){
		if(vars_i.at(which_file)[3] == 13 || vars_i.at(which_file)[3] == -13  ){
			ans = true;
		}


	}

	return ans;
}

int SBNcovar::fillHistograms(int file, int uni, double wei){
	//Fill the histograms
	//
	/*	TRandom3 rangen(0);
		double sigma = 0;
		if(file==0){
		sigma=;
		}else if(file==1){
		sigma=;
		}	
		double en = rangen.Gaus( vars_d.at(file), sigma );
	 */	
	double en = vars_d.at(file)[0];
	//multi_sbnspec.at(uni).hist.at(file).Fill(en, wei);

	return 0;
}

/***************************************************************
 *		And then form a covariance matrix (well 3)
 * ************************************************************/


int SBNcovar::formCovarianceMatrix(){

	int num_bins_total2D = vecspec2DCV.size();

	full_covariance2D.ResizeTo(vecspec2DCV.size(), vecspec2DCV.size());
	frac_covariance2D.ResizeTo(vecspec2DCV.size(), vecspec2DCV.size());
	full_correlation2D.ResizeTo(vecspec2DCV.size(), vecspec2DCV.size());

	full_covariance.ResizeTo(num_bins_total, num_bins_total);
	frac_covariance.ResizeTo(num_bins_total, num_bins_total);
	full_correlation.ResizeTo(num_bins_total, num_bins_total);

	//prepare three TH2D for plotting 
	TH2D * hist_frac_cov = new TH2D("Frac Cov","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
	TH2D * hist_full_cor = new TH2D("Corr","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
	TH2D * hist_full_cov = new TH2D("Full Cov","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);

	TH2D * hist_frac_cov2D = new TH2D("Frac Cov","",vecspec2DCV.size(),1,vecspec2DCV.size(), vecspec2DCV.size(),1,vecspec2DCV.size());
	TH2D * hist_full_cor2D = new TH2D("Corr","",vecspec2DCV.size(),1,vecspec2DCV.size(), vecspec2DCV.size(),1,vecspec2DCV.size());
	TH2D * hist_full_cov2D= new TH2D("Full Cov","",vecspec2DCV.size(),1,vecspec2DCV.size(), vecspec2DCV.size(),1,vecspec2DCV.size());


	for(auto &h: multi_sbnspec){
		h.calcFullVector();
	}


	spec_CV.calcFullVector();	
	std::vector<double> CV = spec_CV.fullVec;

	//for(auto k=0;k<spec_CV.fullVec.size();k++){
	//	std::cout<<" "<<spec_CV.fullVec[k]<<" "<<multi_sbnspec.at(0).fullVec[k]<<std::endl;
	//}
	std::cout<<"multi_sbnspec.size(): "<<multi_vecspec.size()<<" universes_used: "<<universes_used<<std::endl;

	for(int i=0; i<num_bins_total; i++){
		for(int j=0; j<num_bins_total; j++){

			full_covariance(i,j)=0;

			for(int m=0; m < universes_used; m++){

				//full_covariance(i,j) += (CV[i]-multi_sbnspec.at(m).fullVec.at(i))*(CV[j]-multi_sbnspec.at(m).fullVec.at(j));
				full_covariance(i,j) += (CV[i]-multi_vecspec.at(m).at(i))*(CV[j]-multi_vecspec.at(m).at(j));


				if(full_covariance(i,j)!=full_covariance(i,j)){
					std::cout<<"ERROR: nan : at (i,j):  "<<i<<" "<<j<<" fullcov: "<<full_covariance(i,j)<<" multi hist sise "<<multi_vecspec.size()<<" CV: "<<CV[i]<<" "<<CV[j]<<" multihisg "<<multi_vecspec[m][i]<<" "<<multi_vecspec[m][j]<<" on dim : "<<m<<std::endl;
				}



			}
			full_covariance(i,j) = full_covariance(i,j)/( (double)universes_used-1.0);

		}
	}




	for(int i=0; i<num_bins_total; i++){
		for(int j=0; j<num_bins_total; j++){

			frac_covariance(i,j) = full_covariance(i,j)/(spec_CV.fullVec[i]*spec_CV.fullVec[j]) ;
			full_correlation(i,j)= full_covariance(i,j)/(sqrt(full_covariance(i,i))*sqrt(full_covariance(j,j)));
			//	std::cout<<i<<" "<<j<<" "<<full_correlation(i,j)<<" "<<full_covariance(i,j)<<" "<<full_covariance(i,i)<<" "<<full_covariance(j,j)<<" uni-1 "<<universes_used-1.0<<std::endl;
			hist_frac_cov->SetBinContent(i+1,j+1,frac_covariance(i,j));
			hist_full_cor->SetBinContent(i+1,j+1,full_correlation(i,j));
			hist_full_cov->SetBinContent(i+1,j+1,full_covariance(i,j));

		}
	}


	std::cout<<"Starting on 2D covariance construction"<<std::endl;


	for(int i=0; i<num_bins_total2D; i++){
		std::cout<<i<<"/"<<num_bins_total2D<<std::endl;
		for(int j=0; j<num_bins_total2D; j++){

			full_covariance2D(i,j)=0;

			for(int m=0; m < universes_used; m++){
				
				full_covariance2D(i,j) += (vecspec2DCV.at(i)-multi_vecspec2D.at(m).at(i) )*(vecspec2DCV.at(j)-multi_vecspec2D.at(m).at(j));
				
				if(full_covariance2D(i,j) != full_covariance2D(i,j)){

					std::cout<<"ERROR: nan in 2d full cov @ i "<<i<<" j  "<<j<<" m "<<m<<" cvi "<<vecspec2DCV.at(i)<<" mi "<<multi_vecspec2D.at(m).at(i)<<" cvj "<<vecspec2DCV.at(j)<<" mj "<<multi_vecspec2D.at(m).at(j)<<std::endl;

				}
	

			}
			full_covariance2D(i,j) = full_covariance2D(i,j)/( (double)universes_used-1.0);

		}
	}




	for(int i=0; i<num_bins_total2D; i++){
		for(int j=0; j<num_bins_total2D; j++){

			if(vecspec2DCV.at(i)==0 || vecspec2DCV.at(j) == 0){
			//	std::cout<<" frac covariance  is nan! vec2D @i,j "<<i<<","<<j<<" "<<vecspec2DCV.at(i)<<" "<<vecspec2DCV.at(j)<<std::endl;
			}
			
			if(full_covariance2D(i,i)==0 || full_covariance2D(j,j) == 0){
			//	std::cout<<" correlation  is nan! v @ ii,jj "<<i<<","<<j<<" "<<full_covariance2D(i,i)<<" "<<full_covariance(j,j)<<std::endl;
			}

			frac_covariance2D(i,j) = full_covariance2D(i,j)/(vecspec2DCV[i]*vecspec2DCV[j]) ;
			full_correlation2D(i,j)= full_covariance2D(i,j)/(sqrt(full_covariance2D(i,i))*sqrt(full_covariance2D(j,j)));
			//	std::cout<<i<<" "<<j<<" "<<full_correlation(i,j)<<" "<<full_covariance(i,j)<<" "<<full_covariance(i,i)<<" "<<full_covariance(j,j)<<" uni-1 "<<universes_used-1.0<<std::endl;
			hist_frac_cov2D->SetBinContent(i+1,j+1,frac_covariance2D(i,j));
			hist_full_cor2D->SetBinContent(i+1,j+1,full_correlation2D(i,j));
			hist_full_cov2D->SetBinContent(i+1,j+1,full_covariance2D(i,j));

		}
	}





	/************************************************************
	 *			Saving to file				    *
	 * *********************************************************/
	std::string nn = "covariance_matrix_" + parameter_names[0][0]+".root";

	TFile *ftest=new TFile(nn.c_str(),"RECREATE");
	ftest->cd();

	/*
	   TCanvas *cspline =  new TCanvas("Splines");
	   int num_hists = multi_sbnspec.at(0).hist.size();
	   cspline->Divide(num_hists,1);

	   for(int h=0; h<spec_CV.hist.size(); h++){
	   cspline->cd(h+1);
	   spec_CV.hist.at(h).SetLineColor(kBlack);
	   spec_CV.hist.at(h).SetLineWidth(4);
	   spec_CV.hist.at(h).SetTitle( (fullnames[h]+" : "+parameter_names[0][0]).c_str() );
	   spec_CV.hist.at(h).Draw("L SAME");
	   }

	   TRandom3 * rangen = new TRandom3(0);
	   for(int m=0; m< universes_used; m++){
	   for(int h=0; h<multi_sbnspec.at(m).hist.size(); h++){
	   cspline->cd(h+1);
	   TGraph *g = new TGraph( );
	   multi_sbnspec.at(m).hist.at(h).SetLineColor(rangen->Uniform(400,900));

	   multi_sbnspec.at(m).hist.at(h).Draw("L SAME");
	   }
	   }
	   delete rangen;

	   for(int h=0; h<spec_CV.hist.size(); h++){
	   cspline->cd(h+1);
	   spec_CV.hist.at(h).SetLineWidth(4);
	   spec_CV.hist.at(h).Draw("L SAME");
	   }


	   cspline->Write();
	   std::string ppsp = "splines_"+parameter_names[0][0]+".pdf";

	   cspline->SaveAs(ppsp.c_str());
	 */




	//matricies
	TCanvas *c1 =  new TCanvas("Fractional Covariance Matrix");
	c1->cd();

	hist_frac_cov->SetTitle("Fractional Covariance Matrix (sys only)");
	hist_frac_cov->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	hist_frac_cov->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	hist_frac_cov->Draw("COLZ");
	c1->Write();

	TCanvas *c2 =  new TCanvas("Correlation Matrix");
	c2->cd();

	hist_full_cor->SetTitle("Correlation Matrix (sys only)");
	hist_full_cor->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cor->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cor->Draw("COLZ");
	c2->Write();	

	TCanvas *c3 =  new TCanvas("Covariance Matrix");
	c3->cd();

	hist_full_cov->SetTitle("Covariance Matrix (sys only)");
	hist_full_cov->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cov->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cov->Draw("COLZ");
	c3->Write();

	std::string pp = "Fractional Covarariance and Correlation: "+parameter_names[0][0];
	TCanvas *cboth = new TCanvas(pp.c_str());
	cboth->SetCanvasSize(1800,600);

	hist_frac_cov->SetTitle( ("Fractional Covariance: "+ parameter_names[0][0]).c_str());
	hist_full_cor->SetTitle( ("Correlation: "+ parameter_names[0][0]).c_str());
	cboth->Divide(2,1);

	//cboth->SetFixedAspectRatio();
	cboth->cd(1);
	//	cboth->SetBorderSize(20);

	hist_frac_cov->Draw("COLZ");

	//cboth->SetRightMargin(0.30);
	cboth->cd(2);


	hist_full_cor->Draw("COLZ");
	//cboth->SetRightMargin(0.30);
	//cboth->Update();
	cboth->Write();

	hist_full_cov->Write();
	hist_frac_cov->Write();
	hist_full_cor->Write();

	std::string ppdf = "covar_plots_"+parameter_names[0][0]+".pdf";
	cboth->SaveAs(ppdf.c_str());


	frac_covariance.Write();
	full_covariance.Write();
	full_correlation.Write();


	spec_CV.writeOut("CV.root");
	spec_sig.writeOut("SIG.root");

	std::cout<<"Starting to print 2d:"<<std::endl;
	TCanvas * c2d = new TCanvas();
	c2d->Divide(2,2);
	c2d->cd(1);
	h_spec2d_CV.at(0).Draw("COLZ");
	c2d->cd(2);
	h_spec2d_CV.at(1).Draw("COLZ");
	
	c2d->SaveAs("2d_cv.pdf","pdf");

	TCanvas * clong = new TCanvas();	
	TPad * p = (TPad*)clong->cd();
	p->SetLogy();
	std::vector<double> tmpx;
	for(int i=0; i< vecspec2DCV.size();i++) {
		tmpx.push_back(i);
	}

	TGraph * g2d = new TGraph(vecspec2DCV.size(), &tmpx[0],&vecspec2DCV[0]);
	g2d->Draw("ACP");

	clong->SaveAs("2d_long.pdf","pdf");



	//matricies
	TCanvas *c12D =  new TCanvas("Fractional Covariance Matrix 2D");
	c12D->cd();

	hist_frac_cov2D->SetTitle("Fractional Covariance Matrix (sys only)");
	hist_frac_cov2D->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	hist_frac_cov2D->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	hist_frac_cov2D->Draw("COLZ");

	TCanvas *c22D =  new TCanvas("Correlation Matrix 2D");
	c22D->cd();

	hist_full_cor2D->SetTitle("Correlation Matrix (sys only)");
	hist_full_cor2D->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cor2D->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cor2D->Draw("COLZ");

	TCanvas *c32D =  new TCanvas("Covariance Matrix 2D");
	c32D->cd();

	hist_full_cov2D->SetTitle("Covariance Matrix (sys only) 2D");
	hist_full_cov2D->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cov2D->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	hist_full_cov2D->Draw("COLZ");
	
	std::string pp2D = "Fractional Covarariance and Correlation 2D: "+parameter_names[0][0];
	TCanvas *cboth2D = new TCanvas(pp2D.c_str());
	cboth2D->SetCanvasSize(1800,600);

	hist_frac_cov2D->SetTitle( ("Fractional Covariance 2D: "+ parameter_names[0][0]).c_str());
	hist_full_cor2D->SetTitle( ("Correlation 2D: "+ parameter_names[0][0]).c_str());
	cboth2D->Divide(2,1);

	//cboth->SetFixedAspectRatio();
	cboth2D->cd(1);
	//	cboth->SetBorderSize(20);

	hist_frac_cov2D->Draw("COLZ");

	//cboth->SetRightMargin(0.30);
	cboth2D->cd(2);



	hist_full_cor2D->Draw("COLZ");
	std::string ppdf2D = "covar_plots_2D_"+parameter_names[0][0]+".pdf";
	cboth2D->SaveAs(ppdf2D.c_str());

	ftest->cd();
	h_spec2d_CV.at(0).Write();
	h_spec2d_CV.at(1).Write();

	
	h_spec2d_sig.at(0).Write();
	h_spec2d_sig.at(1).Write();
	

	cboth2D->Write();

	c12D->Write();
	c22D->Write();
	c32D->Write();

	hist_full_cov2D->Write();
	hist_frac_cov2D->Write();
	hist_full_cor2D->Write();

	frac_covariance2D.Write();
	full_covariance2D.Write();
	full_correlation2D.Write();



	ftest->Close();
	/************************************************************
	 *		Quality Testing Suite			    *
	 * *********************************************************/


	if(full_covariance.IsSymmetric()){
		std::cout<<"Generated covariance matrix is symmetric"<<std::endl;
	}else{
		std::cerr<<"ERROR: SBNcovar::formCovarianceMatrix, result is not symmetric!"<<std::endl;
		exit(EXIT_FAILURE);
	}


	//if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
	TMatrixDEigen eigen (full_covariance);
	TVectorD eigen_values = eigen.GetEigenValuesRe();


	for(int i=0; i< eigen_values.GetNoElements(); i++){
		if(eigen_values(i)<0){
			is_small_negative_eigenvalue = true;
			if(fabs(eigen_values(i))> tolerence_positivesemi ){
				std::cerr<<"ERROR: SBNcovar::formCovarianceMatrix, contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}


	if(is_small_negative_eigenvalue){	
		std::cout<<"Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
	}else{
		std::cout<<"Generated covariance matrix is also positive semi-definite."<<std::endl;
	}


	/*
	   for (int i=1;i<=11 ;i++){
	   hist_frac_cov->GetYaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
	   hist_frac_cov->GetXaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
	   hist_full_cov->GetYaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
	   hist_full_cov->GetXaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
	   hist_full_cor->GetYaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
	   hist_full_cor->GetXaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
	   }
	 */

	return 0;
}


