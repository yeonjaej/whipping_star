#include "SBNgenerate.h"
//#include "MCEventWeight.h"

using namespace sbn;

SBNgenerate::SBNgenerate(std::string xmlname) {
	neutrinoModel nullModel;	
	SBNgenerate(xmlname, nullModel);
}

SBNgenerate::SBNgenerate(std::string xmlname, neutrinoModel inModel ) : SBNconfig(xmlname), nu_model(inModel) {


	std::map<std::string, int> parameter_sims;

	//Initialise the central value SBNspec.
	SBNspec tm(xmlname,-1,false);
	spec_CV = tm;
	spec_OSC_sin  = tm;
	spec_OSC_sinsq  = tm;

	int Nfiles = multisim_file.size();

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

				std::cout<<"SBNmultisim::SBNmultisim\t|| Adding a friend tree  "<< treefriendfile<<" to file "<<multisim_file.at(i)<<std::endl;

                         	 trees.at(i)->AddFriend( treefriendname.c_str()   ,  treefriendfile.c_str()   );
			}
                  }
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
			if(branch_variables.at(i).at(k)->getOscillate()){
				trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->true_param_name.c_str(), branch_variables.at(i).at(k)->getTrueValue() );
				trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->true_L_name.c_str(), branch_variables.at(i).at(k)->getTrueL() );
			}
		}
	}


	std::cout<<"SBNgenerate::SBNgenerate\t|| -------------------------------------------------------------\n";
	std::cout<<"SBNgenerate::SBNgenerate\t|| -------------------------------------------------------------\n";
	std::vector<double> base_vec (spec_CV.num_bins_total,0.0);

	for(int j=0;j<Nfiles;j++){

		delete fWeights->at(j);
		fWeights->at(j)=0;

		for(int i=0; i< std::min(  multisim_maxevents.at(j)  ,nentries.at(j)); i++){
			trees.at(j)->GetEntry(i);
			std::map<std::string, std::vector<double>> * thisfWeight = fWeights->at(j);

			if(i%2500==0) std::cout<<"SBNgenerate::SBNgenerate\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<multisim_file[j]<<std::endl;

			double global_weight = 1;
			global_weight = global_weight*multisim_scale.at(j);

			if(thisfWeight->count("bnbcorrection_FluxHist")>0){
				global_weight = global_weight*thisfWeight->at("bnbcorrection_FluxHist").front();
			}

			if(std::isinf(global_weight) || global_weight != global_weight){
				std::cout<<"SBNgenerate::SBNgenerate\t|| ERROR  error @ "<<i<<" in File "<<multisim_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
				exit(EXIT_FAILURE);
			}


			if( this->eventSelection(j) ){
					
				for(int t=0; t<branch_variables.at(j).size();t++){
					//Need the histogram index, the value, the global bin...
					int ih = spec_CV.map_hist.at(branch_variables.at(j).at(t)->associated_hist);
					double reco_var = *(static_cast<double*>(branch_variables.at(j).at(t)->getValue()));
					int reco_bin = spec_CV.getGlobalBinNumber(reco_var,ih);

					//Find if this event should be oscillated
					if(branch_variables.at(j).at(t)->getOscillate()){
						double true_var = *(static_cast<double*>(branch_variables.at(j).at(t)->getTrueValue()));		
						double true_L = 0.1;//*(static_cast<double*>(branch_variables.at(j).at(t)->getTrueLength()));		
						double osc_probability_sin = nu_model.oscProbSin(true_var, true_L);
						double osc_probability_sinsq = nu_model.oscProbSinSq(true_var, true_L);

						spec_OSC_sinsq.hist.at(ih).Fill(reco_var, global_weight*osc_probability_sinsq);
						spec_OSC_sin.hist.at(ih).Fill(reco_var, global_weight*osc_probability_sin);
						
						spec_CV.hist.at(ih).Fill(reco_var,global_weight);
					}else{
						spec_CV.hist.at(ih).Fill(reco_var,global_weight);
						spec_OSC_sinsq.hist.at(ih).Fill(reco_var, global_weight);
						spec_OSC_sin.hist.at(ih).Fill(reco_var, global_weight);
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

int SBNgenerate::writePrecomputedOscSpecs(std::string tag){
	spec_OSC_sinsq.writeOut(tag+"_SINSQ_dm_"+nu_model.mass_tag);
	spec_OSC_sin.writeOut(tag+"_SIN_dm_"+nu_model.mass_tag);

	return 0;
}


bool SBNgenerate::eventSelection(int which_file){
	return true;
}

int SBNgenerate::fillHistograms(int file, int uni, double wei){
	return 0;
}

