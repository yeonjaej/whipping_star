#include "SBNspec.h"
using namespace sbn;





SBNspec::SBNspec(std::string whichxml, int which_universe, bool isverbose) : SBNconfig(whichxml,isverbose){

	//Initialise all the things
	//for every multisim, create a vector of histograms, one for every subchannel we want 
	int ctr=0;
	for(auto fn: fullnames){
		for(int c=0; c<channel_names.size(); c++){
			if(fn.find(channel_names[c])!=std::string::npos ){
				double * tbins =&bin_edges[c][0];
				std::string thisname;
				if(which_universe<0){
					thisname = fn;
				}else{
					thisname = fn+"_MS"+std::to_string(which_universe);
				}
				TH1D thischan(thisname.c_str(),"",num_bins[c], tbins );
				hist.push_back(thischan);
				//auto it = hist.begin()+ctr;
				//map_hist[fn] = &(*it);
				map_hist[fn] = ctr;

				ctr++;
			}

		}
	}

	has_been_scaled = false;

	this->collapseVector();	

}

SBNspec::SBNspec(std::string whichxml): SBNspec(whichxml,-1,true){}
SBNspec::SBNspec(std::string whichxml, int which_universe): SBNspec(whichxml,which_universe, true){}

SBNspec::SBNspec(std::string rootfile, std::string whichxml, bool isverbose) : SBNconfig(whichxml, isverbose) {
         //Contruct from a prexisting histograms that exist in a rootfile
         TFile *f = new TFile(rootfile.c_str(),"read");
 
         //Loop over all filenames that should be there, and load up the histograms.
         for(auto fn: fullnames){
                 hist.push_back(*((TH1D*)f->Get(fn.c_str()))); 
         }
         
	has_been_scaled=false;
 

         f->Close();
 
 
 }

SBNspec::SBNspec(std::string rootfile, std::string whichxml) : SBNspec(rootfile, whichxml, true){ };




 
int SBNspec::Add(SBNspec *in){
	//Addes all hists together
	if(xmlname != in->xmlname){ std::cout<<"ERROR: SBNspec::Add, trying to add differently configured SBNspecs!"<<std::endl; exit(EXIT_FAILURE);}

	for(int i=0; i< hist.size(); i++){
		hist[i].Add( &(in->hist[i]));
	}	

	this->collapseVector();
	return 0;
}

int SBNspec::setAsGaussian(double mean, double sigma, int ngen){
	TRandom3 *seedgetter = new TRandom3(0);
	int seed = seedgetter->Integer(1000000);

	for(auto &h: hist){
		TRandom3 *rangen = new TRandom3(seed);
		h.Reset();
		for(int i=0; i<ngen; i++){
			double eve = rangen->Gaus(mean,sigma); 
			h.Fill( eve ); 
		}	
	}

	return 0;

}

int SBNspec::setAsFlat(double val){
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, val );
		}	
	}
}



//All scaling functions are quite self explanatory
int SBNspec::poissonScale(){
	TRandom3 *rangen = new TRandom3(0);
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
		}	
	}
	return 0;
}

int SBNspec::poissonScale(TRandom3* rangen){
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
		}	
	}
	return 0;
}



int SBNspec::randomScale(){
	TRandom3 *rangen    = new TRandom3(0);

	for(auto& h: hist){
		h.Scale(rangen->Uniform(0,2));

	}
	return 0;
}


int SBNspec::Scale(std::string name, TF1 * func){
	for(auto& h: hist){
		std::string test = h.GetName();
		if(test.find(name)!=std::string::npos ){
			for(int b=0; b<=h.GetNbinsX(); b++){
				//std::cout<<h.GetBinContent(b)<<" "<<h.GetBinCenter(b)<<" "<<func->Eval(h.GetBinCenter(b) )<<std::endl; 
				h.SetBinContent(b, h.GetBinContent(b)*func->Eval(h.GetBinCenter(b) ) );
			}
		}

	}

	return 0;
}


int SBNspec::ScaleAll(double sc){
	for(auto& h: hist){
		h.Scale(sc);
	}
	this->collapseVector();

	return 0;
}

int SBNspec::Scale(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();

		if(test.find(name)!=std::string::npos ){
			//	std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
			h.Scale(val);
		}

	}

	has_been_scaled = true;
	scale_hist_name =name;
	scale_hist_val = val;

	this->collapseVector();
	return 0;
}

int SBNspec::NormAll(double n){

	for(auto& h: hist) {
		h.Scale(n/h.GetSumOfWeights());
	}
	return 0;
}

int SBNspec::Norm(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();

		if(test.find(name)!=std::string::npos ){
			//std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
			h.Scale(val/h.GetSumOfWeights());
		}

	}
	return 0;
}

int SBNspec::calcFullVector(){
	fullVec.clear();

	for(auto& h: hist){
		//std::cout<<"Hist size: "<<h.GetSize()-2<<std::endl;
		for(int i = 1; i <= h.GetSize()-2; i++){
			//std::cout<<h.GetBinContent(i)<<" ";
			fullVec.push_back(h.GetBinContent(i));
		}	
	}

	return 0;
}

int SBNspec::collapseVector(){

	compVec.clear();
	calcFullVector();

	for(int im = 0; im < num_modes; im++){
		for(int id =0; id < num_detectors; id++){
			int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector

			for(int ic = 0; ic < num_channels; ic++){
				int corner=edge;

				for(int j=0; j< num_bins.at(ic); j++){

					double tempval=0;

					for(int sc = 0; sc < num_subchannels.at(ic); sc++){

						//std::cout<<im<<"/"<<num_modes<<" "<<id<<"/"<<num_detectors<<" "<<ic<<"/"<<num_channels<<" "<<j<<"/"<<num_bins[ic]<<" "<<sc<<"/"<<num_subchannels[ic]<<std::endl;
						tempval += fullVec.at(j+sc*num_bins.at(ic)+corner);
						edge +=1;	//when your done with a channel, add on every bin you just summed
					}
					compVec.push_back(tempval);
				}






			}
		}
	}
	return 0;
}

double SBNspec::getTotalEvents(){
	double ans =0;
	this->calcFullVector();
	
	for(double d: fullVec){
		ans+=d;
	}

	return ans;


}

int SBNspec::printFullVec(){
	for(double d: fullVec){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
	return 0;
}

int SBNspec::printCompVec(){ 
	for(double d: compVec){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
	return 0;
}


int SBNspec::writeOut(std::string tag){
	//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
	//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
	//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

	std::vector<int> mycol = {kGreen+1, kRed-7, kBlue-4, kOrange-3, kMagenta+1, kCyan-3,kYellow, kGreen-3 }; 				
	int colindex;
	TFile *f2 = new TFile((tag+".SBNspec.root").c_str(),"recreate"); 

	for(auto& h: hist){
		h.Write();
	}
	f2->Close();	


	TFile *f = new TFile(("SBNfit_spectrum_plots_"+tag+".root").c_str(),"RECREATE"); 
	f->cd();

	std::vector<TH1D> temp_hists = hist;


	for(int im = 0; im <mode_names.size(); im++){
	for(int id = 0; id <detector_names.size(); id++){
	for(int ic = 0; ic <channel_names.size(); ic++){

			
				std::string canvas_name = mode_names.at(im)+"_"+detector_names.at(id)+"_"+channel_names.at(ic);

				bool this_run = false;

				TCanvas* Cstack= new TCanvas(canvas_name.c_str(),canvas_name.c_str());
				Cstack->cd();
				THStack * hs 	   = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend legStack(0.59,0.59,0.89,0.89);
				int n=0;

				for(auto &h : temp_hists){
					std::string test = h.GetName();
					if(test.find(canvas_name)!=std::string::npos ){
						colindex++;
						if(colindex ==mycol.size()-1) colindex=0;

						std::ostringstream numberofevents;
						numberofevents << std::setprecision(6) << h.GetSumOfWeights();


						h.Sumw2(false);
						h.Scale(1,"width,nosw2");
											h.SetMarkerStyle(20);
						h.SetMarkerColor(mycol.at( colindex ) );
						h.SetFillColor(mycol.at(colindex));
						h.SetLineColor(kBlack);
						h.SetTitle(h.GetName());
						//h.Write();

						std::string legend_name = h.GetName()  + numberofevents.str();
						legStack.AddEntry(&h, legend_name.c_str() , "f");

						hs->Add(&h);
						n++;

						this_run=true;

					}
				}
				
				if(this_run){
					hs->Draw();
					hs->GetYaxis()->SetTitle(("Events/"+channel_units.at(ic)).c_str());
					hs->GetXaxis()->SetTitle(channel_units.at(ic).c_str());
					Cstack->Update();
					legStack.Draw();	
					Cstack->Write();
				}

			}
		}
	}

	f->Close();

	return 0;
}



int SBNspec::getLocalBinNumber(double invar, int which_hist)
{
	int localbin = hist.at(which_hist).GetXaxis()->FindBin(invar);
	double bin = localbin-1;

	if(localbin==0 || localbin > hist.at(which_hist).GetNbinsX() ){bin = -99;}
	return bin;
}


int SBNspec::getGlobalBinNumber(double invar, int which_hist)
{
	int localbin = hist.at(which_hist).GetXaxis()->FindBin(invar);
	double bin = localbin-1;

	for(int i=0; i<which_hist; i++){
		bin += hist.at(i).GetNbinsX();	
	}

	if(localbin==0 || localbin > hist.at(which_hist).GetNbinsX() ){bin = -99;}
	return bin;
}

