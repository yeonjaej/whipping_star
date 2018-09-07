#include "SBNcls.h"
using namespace sbn;

int SBNcls::setSamplePoisson(){
    which_sample = 0;
   return which_sample;
}

int SBNcls::setSampleCovariance(){
    which_sample = 1;
   return which_sample;
}

int SBNcls::calcCLS(int numMC, std::string tag){

	if(which_sample == 0){
		std::cout<<"SBNcls::calcCLS\t|| Running in Poission sampling mode!"<<std::endl;
	}
	 else if(which_sample ==1){ 
		std::cout<<"SBNcls::calcCLS\t|| Running in Covariance sampling mode!"<<std::endl;
	}
	std::vector<double> ven={0};
	
	std::vector<double>  vec_CLs;

	
	double N_H1 = H1->getTotalEvents();
	double N_H0 = H0->getTotalEvents();

	//step one, find median H1 
	TH1D H1_pdf;
	if(which_sample == 0) H1_pdf = chi.samplePoisson_varyInput(H1, numMC);
	else if(which_sample ==1) H1_pdf = chi.sampleCovariance_varyInput(H1, numMC);
	
	double sig1 = 0.5-(0.6827)/2.0;
	double sig2 = 0.5-(0.9545)/2.0;

	std::vector<double> prob_values = {1-sig2, 1-sig1, 0.5, sig1, sig2};
	std::vector<double> quantiles(prob_values.size());	
	H1_pdf.ComputeIntegral(); 
	H1_pdf.GetQuantiles(prob_values.size(), &quantiles[0], &prob_values[0]);
	std::cout<<" BKG Quantile: "<<quantiles[0]<<" "<<quantiles[1]<<" "<<quantiles[2]<<std::endl;

	//Gotten medians ..etc..

	//Now calculate the pvalues associated with those H1 variations. 
	TH1D H0_pdf;
	std::vector<double> pval = quantiles; 
	if(which_sample == 0){
	 	H0_pdf = chi.samplePoisson_varyInput(H0, numMC, &pval);
	}
	 else if(which_sample ==1){ 
		H0_pdf = chi.sampleCovariance_varyInput(H0, numMC, &pval);
	}


	//lets do CLs
	for(int p=0; p<pval.size();p++){
		vec_CLs.push_back(pval.at(p)/(1-prob_values.at(p)) );
	}

	std::cout<<"NUMBER @ scale : "<<"\t"<< pval.at(0)/(1-prob_values.at(0))<<"\tsig: "<<N_H0<<"\tbkg: "<<N_H1<<"\tsig_only "<<N_H0-N_H1<<"\t\ts/sqrt(s+b): "<<(N_H0-N_H1)/sqrt(N_H0)<<std::endl;

	TFile * fp = new TFile(("SBNfit_CLs_"+tag+".root").c_str(),"recreate");
	fp->cd();
	TCanvas *cp=new TCanvas();


	H0_pdf.SetStats(false);
	H1_pdf.SetStats(false);

	H0_pdf.Scale(1/H0_pdf.GetSumOfWeights());
	H1_pdf.Scale(1/H1_pdf.GetSumOfWeights());

	H0_pdf.SetLineColor(kRed-7);
	H1_pdf.SetLineColor(kBlue-4);
	H0_pdf.SetFillColor(kRed-7);
	H1_pdf.SetFillColor(kBlue-4);
	H0_pdf.SetFillStyle(3445);
	H1_pdf.SetFillStyle(3454);

	H0_pdf.Draw("hist");
		
	double maxval =std::max(H0_pdf.GetMaximum(),H1_pdf.GetMaximum());
	double minval = std::min( H0_pdf.GetBinContent(H0_pdf.FindFirstBinAbove(0)), H1_pdf.GetBinContent(H1_pdf.FindFirstBinAbove(0)));
	std::cout<<"SBNcls::calcCLS() || Minimum value: "<<minval<<" Maximum value: "<<maxval<<std::endl;
	H0_pdf.SetMinimum(minval);
	H0_pdf.SetMaximum(maxval*1.35);

	double minbin = std::min(H0_pdf.GetBinLowEdge(H0_pdf.FindFirstBinAbove(0))+H0_pdf.GetBinWidth(H0_pdf.FindFirstBinAbove(0)), H1_pdf.GetBinLowEdge(H1_pdf.FindFirstBinAbove(0))+H1_pdf.GetBinWidth(H1_pdf.FindFirstBinAbove(0)));
	double maxbin = std::max(H0_pdf.GetBinLowEdge(H0_pdf.FindLastBinAbove(0))+H0_pdf.GetBinWidth(H0_pdf.FindLastBinAbove(0)), H1_pdf.GetBinLowEdge(H1_pdf.FindLastBinAbove(0))+H1_pdf.GetBinWidth(H1_pdf.FindLastBinAbove(0)));

	H0_pdf.GetXaxis()->SetRangeUser(minbin,maxbin);

	std::vector<std::string> quantile_names = {"-2#sigma","-1#sigma","Median","+1#sigma","+2#sigma"};
	std::vector<int> cols ={kYellow-7, kGreen+1, kBlack,kGreen+1, kYellow-7};	
	
	for(int i=0; i< quantiles.size(); i++){	
		TLine *l = new TLine(quantiles.at(i),minval, quantiles.at(i),maxval*1.05);
		l->SetLineColor(cols.at(i));
		l->SetLineWidth(2);
	        TLatex * qnam = new TLatex();
   	        qnam->SetTextSize(0.045);
   		qnam->SetTextAlign(12);  //align at top
		qnam->SetTextAngle(-90);
 		qnam->DrawLatex(quantiles.at(i), maxval*1.3 ,quantile_names.at(i).c_str());
		l->Draw("same");
	
		TLatex * qvals = new TLatex();
		qvals->SetTextSize(0.03);
		qvals->SetTextAlign(32);
		std::string details =  ("#splitline{"+quantile_names.at(i)+"}{1-#beta(" +to_string_prec(1-prob_values.at(i),3) + ") #alpha("+ to_string_prec(pval.at(i),3) +") CL_{s}("+to_string_prec(vec_CLs.at(i),3)+")}");
		std::cout<<details<<std::endl;
		qvals->DrawLatexNDC(0.875, 0.2+i*0.1,details.c_str()  );
	}
	
	TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
	leg->SetLineWidth(0);
	leg->SetFillStyle(0);
	leg->AddEntry(&H0_pdf,"H_{0}","lf");
	leg->AddEntry(&H1_pdf,"H_{1}","lf");
	leg->Draw();

	H0_pdf.GetXaxis()->SetTitle("#chi^{2}");
	H0_pdf.GetYaxis()->SetTitle("PDF");

	H1_pdf.Draw("hist same");	
	cp->Write();	
	cp->SaveAs(("SBNfit_Cls_"+tag+".pdf").c_str(),"pdf");	
	fp->Close();
			
	return 0;
}
