#include "SBNcls.h"
using namespace sbn;




int SBNcls::calcCLS( int numMC){


	std::vector<double> ven={0};
	
	std::vector<std::vector<double>> vec_CLs;

	for(int i=0; i<5;i++){
		vec_CLs.push_back(ven);
	}
	
	double N_H1 = H1->getTotalEvents();
	double N_H0 = H0->getTotalEvents();

	//step one, find median H1 
	TH1D H1_pdf = chi.toyMC_varyInput(H1, numMC);

	std::vector<double> prob_values = {0.05, 1-0.68, 0.5, 0.68, 0.95};
	std::vector<double> quantiles(prob_values.size());	
	H1_pdf.ComputeIntegral(); 
	H1_pdf.GetQuantiles(prob_values.size(), &quantiles[0], &prob_values[0]);
	std::cout<<" BKG Quantile: "<<quantiles[0]<<" "<<quantiles[1]<<" "<<quantiles[2]<<std::endl;

	//Gotten medians ..etc..

	//Now calculate the pvalues associated with those H1 variations. Definitely some double calculating here can optimize later.
	std::vector<double> pval = chi.toyMC_varyInput_getpval(H0, numMC, quantiles);
		
	//lets do CLs
	for(int p=0; p<pval.size();p++){
		vec_CLs.at(p).push_back(pval.at(p)/(1-prob_values.at(p)) );
	}

	std::cout<<"NUMBER @ scale : "<<"\t"<< pval.at(0)/(1-prob_values.at(0))<<"\tsig: "<<N_H0<<"\tbkg: "<<N_H1<<"\tsig_only "<<N_H0-N_H1<<"\t\ts/sqrt(s+b): "<<(N_H0-N_H1)/sqrt(N_H0)<<std::endl;

	TFile * fp = new TFile("printed2.root","recreate");
	fp->cd();
	TCanvas *cp=new TCanvas();

	H1_pdf.SetLineColor(kBlue-6);
	TH1D H0_pdf=chi.toyMC_varyInput(H0,numMC);
	
	H0_pdf.SetStats(false);
	H1_pdf.SetStats(false);

	H0_pdf.Scale(1/H0_pdf.GetSumOfWeights());
	H1_pdf.Scale(1/H1_pdf.GetSumOfWeights());

	H0_pdf.SetLineColor(kRed-6);
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
	
	
	for(int i=0; i< quantiles.size(); i++){	
		TLine *l = new TLine(quantiles.at(i),minval, quantiles.at(i),maxval*1.05);
	        TLatex * qnam = new TLatex();
   	        qnam->SetTextSize(0.05);
   		qnam->SetTextAlign(13);  //align at top
		qnam->SetTextAngle(-90);
 		qnam->DrawLatex(quantiles.at(i), maxval*1.3 ,quantile_names.at(i).c_str());
		l->Draw("same");
	}
	H0_pdf.GetXaxis()->SetTitle("#chi^{2}");
	H0_pdf.GetYaxis()->SetTitle("PDF");


	H1_pdf.Draw("hist same");	
	cp->Write();	
	fp->Close();
			
//			std::cout<<"PRINTED: Chi2: \t\t"<<x.at(0)<<"\t\t"<<x.at(1)<<"\t\t"<<x.at(2)<<"\t\t"<<x.at(3)<<"\t\t"<<x.at(4)<<std::endl;
//			std::cout<<"PRINTED: Power: \t\t"<<q.at(0)<<"\t\t"<<q.at(1)<<"\t\t"<<q.at(2)<<"\t\t"<<q.at(3)<<"\t\t"<<q.at(4)<<std::endl;
//			std::cout<<"PRINTED: Pval: \t\t"<<pval.at(0)<<"\t\t"<<pval.at(1)<<"\t\t"<<pval.at(2)<<"\t\t"<<pval.at(3)<<"\t\t"<<pval.at(4)<<std::endl;

	return 0;
}
