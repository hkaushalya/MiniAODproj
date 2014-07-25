#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>
#include "TStyle.h"
#include "assert.h"
#include "TF1.h"
#include "TLegend.h"
#include <iomanip>
#include <algorithm>

using namespace std;

/* Global variables defining input histograms */
const static string s8TeV_FILE_NAME   = "data/TTbar_8TeV.root"; 
const static string s13TeV_FILE_NAME  = "data/TTbar_13TeV.root"; 


/**************************************
 * Can be used look into a histogram
 *************************************/
void DumpHist(const TH1* hist)
{
	assert (hist != NULL && "DumpHist:: hist passed is null!"); 
	cout << std::setw(5) << "bin " << std::setw(15) << "[edges]" << std::setw(10) << "content" << std::setw(10) << "error" << endl;
	double sum = 0;
	for (int bin = 0; bin <= hist->GetNbinsX()+1; ++bin)
	{
		if (hist->GetBinContent(bin)>0)
		{
			cout  << std::setprecision(4) << std::setw(5) << bin << std::setw(15) << "[" 
			<< hist->GetBinLowEdge(bin) << ", " << hist->GetXaxis()->GetBinUpEdge(bin)<< "]" 
			<< std::setw(10) << hist->GetBinContent(bin) 
			<< std::setw(10) << hist->GetBinError(bin) << endl;
		}
	}
}

/**************************************
 * Scale and adjust bin sizes of hist
 *************************************/
void ScaleHist(TH1* h, const int rebin) {
	h->Sumw2();
	h->SetLineWidth(2);
	h->Rebin(rebin);
	h->Scale(1.0/h->Integral());
}

/**************************************
 * Main body. 
 * Open files, read histograms, 
 * set aesthetics, draw, and print
 *************************************/
void miniAodOverlay(const string histname, const string title, const int rebin=1, const string epsname="", const bool logScale =1) 
{
	TFile* rootFile8TeV = new TFile (s8TeV_FILE_NAME.c_str());
	if (rootFile8TeV->IsZombie())
	{
		cout << "8 TeV ROOT file not found!" <<  endl;
		assert (false);
	}

	TFile* rootFile13TeV = new TFile (s13TeV_FILE_NAME.c_str());
	if (rootFile13TeV->IsZombie())
	{
		cout << "13 TeV ROOT file not found!" <<  endl;
		assert (false);
	}

	/* hack to pick up the hist with difference names */
	stringstream s8TeVHistName, s13TeVHistName, s8TeVHistPath, s13TeVHistPath;
	if (histname == "Jet1Pt") {
		s8TeVHistName  << "jet1_pt";
		s13TeVHistName << "nLeadPt";
	} else if (histname == "Jet2Pt") {
		s8TeVHistName  << "jet2_pt";
		s13TeVHistName << "nSecondPt";
	} else if (histname == "Jet3Pt") {
		s8TeVHistName  << "jet3_pt";
		s13TeVHistName << "nThirdPt";
	} else  {
		cout << "ERR: Cannot find matching hist for your request of " << histname << "!"<< endl;
		assert(false);
	}

	s13TeVHistPath << "/CfiFile/" << s13TeVHistName.str() ;
	s8TeVHistPath  << "/Factorization/HT0to5000/" << s8TeVHistName.str();

	TH1* hist8TeV = dynamic_cast<TH1*> (rootFile8TeV->Get(s8TeVHistPath.str().c_str()));
	if (hist8TeV == NULL) { cout << "ERR: 8 TeV hist " << s8TeVHistPath.str() << " not found!" << endl; assert (false); }

	TH1* hist13TeV = dynamic_cast<TH1*> (rootFile13TeV->Get(s13TeVHistPath.str().c_str()));
	if (hist13TeV == NULL) { cout << "ERR: 13 TeV hist " << s13TeVHistPath.str() << " not found!" << endl; assert (false); }

	/* now normalize hist to 1.0 */
	ScaleHist(hist8TeV, rebin);
	ScaleHist(hist13TeV, rebin);

	/* make them pretty */
	const int i8TeVColor  = 4;
	const int i13TeVColor = 6;

	hist8TeV->SetLineColor(i8TeVColor);
	hist13TeV->SetLineColor(i13TeVColor);
	hist8TeV->SetTitle(title.c_str());

	/* add a legend */
	TLegend *leg1  = new TLegend(0.6,0.75,0.9,0.9);
	leg1->AddEntry(hist8TeV,  "t#bar{t} (8 TeV)");
	leg1->AddEntry(hist13TeV, "t#bar{t} (13 TeV)");

	/* create new canvas and draw */
	new TCanvas();
	gStyle->SetOptStat(0);
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();
	hist8TeV->Draw();
	hist13TeV->Draw("same");
	leg1->Draw();
	gPad->Print(epsname.c_str());

}

/**************************************
 * This is the entry point to the code.
 * This method is overloaded
 *************************************/
void miniAodOverlay()
{
	const bool logScale = 0;
	const int rebin = 4;
	miniAodOverlay("Jet1Pt", "Lead Jet (P_{T}>50 & |#eta | < 2.5);P_{T} [GeV];Fraction of Events;", rebin, "Jet1Pt.eps", logScale);
	miniAodOverlay("Jet2Pt", "2nd Lead Jet (P_{T}>50 && |#eta | < 2.5));P_{T} [GeV];Fraction of Events;", rebin, "Jet2Pt.eps", logScale);
	miniAodOverlay("Jet3Pt", "3rd Lead Jet (P_{T}>50 && |#eta | < 2.5));P_{T} [GeV];Fraction of Events;", rebin, "Jet3Pt.eps", logScale);
}
