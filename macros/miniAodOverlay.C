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
#include <utility>
#include <algorithm>

using namespace std;

/* Global variables defining input histograms */
const static string s8TeV_FILE_NAME   = "data/TTbar_8TeV.root"; 
const static string s13TeV_FILE_NAME  = "data/TTbar_13TeV.root"; 

class HistInfo {
	public:
		HistInfo(const string histname_, const int rebin_=1, const bool logscale_ = true) {
			name = histname_;
			rebin = rebin_;
			logScale = logscale_;
		};
		HistInfo(vector<TFile*> files_,
				const string fullPath_,
				const string titleText_,
				const vector<string> legTexts_,
				const int rebin_,
				const int logScale_) {
			files = files_;
			fullPath = fullPath_;
			titleText = titleText_;
			legTexts = legTexts_;
			rebin = rebin_;
			logScale = logScale_;
		};

		vector<TFile*> files;
		string name, fullPath, titleText, epsName; 
		vector<string> legTexts;
		int rebin; 
		vector<int> lineColors, markerColors, markerStyles;
		bool logScale;
		vector<TH1*> hists;
	
};



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
	//h->Sumw2();
	h->SetLineWidth(2);
	h->Rebin(rebin);
	h->Scale(1.0/h->Integral());
}

/**************************************
 * Main body. 
 * Open files, read histograms, 
 * set aesthetics, draw, and print
 *************************************/
void miniAodOverlay(HistInfo histInfo) 
{


	/* create new canvas and draw */
	TCanvas *c1 = new TCanvas();
	gStyle->SetOptStat(0);
	gPad->SetTickx();
	gPad->SetTicky();


	TLegend *leg  = new TLegend(0.6,0.75,0.9,0.9);
					cout << __LINE__ << endl;
	for (int i =0; i< (histInfo.files).size(); ++i)
	{
		TH1 *h = dynamic_cast<TH1*> ((histInfo.files.at(i))->Get( (histInfo.fullPath).c_str()));
		if (h == NULL) { cout << "Hist " << histInfo.fullPath << " not in " << (histInfo.files.at(i))->GetName() << endl; assert(false); }
		h->Print();

		h->Rebin(histInfo.rebin);
		h->SetLineWidth(2);
		h->SetLineColor(histInfo.lineColors.at(i));
		h->SetMarkerColor(histInfo.markerColors.at(i));
		h->SetTitle(histInfo.titleText.c_str());
		h->Scale(1.0/(1.0* h->Integral()));
		leg->AddEntry(h, histInfo.legTexts.at(i).c_str());
		
		histInfo.hists.push_back(h);

		if (i==0) h->Draw();
		else h->Draw("same");
	}
	gPad->SetLogy(histInfo.logScale);
	leg->Draw();
	gPad->Print(histInfo.epsName.c_str());

}

/**************************************
 * This is the entry point to the code.
 * This method is overloaded
 *************************************/
int miniAodOverlay()
{

	std::vector<std::pair<unsigned, unsigned> > vHtBins, vJetBins;
	const char* metTypeArr[] = {"metall","methad","metlep"};
	const vector<string> svMetType(metTypeArr, metTypeArr+3);       //all met/ all hadmet(fake-met)/ leptonic met (true met)
	vHtBins.push_back(make_pair(0,10000));
	//vHtBins.push_back(make_pair(500,1000));
	//vHtBins.push_back(make_pair(1000,2000));
	//vHtBins.push_back(make_pair(2000,10000));
	//vMetBins.push_back(make_pair(0,10));
	vJetBins.push_back(make_pair(0,40));
	//vJetBins.push_back(make_pair(3,5));
	//vJetBins.push_back(make_pair(6,7));
	//vJetBins.push_back(make_pair(8,40));

	vector<TFile*> vFiles;
	TFile *f8tev = new TFile("data/TTbar_8TeV.root");
	TFile *f13tev50bx = new TFile("data/TTbar_13TeV_PUS13_50bx.root");
	TFile *f13tev25bx = new TFile("data/TTbar_13TeV_PU20_25bx.root");

	if (! f8tev->IsOpen()) { cout << "8tev file not found!" << endl; return 1; }
	if (! f13tev50bx->IsOpen()) { cout << "13tev50bx file not found!" << endl; return 1; }
	if (! f13tev25bx->IsOpen()) { cout << "13tev25bx file not found!" << endl; return 1; }
	
	vFiles.push_back(f8tev);
	vFiles.push_back(f13tev50bx);
	vFiles.push_back(f13tev25bx);

	vector<string> vLegTexts;
	vLegTexts.push_back("8 TeV");
	vLegTexts.push_back("13 TeV 50bx");
	vLegTexts.push_back("13 TeV 25bx");
	

	vector<int> lineCols, markerCols ,markerStls;
	lineCols.push_back(2);
	lineCols.push_back(3);
	lineCols.push_back(4);
	markerCols = lineCols;
	markerStls.push_back(12);
	markerStls.push_back(13);
	markerStls.push_back(14);


	vector<HistInfo> vHist;
	vHist.push_back(HistInfo("ht", 2, 1));
	vHist.push_back(HistInfo("met", 2));
	vHist.push_back(HistInfo("metphi", 2));
	//vHist.push_back(HistInfo("mht", 1));
	//vHist.push_back(HistInfo("mhtphi", 1));
	vHist.push_back(HistInfo("jet1_pt", 2));
	vHist.push_back(HistInfo("jet1_eta", 1));
	vHist.push_back(HistInfo("jet1_phi", 2));
	vHist.push_back(HistInfo("jet1_dphimet", 2));

	for (unsigned mt =0;  mt< svMetType.size(); ++mt)
	{
		for (unsigned jetbin=0; jetbin< vJetBins.size(); ++jetbin)
		{
			const int jlo = vJetBins.at(jetbin).first;
			const int jhi = vJetBins.at(jetbin).second;

			for (unsigned htbin=0; htbin< vHtBins.size(); ++htbin)
			{
				const int htlo = vHtBins.at(htbin).first;
				const int hthi = vHtBins.at(htbin).second;

				for (int ihist = 0; ihist < vHist.size(); ++ihist)
				{
					string metType = "", forEPSname="";
					if (mt == 0) { metType += "All MET Types"; forEPSname += "all";}
					else if (mt == 1) { metType += "Had MET Only"; forEPSname += "had";}
					else if (mt == 2) { metType += "Leptonic MET Only"; forEPSname += "lep";}

					stringstream path, title, epsname;
					path << svMetType.at(mt) <<"/Jets"<< jlo <<"to" << jhi << "HT" << htlo << "to"  << hthi << "/" << vHist.at(ihist).name;
					title << metType << ", Jets[" << jlo <<"," << jhi << "], HT[" << htlo << "," << hthi << "]";
					epsname << vHist.at(ihist).name << "_" << forEPSname << "_jets" << jlo <<"to" << jhi << "ht" << htlo << "to" << hthi << ".eps";

					vHist.at(ihist).files = vFiles;
					vHist.at(ihist).epsName = epsname.str();
					vHist.at(ihist).fullPath = path.str();
					vHist.at(ihist).titleText = title.str();
					vHist.at(ihist).legTexts = vLegTexts;
					vHist.at(ihist).lineColors = lineCols;
					vHist.at(ihist).markerColors = markerCols;
					vHist.at(ihist).markerStyles = markerStls;
					miniAodOverlay(vHist.at(ihist));

				}
			}
		}
	}

	bool logScale = 1; int rebin = 2;
	/*miniAodOverlay("jet1_pt", "nLeadPt"  , "P_{T} of Lead Jet^{P_{T}>50 & |#eta | < 2.5};P_{T} [GeV];Fraction of Events;"    , rebin, "Jet1Pt.eps", logScale);
	miniAodOverlay("jet2_pt", "nSecondPt", "P_{T} of 2nd Lead Jet^{P_{T}>50 & |#eta |<2.5};P_{T} [GeV];Fraction of Events;", rebin, "Jet2Pt.eps", logScale);
	miniAodOverlay("jet3_pt", "nThirdPt" , "P_{T} of 3rd Lead Jet^{P_{T}>50 & |#eta |<2.5};P_{T} [GeV];Fraction of Events;", rebin, "Jet3Pt.eps", logScale);
	logScale = 0; rebin = 2;
	miniAodOverlay("jet1_phi", "nLeadPhi"  , "#Phi of Lead Jet^{P_{T}>50 & |#eta |<2.5};#Phi;Fraction of Events;"    , rebin, "Jet1Phi.eps", logScale);
	miniAodOverlay("jet2_phi", "nSecondPhi", "#Phi of 2nd Lead Jet^{P_{T}>50 & |#eta |<2.5};#Phi;Fraction of Events;", rebin, "Jet2Phi.eps", logScale);
	miniAodOverlay("jet3_phi", "nThirdPhi" , "#Phi of 3rd Lead Jet^{P_{T}>50 & |#eta |<2.5};#Phi;Fraction of Events;", rebin, "Jet3Phi.eps", logScale);
	rebin = 2;
	miniAodOverlay("jet1_eta", "nLeadEta"  , "#eta of Lead Jet^{P_{T}>50 & |#eta | < 2.5};#eta;Fraction of Events;"    , rebin, "Jet1Eta.eps", logScale);
	miniAodOverlay("jet2_eta", "nSecondEta", "#eta of 2nd Lead Jet^{P_{T}>50 & |#eta | < 2.5};#eta;Fraction of Events;", rebin, "Jet2Eta.eps", logScale);
	miniAodOverlay("jet3_eta", "nThirdEta" , "#eta of 3rd Lead Jet^{P_{T}>50 & |#eta | < 2.5};#eta;Fraction of Events;", rebin, "Jet3Eta.eps", logScale);
	
	logScale = 1; rebin = 1;
	miniAodOverlay("pfmet", "MetPt" , "pfMet (8TeV)vs slimmedMETs (13 TeV);MET [GeV];Fraction of Events;", rebin, "met.eps", logScale);
 	logScale = 0; rebin = 1;
//	miniAodOverlay("pfmet", "MetPhi" , "pfMet (8TeV)vs slimmedMETs (13 TeV);MET #Phi [GeV];Fraction of Events;", rebin, "metphi.eps", logScale);
//	miniAodOverlay("pfht", "ht" , "pfHt (8TeV)vs slimmedHTs (13 TeV);HT [GeV];Fraction of Events;", rebin, "ht.eps", logScale);
	logScale = 0, rebin = 0;
//	miniAodOverlay("njet50", "njet50" , ";Jet Multiplicity;# of Jets [P_{T}>50 & |#eta | < 2.5];Fraction of Events;", rebin, "njet50.eps", logScale);
//	miniAodOverlay("njet30", "njet30" , ";Jet Multiplicity;# of Jets [P_{T}>30 & |#eta | < 5.0];Fraction of Events;", rebin, "njet30.eps", logScale);
//	miniAodOverlay("pfht", "ht" , ";;Fraction of Events;", rebin, ".eps", logScale);

*/
	return 0;
}
