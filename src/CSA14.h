#ifndef __CSA14_H__
#define __CSA14_H__

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TChain.h"


class CSA14 {

	public:
		CSA14();
		~CSA14();


		struct Hist_t{
			TH1D* h1_nVtx;
			TH1D* h1_ht;
			TH1D* h1_mht;
			TH1D* h1_mhtphi;
			TH1D* h1_met;
			TH1D* h1_metphi;
			TH1D* h1_uncorrmetphi;
			TH1D* h1_nJets_pt30eta50;
			TH1D* h1_nJets_pt50eta25;
			TH1D* h1_jet1_pt;
			TH1D* h1_jet1_eta;
			TH1D* h1_jet1_phi;
			TH1D* h1_jet1_dphimet;
			TH1D* h1_jet1_dphimht;
			TH1D* h1_jet2_pt;
			TH1D* h1_jet2_eta;
			TH1D* h1_jet2_phi;
			TH1D* h1_jet2_dphimet;
			TH1D* h1_jet2_dphimht;
			TH1D* h1_jet3_pt;
			TH1D* h1_jet3_eta;
			TH1D* h1_jet3_phi;
			TH1D* h1_jet3_dphimet;
			TH1D* h1_jet3_dphimht;
		};
		typedef std::vector<std::vector<std::vector<Hist_t> > > nestHistVec;

		int Process(const std::string inputFileList, const std::string outFileName, const int evts);
		int countJets(const std::vector<TLorentzVector> &inputJets, const double minAbsEta, 
				const double maxAbsEta, const double minPt, const double maxPt);
		int countJets(const std::vector<TLorentzVector> &inputJets, const double *jetCutsArr);
		std::vector<double> calcDPhi(const std::vector<TLorentzVector> &inputJets, const double metphi, 
				const int nDPhi, const double minAbsEta, const double maxAbsEta, 
				const double minPt, const double maxPt);
		std::vector<double> calcDPhi(const std::vector<TLorentzVector> &inputJets, 
				const double metphi, const int nDPhi, const double *jetCutsArr);
		void BookHistograms(const std::string outFileName, nestHistVec& vHist);
		unsigned GetVectorIndex(const std::vector< std::pair<unsigned, unsigned> >& binEdges, const unsigned& val);


	private:

		TFile *outFile;
		TChain *chain;
		std::vector<std::pair<unsigned, unsigned> > vHtBins, vJetBins;
		nestHistVec vHist;

};
#endif
