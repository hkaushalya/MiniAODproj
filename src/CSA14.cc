/* directive for the compiler so we can retrive high level objects stuffed inside
 * std::vector s in our tree 
 * You need these directives if you are running interactively.
 * If you are building the executable and running that this section
 * (from ifdef to endif)
 * can be safely removed (because this is already in mydict.cxx)
 */
#ifdef __MAKECINT__
#pragma link C++ class vector<int>;
#pragma link C++ class vector<vector<int> >;
#pragma link C++ class vector<vector<vector<int> > >;
#pragma link C++ class vector<TLorentzVector>;
#pragma link C++ class vector<vector<TLorentzVector> >;
#endif

/* all header files goes here */

#include <assert.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TFile.h"


/* define some gloable variables
 * needs cleaning up. this level of complexity
 * is not needed for this
 */

using namespace std;

static const bool verbose = true;
static const double      pt30Arr[] = {   -1,        -1,      30,    -1  };
static const double pt30Eta24Arr[] = {   -1,       2.4,      30,    -1  };
static const double pt50Eta24Arr[] = {   -1,       2.4,      50,    -1  };
static const double pt70Eta24Arr[] = {   -1,       2.4,      70,    -1  };
static const double      dphiArr[] = {   -1,       4.7,      30,    -1  };
static const double      bTagArr[] = {   -1,       2.4,      30,    -1  };
static const int nJetsSel = 5, nJetsSelPt30Eta24 = 5, nJetsSelPt50Eta24 = 4, nJetsSelPt70Eta24 = 2;
const double minMETcut = 200, maxMETcut = -1;



/* DEFINE all funtion prototypes here and implement anywhere in the file
 */

int countJets(const vector<TLorentzVector> &inputJets, const double minAbsEta, 
				const double maxAbsEta, const double minPt, const double maxPt);
int countJets(const vector<TLorentzVector> &inputJets, const double *jetCutsArr);
vector<double> calcDPhi(const vector<TLorentzVector> &inputJets, const double metphi, 
		const int nDPhi, const double minAbsEta, const double maxAbsEta, 
		const double minPt, const double maxPt);
vector<double> calcDPhi(const vector<TLorentzVector> &inputJets, 
		const double metphi, const int nDPhi, const double *jetCutsArr);



vector<TH1D*> h1_metVec, h1_metphiVec;
vector<TH1D*> h1_met_allhadVec, h1_metphi_allhadVec;
vector<TH1D*> h1_met_leptonicVec, h1_metphi_leptonicVec;

vector<TH1D*> h1_nJetsVec, h1_nJets_allhadVec, h1_nJets_leptonicVec;
vector<TH1D*> h1_vtxSizeVec;


int countJets(const vector<TLorentzVector> &inputJets, const double minAbsEta, const double maxAbsEta, const double minPt, const double maxPt){

	int cntNJets =0;
	for(unsigned int ij=0; ij<inputJets.size(); ij++){
		double perjetpt = inputJets[ij].Pt(), perjeteta = inputJets[ij].Eta();
		if(   ( minAbsEta == -1 || fabs(perjeteta) >= minAbsEta )
				&& ( maxAbsEta == -1 || fabs(perjeteta) < maxAbsEta )
				&& (     minPt == -1 || perjetpt >= minPt )
				&& (     maxPt == -1 || perjetpt < maxPt ) ){
			cntNJets ++;
		}
	}
	return cntNJets;
}

int countJets(const vector<TLorentzVector> &inputJets, const double *jetCutsArr){
	return countJets(inputJets, jetCutsArr[0], jetCutsArr[1], jetCutsArr[2], jetCutsArr[3]);
}
vector<double> calcDPhi(const vector<TLorentzVector> &inputJets, const double metphi, 
		const int nDPhi, const double minAbsEta, const double maxAbsEta, 
		const double minPt, const double maxPt){

	int cntNJets =0;
	vector<double> outDPhiVec(nDPhi, 999);
	for(unsigned int ij=0; ij<inputJets.size(); ij++){
		double perjetpt = inputJets[ij].Pt(), perjeteta = inputJets[ij].Eta();
		if(   ( minAbsEta == -1 || fabs(perjeteta) >= minAbsEta )
				&& ( maxAbsEta == -1 || fabs(perjeteta) < maxAbsEta )
				&& (     minPt == -1 || perjetpt >= minPt )
				&& (     maxPt == -1 || perjetpt < maxPt ) ){

			if( cntNJets < nDPhi ){
				double perDPhi = fabs(TVector2::Phi_mpi_pi( inputJets[ij].Phi() - metphi ));
				outDPhiVec[cntNJets] = perDPhi;
			}
			cntNJets ++;
		}
	}

	return outDPhiVec;
}

vector<double> calcDPhi(const vector<TLorentzVector> &inputJets, 
		const double metphi, const int nDPhi, const double *jetCutsArr){
	return calcDPhi(inputJets, metphi, nDPhi, jetCutsArr[0], jetCutsArr[1], jetCutsArr[2], jetCutsArr[3]);
}

// need to check if trees with different cycles gets attached or not!!
//Long64_t LoadTree(Long64_t entry);


int main(int argc, char* argv[])
{

	/* read user inputs and validate them */
	if (argc <= 3) {
		cerr <<"Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "# of evts  nJet50Min  nJet50Max" << endl;
		cerr <<" Valid configurations are " << std::endl;
		cerr << " ./optimize filelist.txt out.root nevts2process[-1 = all] smearingSyst[0=mean, 1-6 systs]" << std::endl;
		cerr << "Eg:  ./optimize filelist.txt out.root 100" << std::endl;
		return -1;
	}

	const char *inputFileList = argv[1];
	const char *outFileName   = argv[2];
	const char *evt2Proc      = argv[3];

	int evts    = -1;
	int systematic_var = 0;

	if (isdigit(evt2Proc[0]))
	{
		int num = atoi(argv[3]);
		if (num>=1) evts = num;
		if (verbose) cout << __FUNCTION__ << ": evts = " << evts << endl;
	} else {
		cout << "argument 4 is not a number. using default value for evts = " << evts << endl;
	}

	if (argc>4)
	{
		const char *g4 = argv[4];
		if (isdigit(g4[0])) systematic_var = atoi(g4);
		else {
			cout << "argument 4 is not a number. using default value for systematic_var = " << systematic_var << endl;
		}
	}


	/* check the input file and read in the contents and attach all the files to TChain (almost like TTree) */
	ifstream infile(inputFileList, ifstream::in);
	std::string buffer;

	if(!infile.is_open()) {
		std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
		return kFALSE;
	}

	/* this is where you tell ROOT which tree to retrive (if there is 
	 * more than one with the same name in different folders. Like in Hongxuan's
	 * files. Then loop over the inputFile list and attach each to TChain
	 * which you can use as you would interact with any TTree object
	 * (of course TChain has its own features)
	 */

	TChain* chain = new TChain("stdStop_histAndTree/AUX");
	while(1) {
		infile >> buffer;
		if(!infile.good()) break;
		if (verbose) std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
		TFile *f = new TFile(buffer.c_str());
		if (f->IsZombie()) { cout << buffer << " file not found!" << endl; }
		delete f;
		chain->Add(buffer.c_str());
	}
	if (verbose) chain->Print();
	std::cout << "No. of Entries in this tree = " << chain->GetEntries() << std::endl;


	/* here you define and allocate memory to store the retrived information (from tree)
	 */

	int run, event, lumi, nm1, n0, np1, vtxSize;
	double avg_npv, tru_npv;
	int nJets;
	double evtWeight;
	double met, metphi;

	int nMuons, nElectrons;

	vector<TLorentzVector> *oriJetsVec = new vector<TLorentzVector>(); vector<double> *recoJetsBtagCSVS = new vector<double>();

	vector<TLorentzVector> *genDecayLVec =0;
	vector<int> *genDecayPdgIdVec =0, *genDecayIdxVec =0, *genDecayMomIdxVec =0;
	vector<string> *genDecayStrVec =0, *smsModelFileNameStrVec =0, *smsModelStrVec =0;
	double smsModelMotherMass, smsModelDaughterMass;

	vector<TLorentzVector> *genJetsLVec_myak5GenJetsNoNu =0, *genJetsLVec_myak5GenJetsNoNuNoStopDecays =0, *genJetsLVec_myak5GenJetsNoNuOnlyStopDecays =0;


	/*now tell ROOT about those memory space */

	//chain->SetBranchStatus("*", 0); 	//not necessay, unless you want to read 
													//	hanful of branches and speed up the
													//	processing time.

	chain->SetBranchStatus("run", 1); chain->SetBranchAddress("run", &run);
	chain->SetBranchStatus("event", 1); chain->SetBranchAddress("event", &event);
	chain->SetBranchStatus("lumi", 1); chain->SetBranchAddress("lumi", &lumi);
	chain->SetBranchStatus("avg_npv", 1); chain->SetBranchAddress("avg_npv", &avg_npv);
	chain->SetBranchStatus("nm1", 1); chain->SetBranchAddress("nm1", &nm1);
	chain->SetBranchStatus("n0", 1); chain->SetBranchAddress("n0", &n0);
	chain->SetBranchStatus("np1", 1); chain->SetBranchAddress("np1", &np1);
	chain->SetBranchStatus("tru_npv", 1); chain->SetBranchAddress("tru_npv", &tru_npv);
	chain->SetBranchStatus("vtxSize", 1); chain->SetBranchAddress("vtxSize", &vtxSize);
	chain->SetBranchStatus("nJets", 1); chain->SetBranchAddress("nJets", &nJets);
	chain->SetBranchStatus("jetsLVec", 1); chain->SetBranchAddress("jetsLVec", &oriJetsVec);
	chain->SetBranchStatus("recoJetsBtag_0", 1); chain->SetBranchAddress("recoJetsBtag_0", &recoJetsBtagCSVS);
	chain->SetBranchStatus("evtWeight", 1); chain->SetBranchAddress("evtWeight", &evtWeight);
	chain->SetBranchStatus("met", 1); chain->SetBranchAddress("met", &met);
	chain->SetBranchStatus("metphi", 1); chain->SetBranchAddress("metphi", &metphi);

	chain->SetBranchStatus("nMuons", 1); chain->SetBranchAddress("nMuons", &nMuons);
	chain->SetBranchStatus("nElectrons", 1); chain->SetBranchAddress("nElectrons", &nElectrons);

	chain->SetBranchStatus("genDecayLVec", 1); chain->SetBranchAddress("genDecayLVec", &genDecayLVec);
	chain->SetBranchStatus("genDecayPdgIdVec", 1); chain->SetBranchAddress("genDecayPdgIdVec", &genDecayPdgIdVec);
	chain->SetBranchStatus("genDecayIdxVec", 1); chain->SetBranchAddress("genDecayIdxVec", &genDecayIdxVec);
	chain->SetBranchStatus("genDecayMomIdxVec", 1); chain->SetBranchAddress("genDecayMomIdxVec", &genDecayMomIdxVec);
	chain->SetBranchStatus("genDecayStrVec", 1); chain->SetBranchAddress("genDecayStrVec", &genDecayStrVec);

	/*quick sanity check to count how many entries are in ALL the trees in the TChain 
	 */
	int Entries = chain->GetEntries();
	std::cout<<"\n\n"<<"Entries : "<<Entries<<std::endl;


	/* create output file and add histograms to it */
	TFile *outFile = new TFile(outFileName, "RECREATE");

	TH1D *h1_met = new TH1D("h1_met", ":  met; met", 100, 0, 1000); h1_met->Sumw2();
	TH1D *h1_metphi = new TH1D("h1_metphi", ":  metphi; metphi", 100, -3.2, 3.2); h1_metphi->Sumw2();

	TH1D *h1_met_allhad = new TH1D("h1_met_allhad", ":  met_allhad; met_allhad", 100, 0, 1000); h1_met_allhad->Sumw2();
	TH1D *h1_metphi_allhad = new TH1D("h1_metphi_allhad", ":  metphi_allhad; metphi_allhad", 100, -3.2, 3.2); h1_metphi_allhad->Sumw2();

	TH1D *h1_met_leptonic = new TH1D("h1_met_leptonic", ":  met_leptonic; met_leptonic", 100, 0, 1000); h1_met_leptonic->Sumw2();
	TH1D *h1_metphi_leptonic = new TH1D("h1_metphi_leptonic", ":  metphi_leptonic; metphi_leptonic", 100, -3.2, 3.2); h1_metphi_leptonic->Sumw2();

	TH1D *h1_nJets = new TH1D("h1_nJets", ":  number of jets; nJets", 20, 0, 20); h1_nJets->Sumw2();
	TH1D *h1_nJets_allhad = new TH1D("h1_nJets_allhad", ":  number of jets; nJets_allhad", 20, 0, 20); h1_nJets_allhad->Sumw2();
	TH1D *h1_nJets_leptonic = new TH1D("h1_nJets_leptonic", ":  number of jets; nJets_leptonic", 20, 0, 20); h1_nJets_leptonic->Sumw2();
	TH1D *h1_vtxSize = new TH1D("h1_vtxSize", ":  number of vertices; vtxSize", 100, 0, 100); h1_vtxSize->Sumw2();



	for(int ie=0; ie<Entries; ie++){

		chain->GetEntry(ie);

		if( ie==0 || ie == Entries-1 || ie%(Entries/10) == 0 ) std::cout<<"\n   Processing the "<<ie<<"th event ..."<<std::endl;

		if( oriJetsVec->size() != recoJetsBtagCSVS->size() ){
			std::cout<<"oriJetsVec->size : "<<oriJetsVec->size()<<"  recoJetsBtagCSVS : "<<recoJetsBtagCSVS->size()<<std::endl;
		}

		TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);

		int cntNJetsPt30 = countJets((*oriJetsVec), pt30Arr);
		int cntNJetsPt30Eta24 = countJets((*oriJetsVec), pt30Eta24Arr);
		int cntNJetsPt50Eta24 = countJets((*oriJetsVec), pt50Eta24Arr);
		int cntNJetsPt70Eta24 = countJets((*oriJetsVec), pt70Eta24Arr);
		vector<double> dPhiVec = calcDPhi((*oriJetsVec), metphi, 3, dphiArr);

		double dPhi0 = dPhiVec[0], dPhi1 = dPhiVec[1], dPhi2 = dPhiVec[2];

		bool passExtraCuts = true;
		bool passnJets = true, passdPhis = true, passBJets = true, passMET = true;

		if( cntNJetsPt70Eta24 < nJetsSelPt70Eta24 ){ passExtraCuts = false; passnJets = false; }
		if( cntNJetsPt50Eta24 < nJetsSelPt50Eta24 ){ passExtraCuts = false; passnJets = false; }
		if( cntNJetsPt30Eta24 < nJetsSelPt30Eta24 ){ passExtraCuts = false; passnJets = false; }

		if( dPhi0 < 0.5 || dPhi1 < 0.5 || dPhi2 < 0.3 ){ passExtraCuts = false; passdPhis = false; }

		if( !( (minMETcut == -1 || met >= minMETcut ) && (maxMETcut == -1 || met < maxMETcut) ) ){ passExtraCuts = false; passMET = false; }

		// Parsing the gen information ...
		int cntgenTop = 0, cntleptons =0;

		for(int iv=0; iv<(int)genDecayLVec->size(); iv++){
			int pdgId = genDecayPdgIdVec->at(iv);
			if( abs(pdgId) == 6 ) cntgenTop ++;
			if( abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15 ) cntleptons++;
		}

		//h1_vtxSize->Fill(vtxSize);
		h1_nJets->Fill(cntNJetsPt30Eta24);
		h1_met->Fill(met);
		h1_metphi->Fill(metphi);

		if( cntleptons ==0 ){
			h1_nJets_allhad->Fill(cntNJetsPt30Eta24);
			h1_met_allhad->Fill(met);
			h1_metphi_allhad->Fill(metphi);
		}else{
			h1_nJets_leptonic->Fill(cntNJetsPt30Eta24);
			h1_met_leptonic->Fill(met);
			h1_metphi_leptonic->Fill(metphi);
		}

		if( verbose >=1 ){
			std::cout<<"\nie : "<<ie<<std::endl; 
			std::cout<<"genDecayStr : "<<genDecayStrVec->front().c_str()<<std::endl;
			for(int iv=0; iv<(int)genDecayLVec->size(); iv++){
				int pdgId = genDecayPdgIdVec->at(iv);
				printf("((%d,%d/%d):(%6.2f/%6.2f))  ", pdgId, genDecayIdxVec->at(iv), genDecayMomIdxVec->at(iv), genDecayLVec->at(iv).E(), genDecayLVec->at(iv).Pt());
			}
		}
	}



	/* clean up */
	outFile->Write();
	outFile->Close();
	delete chain;

	return 0;
}


