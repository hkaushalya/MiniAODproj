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

#include "CSA14.h"
#include <assert.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <utility>
#include "TDirectory.h"
#include <sstream>
#include <utility>

/* define some gloable variables
 * needs cleaning up. this level of complexity
 * is not needed for this
 */

using namespace std;

static const bool verbose = false;
//static const double      pt30Arr[] = {   -1,        -1,      30,    -1  };
static const double pt30Eta50Arr[] = {   -1,       5.0,      30,    -1  };
static const double pt50Eta25Arr[] = {   -1,       2.5,      50,    -1  };
//static const double pt70Eta24Arr[] = {   -1,       2.4,      70,    -1  };
static const double      dphiArr[] = {   -1,       5.0,      30,    -1  };
//static const double      bTagArr[] = {   -1,       2.4,      30,    -1  };
static const int nJetsSel = 0, nJetsSelPt30Eta24 = 0, nJetsSelPt50Eta24 = 3, nJetsSelPt70Eta24 = 0;
const double minMETcut = 0, maxMETcut = -1;



/* DEFINE all funtion prototypes here and implement anywhere in the file
 */
/* be caucious when changing mettype! I am using absolute indices
 * when filling histograms and there is not sanity checks for vector size!!
 */
const char* metTypeArr[] = {"metall","methad","metlep"};
const vector<string> svMetType(metTypeArr, metTypeArr+3);       //all met/ all hadmet(fake-met)/ leptonic met (true met)

CSA14::CSA14()
{
	outFile = 0;
	chain = 0;
	vHtBins.push_back(make_pair(0,10000));
	vHtBins.push_back(make_pair(500,1000));
	vHtBins.push_back(make_pair(1000,2000));
	vHtBins.push_back(make_pair(2000,10000));
	//vMetBins.push_back(make_pair(0,10));
	vJetBins.push_back(make_pair(0,40));
	vJetBins.push_back(make_pair(3,5));
	vJetBins.push_back(make_pair(6,7));
	vJetBins.push_back(make_pair(8,40));

}

CSA14::~CSA14()
{
	/* clean up */
	outFile->Write();
	outFile->Close();
	delete chain;
	
}


int CSA14::countJets(const vector<TLorentzVector> &inputJets, const double minAbsEta, 
				const double maxAbsEta, const double minPt, const double maxPt){

	int cntNJets =0;
	for(unsigned int ij=0; ij<inputJets.size(); ij++){
		double perjetpt = inputJets[ij].Pt(), perjeteta = inputJets[ij].Eta();
		//cout << "pt/eta" << perjetpt << "/" << perjeteta << endl;
		if(   ( minAbsEta == -1 || fabs(perjeteta) >= minAbsEta )
				&& ( maxAbsEta == -1 || fabs(perjeteta) < maxAbsEta )
				&& (     minPt == -1 || perjetpt >= minPt )
				&& (     maxPt == -1 || perjetpt < maxPt ) ){
			cntNJets ++;
		}
	}
	return cntNJets;
}

int CSA14::countJets(const vector<TLorentzVector> &inputJets, const double *jetCutsArr){
	return countJets(inputJets, jetCutsArr[0], jetCutsArr[1], jetCutsArr[2], jetCutsArr[3]);
}
vector<double> CSA14::calcDPhi(const vector<TLorentzVector> &inputJets, const double metphi, 
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

vector<double> CSA14::calcDPhi(const vector<TLorentzVector> &inputJets, 
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
		cerr << " ./runcsa filelist.txt out.root nevts2process[-1 == all]" << std::endl;
		cerr << "Eg:  ./runcsa filelist.txt out.root 100" << std::endl;
		return -1;
	}

	const char *inputFileList = argv[1];
	const char *outFileName   = argv[2];
	const char *evt2Proc      = argv[3];


	if ((string(inputFileList)).length() < 1)
	{
		cout << "ERR! Argument 1 (inputFileList) is not a valid string!!" << endl;
		return 1;
	}

	if ((string(outFileName)).length() < 1)
	{
		cout << "ERR! Argument 2 (outputFileName) is not a valid string!!" << endl;
		return 1;
	}


	int evts    = -1;
	if (isdigit(evt2Proc[0]))
	{
		int num = atoi(argv[3]);
		if (num>=1) evts = num;
		if (verbose) cout << __FUNCTION__ << ": evts = " << evts << endl;
	} else {
		cout << "argument 3 is not a number. Processing all events." << endl;
	}

	CSA14 csa = CSA14();
	csa.Process(inputFileList, outFileName, evts);
}


int CSA14::Process(const string inputFileList, const string outFileName, const int evts)
{
	/* check the input file and read in the contents and attach all the files to TChain (almost like TTree) */
	ifstream infile(inputFileList.c_str(), ifstream::in);
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

	chain = new TChain("stdStop_histAndTree/AUX");
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


	/* here you define and allocate memory to store the retrived information (from tree)
	*/

	int run, event, lumi, nm1, n0, np1, vtxSize;
	double avg_npv, tru_npv;
	int nJets;
	double evtWeight;
	double ht, met, metphi, type1metphi, mht, mhtphi;
	int nMuons, nElectrons;
	vector<TLorentzVector> *oriJetsVec = new vector<TLorentzVector>(); 
	vector<double> *recoJetsBtagCSVS = new vector<double>();
	vector<TLorentzVector> *genDecayLVec =0;
	vector<int> *genDecayPdgIdVec =0, *genDecayIdxVec =0, *genDecayMomIdxVec =0;
	vector<string> *genDecayStrVec =0, *smsModelFileNameStrVec =0, *smsModelStrVec =0;
	double smsModelMotherMass, smsModelDaughterMass;
	vector<TLorentzVector> *genJetsLVec_myak5GenJetsNoNu =0, *genJetsLVec_myak5GenJetsNoNuNoStopDecays =0, *genJetsLVec_myak5GenJetsNoNuOnlyStopDecays =0;


	/*now tell ROOT about those memory space */

	chain->SetBranchStatus("*", 0); 	//not necessay, unless you want to read handful 
	//of branches and speed up the processing time.

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
	chain->SetBranchStatus("ht", 1); chain->SetBranchAddress("ht", &ht);
	chain->SetBranchStatus("met", 1); chain->SetBranchAddress("met", &met);
	chain->SetBranchStatus("metphi", 1); chain->SetBranchAddress("metphi", &metphi);
	chain->SetBranchStatus("type1metphi", 1); chain->SetBranchAddress("type1metphi", &type1metphi);  //uncorrected metphi
	chain->SetBranchStatus("mht", 1); chain->SetBranchAddress("mht", &mht);
	chain->SetBranchStatus("mhtphi", 1); chain->SetBranchAddress("mhtphi", &mhtphi);
	chain->SetBranchStatus("nMuons", 1); chain->SetBranchAddress("nMuons", &nMuons);
	chain->SetBranchStatus("nElectrons", 1); chain->SetBranchAddress("nElectrons", &nElectrons);
	chain->SetBranchStatus("genDecayLVec", 1); chain->SetBranchAddress("genDecayLVec", &genDecayLVec);
	chain->SetBranchStatus("genDecayPdgIdVec", 1); chain->SetBranchAddress("genDecayPdgIdVec", &genDecayPdgIdVec);
	chain->SetBranchStatus("genDecayIdxVec", 1); chain->SetBranchAddress("genDecayIdxVec", &genDecayIdxVec);
	chain->SetBranchStatus("genDecayMomIdxVec", 1); chain->SetBranchAddress("genDecayMomIdxVec", &genDecayMomIdxVec);
	chain->SetBranchStatus("genDecayStrVec", 1); chain->SetBranchAddress("genDecayStrVec", &genDecayStrVec);

	/*quick sanity check to count how many entries are in ALL the trees in the TChain 
	*/
	const int Entries = chain->GetEntries();
	std::cout << "No. of Entries in this tree = " << Entries << std::endl;
	const int iEVTS_TO_PROCESS = evts == -1 ? Entries: evts;

	BookHistograms(outFileName, vHist);

	for (int ie=0; ie< min(Entries, iEVTS_TO_PROCESS); ie++){
		chain->GetEntry(ie);
		if( ie==0 || ie == Entries-1 || ie%(max(iEVTS_TO_PROCESS,10)/10) == 0 ) 
			std::cout<<"\tProcessing the "<< setw(10) << ie<<"th event ...[~" << ceil(100*ie/(float)iEVTS_TO_PROCESS) << "%]"<<std::endl;

		if( oriJetsVec->size() != recoJetsBtagCSVS->size() ){
			std::cout<<"oriJetsVec->size : "<<oriJetsVec->size()<<"  recoJetsBtagCSVS : "<<recoJetsBtagCSVS->size()<<std::endl;
		}

		TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
		TLorentzVector mhtLVec; mhtLVec.SetPtEtaPhiM(mht, 0, mhtphi, 0);

		//cout << "\tentry " << ie << endl;
		const int cntNJetsPt30Eta50 = countJets((*oriJetsVec), pt30Eta50Arr);
		const int cntNJetsPt50Eta25 = countJets((*oriJetsVec), pt50Eta25Arr);
		const vector<double> dPhiVec_met = calcDPhi((*oriJetsVec), metphi, 3, dphiArr);
		const vector<double> dPhiVec_mht = calcDPhi((*oriJetsVec), mhtphi, 3, dphiArr);

		const double dPhi0_met = dPhiVec_met[0], dPhi1_met = dPhiVec_met[1], dPhi2_met = dPhiVec_met[2];
		const double dPhi0_mht = dPhiVec_mht[0], dPhi1_mht = dPhiVec_mht[1], dPhi2_mht = dPhiVec_mht[2];

		bool passExtraCuts = true;
		bool passnJets = true, passdPhis = true, passBJets = true, passMET = true;

		// Parsing the gen information ...
		int cntgenTop = 0, cntleptons =0;
		for(int iv=0; iv<(int)genDecayLVec->size(); iv++){
			int pdgId = genDecayPdgIdVec->at(iv);
			if( abs(pdgId) == 6 ) cntgenTop ++;
			if( abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15 ) cntleptons++;
		}

		//cout << "\ncntleptons/cntNJetsPt50Eta25/ht = " << cntleptons << " / " << cntNJetsPt50Eta25 << "/" << ht << endl;
		//cout << "njets/ht = " << cntNJetsPt50Eta25 << "/" << ht <<endl;
		for (int mt =0;  mt< svMetType.size(); ++mt)
		{
			if (mt==1 && cntleptons !=0) continue; //do not fill had-met with leptonic events
			if (mt==2 && cntleptons ==0) continue;	//do not fill lep-met with hadronic events

			for (int jetbin=0; jetbin< vJetBins.size(); ++jetbin)
			{
				const int jlo = vJetBins.at(jetbin).first;
				const int jhi = vJetBins.at(jetbin).second;
				if (cntNJetsPt50Eta25 < jlo || cntNJetsPt50Eta25 > jhi) continue;

				for (int htbin=0; htbin< vHtBins.size(); ++htbin)
				{
					const int htlo = vHtBins.at(htbin).first;
					const int hthi = vHtBins.at(htbin).second;
					if (ht<htlo || ht>hthi) continue;

					//cout << "\tFilling mtbin,jetbin,htbin:" << mt << ", " << jetbin << ", " << htbin << endl; 
					const int idx = mt+jetbin+htbin;
					vHist.at(mt).at(jetbin).at(htbin).h1_nVtx->Fill(vtxSize);
					vHist.at(mt).at(jetbin).at(htbin).h1_ht->Fill(ht);
					vHist.at(mt).at(jetbin).at(htbin).h1_nJets_pt30eta50->Fill(cntNJetsPt30Eta50);
					vHist.at(mt).at(jetbin).at(htbin).h1_nJets_pt50eta25->Fill(cntNJetsPt50Eta25);
					vHist.at(mt).at(jetbin).at(htbin).h1_met->Fill(met);
					vHist.at(mt).at(jetbin).at(htbin).h1_metphi->Fill(metphi);
					vHist.at(mt).at(jetbin).at(htbin).h1_mht->Fill(mht);
					vHist.at(mt).at(jetbin).at(htbin).h1_mhtphi->Fill(mhtphi);
					vHist.at(mt).at(jetbin).at(htbin).h1_uncorrmetphi->Fill(type1metphi);
					//Fill jet histograms
					if (oriJetsVec->size()>=1) {
						vHist.at(mt).at(jetbin).at(htbin).h1_jet1_pt->Fill(oriJetsVec->at(0).Pt());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet1_eta->Fill(oriJetsVec->at(0).Eta());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet1_phi->Fill(oriJetsVec->at(0).Phi());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet1_dphimet->Fill(dPhi0_met);
						vHist.at(mt).at(jetbin).at(htbin).h1_jet1_dphimht->Fill(dPhi0_mht);
					}
					if (oriJetsVec->size()>=2) {
						vHist.at(mt).at(jetbin).at(htbin).h1_jet2_pt->Fill(oriJetsVec->at(1).Pt());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet2_eta->Fill(oriJetsVec->at(1).Eta());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet2_phi->Fill(oriJetsVec->at(1).Phi());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet2_dphimet->Fill(dPhi1_met);
						vHist.at(mt).at(jetbin).at(htbin).h1_jet2_dphimht->Fill(dPhi1_mht);
					}
					if (oriJetsVec->size()>=3) {
						vHist.at(mt).at(jetbin).at(htbin).h1_jet3_pt->Fill(oriJetsVec->at(2).Pt());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet3_eta->Fill(oriJetsVec->at(2).Eta());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet3_phi->Fill(oriJetsVec->at(2).Phi());
						vHist.at(mt).at(jetbin).at(htbin).h1_jet3_dphimet->Fill(dPhi2_met);
						vHist.at(mt).at(jetbin).at(htbin).h1_jet3_dphimht->Fill(dPhi2_mht);
					}

				}
			}
		}




		/*		if( verbose >=1 ){
				std::cout<<"\nie : "<<ie<<std::endl; 
				std::cout<<"genDecayStr : "<<genDecayStrVec->front().c_str()<<std::endl;
				for(int iv=0; iv<(int)genDecayLVec->size(); iv++){
				int pdgId = genDecayPdgIdVec->at(iv);
				printf("((%d,%d/%d):(%6.2f/%6.2f))  ", pdgId, genDecayIdxVec->at(iv), genDecayMomIdxVec->at(iv), genDecayLVec->at(iv).E(), genDecayLVec->at(iv).Pt());
				}
				}
				*/
	}

		/********************** JOB SUMMARY *******************************/
		cout << "\n -------------- BEGIN SUMMARY -------------" << endl;
		cout << "  Processed file list = " << inputFileList << endl;
		cout << "  Output written to   = " << outFileName << endl;
		cout << "  Processed events    = " << iEVTS_TO_PROCESS <<  " ( out of " << Entries << ")" << endl;
		cout << " -------------- END  SUMMARY -------------\n\n" << endl;




		for (int mt =0;  mt< svMetType.size(); ++mt)
		{
			for (int jetbin=0; jetbin< vJetBins.size(); ++jetbin)
			{
				for (int htbin=0; htbin< vHtBins.size(); ++htbin)
				{
					cout << "\tmtbin,jetbin,htbin:" << mt << ", " << jetbin << ", " << htbin << endl; 
					vHist.at(mt).at(jetbin).at(htbin).h1_nVtx->Print();
					vHist.at(mt).at(jetbin).at(htbin).h1_met->Print();
				}
			}
		}
		return 0;
}


void CSA14::BookHistograms(const string outFileName, nestHistVec& vHist)
{
	vHist.clear();
	/* create output file and add histograms to it */
	outFile = new TFile(outFileName.c_str(), "RECREATE");
	for (int mt =0;  mt< svMetType.size(); ++mt)
	{
		outFile->cd();
		TDirectory *dir_0 = outFile->mkdir(svMetType.at(mt).c_str()); 
		dir_0->cd();

		vector <vector<Hist_t> > vJetBinnedHist; 
		for (int jetbin=0; jetbin< vJetBins.size(); ++jetbin)
		{
			const int jlo = vJetBins.at(jetbin).first;
			const int jhi = vJetBins.at(jetbin).second;

			vector<Hist_t> vHtBinnedHist; 
			for (int htbin=0; htbin< vHtBins.size(); ++htbin)
			{
				const int htlo = vHtBins.at(htbin).first;
				const int hthi = vHtBins.at(htbin).second;

				stringstream dirname;
				dirname <<"Jets"<< jlo <<"to" << jhi << "HT" << htlo << "to" << hthi;
				TDirectory *dir_1 = dir_0->mkdir(dirname.str().c_str()); 
				dir_1->cd();

				stringstream title_suffix;
				title_suffix << svMetType.at(mt) << ", Jets["<< jlo <<"," << jhi << "], HT[" << htlo << "," << hthi << "]";

				stringstream tit_ht, tit_mht, tit_mhtphi, tit_met, tit_metphi, tit_uncorrmetphi, tit_njetpt30, tit_njetpt50, tit_nvtx;
				stringstream tit_jet1_pt, tit_jet2_pt, tit_jet3_pt;
				stringstream tit_jet1_eta, tit_jet2_eta, tit_jet3_eta;
				stringstream tit_jet1_phi, tit_jet2_phi, tit_jet3_phi;
				stringstream tit_jet1_dphimet, tit_jet2_dphimet, tit_jet3_dphimet;
				stringstream tit_jet1_dphimht, tit_jet2_dphimht, tit_jet3_dphimht;

				tit_ht << title_suffix.str() << ";HT;Events;"; 
				tit_mht << title_suffix.str() << ";MHT;Events;"; 
				tit_mhtphi << title_suffix.str() << ";Corrected MHT-#Phi;Events;"; 
				tit_met << title_suffix.str() << ";MET;Events;"; 
				tit_metphi << title_suffix.str() << ";Corrected MET-#Phi;Events;"; 
				tit_uncorrmetphi << title_suffix.str() << ";Uncorrected MET-#Phi;Events;"; 
				tit_njetpt30 << title_suffix.str() << ";Jets [Pt>30, |#eta |< 5.0];Events;"; 
				tit_njetpt50 << title_suffix.str() << ";Jets [Pt>50, |#eta |< 2.5];Events;"; 
				tit_nvtx << title_suffix.str() << ";Primary Vertices;Events;"; 

				tit_jet1_pt << title_suffix.str() << ";Pt of Lead Jet^{Pt>30, |#eta |<5.0};Events;"; 
				tit_jet2_pt << title_suffix.str() << ";Pt of 2nd Lead Jet^{Pt>30, |#eta |<5.0};Events;"; 
				tit_jet3_pt << title_suffix.str() << ";Pt of 3rd Lead Jet^{Pt>30, |#eta |<5.0};Events;"; 

				tit_jet1_eta << title_suffix.str() << ";#eta Lead Jet^{Pt >30, |#eta |<5.0};Events;"; 
				tit_jet2_eta << title_suffix.str() << ";#eta 2nd Lead Jet^{Pt >30, |#eta |<5.0};Events;"; 
				tit_jet3_eta << title_suffix.str() << ";#eta 3rd Lead Jet^{Pt >30, |#eta |<5.0};Events;"; 

				tit_jet1_phi << title_suffix.str() << ";#Phi of Lead Jet^{Pt>30, |#eta |<5.0};Events;"; 
				tit_jet2_phi << title_suffix.str() << ";#Phi of 2nd Lead Jet^{Pt>30, |#eta |<5.0};Events;"; 
				tit_jet3_phi << title_suffix.str() << ";#Phi of 3rd Lead Jet^{Pt>30, |#eta |<5.0};Events;"; 

				tit_jet1_dphimet << title_suffix.str() << ";#delta#Phi (#slash{E}_{T}, Lead Jet^{Pt>30, |#eta |<5.0});Events;"; 
				tit_jet2_dphimet << title_suffix.str() << ";#delta#Phi (#slash{E}_{T}, 2nd Lead Jet^{Pt>30, |#eta |<5.0});Events;"; 
				tit_jet3_dphimet << title_suffix.str() << ";#delta#Phi (#slash{E}_{T}, 3rd Lead Jet^{Pt>30, |#eta |<5.0});Events;"; 

				tit_jet1_dphimht << title_suffix.str() << ";#delta#Phi (MHT, Lead Jet^{Pt>30, |#eta |<5.0});Events;"; 
				tit_jet2_dphimht << title_suffix.str() << ";#delta#Phi (MHT, 2nd Lead Jet^{Pt>30, |#eta |<5.0});Events;"; 
				tit_jet3_dphimht << title_suffix.str() << ";#delta#Phi (MHT, 3rd Lead Jet^{Pt>30, |#eta |<5.0});Events;"; 


				Hist_t hist;
				hist.h1_nVtx 	= new TH1D("nVtx", tit_nvtx.str().c_str(), 100, 0, 100); hist.h1_nVtx->Sumw2();
				hist.h1_ht 	= new TH1D("ht", tit_ht.str().c_str(), 160, 0, 8000); hist.h1_ht->Sumw2();
				hist.h1_met 	= new TH1D("met", tit_met.str().c_str(), 150, 0, 1500); hist.h1_met->Sumw2();
				hist.h1_metphi = new TH1D("metphi", tit_metphi.str().c_str(), 140, -3.5, 3.5); hist.h1_metphi->Sumw2();
				hist.h1_mht 	= new TH1D("mht", tit_mht.str().c_str(), 150, 0, 1500); hist.h1_mht->Sumw2();
				hist.h1_mhtphi = new TH1D("mhtphi", tit_mhtphi.str().c_str(), 140, -3.5, 3.5); hist.h1_mhtphi->Sumw2();
				hist.h1_uncorrmetphi 	= new TH1D("uncorrmetphi", tit_uncorrmetphi.str().c_str(), 140, -3.5, 3.5); hist.h1_uncorrmetphi->Sumw2();
				hist.h1_nJets_pt30eta50 = new TH1D("nJets_pt30eta50", tit_njetpt30.str().c_str(), 20, 0, 20); hist.h1_nJets_pt30eta50->Sumw2();
				hist.h1_nJets_pt50eta25 = new TH1D("nJets_pt50eta25", tit_njetpt50.str().c_str(), 20, 0, 20); hist.h1_nJets_pt50eta25->Sumw2();

				hist.h1_jet1_pt = new TH1D("jet1_pt", tit_jet1_pt.str().c_str(), 100, 0, 1000); hist.h1_jet1_pt->Sumw2();
				hist.h1_jet2_pt = new TH1D("jet2_pt", tit_jet2_pt.str().c_str(), 100, 0, 1000); hist.h1_jet2_pt->Sumw2();
				hist.h1_jet3_pt = new TH1D("jet3_pt", tit_jet3_pt.str().c_str(), 100, 0, 1000); hist.h1_jet3_pt->Sumw2();

				hist.h1_jet1_eta = new TH1D("jet1_eta", tit_jet1_eta.str().c_str(), 24, -6, 6); hist.h1_jet1_eta->Sumw2();
				hist.h1_jet2_eta = new TH1D("jet2_eta", tit_jet2_eta.str().c_str(), 24, -6, 6); hist.h1_jet2_eta->Sumw2();
				hist.h1_jet3_eta = new TH1D("jet3_eta", tit_jet3_eta.str().c_str(), 24, -6, 6); hist.h1_jet3_eta->Sumw2();

				hist.h1_jet1_phi = new TH1D("jet1_phi", tit_jet1_phi.str().c_str(), 140, -3.5, 3.5); hist.h1_jet1_phi->Sumw2();
				hist.h1_jet2_phi = new TH1D("jet2_phi", tit_jet2_phi.str().c_str(), 140, -3.5, 3.5); hist.h1_jet2_phi->Sumw2();
				hist.h1_jet3_phi = new TH1D("jet3_phi", tit_jet3_phi.str().c_str(), 140, -3.5, 3.5); hist.h1_jet3_phi->Sumw2();

				hist.h1_jet1_dphimet = new TH1D("jet1_dphimet", tit_jet1_dphimet.str().c_str(), 70, 0, 3.5); hist.h1_jet1_dphimet->Sumw2();
				hist.h1_jet2_dphimet = new TH1D("jet2_dphimet", tit_jet2_dphimet.str().c_str(), 70, 0, 3.5); hist.h1_jet2_dphimet->Sumw2();
				hist.h1_jet3_dphimet = new TH1D("jet3_dphimet", tit_jet3_dphimet.str().c_str(), 70, 0, 3.5); hist.h1_jet3_dphimet->Sumw2();

				hist.h1_jet1_dphimht = new TH1D("jet1_dphimht", tit_jet1_dphimht.str().c_str(), 70, 0, 3.5); hist.h1_jet1_dphimht->Sumw2();
				hist.h1_jet2_dphimht = new TH1D("jet2_dphimht", tit_jet2_dphimht.str().c_str(), 70, 0, 3.5); hist.h1_jet2_dphimht->Sumw2();
				hist.h1_jet3_dphimht = new TH1D("jet3_dphimht", tit_jet3_dphimht.str().c_str(), 70, 0, 3.5); hist.h1_jet3_dphimht->Sumw2();

				vHtBinnedHist.push_back(hist);	

			} /*end ht bin loop*/
			vJetBinnedHist.push_back(vHtBinnedHist);	
		} /*end jetbin loop*/
		vHist.push_back(vJetBinnedHist);	
	} /* end mt loop */
}

unsigned CSA14::GetVectorIndex(const vector< pair<unsigned, unsigned> >& binEdges, const unsigned& val)
{
	for (unsigned bin=0; bin < binEdges.size(); ++bin)
	{
		const double min = binEdges.at(bin).first;
		const double max = binEdges.at(bin).second;
		if (val>=min && val<=max) return bin;
	}
	return 99999;
}

