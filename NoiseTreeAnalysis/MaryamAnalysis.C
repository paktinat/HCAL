//////////////////////////////////////////////////////////
// Mon Dec 29 05:22:57 2014 by ROOT version 5.34/22
//////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
#include <memory>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <TH2.h>
//---------------------------------------------------------------------------
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TFileCollection.h"
#include "TChain.h"
//---------------------------------------------------------------------------
//#include "TreeReader_231228.h"
#include "treeBase.h"
#include "HCALmap.C"





using namespace std;
//---------------------------------------------------------------------------

double ieta2eta(int ieta) {
    static const double theHBHEEtaBounds[] = {0.000, 0.087, 0.087 * 2, 0.087 * 3, 0.087 * 4,
        0.087 * 5, 0.087 * 6, 0.087 * 7, 0.087 * 8, 0.087 * 9,
        0.087 * 10, 0.087 * 11, 0.087 * 12, 0.087 * 13, 0.087 * 14,
        0.087 * 15, 0.087 * 16, 0.087 * 17, 0.087 * 18, 0.087 * 19,
        1.74, 1.83, 1.93, 2.043, 2.172,
        2.332, 2.5, 2.65, 2.868, 3.000};
    double eta = 0.5 * (theHBHEEtaBounds[abs(ieta)] + theHBHEEtaBounds[abs(ieta) + 1]);
    if (ieta < 0) eta *= -1;
    return eta;
}

double iphi2phi(int iphi) {
    return 0.0875 * (iphi - 0.5);
}

static const double theHFEtaBounds[] = {2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839,
    4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

double UpperLimit(double energy) {

    double TS4TS5UpperThreshold[5] = {70, 90, 100, 400, 4000};
    double TS4TS5UpperCut[5] = {1, 0.8, 0.75, 0.72, 0.72};

    Int_t UpperBin = -1;

    for (int j = 0; j < 4; j++) {
        if (TS4TS5UpperThreshold[j] < energy)
            UpperBin = j;
    }

    if (UpperBin == -1)
        return 10000.0;

    else
        return (energy - TS4TS5UpperThreshold[UpperBin]) / (TS4TS5UpperThreshold[UpperBin + 1] - TS4TS5UpperThreshold[UpperBin]) * (TS4TS5UpperCut[UpperBin + 1] - TS4TS5UpperCut[UpperBin]) + TS4TS5UpperCut[UpperBin];
}

double LowerLimitRecHit(double energy) {
    //R45 per RecHit RecoLocalCalo/HcalRecProducers/python/HcalHitReconstructor_hbhe_cfi.py:
    double TS4TS5LowerThreshold[7] = {100, 120, 160, 200, 300, 500, 4000};
    double TS4TS5LowerCut[7] = {-1, -0.7, -0.5, -0.4, -0.3, 0.1, 0.1};

    Int_t LowerBin = -1;

    for (int j = 0; j < 7; j++) {
        if (TS4TS5LowerThreshold[j] < energy)
            LowerBin = j;
    }

    if (LowerBin == -1)
        return -10000.0;

    else{
	double outvalue = (energy - TS4TS5LowerThreshold[LowerBin]) / (TS4TS5LowerThreshold[LowerBin + 1] - TS4TS5LowerThreshold[LowerBin]) * (TS4TS5LowerCut[LowerBin + 1] - TS4TS5LowerCut[LowerBin]) + TS4TS5LowerCut[LowerBin];
	cout<<" => the output value: "<<outvalue<<endl;
        return (energy - TS4TS5LowerThreshold[LowerBin]) / (TS4TS5LowerThreshold[LowerBin + 1] - TS4TS5LowerThreshold[LowerBin]) * (TS4TS5LowerCut[LowerBin + 1] - TS4TS5LowerCut[LowerBin]) + TS4TS5LowerCut[LowerBin];
	}
}

double LowerLimitRBX(double energy) {
    //R45 per RBX RecoMET/METProducers/python/hcalnoiseinfoproducer_cfi.py:
    double TS4TS5LowerThreshold[8] = {100, 120, 150, 200, 300, 400, 500, 4000};
    double TS4TS5LowerCut[8] = {-1, -0.7, -0.4, -0.2, -0.08, 0, 0.1, 0.1};

    Int_t LowerBin = -1;

    for (int j = 0; j < 7; j++) {
        if (TS4TS5LowerThreshold[j] < energy)
            LowerBin = j;
    }

    if (LowerBin == -1)
        return -10000.0;

    else
        return (energy - TS4TS5LowerThreshold[LowerBin]) / (TS4TS5LowerThreshold[LowerBin + 1] - TS4TS5LowerThreshold[LowerBin]) * (TS4TS5LowerCut[LowerBin + 1] - TS4TS5LowerCut[LowerBin]) + TS4TS5LowerCut[LowerBin];
}

void prn(int n, int N, int d, TString& cf, TTree* ch) {

    if (cf != (ch->GetCurrentFile()->GetName())) {
        cout << "\n >> Current File: " << (ch->GetCurrentFile()->GetName()) << endl;
        cf = (ch->GetCurrentFile()->GetName());
    }

    if (n % d == 0) {
        //        fprintf(stdout, "\rProcessed events: \033[1;36;40m%6d of %6d\033[0m", n, N);
        fprintf(stdout, "\rProcessed events: %6d of %6d", n, N);
        fflush(stdout);
    }
    if (n == (N - 1)) {

        fprintf(stdout, "\n");
        fflush(stdout);
    }
}

double DeltaPhi(double v1, double v2) {
    double diff = fabs(v2 - v1);
    double corr = 2 * acos(-1.) - diff;
    if (diff < acos(-1.)) {
        return diff;
    } else {

        return corr;
    }
}

double GetDeltaR(double eta1, double eta2, double phi1, double phi2) {

    return sqrt((eta1 - eta2)*(eta1 - eta2) + DeltaPhi(phi1, phi2) * DeltaPhi(phi1, phi2));
}

Int_t isRBXNoisy(double c4, double c5, TH2D* h) {
    Int_t isRBXNsy = 0;
    double P45 = (c4 + c5);
    if(P45 < 1) P45 = 1; // this is put because http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0022
    double R45 = (c4 - c5) / P45;
    h->Fill(P45, R45);
    if ((R45 < LowerLimitRBX(P45) || R45 > UpperLimit(P45))) isRBXNsy = 1;
	//if(isRBXNsy==1) cout<<"\n			R45: "<<R45<<endl;
	//if(isRBXNsy==1 && R45 < LowerLimitRBX(P45)) cout<<"			Lower Limit Fired ------- "<<LowerLimitRBX(P45)<<endl;
	//if(isRBXNsy==1 && R45 > UpperLimit(P45)) cout<<"			Upper Limit Fired ------- "<<UpperLimit(P45)<<endl;
    return isRBXNsy;
}

Int_t isRecHitNoisy(double c4, double c5, TH2D* h) {
    Int_t isRHNsy = 0;
    double P45 = (c4 + c5);
    double R45 = (c4 - c5) / P45;
    h->Fill(P45, R45);
    if ((R45 < LowerLimitRecHit(P45) || R45 > UpperLimit(P45))) isRHNsy = 1;
	if(isRHNsy==1) cout<<"\n		 ========= bad rechits ========="<<endl; 
	cout<<"\nR45: "<<R45<<", lower: "<<LowerLimitRecHit(P45)<<", upper: "<<UpperLimit(P45)<<endl;
	if(isRHNsy==1 && R45 < LowerLimitRecHit(P45)) cout<<"			Lower Limit Fired ------- "<<LowerLimitRecHit(P45)<<endl;
	if(isRHNsy==1 && R45 > UpperLimit(P45)) cout<<"			Upper Limit Fired ------- "<<UpperLimit(P45)<<endl;
    return isRHNsy;
}

Int_t RBX_X(Int_t id) {
    int x2x[] = {2, 1, 3, 0};
    return x2x[(int) ((id - fmod(id, 18)) / 18)];
}

Int_t RBX_Y(Int_t id) {
    return fmod(id, 18);
}

/*
 * 
 */


//To run "root MaryamAnalysis.C'("Test")'   "
//---------------------------------------------------------------------------

int MaryamAnalysis(TString filelist = "TEST", bool addPU=false) {

    TString sFOut = "res_";
    sFOut.Append(filelist);

    if(addPU) 
      sFOut.Append("_ManPU");
    
    sFOut.Append(".root");
    TFile fout(sFOut, "RECREATE");

    TH2D* h_RBX_R45_P45 = new TH2D("RBX_R45_P45", "RBX_R45_P45", 100, 0, 1000, 30, -1.5, 1.5);
    TH2D* h_RecHit_R45_P45 = new TH2D("RecHit_R45_P45", "RecHit_R45_P45", 100, 0, 1000, 30, -1.5, 1.5);

    TH2D* h_RBXNoisy_NoisyNmFr_NoisyEnFr = new TH2D("RBXNoisy_NoisyNmFr_NoisyEnFr", "RBXNoisy_NoisyNmFr_NoisyEnFr", 21, 0, 1.05, 21, 0, 1.05);
    TH2D* h_RBXHealthy_NoisyNmFr_NoisyEnFr = new TH2D("RBXHealthy_NoisyNmFr_NoisyEnFr", "RBXHealthy_NoisyNmFr_NoisyEnFr", 21, 0, 1.05, 21, 0, 1.05);

    TH2D* h_RBXNoisy_RBXPhi_RBXEta = new TH2D("RBXNoisy_RBXPhi_RBXEta", "RBXNoisy_RBXPhi_RBXEta", 4, 0, 4, 18, 0, 18);
    TH2D* h_RBXHealthy_RBXPhi_RBXEta = new TH2D("RBXHealthy_RBXPhi_RBXEta", "RBXHealthy_RBXPhi_RBXEta", 4, 0, 4, 18, 0, 18);

    TH1D* h_RBXNoisy_En = new TH1D("RBXNoisy_En", "RBXNoisy_En", 100,0,1000);

    TH1D* h_RBXHealthy_En = new TH1D("RBXHealthy_En", "RBXHealthy_En", 100,0,1000);

    TH1D* h_TEST_TS = new TH1D("TEST_TS", "TEST_TS", 10, 0, 10);

    TH1D * h_PS_TS[10];

    for (int j = 0; j < 10; j++) {
        TString s_ps = "PS_TS_";
        stringstream ss;
        ss << j;
        s_ps.Append(ss.str());
        h_PS_TS[j] = new TH1D(s_ps, s_ps, 10, 0, 10);
    }

/*
    TH1F* h_RecHitEnergy = new TH1F("RecHitEnergy", "RecHitEnergy", 50,0,50);
    TH1F* h_RBXEnergy = new TH1F("RBXEnergy", "RBXEnergy", 50,0,50);

    TH2I* h_RBX_X_Y = new TH2I("RBX_X_Y", "RBX_X_Y", 4, 0, 4, 17, 0, 17);
    h_RBX_X_Y->GetXaxis()->SetBinLabel(1,"HE-");
    h_RBX_X_Y->GetXaxis()->SetBinLabel(2,"HB-");
    h_RBX_X_Y->GetXaxis()->SetBinLabel(3,"HB+");
    h_RBX_X_Y->GetXaxis()->SetBinLabel(4,"HE+");

    TH1D* h_Diff_RBXEn_RHSumEn = new TH1D("Diff_RBXEn_RHSumEn", "Diff_RBXEn_RHSumEn", 40, -10, 10);
    h_Diff_RBXEn_RHSumEn->GetXaxis()->SetTitle("#varepsilon_{RBX}-#varepsilon_{#Sigma RH}");
    h_Diff_RBXEn_RHSumEn->GetYaxis()->SetTitle("Energy Difference per RBX per Event");
    h_Diff_RBXEn_RHSumEn->GetXaxis()->SetTitleSize(0.07);
    h_Diff_RBXEn_RHSumEn->GetXaxis()->SetTitleOffset(0.6);
    h_Diff_RBXEn_RHSumEn->GetYaxis()->SetTitleSize(0.05);
    h_Diff_RBXEn_RHSumEn->GetYaxis()->SetTitleOffset(1);
    TProfile* h_EnDiff_Vs_RBXIndex = new TProfile("EnDiff_Vs_RBXIndex", "EnDiff_Vs_RBXIndex", 72, 0, 72);
    h_EnDiff_Vs_RBXIndex->GetXaxis()->SetTitle("RBX Index");
    h_EnDiff_Vs_RBXIndex->GetYaxis()->SetTitle("|#varepsilon_{RBX}-#varepsilon_{#Sigma RH}|");
    h_EnDiff_Vs_RBXIndex->GetXaxis()->SetTitleSize(0.07);
    h_EnDiff_Vs_RBXIndex->GetXaxis()->SetTitleOffset(0.6);
    h_EnDiff_Vs_RBXIndex->GetYaxis()->SetTitleSize(0.07);
    h_EnDiff_Vs_RBXIndex->GetYaxis()->SetTitleOffset(0.6);

    TH1I* h_counting = new TH1I("EnDiff", "#varepsilon_{RBX}-#varepsilon_{#Sigma RH}", 2,0,2);
    h_counting->GetXaxis()->SetBinLabel(1,"|#varepsilon_{RBX}-#varepsilon_{#Sigma RH}| <= 0.001");
    h_counting->GetXaxis()->SetBinLabel(2,"|#varepsilon_{RBX}-#varepsilon_{#Sigma RH}| > 0.001");
    h_counting->GetXaxis()->SetLabelSize(0.11);
*/
    TH2I* h_RBXR45Noise_vs_HasBadRBXR45 = new TH2I("h_RBXR45Noise_vs_HasBadRBXR45", "h_RBXR45Noise_vs_HasBadRBXR45", 2, 0, 2, 2, 0, 2);
    TH2I* h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose = new TH2I("h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose", "h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose", 2, 0, 2, 2, 0, 2);
    TH2I* h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight = new TH2I("h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight", "h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight", 2, 0, 2, 2, 0, 2);

    TH1D* MET_RBXR45Noise = new TH1D("MET_RBXR45Noise","MET_RBXR45Noise", 100, 0, 500);
    TH1D* MET_healthyEvents = new TH1D("MET_healthyEvents","MET_healthyEvents", 100, 0, 500);

    double recHitEnThreshold[5] = {1.5,2.5,5.,7.5,10.};	
    int arraysize = sizeof(recHitEnThreshold)/sizeof(recHitEnThreshold[0]);	

    TH1D* ntupleMET = new TH1D("ntupleMET","ntupleMET", 100, 0, 500);
    TH1D* hMETloose[5];
    TH1D* hMETtight[5];

    TString metThresholds[5] = {"1p5","2p5","5","7p5","10"};

    for (int j = 0; j < 5; j++) {
        TString histonamel = "MET_loose";
        TString histonamet = "MET_tight";
        histonamel.Append(metThresholds[j]);
        histonamet.Append(metThresholds[j]);
        hMETloose[j] = new TH1D(histonamel,histonamel, 100, 0, 500);
        hMETtight[j] = new TH1D(histonamet,histonamet, 100, 0, 500);
    }

    TString histos[6] = {"0","1","2","3","4","5"};

    TH1D* hpdHits[6];
    TH1D* hpdNoHits[6];

    for (int j = 0; j < 6; j++) {
        TString histoname = "HPDHits_";
        TString histoname1 = "HPDNoOtherHits_";
        histoname.Append(histos[j]);
        histoname1.Append(histos[j]);
        hpdHits[j] = new TH1D(histoname,histoname, 20, 0, 20);
        hpdNoHits[j] = new TH1D(histoname1,histoname1, 20, 0, 20);
    }

	vector<std::string> inputfile;

/*
	// Commissioning2014
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_2402393C_D869_E411_89A7_02163E010EC9.root");
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_70466B54_DA69_E411_91CF_02163E010DD6.root");
*/
	// MET 
        //inputfile.push_back("/home/HCAL/data-Run2015A-MET-RAW-v1-248-038-00000-DA2F1ED8-0513-E511-BBE5-02163E01451E.root");
        //inputfile.push_back("/home/HCAL/data-Run2015B-MET-RAW-v1-251-252-00000-6E03DF13-E925-E511-8DD5-02163E0134D6.root");
        //inputfile.push_back("/home/HCAL/data-Run2015B-MET-RAW-v1-251-562-00000-C068C85C-6B28-E511-A828-02163E012787.root");
        inputfile.push_back("/home/HCAL/data-Run2015B-MET-RAW-v1-251-168-00000-12E48CA8-2C25-E511-9631-02163E0136B4.root");
/*
	// Jet
        inputfile.push_back("/home/HCAL/data-Run2015A-Jet-RAW-v1-246-963-00000-E24810A5-5B0A-E511-9906-02163E014646.root");
        inputfile.push_back("/home/HCAL/data-Run2015B-Jet-RAW-v1-251-160-00000-460A8FEE-2525-E511-8453-02163E01410A.root");
        inputfile.push_back("/home/HCAL/data-Run2015B-Jet-RAW-v1-251-028-00000-027F6477-7924-E511-8B35-02163E012543.root");
	// JetHT
        inputfile.push_back("/home/HCAL/data-Run2015A-JetHT-RAW-v1-248-038-00000-7081F977-C612-E511-80E3-02163E014614.root");
        inputfile.push_back("/home/HCAL/data-Run2015B-JetHT-RAW-v1-1-251-168-00000-50A19724-1E25-E511-9552-02163E0135D7.root");
        inputfile.push_back("/home/HCAL/data-Run2015C-JetHT-RAW-v1-1-254-790-00000-FE3D684C-3F48-E511-BFCC-02163E0139AD.root");
        inputfile.push_back("/home/HCAL/data-Run2015C-JetHT-RAW-v1-1-254-833-00000-080BD87D-1649-E511-BBD7-02163E01354F.root");
        inputfile.push_back("/home/HCAL/data-Run2015C-JetHT-RAW-v1-1-254-833-00000-A676631F-1749-E511-B539-02163E014515.root");
*/

	int badEvt_tightTag = 0;
	int goodEvt_notTag = 0;
	for (unsigned int file=0;file<inputfile.size();file++) {

	TFile *f = TFile::Open(inputfile[file].c_str()); 
	TTree *t; f->GetObject("hcalTupleTree/tree",t);


    // ========== B E G I N ==========
    Long64_t nentries = t->GetEntries();
    //Long64_t nentries = 100;

    unsigned int event = 0;
    unsigned int run = 0;
    std::vector<float> * HBHERecHitEnergy = 0;
    std::vector<float> * HBHERecHitEnergyMethod0 = 0;
    std::vector<int> * HBHERecHitRBXid = 0;
    std::vector<int> * HBHERecHitIEta = 0;
    std::vector<int> * HBHERecHitIPhi = 0;
    std::vector<float> * HBHERecHitEta = 0;
    std::vector<float> * HBHERecHitPhi = 0;
    std::vector<float> * HBHERecHitDepth = 0;
    std::vector< std::vector<float> > *HBHEDigiFC = 0;
    std::vector< std::vector<float> > *HBHEDigiAllFC = 0;
    std::vector< std::vector<float> > *HBHEDigiEnergy = 0;
    std::vector<float> *HBHEDigiRecEnergy = 0;
    std::vector<float> *HBHEDigiEta = 0;
    std::vector<double> *RBXEnergy = 0;
    std::vector<double> *RBXEnergy15Method0 = 0;
    std::vector< std::vector<double> > *RBXCharge = 0;
    std::vector<int>     *HPDHits = 0;
    std::vector<int>     *HPDNoOtherHits = 0;
    std::vector<int>     *HasBadRBXR45Method0 = 0;
    std::vector<int>     *HasBadRBXRechitR45LooseMethod0 = 0;
    std::vector<int>     *HasBadRBXRechitR45TightMethod0 = 0;
    std::vector<double> *HBET = 0;
    std::vector<double> *HEET = 0;

    t->SetBranchAddress("event",&event);
    t->SetBranchAddress("run",&run);
    t->SetBranchAddress("HBHERecHitEnergy",&HBHERecHitEnergy);
    t->SetBranchAddress("HBHERecHitEnergyMethod0",&HBHERecHitEnergyMethod0);
    t->SetBranchAddress("HBHERecHitRBXid",&HBHERecHitRBXid);
    t->SetBranchAddress("HBHERecHitIEta",&HBHERecHitIEta);
    t->SetBranchAddress("HBHERecHitIPhi",&HBHERecHitIPhi);
    t->SetBranchAddress("HBHERecHitEta",&HBHERecHitEta);
    t->SetBranchAddress("HBHERecHitPhi",&HBHERecHitPhi);
    t->SetBranchAddress("HBHERecHitDepth",&HBHERecHitDepth);
    t->SetBranchAddress("HBHEDigiFC",&HBHEDigiFC);
    t->SetBranchAddress("HBHEDigiAllFC",&HBHEDigiAllFC);
    t->SetBranchAddress("HBHEDigiEnergy",&HBHEDigiEnergy);
    t->SetBranchAddress("HBHEDigiRecEnergy",&HBHEDigiRecEnergy);
    t->SetBranchAddress("HBHEDigiEta",&HBHEDigiEta);
    t->SetBranchAddress("RBXEnergy",&RBXEnergy);
    t->SetBranchAddress("RBXEnergy15Method0",&RBXEnergy15Method0);
    t->SetBranchAddress("RBXCharge",&RBXCharge);
    t->SetBranchAddress("HPDHits",&HPDHits);
    t->SetBranchAddress("HPDNoOtherHits",&HPDNoOtherHits);
    t->SetBranchAddress("HasBadRBXR45Method0",&HasBadRBXR45Method0);
    t->SetBranchAddress("HasBadRBXRechitR45LooseMethod0",&HasBadRBXRechitR45LooseMethod0);
    t->SetBranchAddress("HasBadRBXRechitR45TightMethod0",&HasBadRBXRechitR45TightMethod0);
    t->SetBranchAddress("HBET",&HBET);
    t->SetBranchAddress("HEET",&HEET);

    int i_ps(0);

    TString CurrentFile = "";
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        t->GetEntry(jentry);

        prn(jentry, nentries, 1, CurrentFile, t);

        //if (!(tr.OfficialDecision == 1)) continue;//1 no error, 0 error

        double NRH[72][5] = {0};
        double RecHitNoisy_NRH[72][5] = {0};
        double ERH[72][5] = {0};
        double RecHitNoisy_ERH[72][5] = {0};

        double TS4RBX[72] = {};
        double TS5RBX[72] = {};

	double TS4RBXall[72] = {};
        double TS5RBXall[72] = {};

	double MET[4] = {0};
	double METx[4] = {0};
	double METy[4] = {0};

        // ========== LOOP OVER Channels ==========

	if(jentry!=735 && jentry!=1141) continue;

	for (int k = 0; k < arraysize; k++) {
	if (k!=2) continue;

        for (unsigned int i = 0; i < HBHERecHitRBXid->size(); i++) {

	    if(!(HBHERecHitEnergyMethod0->at(i) > recHitEnThreshold[k])) continue;

	    int RBXIndex = HBHERecHitRBXid->at(i);
/*
	    METx[0] += sin(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));
	    METy[0] += cos(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));

	    if(HBHERecHitEnergy->at(i) > 1.5) {
		METx[1] += sin(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));	    
		METy[1] += cos(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));
            }
	    if(HBHERecHitEnergy->at(i) > 3) {
		METx[2] += sin(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));	    
		METy[2] += cos(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));
            }
	    if(HBHERecHitEnergy->at(i) > 5) {
		METx[3] += sin(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));	    
		METy[3] += cos(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));
            }

	    //int RBXIndex = -10;		
	    //RBXIndex = halilCode->EtaPhitoRBX(HBHERecHitIEta->at(i),HBHERecHitIPhi->at(i),HBHERecHitDepth->at(i));
*/
            double C4 = (*HBHEDigiFC)[i][4];
            double C5 = (*HBHEDigiFC)[i][5];

	    TS4RBX[RBXIndex] += C4;
	    TS5RBX[RBXIndex] += C5;
/*
            double C4all = (*HBHEDigiAllFC)[i][4];
            double C5all = (*HBHEDigiAllFC)[i][5];

            TS4RBXall[RBXIndex] += C4all;
            TS5RBXall[RBXIndex] += C5all;
*/
	    //if(!(HBHERecHitEnergyMethod0->at(i) > 1.5)) continue;
/*
	    int ii_ps(0);
            for (int ii = 0; ii < 10; ii++) {
                h_TEST_TS->Fill(ii, (*HBHEDigiFC)[i][ii]);
                if (i_ps < 10 && (C4 + C5) > 50) {
                    h_PS_TS[i_ps]->Fill(ii, (*HBHEDigiFC)[i][ii]);
                    ii_ps = 1;
                }
            }

            if (ii_ps == 1) i_ps++;
*/

            NRH[RBXIndex][k] += 1.0;
            ERH[RBXIndex][k] += HBHERecHitEnergyMethod0->at(i);

	    //cout<<"\npass threshold, RBX: "<<HBHERecHitRBXid->at(i)<<", HBHERecHitEnergyMethod0: "<<HBHERecHitEnergyMethod0->at(i)<<endl;
            if (isRecHitNoisy(C4, C5, h_RecHit_R45_P45) && (fabs(HBHERecHitIEta->at(i))!=28 && fabs(HBHERecHitIEta->at(i))!=29)) {
	        cout<<"in the first loop over rechits, !!! bad rechit fired !!!"<<", k = "<<k<<endl;
		//cout<<"RBXIndex: "<<RBXIndex<<endl;
                RecHitNoisy_NRH[RBXIndex][k] += 1.0;
                RecHitNoisy_ERH[RBXIndex][k] += HBHERecHitEnergyMethod0->at(i);
            }
	}
        }// i over channels

        if(jentry == 735) cout<<"\nNm: "<<NRH[68][2]<<", Noisy Nm: "<<RecHitNoisy_NRH[68][2]<<", En: "<<ERH[68][2]<<", Noisy En: "<<RecHitNoisy_ERH[68][2]<<endl;
        if(jentry == 1141) cout<<"\nNm: "<<NRH[48][2]<<", Noisy Nm: "<<RecHitNoisy_NRH[48][2]<<", En: "<<ERH[48][2]<<", Noisy En: "<<RecHitNoisy_ERH[48][2]<<endl;
/*
	// calculate MET per event
	for(unsigned int l = 0; l < 4; l++) {
	if(!(METx[l]!=0 || METy[l]!=0)) continue;
	MET[l] = sqrt(pow(METx[l],2)+pow(METy[l],2));
	//hMET[l]->Fill(MET[l]);
	}
*/
	// compare the calculated MET with the one saved into the ntuples
	ntupleMET->Fill(sqrt(pow(HBET->at(0)+HEET->at(0),2)+pow(HBET->at(1)+HEET->at(1),2)));

	bool RBXR45Noise_EnergyGt50 = false;
	bool RBXRechitR45Tight[5] = {false,false,false,false,false};
	bool RBXRechitR45Loose[5] = {false,false,false,false,false};
        for (int i = 0; i < 72; i++) {

            double EnFr[5] = {0.};
            double NmFr[5] = {0.};

	for (int k = 0; k < arraysize; k++) {

            EnFr[k] = RecHitNoisy_ERH[i][k] / ERH[i][k];
            NmFr[k] = RecHitNoisy_NRH[i][k] / NRH[i][k];

		if (EnFr[k] > 0.2 || NmFr[k] > 0.2)  
			if (ERH[i][k] != 0 && NRH[i][k] != 0) 
				RBXRechitR45Tight[k] = true; 

		if (EnFr[k] > 0.5 || NmFr[k] > 0.5)  
			if (ERH[i][k] != 0 && NRH[i][k] != 0)
				RBXRechitR45Loose[k] = true; 

	}
            if (isRBXNoisy((*RBXCharge)[i][4], (*RBXCharge)[i][5], h_RBX_R45_P45) && RBXEnergy15Method0->at(i) > 50) {
	        //cout<<"\nRBXEnergy15Method0->at("<<i<<"): "<<RBXEnergy15Method0->at(i)<<", output of ntuple?? "<<HasBadRBXR45Method0->at(0)<<endl;
	    //if (isRBXNoisy(TS4RBXall[i], TS5RBXall[i], h_RBX_R45_P45))  // According to http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0020, it seems that instead of HBHEDigiFC, HBHEDigiAllFC is used.

		RBXR45Noise_EnergyGt50 = true; // there is a cut on RBXEnergy according to http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0018

                h_RBXNoisy_NoisyNmFr_NoisyEnFr->Fill(EnFr[2], NmFr[2]);
                h_RBXNoisy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXNoisy_En->Fill(RBXEnergy15Method0->at(i));
	    
	    } else {

                h_RBXHealthy_NoisyNmFr_NoisyEnFr->Fill(EnFr[2], NmFr[2]);
                h_RBXHealthy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXHealthy_En->Fill(RBXEnergy15Method0->at(i));
	    
	    }
        }// i over RBX

	if (RBXR45Noise_EnergyGt50 && RBXRechitR45Tight[2]) badEvt_tightTag++;
	if (!RBXR45Noise_EnergyGt50 && !RBXRechitR45Tight[2]) goodEvt_notTag++;
 	 
            h_RBXR45Noise_vs_HasBadRBXR45->Fill(HasBadRBXR45Method0->at(0),RBXR45Noise_EnergyGt50);
            h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose->Fill(HasBadRBXRechitR45LooseMethod0->at(0),RBXRechitR45Loose[2]);
            h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight->Fill(HasBadRBXRechitR45TightMethod0->at(0),RBXRechitR45Tight[2]);

	for (int k = 0; k < arraysize; k++) {

		double met = sqrt(pow(HBET->at(0)+HEET->at(0),2)+pow(HBET->at(1)+HEET->at(1),2));
		if(RBXRechitR45Loose[k] == false) hMETloose[k]->Fill(met);
		if(RBXRechitR45Tight[k] == false) hMETtight[k]->Fill(met);

	}

	double met = sqrt(pow(HBET->at(0)+HEET->at(0),2)+pow(HBET->at(1)+HEET->at(1),2));
	if(RBXR45Noise_EnergyGt50) MET_RBXR45Noise->Fill(met);
	else MET_healthyEvents->Fill(met);


	if((HasBadRBXR45Method0->at(0) && !RBXR45Noise_EnergyGt50) || (!HasBadRBXR45Method0->at(0) && RBXR45Noise_EnergyGt50))
		for (int i = 0; i < 72; i++) {
		if(RBXEnergy15Method0->at(i) > 50) {
	        cout<<"\nRBXEnergy15Method0->at("<<i<<"): "<<RBXEnergy15Method0->at(i)<<endl;
                cout<<"(*RBXCharge)["<<i<<"][4]: "<<(*RBXCharge)[i][4]<<endl;
                cout<<"(*RBXCharge)["<<i<<"][5]: "<<(*RBXCharge)[i][5]<<endl;
		double P45 = (*RBXCharge)[i][4]+(*RBXCharge)[i][5];
		if(P45 < 1) P45 = 1; // this is put because http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0022
		double R45 = ((*RBXCharge)[i][4] - (*RBXCharge)[i][5]) / P45;
		cout<<"\n                       R45: "<<R45<<endl;
		cout<<"                     ^^^^^^^^^ Lower Limit: "<<LowerLimitRBX(P45)<<endl;
		cout<<"                     ^^^^^^^^^ Upper Limit: "<<UpperLimit(P45)<<endl;

		}
		}
	if(HasBadRBXR45Method0->at(0) && !RBXR45Noise_EnergyGt50) cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is NOT recognized as NOISE by our R45 definition !! *** \n\n"<<endl;
	if(!HasBadRBXR45Method0->at(0) && RBXR45Noise_EnergyGt50) cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is recognized as NOISE by our R45 definition !! *** \n\n"<<endl;


	//if((HasBadRBXRechitR45TightMethod0->at(0) && !RBXRechitR45Tight[2]) || (!HasBadRBXRechitR45TightMethod0->at(0) && RBXRechitR45Tight[2]) || (HasBadRBXRechitR45LooseMethod0->at(0) && !RBXRechitR45Loose[2]) || (!HasBadRBXRechitR45LooseMethod0->at(0) && RBXRechitR45Loose[2])) { 
        for (unsigned int i = 0; i < HBHERecHitRBXid->size(); i++) {

	//cout<<"*** conflict *** Now in the last loop over rechits..."<<endl;
            double C4 = (*HBHEDigiFC)[i][4];
            double C5 = (*HBHEDigiFC)[i][5];

	    if(!(HBHERecHitEnergyMethod0->at(i) > 5)) continue;
	    cout<<"\nRecHitIEta: "<<HBHERecHitIEta->at(i)<<endl;
	    cout<<"pass threshold, RBX: "<<HBHERecHitRBXid->at(i)<<", HBHERecHitEnergyMethod0: "<<HBHERecHitEnergyMethod0->at(i)<<endl;

            if (isRecHitNoisy(C4, C5, h_RecHit_R45_P45) && (fabs(HBHERecHitIEta->at(i))!=28 && fabs(HBHERecHitIEta->at(i))!=29)) {
	        //cout<<"*** conflict *** Now in the last loop over rechits, !!! bad rechit fired !!!"<<endl;
	        cout<<"\nRBXEnergy15Method0->at("<<HBHERecHitRBXid->at(i)<<"): "<<RBXEnergy15Method0->at(HBHERecHitRBXid->at(i))<<endl;
                cout<<"(*HBHEDigiFC)[i][4]: "<<C4<<endl;
                cout<<"(*HBHEDigiFC)[i][5]: "<<C5<<endl;
		double P45 = C4+C5;
		//if(P45 < 1) P45 = 1; // this is put because http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0022
		double R45 = (C4-C5) / P45;
		cout<<"\n                       R45: "<<R45<<endl;
		cout<<"                     ^^^^^^^^^ Lower Limit: "<<LowerLimitRecHit(P45)<<endl;
		cout<<"                     ^^^^^^^^^ Upper Limit: "<<UpperLimit(P45)<<endl;
            }
		if(HBHERecHitEnergyMethod0->at(i)>2000){
		double P45 = C4+C5;
		double R45 = (C4-C5) / P45;
		cout<<"\n                       R45: "<<R45<<endl;
		cout<<"                     ^^^^^^^^^ Lower Limit: "<<LowerLimitRecHit(P45)<<endl;
		cout<<"                     ^^^^^^^^^ Upper Limit: "<<UpperLimit(P45)<<endl;
}
        }
        for (int i = 0; i < 72; i++) {

            double EnFr = RecHitNoisy_ERH[i][2] / ERH[i][2];
            double NmFr = RecHitNoisy_NRH[i][2] / NRH[i][2];

        if((NmFr!=0 || EnFr!=0) && (ERH[i][2] != 0 && NRH[i][2] != 0)) cout<<"\nNmFr: "<<NmFr<<", EnFr: "<<EnFr<<endl;

		if(isRBXNoisy((*RBXCharge)[i][4], (*RBXCharge)[i][5], h_RBX_R45_P45) && RBXEnergy15Method0->at(i) > 50) {
		cout<<"\n --------- RBX NOISY ---------"<<endl;
	        cout<<"\nRBXEnergy15Method0->at("<<i<<"): "<<RBXEnergy15Method0->at(i)<<endl;
                cout<<"(*RBXCharge)["<<i<<"][4]: "<<(*RBXCharge)[i][4]<<endl;
                cout<<"(*RBXCharge)["<<i<<"][5]: "<<(*RBXCharge)[i][5]<<endl;
		double P45 = (*RBXCharge)[i][4]+(*RBXCharge)[i][5];
		if(P45 < 1) P45 = 1; // this is put because http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0022
		double R45 = ((*RBXCharge)[i][4] - (*RBXCharge)[i][5]) / P45;
		cout<<"\n                       R45: "<<R45<<endl;
		cout<<"                     ^^^^^^^^^ Lower Limit: "<<LowerLimitRBX(P45)<<endl;
		cout<<"                     ^^^^^^^^^ Upper Limit: "<<UpperLimit(P45)<<endl;

		}
	}

	//}
	if(HasBadRBXRechitR45TightMethod0->at(0) && RBXRechitR45Tight[2]) {
		hpdHits[0]->Fill(HPDHits->at(0)); hpdNoHits[0]->Fill(HPDNoOtherHits->at(0));
		//cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ***   "<<HPDHits->at(0)<<"\n\n"<<endl;
	}
	if(!HasBadRBXRechitR45TightMethod0->at(0) && !RBXRechitR45Tight[2]) {
		hpdHits[1]->Fill(HPDHits->at(0)); hpdNoHits[1]->Fill(HPDNoOtherHits->at(0));
		//cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ***   "<<HPDHits->at(0)<<"\n\n"<<endl;
	}
	if(HasBadRBXRechitR45TightMethod0->at(0) && !RBXRechitR45Tight[2]) {
		hpdHits[2]->Fill(HPDHits->at(0)); hpdNoHits[2]->Fill(HPDNoOtherHits->at(0));
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is NOT recognized as TIGHT NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	}
	if(HasBadRBXRechitR45LooseMethod0->at(0) && !RBXRechitR45Loose[2]) {
		hpdHits[3]->Fill(HPDHits->at(0)); hpdNoHits[3]->Fill(HPDNoOtherHits->at(0));
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is NOT recognized as LOOSE NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	}
	if(!HasBadRBXRechitR45TightMethod0->at(0) && RBXRechitR45Tight[2]) {
		hpdHits[4]->Fill(HPDHits->at(0)); hpdNoHits[4]->Fill(HPDNoOtherHits->at(0));
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is recognized as TIGHT NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	}
	if(!HasBadRBXRechitR45LooseMethod0->at(0) && RBXRechitR45Loose[2]) {
		hpdHits[5]->Fill(HPDHits->at(0)); hpdNoHits[5]->Fill(HPDNoOtherHits->at(0));
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is recognized as LOOSE NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	}

    }//jentry	
}
	cout<<"nr of badEvt_tightTag: "<<badEvt_tightTag<<", and eff: "<<(double)badEvt_tightTag/MET_RBXR45Noise->GetEntries()<<endl;
	cout<<"nr of goodEvt_notTag: "<<goodEvt_notTag<<", and eff: "<<(double)goodEvt_notTag/MET_healthyEvents->GetEntries()<<endl;


    // ======= E N D =======
    fout.Write();
    fout.Close();


    return 0;
}
