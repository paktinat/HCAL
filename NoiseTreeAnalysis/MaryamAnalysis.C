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

    else
        return (energy - TS4TS5LowerThreshold[LowerBin]) / (TS4TS5LowerThreshold[LowerBin + 1] - TS4TS5LowerThreshold[LowerBin]) * (TS4TS5LowerCut[LowerBin + 1] - TS4TS5LowerCut[LowerBin]) + TS4TS5LowerCut[LowerBin];
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

    TH1D* h_RBXNoisy_En = new TH1D("RBXNoisy_En", "RBXNoisy_En", 50,0,50);

    TH1D* h_RBXHealthy_En = new TH1D("RBXHealthy_En", "RBXHealthy_En", 50,0,50);

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

    TH1D* ntupleMET = new TH1D("ntupleMET","ntupleMET", 100, 0, 100);
    TH1D* hMET[4];

    TString metThresholds[4] = {"noCut","1p5","3","5"};

    for (int j = 0; j < 4; j++) {
        TString histoname = "MET_";
        histoname.Append(metThresholds[j]);
        hMET[j] = new TH1D(histoname,histoname, 100, 0, 100);
    }

	vector<std::string> inputfile;

	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015A_MET_RAW_v1_248_038_DA2F1ED8-0513-E511-BBE5-02163E01451E.root");
/*
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015B_MET_RAW_v1_251_168_12E48CA8-2C25-E511-9631-02163E0136B4.root");
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015B_MET_RAW_v1_251_244_FE05BEF9-B425-E511-BD59-02163E0121CC.root");
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015B_MET_RAW_v1_251_251_16CD9951-CA25-E511-AE31-02163E013823.root");
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015B_MET_RAW_v1_251_252_6E03DF13-E925-E511-8DD5-02163E0134D6.root");
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015B_MET_RAW_v1_251_562_C068C85C-6B28-E511-A828-02163E012787.root");
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015B_MET_RAW_v1_251_562_CA73FD38-D128-E511-9223-02163E0146A4.root");
	inputfile.push_back("/dataLOCAL/dataLOCAL/HCAL/RootFiles/Schenara750Pre5/data_Run2015B_MET_RAW_v1_251_562_F8665038-D128-E511-9101-02163E014275.root");
*/

	for (unsigned int file=0;file<inputfile.size();file++) {

	TFile *f = TFile::Open(inputfile[file].c_str()); 
	TTree *t; f->GetObject("hcalTupleTree/tree",t);


    // ========== B E G I N ==========
    //Long64_t nentries = t->GetEntries();
    Long64_t nentries = 100;

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
    std::vector<int>     *HasBadRBXR45 = 0;
    std::vector<int>     *HasBadRBXRechitR45Loose = 0;
    std::vector<int>     *HasBadRBXRechitR45Tight = 0;
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
    t->SetBranchAddress("HasBadRBXR45",&HasBadRBXR45);
    t->SetBranchAddress("HasBadRBXRechitR45Loose",&HasBadRBXRechitR45Loose);
    t->SetBranchAddress("HasBadRBXRechitR45Tight",&HasBadRBXRechitR45Tight);
    t->SetBranchAddress("HBET",&HBET);
    t->SetBranchAddress("HEET",&HEET);

    int i_ps(0);

    TString CurrentFile = "";
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        t->GetEntry(jentry);

        prn(jentry, nentries, 1, CurrentFile, t);

        //if (!(tr.OfficialDecision == 1)) continue;//1 no error, 0 error

        double NRH[72] = {};
        double RecHitNoisy_NRH[72] = {};
        double ERH[72] = {};
        double RecHitNoisy_ERH[72] = {};

        double TS4RBX[72] = {};
        double TS5RBX[72] = {};

	double TS4RBXall[72] = {};
        double TS5RBXall[72] = {};

	double MET[4] = {0};
	double METx[4] = {0};
	double METy[4] = {0};

        // ========== LOOP OVER Channels ==========

        for (unsigned int i = 0; i < HBHERecHitRBXid->size(); i++) {

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

	    int RBXIndex = HBHERecHitRBXid->at(i);

            float C4 = (*HBHEDigiFC)[i][4];
            float C5 = (*HBHEDigiFC)[i][5];

	    TS4RBX[RBXIndex] += C4;
	    TS5RBX[RBXIndex] += C5;

            float C4all = (*HBHEDigiAllFC)[i][4];
            float C5all = (*HBHEDigiAllFC)[i][5];

            TS4RBXall[RBXIndex] += C4all;
            TS5RBXall[RBXIndex] += C5all;

	    //if(!(HBHERecHitEnergyMethod0->at(i) > 1.5)) continue;
	    if(!(HBHERecHitEnergyMethod0->at(i) > 5)) continue;

            NRH[RBXIndex] += 1.0;
            ERH[RBXIndex] += HBHERecHitEnergyMethod0->at(i);

	    int ii_ps(0);
            for (int ii = 0; ii < 10; ii++) {
                h_TEST_TS->Fill(ii, (*HBHEDigiFC)[i][ii]);
                if (i_ps < 10 && (C4 + C5) > 50) {
                    h_PS_TS[i_ps]->Fill(ii, (*HBHEDigiFC)[i][ii]);
                    ii_ps = 1;
                }
            }
            if (ii_ps == 1) i_ps++;

            if (isRecHitNoisy(C4, C5, h_RecHit_R45_P45)) {
                RecHitNoisy_NRH[RBXIndex] += 1.0;
                RecHitNoisy_ERH[RBXIndex] += HBHERecHitEnergyMethod0->at(i);
            }
        }// i over channels


	// calculate MET per event
	for(unsigned int l = 0; l < 4; l++) {

	if(!(METx[l]!=0 || METy[l]!=0)) continue;
	MET[l] = sqrt(pow(METx[l],2)+pow(METy[l],2));
	hMET[l]->Fill(MET[l]);

	}
	// compare the calculated MET with the one saved into the ntuples
	ntupleMET->Fill(sqrt(pow(HBET->at(0)+HEET->at(0),2)+pow(HBET->at(1)+HEET->at(1),2)));

	bool RBXR45Noise_EnergyGt50 = false;
	bool RBXRechitR45Tight = false;
	bool RBXRechitR45Loose = false;

        for (int i = 0; i < 72; i++) {

            double EnFr = RecHitNoisy_ERH[i] / ERH[i];
            double NmFr = RecHitNoisy_NRH[i] / NRH[i];

		if (EnFr > 0.2 || NmFr > 0.2)  
			if (ERH[i] != 0 && NRH[i] != 0) 
				RBXRechitR45Tight = true; 

		if (EnFr > 0.5 || NmFr > 0.5)  
			if (ERH[i] != 0 && NRH[i] != 0)
				RBXRechitR45Loose = true; 


            if (isRBXNoisy((*RBXCharge)[i][4], (*RBXCharge)[i][5], h_RBX_R45_P45)) {
	    //if (isRBXNoisy(TS4RBXall[i], TS5RBXall[i], h_RBX_R45_P45)) { // According to http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0020, it seems that instead of HBHEDigiFC, HBHEDigiAllFC is used.

		if(RBXEnergy15Method0->at(i) > 50) RBXR45Noise_EnergyGt50 = true; // there is a cut on RBXEnergy according to http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0018

                h_RBXNoisy_NoisyNmFr_NoisyEnFr->Fill(EnFr, NmFr);
                h_RBXNoisy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXNoisy_En->Fill(ERH[i]);
	    
	    } else {

                h_RBXHealthy_NoisyNmFr_NoisyEnFr->Fill(EnFr, NmFr);
                h_RBXHealthy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXHealthy_En->Fill(ERH[i]);
	    
	    }
        }// i over RBX

            h_RBXR45Noise_vs_HasBadRBXR45->Fill(HasBadRBXR45->at(0),RBXR45Noise_EnergyGt50);
            h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose->Fill(HasBadRBXRechitR45Loose->at(0),RBXRechitR45Loose);
            h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight->Fill(HasBadRBXRechitR45Tight->at(0),RBXRechitR45Tight);

	if((HasBadRBXR45->at(0) && !RBXR45Noise_EnergyGt50) || (!HasBadRBXR45->at(0) && RBXR45Noise_EnergyGt50))
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
	if(HasBadRBXR45->at(0) && !RBXR45Noise_EnergyGt50) cout<<"\n*** !! event "<<jentry<<" is NOT recognized as NOISE by our R45 definition !! *** \n\n"<<endl;
	if(!HasBadRBXR45->at(0) && RBXR45Noise_EnergyGt50) cout<<"\n*** !! event "<<jentry<<" is recognized as NOISE by our R45 definition !! *** \n\n"<<endl;

    }//jentry	
}


    // ======= E N D =======
    fout.Write();
    fout.Close();


    return 0;
}
