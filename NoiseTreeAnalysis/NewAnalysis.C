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

    for (int j = 0; j < 6; j++) {
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

Int_t isRecHitNoisy2(double c4, double c5, TH2D* h) {
    Int_t isRHNsy = 0;
    double P45 = (c4 + c5);
    double R45 = (c4 - c5) / P45;
    h->Fill(P45, R45);
    if ((R45 < LowerLimitRecHit(P45) || R45 > UpperLimit(P45))) isRHNsy = 1;
    return isRHNsy;
}

/*
Int_t isRecHitNoisy3(double c5, double c6, TH2D* h) {
    Int_t isRHNsy = 0;
    double P56 = (c5 + c6);
    double R56 = (c5 - c6) / P56;
    h->Fill(P56, R56);
    if ((R56 < LowerLimitRecHit(P56) || R56 > UpperLimit(P56))) isRHNsy = 1;
    return isRHNsy;
}
*/




Int_t RBX_X(Int_t id) {
    int x2x[] = {2, 1, 3, 0};
    return x2x[(int) ((id - fmod(id, 18)) / 18)];
}

Int_t RBX_Y(Int_t id) {
    return fmod(id, 18);
}

//To run "root NewAnalysis.C'("Test")'   "
//---------------------------------------------------------------------------

int NewAnalysis(TString filelist = "TEST", bool addPU = false) {

    TString sFOut = "res_";
    sFOut.Append(filelist);

    if(addPU) 
      sFOut.Append("_ManPU");
    
    sFOut.Append(".root");
    TFile fout(sFOut, "RECREATE");

    TH2D* h_RBX_R45_P45 = new TH2D("RBX_R45_P45", "RBX_R45_P45", 100, 0, 1000, 30, -1.5, 1.5);
   


    TH2D* h_RecHit_R45_P45 = new TH2D("RecHit_R45_P45", "RecHit_R45_P45", 100, 0, 1000, 30, -1.5, 1.5);
    TH2D* h_RecHit_R45_P45modified = new TH2D("RecHit_R45_P45modified", "RecHit_R45_P45modified", 100, 0, 1000, 30, -1.5, 1.5);
    TH2D* h_RecHit_R34_P34 = new TH2D("RecHit_R34_P34", "RecHit_R34_P34", 100, 0, 1000, 30, -1.5, 1.5);
    TH2D* h_RecHit_R56_P56 = new TH2D("RecHit_R56_P56", "RecHit_R56_P56", 100, 0, 1000, 30, -1.5, 1.5);


    TH2D* h_RBXNoisy_NoisyNmFr_NoisyEnFr = new TH2D("RBXNoisy_NoisyNmFr_NoisyEnFr", "RBXNoisy_NoisyNmFr_NoisyEnFr", 21, 0, 1.05, 21, 0, 1.05);
    TH2D* h_RBXHealthy_NoisyNmFr_NoisyEnFr = new TH2D("RBXHealthy_NoisyNmFr_NoisyEnFr", "RBXHealthy_NoisyNmFr_NoisyEnFr", 21, 0, 1.05, 21, 0, 1.05);

    TH2D* h_RBXNoisy_RBXPhi_RBXEta = new TH2D("RBXNoisy_RBXPhi_RBXEta", "RBXNoisy_RBXPhi_RBXEta", 4, 0, 4, 18, 0, 18);
    TH2D* h_RBXHealthy_RBXPhi_RBXEta = new TH2D("RBXHealthy_RBXPhi_RBXEta", "RBXHealthy_RBXPhi_RBXEta", 4, 0, 4, 18, 0, 18);

    TH1D* h_RBXNoisy_En = new TH1D("RBXNoisy_En", "RBXNoisy_En", 100,0,1000);

    TH1D* h_RBXHealthy_En = new TH1D("RBXHealthy_En", "RBXHealthy_En", 100,0,1000);

    TH1D* h_RechitEnergy_modifierdrechit = new TH1D("h_RechitEnergy_modifierdrechit","h_RechitEnergy_modifierdrechit", 100,0,1000);
    TH1D* h_RechitEnergyRAW_modifierdrechit = new TH1D("h_RechitEnergyRAW_modifierdrechit","h_RechitEnergyRAW_modifierdrechit", 100,0,1000);
   
 
//------------------------Shirin------------------------------------------
   //Rechit
    TH1D* h_NrechitRBXid = new TH1D("h_NrechitRBXid", "h_NrechitRBXid", 100,0,100);
    
    TH1D* h_TS3R45down = new TH1D("h_TS3R45down", "h_TS3R45down", 40,0,400);
    TH1D* h_TS4R45down= new TH1D("h_TS4R45down", "h_TS4R45down",  40,0,400);
    TH1D* h_TS5R45down = new TH1D("h_TS5R45down", "h_TS5R45down", 40,0,400);
    TH1D* h_TS6R45down = new TH1D("h_TS6R45down", "h_TS6R45down", 40,0,400);

    TH1D* h_TS3R45up = new TH1D("h_TS3R45up", "h_TS3R45up", 40,0,400);
    TH1D* h_TS4R45up= new TH1D("h_TS4R45up", "h_TS4R45up", 40,0,400);
    TH1D* h_TS5R45up = new TH1D("h_TS5R45up", "h_TS5R45up", 40,0,400);
    TH1D* h_TS6R45up = new TH1D("h_TS6R45up", "h_TS6R45up", 40,0,400);


    TH1D* h_TS3TS4up= new TH1D("h_TS3TS4up", "h_TS3TS4up", 80,0,400);
    TH1D* h_TS5TS6up= new TH1D("h_TS5TS6up", "h_TS5TS6up", 80,0,400);

    TH1D* h_TS3TS4down= new TH1D("h_TS3TS4down", "h_TS3TS4down", 80,0,400);
    TH1D* h_TS5TS6down= new TH1D("h_TS5TS6down", "h_TS5TS6down", 80,0,400);

   /* TH1D* h_R56new = new TH1D("h_R56new","h_R56new",100,0,2);
    TH1D* h_R45new = new TH1D("h_R45new","h_R45new",100,0,2);
    TH1D* h_R34new = new TH1D("h_R34new","h_R34new",100,0,2);*/
    TH2D* h_R45RBX_Rechit = new TH2D("h_R45RBX_Rechit", "h_R45RBX_Rechit", 100, 0, 1000, 30, -1.5, 1.5); 
    TH2D* h_R34RBX_Rechitup= new TH2D("h_R34RBX_Rechitup", "h_R34RBX_Rechitup", 100, 0, 1000, 30, -1.5, 1.5);
    TH2D* h_R56RBX_Rechitup = new TH2D("h_R56RBX_Rechitup", "h_R56RBX_Rechitup", 100, 0, 1000, 30, -1.5, 1.5);

    TH2D* h_R34RBX_Rechitdown = new TH2D("h_R34RBX_Rechitdown", "h_R34RBX_Rechitdown", 100, 0, 1000, 30, -1.5, 1.5);
    TH2D* h_R56RBX_Rechitdown = new TH2D("h_R56RBX_Rechitdown", "h_R56RBX_Rechitdown", 100, 0, 1000, 30, -1.5, 1.5);


    TH2D* h_RBX_NmFr_EnFr_Rechit = new TH2D("h_RBX_NmFr_EnFr_Rechit", "h_RBX_NmFr_EnFr_Rechit", 21, 0, 1.05, 21, 0, 1.05);
    TH2D* h_RBX_RBXPhi_RBXEta_RecHit = new TH2D("h_RBX_RBXPhi_RBXEta_RecHit", "h_RBX_RBXPhi_RBXEta_RecHit", 4, 0, 4, 18, 0, 18);
  
    TH1D* hAllRBXs_En = new TH1D("hAllRBXs_En","hAllRBXs_En", 100,0,1000);
    TH1D* hloose_En= new TH1D("hloose_En","hloose_En",100,0,1000);
    TH1D* htight_En= new TH1D("htight_En","htight_En",100,0,1000);
    TH1D* hmodified_tightEn = new TH1D("hmodified_tightEn","hmodified_tightEn", 100,0,1000);
    TH1D* hmodified_looseEn = new TH1D("hmodified_looseEn","hmodified_looseEn", 100,0,1000);
    
     TH1D* hrechit_r45RBX1 = new TH1D("hrechit_r45RBX1","hrechit_r45RBX1",50,0,50);
    
    
   
   
//---------------------------------------------------------------------------  
    TH2I* h_RBXR45Noise_vs_HasBadRBXR45 = new TH2I("h_RBXR45Noise_vs_HasBadRBXR45", "h_RBXR45Noise_vs_HasBadRBXR45", 2, 0, 2, 2, 0, 2);
    TH2I* h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose = new TH2I("h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose", "h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose", 2, 0, 2, 2, 0, 2);
    TH2I* h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight = new TH2I("h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight", "h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight", 2, 0, 2, 2, 0, 2);

TH2I* holdrechitR45loose_vs_modifiedrechitR45loose = new TH2I("holdrechitR45loose_vs_modifiedrechitR45loose","holdrechitR45loose_vs_modifiedrechitR45loose",
 2, 0, 2, 2, 0, 2);
TH2I* holdrechitR45tight_vs_modifiedrechitR45tight = new TH2I("holdrechitR45tight_vs_modifiedrechitR45tight","holdrechitR45tight_vs_modifiedrechitR45tight",
 2, 0, 2, 2, 0, 2);
TH2I* holdrechitR45_vs_modifiedrechitR45 = new TH2I("holdrechitR45_vs_modifiedrechitR45","holdrechitR45_vs_modifiedrechitR45",
 2, 0, 2, 2, 0, 2);


TH2D*  RechitEnergy_vs_RechitEnergyRAW = new TH2D("RechitEnergy_vs_RechitEnergyRAW","RechitEnergy_vs_RechitEnergyRAW",
 20, 0, 200, 20, 0, 200);

TH2I* RechitEnergy_vs_charge45 = new TH2I("RechitEnergy_vs_charge45","RechitEnergy_vs_charge45",
 40, 0, 400, 100, 0, 1000);

TH2I* RechitEnergy_vs_TS4 = new TH2I("RechitEnergy_vs_TS4","RechitEnergy_vs_TS4",
 40, 0, 400, 20, 0, 200);

TH2I* RechitEnergy_vs_TS5 = new TH2I("RechitEnergy_vs_TS5","RechitEnergy_vs_TS5",
 40, 0, 400, 100, 0,1000);
   TH2I* RechitEnergy_vs_TS6 = new TH2I("RechitEnergy_vs_TS6","RechitEnergy_vs_TS6",
 40, 0, 400, 100, 0,1000); 

    TH1D* MET_RBXR45Noise = new TH1D("MET_RBXR45Noise","MET_RBXR45Noise", 100, 0, 500);
    TH1D* MET_healthyEvents = new TH1D("MET_healthyEvents","MET_healthyEvents", 100, 0, 500);

    TH1D* ntupleMET = new TH1D("ntupleMET","ntupleMET", 100, 0, 500);
    TH1D* hMET = new TH1D("MET","MET", 100, 0, 500);
    TH1D* hMETloose = new TH1D("looseMET","looseMET", 100, 0, 500);
    TH1D* hMETtight = new TH1D("tightMET","tightMET", 100, 0, 500);

	vector<std::string> inputfile;


      
      /*  inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-MET-RECO-PromptReco-v3-256-587-AE7486CF-925D-E511-BF47-02163E011B48.root");
        inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-MET-RECO-PromptReco-v3-256-675-2099756C-145F-E511-923F-02163E0146EC.root");
        inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-MET-RECO-PromptReco-v3-256-672-B0B0BE29-075F-E511-B461-02163E0142DD.root");
        inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-MET-RECO-PromptReco-v4-258-159-08C67925-B36B-E511-AC3B-02163E0134A6.root");
        inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-MET-RECO-PromptReco-v4-258-159-043E1FA7-BD6B-E511-B7A5-02163E0141BE.root");
        inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-MET-RECO-PromptReco-v3-256-584-AAD4E42C-855D-E511-B4F5-02163E01206B.root");*/

      //SingleMu
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-705-DE427009-0871-E511-B7EA-02163E012288.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-656-EABEF9F6-4570-E511-A0BB-02163E011DD0.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-215-22F4433A-0D6D-E511-83C6-02163E0134C9.root");
      inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-257-613-EE1BE058-CE66-E511-9253-02163E01368B.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-046D8FD3-125F-E511-B25F-02163E014504.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-630-343A297D-2C5F-E511-81D3-02163E011931.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-926-8A703CBD-8661-E511-B641-02163E014318.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-257-599-0E89AB21-4966-E511-87F7-02163E0144AA.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-926-FC002608-8761-E511-880E-02163E0138A7.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-868-DEC34B16-7261-E511-AF16-02163E0146EB.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-403-FC1DFDD8-FB6D-E511-8A3B-02163E0139B5.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-403-4D696D1F4-F96D-E511-9DCC-02163E0134A6.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-159-FC0426EA-D76B-E511-90EB-02163E0141FB.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-434-EA03C1D8-5E6E-E511-A4E8-02163E011EB8.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-129-BC4871F2-B36A-E511-BC32-02163E01462C.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-136-C024C6AC-A66A-E511-BD94-02163E011FE7.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-158-DE13235B-796B-E511-B10C-02163E011B49.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-434-FA0B3E86-586E-E511-8BFA-02163E012078.root");
       inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-158-E2CBE545-836B-E511-9F15-02163E011B49.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-705-3A26A450-1B71-E511-98CE-02163E014417.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-655-F842BBDA-3670-E511-B046-02163E014451.root");
//inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-705-DE427009-0871-E511-B7EA-02163E012288.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-655-28DC9296-2570-E511-B89B-02163E012AC8.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-705-96B5F515-0571-E511-80A9-02163E0139C9.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-157-8A73D63D-D36A-E511-B652-02163E0143B2.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-714-ECF106F1-7B71-E511-BCEB-02163E0118DB.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-714-4A52288B-8E71-E511-80A8-02163E011C22.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-157-963310A5-DE6A-E511-A22A-02163E014310.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-714-2C2F88A2-7B71-E511-B802-02163E01435B.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-157-2293C54B-C46A-E511-B88A-02163E014615.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-157-4215570C-CB6A-E511-BDB6-02163E011F38.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-714-0470A923-8E71-E511-94BF-02163E01438E.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-157-1690EA75-CA6A-E511-B280-02163E01203C.root");

inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-129-7AF5170C-B66A-E511-A0D2-02163E0144B7.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-444-9EB46645-D86E-E511-A6C4-02163E0144B6.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-444-803901A1-C16E-E511-B298-02163E014395.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-444-1E702F56-C86E-E511-B096-02163E0145DF.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-129-0336E6B-A36A-E511-BD10-02163E011B09.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-129-74C7C30F-9F6A-E511-B7DD-02163E013862.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-129-10CF629D-A66A-E511-8945-02163E014261.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-258-129-0498E756-A16A-E511-B4D5-02163E0134AE.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-677-B4B2972C-2A5F-E511-9446-02163E0079EF.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-3241C3A7-355F-E511-841C-02163E0126BC.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v4-258-444-6898605C-DD6E-E511-8AB6-02163E0118F4.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-677-C8E4420E-255F-E511-A0D2-02163E0125AD.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-677-84E2938D-275F-E511-AB6F-02163E013559.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-677-72C18E67-175F-E511-B38B-02163E01389D.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-677-E28EC3C1-3B5F-E511-9C14-02163E014422.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-677-C0B77DD5-205F-E511-AB97-02163E01339A.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-677-76C45591-015F-E511-9417-02163E0124FF.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-CAF0919F-9A5F-E511-BA0D-02163E0145FF.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-F4F326C5-355F-E511-A884-02163E011D52.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-D4F29CB4-735F-E511-82DC-02163E011EA7.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-B2F75947-325F-E511-97BC-02163E01463F.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-0219CFC2-2A5F-E511-891E-02163E01366D.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-0EC2FEF2-385F-E511-B4C5-02163E011A22.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-675-E00C8A5A-385F-E511-B1E6-02163E0128B2.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015D-SingleMuon-RECO-PromptReco-v3-256-673-323E6084-145F-E511-96EF-02163E0146E4.root");





//SingleEle
/*
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-790-840A13CC-D449-E511-A9EB-02163E011856.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-790-C6B681F4-DB49-E511-9316-02163E0134A4.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-852-4C42DB8D-944B-E511-A24E-02163E011E5D.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-790-C83F3102-DA49-E511-9458-02163E0137B7.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-852-7470A18F-944B-E511-8440-02163E0125A4.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-231-0C33AAE5-1C46-E511-8D78-02163E011E9B.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-232-F48447DD-2046-E511-A0F1-02163E012AB5.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-790-284ECDD4-D949-E511-A252-02163E0146E8.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-852-4C42DB8D-944B-E511-A24E-02163E011E5D.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-879-E40A8E2D-A04B-E511-9E4F-02163E0145B5.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-906-00D92778-DB4B-E511-B00A-02163E011DE0.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-906-E25D0179-DB4B-E511-B593-02163E011D21.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-907-C2101CF4-E44B-E511-88E6-02163E01376B.root");
inputfile.push_back("/cmsdata2/HCAL/data-Run2015C-SingleElectron-RECO-PromptReco-v1-254-914-2E9EB21A-ED4B-E511-BFE3-02163E0142B5.root");
*/
//inputfile.push_back("/cmsdata2/HCAL/MET__Run2015D-PromptReco-v4__RECO_99_1_O7F.root");
//inputfile.push_back("/cmsdata2/HCAL/MET__Run2015D-PromptReco-v4__RECO_79_1_2Yw.root");
/*
inputfile.push_back("/cmsdata2/HCAL/NoBPTX__Run2015D-PromptReco-v4__RECO_33_1_jFd.root");
inputfile.push_back("/cmsdata2/HCAL/NoBPTX__Run2015D-PromptReco-v4__RECO_73_1_wAF.root");
inputfile.push_back("/cmsdata2/HCAL/NoBPTX__Run2015D-PromptReco-v4__RECO_78_1_52L.root");*/


        int badEvt_tightTag = 0;
	int goodEvt_notTag = 0;
	int nevts= 0;
        int nevts2= 0;
	for (unsigned int file=0;file<inputfile.size();file++) {

	TFile *f = TFile::Open(inputfile[file].c_str()); 
	TTree *t; f->GetObject("hcalTupleTree/tree",t);


    // ========== B E G I N ==========
    Long64_t nentries = t->GetEntries();
    //Long64_t nentries = 10;
	nevts+=nentries;
    unsigned int event = 0;
    unsigned int run = 0;
    unsigned int ls = 0;
    std::vector<float> * HBHERecHitEnergy = 0;
    std::vector<double> * HBHERecHitEnergyRaw = 0;
    std::vector<int> * HBHERecHitRBXid = 0;
    std::vector<int> * HBHERecHitIEta = 0;
    std::vector<int> * HBHERecHitIPhi = 0;
    std::vector<float> * HBHERecHitEta = 0;
    std::vector<float> * HBHERecHitPhi = 0;
    std::vector<float> * HBHERecHitDepth = 0;
    std::vector<int> * HBHERecHitFlags = 0;
    std::vector< std::vector<double> > *HBHERecHitAuxFC = 0;
    std::vector<double> *RBXEnergy = 0;
    std::vector<double> *RBXEnergy15 = 0;
    std::vector< std::vector<double> > *RBXCharge = 0;
    std::vector< std::vector<double> > *RBXCharge15 = 0;
    std::vector<int>     *HPDHits = 0;
    std::vector<int>     *HPDNoOtherHits = 0;
    std::vector<int>     *HasBadRBXR45 = 0;
    std::vector<int>     *HasBadRBXRechitR45Loose = 0;
    std::vector<int>     *HasBadRBXRechitR45Tight = 0;
    std::vector<double> *HBET = 0;
    std::vector<double> *HEET = 0;

    t->SetBranchAddress("event",&event);
    t->SetBranchAddress("run",&run);
    t->SetBranchAddress("ls",&ls);
    t->SetBranchAddress("HBHERecHitEnergy",&HBHERecHitEnergy);
    t->SetBranchAddress("HBHERecHitEnergyRaw",&HBHERecHitEnergyRaw);
    t->SetBranchAddress("HBHERecHitRBXid",&HBHERecHitRBXid);
    t->SetBranchAddress("HBHERecHitIEta",&HBHERecHitIEta);
    t->SetBranchAddress("HBHERecHitIPhi",&HBHERecHitIPhi);
    t->SetBranchAddress("HBHERecHitEta",&HBHERecHitEta);
    t->SetBranchAddress("HBHERecHitPhi",&HBHERecHitPhi);
    t->SetBranchAddress("HBHERecHitDepth",&HBHERecHitDepth);
    t->SetBranchAddress("HBHERecHitFlags",&HBHERecHitFlags);
    t->SetBranchAddress("HBHERecHitAuxFC",&HBHERecHitAuxFC);
    t->SetBranchAddress("RBXEnergy",&RBXEnergy);
    t->SetBranchAddress("RBXEnergy15",&RBXEnergy15);
    t->SetBranchAddress("RBXCharge",&RBXCharge);
    t->SetBranchAddress("RBXCharge15",&RBXCharge15);
    t->SetBranchAddress("HPDHits",&HPDHits);
    t->SetBranchAddress("HPDNoOtherHits",&HPDNoOtherHits);
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

        double NRH[72] = {0};
        double RecHitNoisy_NRH[72] = {0};
        double ERH[72] = {0};
        double RecHitNoisy_ERH[72] = {0};

        double NRHNew[72] = {0};
        double RecHitNoisy_NRHNew[72] = {0};
        double ERHNew[72] = {0};
        double RecHitNoisy_ERHNew[72] = {0};



        double TS4RBX[72] = {};
        double TS5RBX[72] = {};

	double MET = {0};
	double METx = {0};
	double METy = {0};
        double MaxRBXEnergy=0;
	

	// apply HPD filters
	if(HPDHits->at(0)>=17 || HPDNoOtherHits->at(0)>=10) continue;
        nevts2+=nentries;
        // apply NEF falg
        bool RecHitFlagsNoisy=false;
        for (unsigned int i = 0; i < HBHERecHitFlags->size(); i++) {
        if(fabs(HBHERecHitIEta->at(i))==28 || fabs(HBHERecHitIEta->at(i))==29) continue;
        if(!(HBHERecHitEnergyRaw->at(i) > 5.)) continue;
        int bitval = HBHERecHitFlags->at(i)>>27;
        int outval = bitval & 1;
        if (outval==1) RecHitFlagsNoisy=true;}



//-----------------------------------------------------------------------
            bool RBX_rechit[72]= {false};
            

            for (unsigned int i = 0; i < HBHERecHitRBXid->size(); i++) {

	    METx += sin(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));
	    METy += cos(HBHERecHitPhi->at(i))*HBHERecHitEnergy->at(i)/cosh(HBHERecHitEta->at(i));

	    if(!(HBHERecHitEnergyRaw->at(i) > 5.)) continue;
            
	    int RBXIndex = HBHERecHitRBXid->at(i);
            double C3 = (*HBHERecHitAuxFC)[i][3];
            double C4 = (*HBHERecHitAuxFC)[i][4];
            double C5 = (*HBHERecHitAuxFC)[i][5];
            double C6 = (*HBHERecHitAuxFC)[i][6];
          
   //      int bitval = HBHERecHitFlags->at(i)>>27;
     //    int outval = bitval & 1;
     //    if (outval==1) continue;
//
            NRH[RBXIndex] += 1.0;
            ERH[RBXIndex] += HBHERecHitEnergyRaw->at(i);
 
            NRHNew[RBXIndex] += 1.0;
            ERHNew[RBXIndex] += HBHERecHitEnergyRaw->at(i);
           
             double P45=(C4+C5);
             double P56=(C5+C6);
             double P34 = (C3 + C4);
             double R34=(C3-C4)/P34;
             double R45=(C4-C5)/P45;
             double R56=(C5-C6)/P56;
            
 
            /* double R45new=C4/P45;
             double R56new=C5/P56;
             double R34new=C4/P34; */
             if    (  -1.1< R45 && R45 <-0.9 && 100 < P45 && RBX_rechit[RBXIndex]==false ) {
             h_NrechitRBXid->Fill(HBHERecHitRBXid->at(i));
             RBX_rechit[RBXIndex] = true;}


              if   (-1.1< R45 && R45 <-0.9  && 100 <P45) {
              h_R34RBX_Rechitdown->Fill(P34, R34);
              h_R56RBX_Rechitdown->Fill(P56, R56);  }  

             if (0.9< R45 && R45 <1.1  && 70 <P45){
             h_R34RBX_Rechitup->Fill(P34, R34);
             h_R56RBX_Rechitup->Fill(P56, R56);
                                                       }

             if ((0.9< R45 && R45 <1.1  && 70 <P45)||(-1.1< R45 && R45 <-0.9  && 100 <P45))
             h_R45RBX_Rechit->Fill(P45, R45);  


               if   (-1.1< R45 && R45 <-0.9  && 100 <P45){
                h_TS3R45down->Fill(C3);
                h_TS4R45down->Fill(C4);
                h_TS5R45down->Fill(C5);
                h_TS6R45down->Fill(C6);
                h_TS3TS4down->Fill(P34);
                h_TS5TS6down->Fill(P56);
               
                                                        }

                if (0.9< R45 && R45 <1.1  && 70 <P45){
                h_TS3R45up->Fill(C3);
                h_TS4R45up->Fill(C4);
                h_TS5R45up->Fill(C5);
                h_TS6R45up->Fill(C6);
                h_TS3TS4up->Fill(P34);
                h_TS5TS6up->Fill(P56);


                                               }

              
              
             

            bool  R56isinsideenvelope=false;              
            bool  R34isinsideenvelope=false;  
           

              if (  -1.1< R45 && R45 <-0.9  && 100 <P45){
              if (!(isRecHitNoisy(C5, C6, h_RecHit_R56_P56))) {
              if (C3+C4 < 20 ){
              R56isinsideenvelope=true;
              }}}


             if (  0.9 < R45 && R45 <1.1  && 70 <P45 ){
              
             if (!(isRecHitNoisy(C3, C4, h_RecHit_R34_P34))) {
              if ( C5+C6 < 20 ) {
             R34isinsideenvelope=true;}}}


             bool   oldrechitisNoisy=false;
             bool modifiedrechitisNoisy=false;              
	     if (isRecHitNoisy(C4, C5, h_RecHit_R45_P45) && (fabs(HBHERecHitIEta->at(i))!=28 && fabs(HBHERecHitIEta->at(i))!=29)) {
              RecHitNoisy_NRH[RBXIndex] += 1.0;
              RecHitNoisy_ERH[RBXIndex] += HBHERecHitEnergyRaw->at(i);
              oldrechitisNoisy=true;


              
              if (R34isinsideenvelope==false && R56isinsideenvelope==false){
	

              RecHitNoisy_NRHNew[RBXIndex] += 1.0;
              RecHitNoisy_ERHNew[RBXIndex] += HBHERecHitEnergyRaw->at(i);
              modifiedrechitisNoisy=true;	

}	
}              
               if((oldrechitisNoisy == false || modifiedrechitisNoisy == true))
               isRecHitNoisy2(C4, C5, h_RecHit_R45_P45modified);
               else
               {
                   if (HBHERecHitEnergy->at(i)>50)
                {int isofilterflag = HBHERecHitFlags->at(i)>>11;
                int IsoFilterFlag = isofilterflag & 1;
                if(IsoFilterFlag) cout<<"this RH fired IsoFilterFlag"<<endl;
                int spike = HBHERecHitFlags->at(i)>>13;
                int Spike = spike & 1;
                if(Spike) cout<<"this RH is a spike !!"<<endl;
                cout<<"TS3"<<C3<<endl;
                cout<<"energy"    <<HBHERecHitEnergy->at(i)   <<endl;
                cout<<"energyRAW"    <<HBHERecHitEnergyRaw->at(i)   <<endl; 
                cout<<"TS4"<<C4<<endl;
                cout<<"TS5"<<C5<<endl;
                cout<<"TS6"<<C6<<endl;
                cout<<"runnumber"<<run<<endl;
                cout<<"event"<<event<<endl;
                cout<<"luminosity"<<ls<<endl;
    
}

          h_RechitEnergy_modifierdrechit->Fill(HBHERecHitEnergy->at(i));
          h_RechitEnergyRAW_modifierdrechit->Fill(HBHERecHitEnergyRaw->at(i));

          RechitEnergy_vs_RechitEnergyRAW->Fill(HBHERecHitEnergyRaw->at(i),HBHERecHitEnergy->at(i));
          RechitEnergy_vs_charge45->Fill(P45,HBHERecHitEnergy->at(i));
          RechitEnergy_vs_TS4->Fill(C4,HBHERecHitEnergy->at(i));
          RechitEnergy_vs_TS5->Fill(C5,HBHERecHitEnergy->at(i));
          RechitEnergy_vs_TS6->Fill(C6,HBHERecHitEnergy->at(i));
          



}

         
              holdrechitR45_vs_modifiedrechitR45->Fill(oldrechitisNoisy,modifiedrechitisNoisy);

           


//}
              




              }// i over channels

	
              // calculate MET per event
	      if(!(METx!=0 || METy!=0)) continue;
	      MET = sqrt(pow(METx,2)+pow(METy,2));
	      hMET->Fill(MET);
              // compare the calculated MET with the one saved into the ntuples
	      double met = sqrt(pow(HBET->at(0)+HEET->at(0),2)+pow(HBET->at(1)+HEET->at(1),2));
	      ntupleMET->Fill(met);

	      bool RBXR45Noise_EnergyGt50 = false;
	      bool RBXRechitR45Tight = false;
	      bool RBXRechitR45Loose = false;
              bool RBXRechitR45TightNew = false;
	      bool RBXRechitR45LooseNew = false;

              for (int i = 0; i < 72; i++) {
              double EnFr = 0.;
              double NmFr = 0.;
              double EnFrNew=0;
              double NmFrNew=0;

              EnFr = RecHitNoisy_ERH[i] / ERH[i];
              NmFr = RecHitNoisy_NRH[i] / NRH[i];

              EnFrNew = RecHitNoisy_ERHNew[i] / ERHNew[i];
              NmFrNew = RecHitNoisy_NRHNew[i] / NRHNew[i];
 

	      if (EnFr > 0.2 || NmFr > 0.2) { 
	      if (ERH[i] != 0 && NRH[i] != 0) 
	      RBXRechitR45Tight = true; }

              if  (EnFrNew > 0.2 || NmFrNew > 0.2) { 
	      if (ERHNew[i] != 0 && NRHNew[i] != 0) 
	      RBXRechitR45TightNew = true; }

	       

               if (EnFr > 0.5 || NmFr > 0.5)  {
	       if (ERH[i] != 0 && NRH[i] != 0)
	       RBXRechitR45Loose = true; }

               if (EnFrNew > 0.5 || NmFrNew > 0.5)  {
	       if (ERHNew[i] != 0 && NRHNew[i] != 0)
	       RBXRechitR45LooseNew = true; }

          
 



               //hAllRBXs_En-> Fill(RBXEnergy15->at(i));
               if (RBX_rechit[i] == true ){
               hrechit_r45RBX1->Fill(NRH[i]);
               h_RBX_RBXPhi_RBXEta_RecHit->Fill(RBX_X(i), RBX_Y(i));
              // hRechitswithR451_En->Fill(RBXEnergy->at(i));
               h_RBX_NmFr_EnFr_Rechit->Fill(EnFr, NmFr);}

               if (isRBXNoisy((*RBXCharge)[i][4], (*RBXCharge)[i][5], h_RBX_R45_P45) && RBXEnergy15->at(i) > 50) {

		RBXR45Noise_EnergyGt50 = true; // there is a cut on RBXEnergy according to http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/HcalNoiseAlgo.cc?v=CMSSW_7_5_0_pre5#0018

                h_RBXNoisy_NoisyNmFr_NoisyEnFr->Fill(EnFr, NmFr);
                h_RBXNoisy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXNoisy_En->Fill(RBXEnergy15->at(i));
	    
	  } else {

                h_RBXHealthy_NoisyNmFr_NoisyEnFr->Fill(EnFr, NmFr);
                h_RBXHealthy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXHealthy_En->Fill(RBXEnergy15->at(i));
	    
	    }
               }



               bool highestRBXEn=false;
               double maxRBXEnergy=0;
              for (int i = 0; i < 72; i++) {
              if (RBXEnergy15->at(i) > maxRBXEnergy){
              maxRBXEnergy=RBXEnergy15->at(i);

               highestRBXEn=true;
          }
            
     }    
                hAllRBXs_En->Fill(maxRBXEnergy);   

                if (RBXRechitR45Loose == true ){
                hloose_En->Fill(maxRBXEnergy);}

                if (RBXRechitR45Tight == true ){
                htight_En->Fill(maxRBXEnergy);}

               if (RBXRechitR45TightNew== true )
 {
              hmodified_tightEn->Fill(maxRBXEnergy);}
               
                                  
                     
              if (RBXRechitR45LooseNew == true )
          
             {
              hmodified_looseEn->Fill(maxRBXEnergy);}
               

            holdrechitR45loose_vs_modifiedrechitR45loose->Fill(RBXRechitR45Loose,RBXRechitR45LooseNew);
            holdrechitR45tight_vs_modifiedrechitR45tight->Fill(RBXRechitR45Tight,RBXRechitR45TightNew);
  

  

         
      //  }// i over RBX
       
            


   //-----------------------------------------------------------------------------------------------------------

//RBX
        /*   for (int j = 0; j < 72; j++) {
            
            double EnFr = 0.;
            double NmFr = 0.;

            EnFr = RecHitNoisy_ERH[j] / ERH[j];
            NmFr = RecHitNoisy_NRH[j] / NRH[j];
           


           hAllRBXs_En-> Fill(RBXEnergy15->at(j));
      

               double C4 = (*RBXCharge)[j][4];
               double C5 = (*RBXCharge)[j][5];
               double P45 = (C4 + C5);
              
               if(P45 < 1) P45 = 1; 
              double R45 = (C4 - C5) / P45;

               double C3 = (*RBXCharge)[j][3];
            
               double P34 = (C3 + C4);
               
              if(P34 < 1) P34 = 1; 
              double R34 = (C3- C4) / P34;
              
              

             if   (0.9< R45 && R45 <1.1 && 100<P45  ) {
            
             hRBXswithR45RBX1_En->Fill(RBXEnergy15->at(j));
          //   h_R45RBX->Fill(P45, R45);
          //   h_R34RBX->Fill(P34, R34);
             
      

            h_RBX_RBXPhi_RBXEta_RBX->Fill(RBX_X(j), RBX_Y(j));
             h_RBX_NmFr_EnFr_RBX->Fill(EnFr, NmFr); 
           }         
 }
*/
//
//--------------------------------------


  

	if (RBXR45Noise_EnergyGt50 && RBXRechitR45Tight) badEvt_tightTag++;
	if (!RBXR45Noise_EnergyGt50 && !RBXRechitR45Tight) goodEvt_notTag++;
 	 
            h_RBXR45Noise_vs_HasBadRBXR45->Fill(HasBadRBXR45->at(0),RBXR45Noise_EnergyGt50);
            h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose->Fill(HasBadRBXRechitR45Loose->at(0),RBXRechitR45Loose);
            h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight->Fill(HasBadRBXRechitR45Tight->at(0),RBXRechitR45Tight);


           
            
          
		if(RBXRechitR45Loose == false) hMETloose->Fill(met);
		if(RBXRechitR45Tight == false) hMETtight->Fill(met);

		if(RBXR45Noise_EnergyGt50) MET_RBXR45Noise->Fill(met);
		else MET_healthyEvents->Fill(met);


	/*if((HasBadRBXR45->at(0) && !RBXR45Noise_EnergyGt50) || (!HasBadRBXR45->at(0) && RBXR45Noise_EnergyGt50))
		for (int i = 0; i < 72; i++) {
		if(RBXEnergy15->at(i) > 50) {
	        cout<<"\nRBXEnergy15->at("<<i<<"): "<<RBXEnergy15->at(i)<<endl;
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
	if(HasBadRBXR45->at(0) && !RBXR45Noise_EnergyGt50) cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is NOT recognized as NOISE by our R45 definition !! *** \n\n"<<endl;
	if(!HasBadRBXR45->at(0) && RBXR45Noise_EnergyGt50) cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is recognized as NOISE by our R45 definition !! *** \n\n"<<endl;


	if((HasBadRBXRechitR45Tight->at(0) && !RBXRechitR45Tight) || (!HasBadRBXRechitR45Tight->at(0) && RBXRechitR45Tight) || (HasBadRBXRechitR45Loose->at(0) && !RBXRechitR45Loose) || (!HasBadRBXRechitR45Loose->at(0) && RBXRechitR45Loose)) { 
        for (unsigned int i = 0; i < HBHERecHitRBXid->size(); i++) {

            double C4 = (*HBHERecHitAuxFC)[i][4];
            double C5 = (*HBHERecHitAuxFC)[i][5];

	    if(!(HBHERecHitEnergyRaw->at(i) > 5)) continue;
	    cout<<"\nRecHitIEta: "<<HBHERecHitIEta->at(i)<<endl;
	    cout<<"pass threshold, RBX: "<<HBHERecHitRBXid->at(i)<<", HBHERecHitEnergyRaw: "<<HBHERecHitEnergyRaw->at(i)<<endl;

            if (isRecHitNoisy(C4, C5, h_RecHit_R45_P45) && (fabs(HBHERecHitIEta->at(i))!=28 && fabs(HBHERecHitIEta->at(i))!=29)) {
	        cout<<"\nRBXEnergy15->at("<<HBHERecHitRBXid->at(i)<<"): "<<RBXEnergy15->at(HBHERecHitRBXid->at(i))<<endl;
                cout<<"(*HBHERecHitAuxFC)[i][4]: "<<C4<<endl;
                cout<<"(*HBHERecHitAuxFC)[i][5]: "<<C5<<endl;
		double P45 = C4+C5;
		double R45 = (C4-C5) / P45;
		cout<<"\n                       R45: "<<R45<<endl;
		cout<<"                     ^^^^^^^^^ Lower Limit: "<<LowerLimitRecHit(P45)<<endl;
		cout<<"                     ^^^^^^^^^ Upper Limit: "<<UpperLimit(P45)<<endl;
            }
        }
        for (int i = 0; i < 72; i++) {

            double EnFr = RecHitNoisy_ERH[i] / ERH[i];
            double NmFr = RecHitNoisy_NRH[i] / NRH[i];

        if((NmFr!=0 || EnFr!=0) && (ERH[i] != 0 && NRH[i] != 0)) cout<<"\nNmFr: "<<NmFr<<", EnFr: "<<EnFr<<endl;

		if(isRBXNoisy((*RBXCharge)[i][4], (*RBXCharge)[i][5], h_RBX_R45_P45) && RBXEnergy15->at(i) > 50) {
		cout<<"\n --------- RBX NOISY ---------"<<endl;
	        cout<<"\nRBXEnergy15->at("<<i<<"): "<<RBXEnergy15->at(i)<<endl;
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

	}
	
	if(HasBadRBXRechitR45Tight->at(0) && !RBXRechitR45Tight) 
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is NOT recognized as TIGHT NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	
	if(HasBadRBXRechitR45Loose->at(0) && !RBXRechitR45Loose) 
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is NOT recognized as LOOSE NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	
	if(!HasBadRBXRechitR45Tight->at(0) && RBXRechitR45Tight) 
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is recognized as TIGHT NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	
	if(!HasBadRBXRechitR45Loose->at(0) && RBXRechitR45Loose) 
		cout<<"\n*** !! event "<<jentry<<", event: "<<event<<", run: "<<run<<" is recognized as LOOSE NOISE by our R45 definition !! ***   "<<HPDHits->at(0)<<", "<<HPDNoOtherHits->at(0)<<"\n\n"<<endl;
	*/


    }//jentry	
}
	cout<<nevts<<endl;
        cout<<nevts2<<endl;
	cout<<hAllRBXs_En->Integral()<<endl;
	cout<<"nr of badEvt_tightTag: "<<badEvt_tightTag<<", and eff: "<<(double)badEvt_tightTag/MET_RBXR45Noise->GetEntries()<<endl;
	cout<<"nr of goodEvt_notTag: "<<goodEvt_notTag<<", and eff: "<<(double)goodEvt_notTag/MET_healthyEvents->GetEntries()<<endl;


    // ======= E N D =======
    fout.Write();
    fout.Close();


    return 0;
}
