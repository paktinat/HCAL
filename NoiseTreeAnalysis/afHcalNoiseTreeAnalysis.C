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
#include "TFileCollection.h"
#include "TChain.h"
//---------------------------------------------------------------------------
#include "TreeReader_231228.h"





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

Float_t UpperLimit(Float_t energy) {

    Float_t TS4TS5UpperThreshold[5] = {70, 90, 100, 400, 4000};
    Float_t TS4TS5UpperCut[5] = {1, 0.8, 0.75, 0.72, 0.72};

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

Float_t LowerLimitRecHit(Float_t energy) {
    //R45 per RecHit RecoLocalCalo/HcalRecProducers/python/HcalHitReconstructor_hbhe_cfi.py:
    Float_t TS4TS5LowerThreshold[7] = {100, 120, 160, 200, 300, 500, 4000};
    Float_t TS4TS5LowerCut[7] = {-1, -0.7, -0.5, -0.4, -0.3, 0.1, 0.1};

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

Float_t LowerLimitRBX(Float_t energy) {
    //R45 per RBX RecoMET/METProducers/python/hcalnoiseinfoproducer_cfi.py:
    Float_t TS4TS5LowerThreshold[8] = {100, 120, 150, 200, 300, 400, 500, 4000};
    Float_t TS4TS5LowerCut[8] = {-1, -0.7, -0.4, -0.2, -0.08, 0, 0.1, 0.1};

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

Float_t DeltaPhi(Float_t v1, Float_t v2) {
    Float_t diff = fabs(v2 - v1);
    Float_t corr = 2 * acos(-1.) - diff;
    if (diff < acos(-1.)) {
        return diff;
    } else {

        return corr;
    }
}

Float_t GetDeltaR(Float_t eta1, Float_t eta2, Float_t phi1, Float_t phi2) {

    return sqrt((eta1 - eta2)*(eta1 - eta2) + DeltaPhi(phi1, phi2) * DeltaPhi(phi1, phi2));
}

Int_t isRBXNoisy(double c4, double c5, TH2D* h) {
    Int_t isRBXNsy = 0;
    Float_t P45 = (c4 + c5);
    Float_t R45 = (c4 - c5) / P45;
    h->Fill(P45, R45);
    if ((R45 < LowerLimitRBX(P45) || R45 > UpperLimit(P45))) isRBXNsy = 1;
    return isRBXNsy;
}

Int_t isRecHitNoisy(double c4, double c5, TH2D* h) {
    Int_t isRHNsy = 0;
    Float_t P45 = (c4 + c5);
    Float_t R45 = (c4 - c5) / P45;
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


//To run "root afHcalNoiseTreeAnalysis.C'("Test")'   "
//---------------------------------------------------------------------------

int afHcalNoiseTreeAnalysis(TString filelist = "TEST", bool addPU=false) {

  float TS5PU[62][73];
  float TS6PU[62][73];

  if(addPU){
    for(int i = 0; i < 62; i++){
      for(int j = 0; j < 73; j++){
	TS5PU[i][j] = -123.45;
	TS6PU[i][j] = -123.45;
      }}

    TChain ch("ExportTree/HcalNoiseTree");

//     TString fIn = "/dataLOCAL/HCAL/RootFiles/Seema/NoiseTree_Commissionig2014_MinimumBias_v3_229684_"+filelist+".root";
    TString fIn = "/dataLOCAL/HCAL/RootFiles/HcalNoiseTree_V00-03-12142_537p6_FT_P_V42D_HBHEYesHFNoHONo_000004/NoiseTree_*.root";

    ch.Add(fIn);

    TreeReader tr(&ch);

    if (tr.fChain == 0) return 0;

    // ========== B E G I N ==========
    //Long64_t nentries = tr.fChain->GetEntries();
    Long64_t nentries = 55000;//This number is sufficient to find the pu for all of the channels.

    TString CurrentFile = "";
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = tr.LoadTree(jentry);
        if (ientry < 0) break;
        tr.fChain->GetEntry(jentry);

        prn(jentry, nentries, 100, CurrentFile, tr.fChain);

//         if (!(tr.OfficialDecision == 1)) continue;//1 no error, 0 error

        // ========== LOOP OVER Channels ==========
        for (int i = 0; i < tr.PulseCount; i++) {
            Int_t id = tr.Hit_RBXIndex[i];
       
	    if (id > 71) continue;
            if (id > 35) continue;
            if (id == 18) continue;
	    
	    if (!(tr.Energy[i] > 2.5))continue;

	    int rediEta = tr.IEta[i] + 30;
	    
	    if(TS5PU[rediEta][tr.IPhi[i]] > -100.0) continue;
	      
// 	    cout<<" iEta[i] "<<tr.IEta[i]<<endl;
// 	    cout<<" rediEta "<<rediEta<<endl;

	    TS5PU[rediEta][tr.IPhi[i]] = tr.Charge[i][5];
	    TS6PU[rediEta][tr.IPhi[i]] = tr.Charge[i][6];
	}//========== LOOP OVER Channels ==========
    }//jentry
  }//if(addPU)
    

    TChain ch("ExportTree/HcalNoiseTree");

    TString fIn = "/dataLOCAL/HCAL/RootFiles/Seema/HcalNoiseTree_HcalHPDNoise_229684_"+filelist+"_*.root";

    ch.Add(fIn);

    TString sFOut = "res_";
    sFOut.Append(filelist);

    if(addPU) 
      sFOut.Append("_ManPU");
    
    sFOut.Append(".root");
    TFile fout(sFOut, "RECREATE");

    TH2D* h_RBX_R45_P45 = new TH2D("RBX_R45_P45", "RBX_R45_P45", 100, 0, 1000, 30, -1.5, 1.5);
    TH2D* h_RecHit_R45_P45 = new TH2D("RecHit_R45_P45", "RecHit_R45_P45", 100, 0, 1000, 30, -1.5, 1.5);

    TH2D* h_RBXNoisy_NoisyNmFr_NoisyEnFr = new TH2D("RBXNoisy_NoisyNmFr_NoisyEnFr", "RBXNoisy_NoisyNmFr_NoisyEnFr", 20, 0, 1.01, 20, 0, 1.01);
    TH2D* h_RBXHealthy_NoisyNmFr_NoisyEnFr = new TH2D("RBXHealthy_NoisyNmFr_NoisyEnFr", "RBXHealthy_NoisyNmFr_NoisyEnFr", 20, 0, 1.1, 20, 0, 1.1);

    TH2D* h_RBXNoisy_RBXPhi_RBXEta = new TH2D("RBXNoisy_RBXPhi_RBXEta", "RBXNoisy_RBXPhi_RBXEta", 4, 0, 4, 18, 0, 18);
    TH2D* h_RBXHealthy_RBXPhi_RBXEta = new TH2D("RBXHealthy_RBXPhi_RBXEta", "RBXHealthy_RBXPhi_RBXEta", 4, 0, 4, 18, 0, 18);

    TH1D* h_RBXNoisy_NRH = new TH1D("RBXNoisy_NRH", "RBXNoisy_NRH", 72, 0, 72);
    TH1D* h_RecHitNoisy_NRH = new TH1D("RecHitNoisy_NRH", "RecHitNoisy_NRH", 72, 0, 72);
    TH1D* h_EpsNoisy_NRH = new TH1D("EpsNoisy_NRH", "EpsNoisy_NRH", 72, 0, 72);

    TH1D* h_RBXHealthy_NRH = new TH1D("RBXHealthy_NRH", "RBXHealthy_NRH", 72, 0, 72);
    TH1D* h_RecHitHealthy_NRH = new TH1D("RecHitHealthy_NRH", "RecHitHealthy_NRH", 72, 0, 72);
    TH1D* h_EpsHealthy_NRH = new TH1D("EpsHealthy_NRH", "EpsHealthy_NRH", 72, 0, 72);

    TH1D* h_TEST_TS = new TH1D("TEST_TS", "TEST_TS", 10, 0, 10);
    //    TH1D* h_NRecHit_TS = new TH1D("TEST_TS", "TEST_TS", 10, 0, 10);

    TH1D* h_RecHitEnergy = new TH1D("RecHitEnergy", "RecHitEnergy", 200,0,200);

    TH1D * h_PS_TS[10];

    for (int i = 0; i < 10; i++) {
        TString s_ps = "PS_TS_";
        stringstream ss;
        ss << i;
        s_ps.Append(ss.str());
        h_PS_TS[i] = new TH1D(s_ps, s_ps, 10, 0, 10);
    }


    TH1D* h_tmp;

    TreeReader tr(&ch);

    if (tr.fChain == 0) return 0;

    // ========== B E G I N ==========
    Long64_t nentries = tr.fChain->GetEntries();
    //Long64_t nentries = 100000;

    int i_ps(0);

    TString CurrentFile = "";
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = tr.LoadTree(jentry);
        if (ientry < 0) break;
        tr.fChain->GetEntry(jentry);

        prn(jentry, nentries, 100, CurrentFile, tr.fChain);

//         if (!(tr.OfficialDecision == 1)) continue;//1 no error, 0 error

        Float_t NRH[72] = {};
        Float_t RecHitNoisy_NRH[72] = {};
        Float_t ERH[72] = {};
        Float_t RecHitNoisy_ERH[72] = {};

        Float_t TS4RBX[72] = {};
        Float_t TS5RBX[72] = {};

        // ========== LOOP OVER Channels ==========
        for (int i = 0; i < tr.PulseCount; i++) {
            Int_t id = tr.Hit_RBXIndex[i];

            if (id > 71) continue;
            if (id > 35) continue;
            if (id == 18) continue;

	    h_RecHitEnergy->Fill(tr.Energy[i]);

            //            if (!(tr.Energy[i] > 1 && tr.Energy[i] < 10))continue;
	    if (!(tr.Energy[i] > 1.5))continue;
	    if (!(tr.Energy[i] > 7.5))continue;

            NRH[id] += 1.0;
            ERH[id] += tr.Energy[i];

            double C4 = tr.Charge[i][4];
            double C5 = tr.Charge[i][5];

	    int ii_ps(0);
            for (int ii = 0; ii < 10; ii++) {
                h_TEST_TS->Fill(ii, tr.Charge[i][ii]);
                if (i_ps < 10 && (C4 + C5) > 50) {
                    h_PS_TS[i_ps]->Fill(ii, tr.Charge[i][ii]);
                    ii_ps = 1;
                }
            }
            if (ii_ps == 1) i_ps++;

	    if(addPU){
	      int rediEta = tr.IEta[i] + 30;
	      C4 += TS5PU[rediEta][tr.IPhi[i]];
	      C5 += TS6PU[rediEta][tr.IPhi[i]];
	    }

	    TS4RBX[id] += C4;
	    TS5RBX[id] += C5;
	    
            if (isRecHitNoisy(C4, C5, h_RecHit_R45_P45)) {
                RecHitNoisy_NRH[id] += 1.0;
                RecHitNoisy_ERH[id] += tr.Energy[i];
            }
        }// i over channels

        for (int i = 0; i < 72; i++) {
            if (ERH[i] == 0 || NRH[i] == 0) continue;

            double EnFr = RecHitNoisy_ERH[i] / ERH[i];
            double NmFr = RecHitNoisy_NRH[i] / NRH[i];
            double C4 = TS4RBX[i];//tr.RBXCharge[i][4];
            double C5 = TS5RBX[i];//tr.RBXCharge[i][5];

            if (isRBXNoisy(C4, C5, h_RBX_R45_P45)) {
                h_RBXNoisy_NoisyNmFr_NoisyEnFr->Fill(EnFr, NmFr);
                h_RBXNoisy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXNoisy_NRH->Fill(NRH[i]);

		if (EnFr >= 0.2 || NmFr >= 0.2) 
		  h_RecHitNoisy_NRH->Fill(NRH[i]);
	    
	    } else {
                h_RBXHealthy_NoisyNmFr_NoisyEnFr->Fill(EnFr, NmFr);
                h_RBXHealthy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXHealthy_NRH->Fill(NRH[i]);

		if (EnFr < 0.2 && NmFr < 0.2) 
		  h_RecHitHealthy_NRH->Fill(NRH[i]);
	    
	    }

	      
//             if (EnFr < 0.2 && NmFr < 0.2) {
//                 h_RecHitHealthy_NRH->Fill(NRH[i]);
//             } else {
//                 h_RecHitNoisy_NRH->Fill(NRH[i]);
//             }

        }// i over RBX

        h_tmp = (TH1D*) h_RBXNoisy_NRH->Clone("tmp");
        h_tmp->Add(h_RecHitNoisy_NRH, -1);
        h_EpsNoisy_NRH->Divide(h_RBXNoisy_NRH, h_tmp);
        delete h_tmp;

        h_EpsHealthy_NRH->Divide(h_RecHitHealthy_NRH, h_RBXHealthy_NRH);

    }//jentry	

    // ======= E N D =======
    fout.Write();
    fout.Close();


    return 0;
}




