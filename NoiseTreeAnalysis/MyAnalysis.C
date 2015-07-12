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


//To run "root MyAnalysis.C'("Test")'   "
//---------------------------------------------------------------------------

int MyAnalysis(TString filelist = "TEST", bool addPU=false) {

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
    TH1D* h_eff = new TH1D("eff", "eff", 72, 0, 72);

    TH1D* h_RBXHealthy_NRH = new TH1D("RBXHealthy_NRH", "RBXHealthy_NRH", 72, 0, 72);
    TH1D* h_RecHitHealthy_NRH = new TH1D("RecHitHealthy_NRH", "RecHitHealthy_NRH", 72, 0, 72);
    TH1D* h_EpsHealthy_NRH = new TH1D("EpsHealthy_NRH", "EpsHealthy_NRH", 72, 0, 72);

    TH1D* h_TEST_TS = new TH1D("TEST_TS", "TEST_TS", 10, 0, 10);
    //    TH1D* h_NRecHit_TS = new TH1D("TEST_TS", "TEST_TS", 10, 0, 10);


    TH1F* h_RecHitEnergy = new TH1F("RecHitEnergy", "RecHitEnergy", 60,-10,50);
    TH1F* h_ERH = new TH1F("RecHitEn_perRBX", "RecHitEn_perRBX", 70,-20,50);
    TH1F* h_RBXEnergy = new TH1F("RBXEnergy", "RBXEnergy", 70,-20,50);

    TH2I* h_RBX_X_Y = new TH2I("RBX_X_Y", "RBX_X_Y", 4, 0, 4, 17, 0, 17);
    h_RBX_X_Y->GetXaxis()->SetBinLabel(1,"HE-");
    h_RBX_X_Y->GetXaxis()->SetBinLabel(2,"HB-");
    h_RBX_X_Y->GetXaxis()->SetBinLabel(3,"HB+");
    h_RBX_X_Y->GetXaxis()->SetBinLabel(4,"HE+");

    TH1D* h_Diff_RBXEn_RHSumEn = new TH1D("Diff_RBXEn_RHSumEn", "Diff_RBXEn_RHSumEn", 100, -10, 10);
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
    h_counting->GetXaxis()->SetBinLabel(1,"#varepsilon_{RBX}-#varepsilon_{#Sigma RH} <= 0.001");
    h_counting->GetXaxis()->SetBinLabel(2,"#varepsilon_{RBX}-#varepsilon_{#Sigma RH} > 0.001");
    h_counting->GetXaxis()->SetLabelSize(0.11);

    TH1D * h_PS_TS[10];

    for (int j = 0; j < 10; j++) {
        TString s_ps = "PS_TS_";
        stringstream ss;
        ss << j;
        s_ps.Append(ss.str());
        h_PS_TS[j] = new TH1D(s_ps, s_ps, 10, 0, 10);
    }

    TH1D* h_tmp;
    TH2I* h2D = new TH2I("official", "R45Noise", 2, 0, 2, 2, 0, 2);
    TH2I* h2D_Tight = new TH2I("official_Tight", "R45Noise", 2, 0, 2, 2, 0, 2);

    TH1D* ntupleMET = new TH1D("savedMET","savedMET", 100, 0, 100);
    TH1D* hMET[4];

    TString metThresholds[4] = {"noCut","1p5","3","5"};

    for (int j = 0; j < 4; j++) {
        TString histoname = "MET_";
        histoname.Append(metThresholds[j]);
        hMET[j] = new TH1D(histoname,histoname, 100, 0, 100);
    }

	vector<std::string> inputfile;

	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_2402393C_D869_E411_89A7_02163E010EC9.root");
	/*inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_3AA6CA83_D669_E411_8BBB_02163E010CAB.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_466D8789_D769_E411_AA9F_02163E00FFC8.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise44404E2B_DA69_E411_9B5B_02163E010F9B.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_30F784B4_D669_E411_BE94_02163E010DEF.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_3685A126_D769_E411_A868_02163E010E36.root");*/

/*
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_2402393Cto3685A126.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_3AA6CA83to466D8789.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_48427151to566D74B5.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_5AAADD05to70F318A1.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_76E91BFEto861B5538.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_8CCAACA7toACB03CF4.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_AEE04716toCC6960E7.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_E0090964toF2CDE4A6.root");
	inputfile.push_back("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_F44047C5-D769to1A5D12E0-D769.root");
*/

	analysisClass * halilCode = new analysisClass();

	for (unsigned int file=0;file<inputfile.size();file++) {

	TFile *f = TFile::Open(inputfile[file].c_str()); 
	TTree *t; f->GetObject("hcalTupleTree/tree",t);


    // ========== B E G I N ==========
    //Long64_t nentries = t->GetEntries();
    Long64_t nentries = 1;

    unsigned int event = 0;
    unsigned int run = 0;
    std::vector<float> * HBHERecHitEnergy = 0;
    std::vector<float> * HBHERecHitEnergyMethod0 = 0;
    std::vector<int> * HBHERecHitIEta = 0;
    std::vector<int> * HBHERecHitIPhi = 0;
    std::vector<float> * HBHERecHitEta = 0;
    std::vector<float> * HBHERecHitPhi = 0;
    std::vector<float> * HBHERecHitDepth = 0;
    std::vector< std::vector<float> > *HBHEDigiFC = 0;
    std::vector< std::vector<float> > *HBHEDigiEnergy = 0;
    std::vector<float> *HBHEDigiRecEnergy = 0;
    std::vector<float> *HBHEDigiEta = 0;
    std::vector<double> *RBXEnergy15 = 0;
    std::vector< std::vector<double> > *RBXCharge = 0;
    std::vector<int>     *HasBadRBXRechitR45Loose = 0;
    std::vector<int>     *HasBadRBXRechitR45Tight = 0;
    std::vector<double> *HBET = 0;
    std::vector<double> *HEET = 0;


    t->SetBranchAddress("event",&event);
    t->SetBranchAddress("run",&run);
    t->SetBranchAddress("HBHERecHitEnergy",&HBHERecHitEnergy);
    t->SetBranchAddress("HBHERecHitEnergyMethod0",&HBHERecHitEnergyMethod0);
    t->SetBranchAddress("HBHERecHitIEta",&HBHERecHitIEta);
    t->SetBranchAddress("HBHERecHitIPhi",&HBHERecHitIPhi);
    t->SetBranchAddress("HBHERecHitEta",&HBHERecHitEta);
    t->SetBranchAddress("HBHERecHitPhi",&HBHERecHitPhi);
    t->SetBranchAddress("HBHERecHitDepth",&HBHERecHitDepth);
    t->SetBranchAddress("HBHEDigiFC",&HBHEDigiFC);
    t->SetBranchAddress("HBHEDigiEnergy",&HBHEDigiEnergy);
    t->SetBranchAddress("HBHEDigiRecEnergy",&HBHEDigiRecEnergy);
    t->SetBranchAddress("HBHEDigiEta",&HBHEDigiEta);
    t->SetBranchAddress("RBXEnergy15",&RBXEnergy15);
    t->SetBranchAddress("RBXCharge",&RBXCharge);
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

        Float_t NRH[72] = {};
        Float_t RecHitNoisy_NRH[72] = {};
        Float_t ERH[72] = {};
        Float_t RecHitNoisy_ERH[72] = {};

        Float_t TS4RBX[72] = {};
        Float_t TS5RBX[72] = {};

	double MET[4] = {0};
	double METx[4] = {0};
	double METy[4] = {0};

        // ========== LOOP OVER Channels ==========

        for (unsigned int i = 0; i < HBHERecHitEnergy->size(); i++) {

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

	    int RBXIndex = -10;		
/*
	    // --- this is the way I calculate the RBXIndex --- is not correct since the depth information are also needed to be used

	    int id = ((HBHERecHitIPhi->at(i)+1)/4)%18;

	    if (HBHERecHitIEta->at(i) >= 1 && HBHERecHitIEta->at(i) <= 16) RBXIndex = id;
	    if (HBHERecHitIEta->at(i) <= -1 && HBHERecHitIEta->at(i) >= -16) RBXIndex = id + 18;
	    if (HBHERecHitIEta->at(i) > 16) RBXIndex = id + 36;
	    if (HBHERecHitIEta->at(i) < -16) RBXIndex = id + 54;
*/

	    RBXIndex = halilCode->EtaPhitoRBX(HBHERecHitIEta->at(i),HBHERecHitIPhi->at(i),HBHERecHitDepth->at(i));

	    h_RecHitEnergy->Fill(HBHERecHitEnergy->at(i));

	    //if(!(HBHERecHitEnergy->at(i) > 1.5)) continue;

            NRH[RBXIndex] += 1.0;
            ERH[RBXIndex] += HBHERecHitEnergy->at(i);

            float C4 = (*HBHEDigiFC)[i][4];
            float C5 = (*HBHEDigiFC)[i][5];

	    int ii_ps(0);
            for (int ii = 0; ii < 10; ii++) {
                h_TEST_TS->Fill(ii, (*HBHEDigiFC)[i][ii]);
                if (i_ps < 10 && (C4 + C5) > 50) {
                    h_PS_TS[i_ps]->Fill(ii, (*HBHEDigiFC)[i][ii]);
                    ii_ps = 1;
                }
            }
            if (ii_ps == 1) i_ps++;

	    TS4RBX[RBXIndex] += C4;
	    TS5RBX[RBXIndex] += C5;
            if (isRecHitNoisy(C4, C5, h_RecHit_R45_P45)) {
                RecHitNoisy_NRH[RBXIndex] += 1.0;
                RecHitNoisy_ERH[RBXIndex] += HBHERecHitEnergy->at(i);
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

	    //h_RBX_X_Y->Fill(RBX_X(RBXIndex),RBX_Y(RBXIndex));

	bool R45Noise = false;

        for (int i = 0; i < 72; i++) {

	    //if (ERH[i] == 0 || NRH[i] == 0) continue;

            h_RBXEnergy->Fill(RBXEnergy15->at(i));
            h_ERH->Fill(ERH[i]);
	    
	    double diff = RBXEnergy15->at(i) - ERH[i];
	    h_Diff_RBXEn_RHSumEn->Fill(diff);		
	    if(fabs(diff) <= 0.001) h_counting->Fill(0); 
	    else if(fabs(diff) > 0.001) h_counting->Fill(1); 
	    h_EnDiff_Vs_RBXIndex->Fill(i,fabs(diff));

            double EnFr = RecHitNoisy_ERH[i] / ERH[i];
            double NmFr = RecHitNoisy_NRH[i] / NRH[i];
            float C4 = TS4RBX[i];
            float C5 = TS5RBX[i];

	    //cout<<"TS4RBX[i]: "<<C4<<" , from (*RBXCharge)[i][4]: "<<(*RBXCharge)[i][4]<<endl;

            if (isRBXNoisy(C4, C5, h_RBX_R45_P45)) {
                h_RBXNoisy_NoisyNmFr_NoisyEnFr->Fill(EnFr, NmFr);
                h_RBXNoisy_RBXPhi_RBXEta->Fill(RBX_X(i), RBX_Y(i));
                h_RBXNoisy_NRH->Fill(NRH[i]);

		if (EnFr >= 0.2 || NmFr >= 0.2) {
		  R45Noise = true;
		  h_RecHitNoisy_NRH->Fill(NRH[i]);
		}
	    
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
            h2D->Fill(HasBadRBXRechitR45Loose->at(0),R45Noise);
            h2D_Tight->Fill(HasBadRBXRechitR45Tight->at(0),R45Noise);

        h_tmp = (TH1D*) h_RBXNoisy_NRH->Clone("tmp");
        h_tmp->Add(h_RecHitNoisy_NRH, -1);
        h_EpsNoisy_NRH->Divide(h_RBXNoisy_NRH, h_tmp);
        delete h_tmp;

        h_eff->Divide(h_RecHitNoisy_NRH, h_RBXNoisy_NRH);
        h_EpsHealthy_NRH->Divide(h_RecHitHealthy_NRH, h_RBXHealthy_NRH);

    }//jentry	
}

    delete halilCode;

    // ======= E N D =======
    fout.Write();
    fout.Close();


    return 0;
}
