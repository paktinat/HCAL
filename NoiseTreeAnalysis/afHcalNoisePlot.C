//To run "root afHcalNoisePlot.C res_Test.root"

#include <TH1.h>
#include <TGraph.h>

std::vector<TString> Args(TString startPoint) {
    std::vector< TString > sin;
    int argc = gApplication->Argc();
    for (int i = 0; i < argc; i++) {
        TString argvi = gApplication->Argv(i);
        if (argvi == startPoint)
            for (int j = i + 1; j < argc; j++) {
                TString argvj = gApplication->Argv(j);
                sin.push_back(argvj);
                i++;
            }
    }
    return sin;
}

void fixTH1D(TH1D* h, TString tt, TString yt, TString xt) {
    h->SetLineWidth(2);
    h->SetTitle(tt);
    h->SetYTitle(yt);
    h->SetXTitle(xt);
}

void fixTH2D(TH2D* h, TString tt, TString yt, TString xt) {
    h->SetTitle(tt);
    h->SetYTitle(yt);
    h->SetXTitle(xt);
}

void afHcalNoisePlot() {
    std::vector<TString> sin = Args("afHcalNoisePlot.C");

    TH1D* h_TS = (TH1D*) TFile::Open(sin[0])->Get("TEST_TS");
    h_TS->Scale(1 / h_TS->Integral());
    fixTH1D(h_TS, "Pulse Shape", "", "Time Slice");

    TH1D* h_EpsHealthy_NRH = (TH1D*) TFile::Open(sin[0])->Get("EpsHealthy_NRH");
    fixTH1D(h_EpsHealthy_NRH, "Epsilon Healthy", "new Healthy / old Healthy", "Number of RecHits");

    TH1D* h_EpsNoisy_NRH = (TH1D*) TFile::Open(sin[0])->Get("EpsNoisy_NRH");
    fixTH1D(h_EpsNoisy_NRH, "Epsilon Noisy", "Rejection Factor", "Number of RecHits");

    TH2D* h_RBXNoisy_RBXPhi_RBXEta = (TH2D*) TFile::Open(sin[0])->Get("RBXNoisy_RBXPhi_RBXEta");
    fixTH2D(h_RBXNoisy_RBXPhi_RBXEta, "RBX Noisy", "Phi", "Eta");

    TH2D* h_RBXHealthy_RBXPhi_RBXEta = (TH2D*) TFile::Open(sin[0])->Get("RBXHealthy_RBXPhi_RBXEta");
    fixTH2D(h_RBXHealthy_RBXPhi_RBXEta, "RBX Healthy", "Phi", "Eta");

    TH2D* h_RBXNoisy_NoisyNmFr_NoisyEnFr = (TH2D*) TFile::Open(sin[0])->Get("RBXNoisy_NoisyNmFr_NoisyEnFr");
    fixTH2D(h_RBXNoisy_NoisyNmFr_NoisyEnFr, "RBX Noisy", "Noisy RecHit Number Fraction", "Noisy RecHit Energy Fraction");

    TH2D* h_RBXHealthy_NoisyNmFr_NoisyEnFr = (TH2D*) TFile::Open(sin[0])->Get("RBXHealthy_NoisyNmFr_NoisyEnFr");
    fixTH2D(h_RBXHealthy_NoisyNmFr_NoisyEnFr, "RBX Healthy", "Noisy RecHit Number Fraction", "Noisy RecHit Energy Fraction");

    TH2D* h_RBX_R45_P45 = (TH2D*) TFile::Open(sin[0])->Get("RBX_R45_P45");
    fixTH2D(h_RBX_R45_P45, "RBX", "R45", "Charge45");

    TH2D* h_RecHit_R45_P45 = (TH2D*) TFile::Open(sin[0])->Get("RecHit_R45_P45");
    fixTH2D(h_RecHit_R45_P45, "RecHit", "R45", "Charge45");

    TH1D* h_PS_TS[10];
    for (int i = 0; i < 10; i++) {
      TString s_ps = "PS_TS_";
      stringstream ss;
      ss << i;
      s_ps.Append(ss.str());
      h_PS_TS[i] = (TH1D*) TFile::Open(sin[0])->Get(s_ps);
      fixTH1D(h_PS_TS[i], "Pulse Shape", "TS", "fC");
    }

    // ===  ===
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    c1->SetGrid();
    c1->cd();
    h_TS->Draw();

    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    c2->SetGrid();
    c2->cd();
    h_EpsHealthy_NRH->Draw();

    TCanvas *c3 = new TCanvas("c3", "", 800, 600);
    c3->SetGrid();
    c3->cd();
    h_EpsNoisy_NRH->Draw();

    TCanvas *c4 = new TCanvas("c4", "", 800, 600);
    c4->SetGrid();
    c4->cd();
    h_RBXNoisy_RBXPhi_RBXEta->Draw("TEXT Color");

    TCanvas *c5 = new TCanvas("c5", "", 800, 600);
    c5->SetGrid();
    c5->cd();
    h_RBXHealthy_RBXPhi_RBXEta->Draw("TEXT Color");

    TCanvas *c6 = new TCanvas("c6", "", 800, 600);
    c6->SetGrid();
    c6->cd();
    h_RBXNoisy_NoisyNmFr_NoisyEnFr->Draw("TEXT Color");

    TCanvas *c7 = new TCanvas("c7", "", 800, 600);
    c7->SetGrid();
    c7->cd();
    h_RBXHealthy_NoisyNmFr_NoisyEnFr->Draw("TEXT Color");

    TCanvas *c8 = new TCanvas("c8", "", 800, 600);
    c8->SetGrid();
    c8->cd();
    
    Float_t TS4TS5UpperThreshold[6] = {70, 90, 100, 400, 500, 1000};
    Float_t TS4TS5UpperCut[6] = {1, 0.8, 0.75, 0.72, 0.72, 0.72};
    TGraph *UpperGraph = new TGraph(6, TS4TS5UpperThreshold, TS4TS5UpperCut);
    Float_t TS4TS5LowerThreshold[7] = {100, 120, 160, 200, 300, 500, 1000};
    Float_t TS4TS5LowerCut[7] = {-1, -0.7, -0.5, -0.4, -0.3, 0.1, 0.1};
    TGraph *LowerGraph = new TGraph(7, TS4TS5LowerThreshold, TS4TS5LowerCut);

    LowerGraph->SetLineWidth(2);
    LowerGraph->SetLineColor(2);
    UpperGraph->SetLineWidth(2);
    UpperGraph->SetLineColor(2);

    h_RBX_R45_P45->Draw("Box");
    LowerGraph->Draw("L");
    UpperGraph->Draw("L");

    TCanvas *c9 = new TCanvas("c9", "", 800, 600);
    c9->SetGrid();
    c9->cd();
    h_RecHit_R45_P45->Draw("Box");
    LowerGraph->Draw("L");
    UpperGraph->Draw("L");

    
    TCanvas *c10 = new TCanvas("c10", "", 800, 600);
    c10->Divide(4,3);
    for (int i = 0; i < 10; i++) {
      c10->cd(i+1);
      h_PS_TS[i]->Draw();
    }

}
