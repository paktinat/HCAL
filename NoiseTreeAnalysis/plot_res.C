{
  //gStyle->SetOptFit(110);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  gStyle->SetStatX(0.99);
  gStyle->SetStatY(0.99);
  gStyle->SetStatH(0.35);
  gStyle->SetStatW(0.25);

  TFile fa("res_test.root");

 //--- 2 ---
 can2 = new TCanvas("RES2");
 gPad->SetLogz();

    Float_t TS4TS5UpperThreshold[6] = {70, 90, 100, 400, 500, 4000};
    Float_t TS4TS5UpperCut[6] = {1, 0.8, 0.75, 0.72, 0.72, 0.72};
    TGraph *UpperGraph = new TGraph(6, TS4TS5UpperThreshold, TS4TS5UpperCut);
    Float_t TS4TS5LowerThreshold[7] = {100, 120, 160, 200, 300, 500, 4000};
    Float_t TS4TS5LowerCut[7] = {-1, -0.7, -0.5, -0.4, -0.3, 0.1, 0.1};
    Float_t LowerThresholdRBX[8] = {100, 120, 150, 200, 300, 400, 500, 4000};
    Float_t LowerCutRBX[8] = {-1, -0.7, -0.4, -0.2, -0.08, 0, 0.1, 0.1};
    TGraph *LowerGraph = new TGraph(7, TS4TS5LowerThreshold, TS4TS5LowerCut);
    TGraph *LowerGraphRBX = new TGraph(8, LowerThresholdRBX, LowerCutRBX);

    LowerGraph->SetLineWidth(2);
    LowerGraph->SetLineColor(2);
    LowerGraphRBX->SetLineWidth(2);
    LowerGraphRBX->SetLineColor(2);
    UpperGraph->SetLineWidth(2);
    UpperGraph->SetLineColor(2);

    RBX_R45_P45->SetTitle("RBX");	
    RBX_R45_P45->GetXaxis()->SetTitle("Charge45");	
    RBX_R45_P45->GetYaxis()->SetTitle("R45");	
    RBX_R45_P45->Draw("Colz");
    LowerGraphRBX->Draw("L");
    UpperGraph->Draw("L");

 //--- 3 ---
 can3 = new TCanvas("RES3");
 gPad->SetLogz();

    RecHit_R45_P45->SetTitle("Rechit");	
    RecHit_R45_P45->GetXaxis()->SetTitle("Charge45");	
    RecHit_R45_P45->GetYaxis()->SetTitle("R45");	
    RecHit_R45_P45->Draw("Colz");
    LowerGraph->Draw("L");
    UpperGraph->Draw("L");

 //--- 4 ---
 can4 = new TCanvas("RES4");
 gPad->SetLogz();

    RBXNoisy_NoisyNmFr_NoisyEnFr->Draw("TEXT Colz");
    RBXNoisy_NoisyNmFr_NoisyEnFr->SetTitle("RBX Noisy");	
    RBXNoisy_NoisyNmFr_NoisyEnFr->GetXaxis()->SetTitle("#frac{En of Noisy RH}{En of RH}");	
    RBXNoisy_NoisyNmFr_NoisyEnFr->GetYaxis()->SetTitle("#frac{Nr of Noisy RH}{Nr of RH}");	

 //--- 5 ---
 can5 = new TCanvas("RES5");
 gPad->SetLogz();

    RBXHealthy_NoisyNmFr_NoisyEnFr->Draw("TEXT Colz");
    RBXHealthy_NoisyNmFr_NoisyEnFr->SetTitle("RBX Healthy");	
    RBXHealthy_NoisyNmFr_NoisyEnFr->GetXaxis()->SetTitle("#frac{En of Noisy RH}{En of RH}");	
    RBXHealthy_NoisyNmFr_NoisyEnFr->GetYaxis()->SetTitle("#frac{Nr of Noisy RH}{Nr of RH}");	

 //--- 6 ---
 can6 = new TCanvas("RES6");
 gPad->SetLogz();

    RBXNoisy_RBXPhi_RBXEta->Draw("TEXT Colz");
    RBXNoisy_RBXPhi_RBXEta->SetTitle("RBX Noisy");	
    RBXNoisy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(1,"HE-");
    RBXNoisy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(2,"HB-");
    RBXNoisy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(3,"HB+");
    RBXNoisy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(4,"HE+");

 //--- 7 ---
 can7 = new TCanvas("RES7");
 gPad->SetLogz();

    RBXHealthy_RBXPhi_RBXEta->Draw("TEXT Colz");
    RBXHealthy_RBXPhi_RBXEta->SetTitle("RBX Healthy");	
    RBXHealthy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(1,"HE-");
    RBXHealthy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(2,"HB-");
    RBXHealthy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(3,"HB+");
    RBXHealthy_RBXPhi_RBXEta->GetXaxis()->SetBinLabel(4,"HE+");

 //--- 8 ---
 can8 = new TCanvas("RES8");

 gPad->SetLogy();
 RBXHealthy_En->SetTitle("RBX Energy");
 RBXHealthy_En->GetXaxis()->SetTitle("RBXEnergy15Method0");
 RBXHealthy_En->GetYaxis()->SetTitle("per RBX pre Event");
 RBXHealthy_En->SetLineColor(4); 
 RBXHealthy_En->SetLineWidth(2); 
 RBXHealthy_En->Draw();
 RBXNoisy_En->SetLineColor(2); 
 RBXNoisy_En->SetLineWidth(2); 
 RBXNoisy_En->SetLineStyle(7); 
 RBXNoisy_En->Draw("same");

 L1 = new TLegend(.4, .4, .89, .89);
 L1->SetBorderSize(0);
 L1->SetFillColor(0);
 L1->SetTextSize(0.04);
 L1->AddEntry(RBXHealthy_En,"Healthy RBXs (from R45)","L");
 L1->AddEntry(RBXNoisy_En,"Noisy RBXs (from R45)","L");
 L1->Draw("same");

 //--- 9 ---
 can9 = new TCanvas("RES9");
 gPad->SetLogz();

    h_RBXR45Noise_vs_HasBadRBXR45->SetTitle("Correlation Plot");	
    h_RBXR45Noise_vs_HasBadRBXR45->GetXaxis()->SetTitle("HasBadRBXR45Method0 (stored in ntuple)");	
    h_RBXR45Noise_vs_HasBadRBXR45->GetYaxis()->SetTitle("RBXR45Noise (RBXEnergy15Method0>50)");	
    h_RBXR45Noise_vs_HasBadRBXR45->Draw("TEXT Colz");

 //--- 10 ---
 can10 = new TCanvas("RES10");
 gPad->SetLogz();

    h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose->SetTitle("Correlation Plot");	
    h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose->GetXaxis()->SetTitle("HasBadRBXRechitR45LooseMethod0 (stored in ntuple)");	
    h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose->GetYaxis()->SetTitle("#frac{Nr of Noisy RH}{Nr of RH}>0.5 || #frac{En of Noisy RH}{En of RH}>0.5");	
    h_RBXRechitR45Loose_vs_HasBadRBXRechitR45Loose->Draw("TEXT Colz");

 //--- 11 ---
 can11 = new TCanvas("RES11");
 gPad->SetLogz();

    h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight->SetTitle("Correlation Plot");	
    h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight->GetXaxis()->SetTitle("HasBadRBXRechitR45TightMethod0 (stored in ntuple)");	
    h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight->GetYaxis()->SetTitle("#frac{Nr of Noisy RH}{Nr of RH}>0.2 || #frac{En of Noisy RH}{En of RH}>0.2");	
    h_RBXRechitR45Tight_vs_HasBadRBXRechitR45Tight->Draw("TEXT Colz");

 //--- 12 ---
 can12 = new TCanvas("RES12");

 gPad->SetLogy();
 MET_healthyEvents->SetTitle("Missing Transverse Energy");
 MET_healthyEvents->GetXaxis()->SetTitle("MET");
 MET_healthyEvents->GetYaxis()->SetTitle("Events");
 MET_healthyEvents->SetLineColor(4); 
 MET_healthyEvents->SetLineWidth(2); 
 MET_healthyEvents->Draw();
 MET_RBXR45Noise->SetLineColor(2); 
 MET_RBXR45Noise->SetLineWidth(2); 
 MET_RBXR45Noise->SetLineStyle(7); 
 MET_RBXR45Noise->Draw("same");

 L1 = new TLegend(.4, .4, .89, .89);
 L1->SetBorderSize(0);
 L1->SetFillColor(0);
 L1->SetTextSize(0.04);
 L1->AddEntry(MET_healthyEvents,"Events pass HasBadRBXR45Method0 filter","L");
 L1->AddEntry(MET_RBXR45Noise,"Events fail HasBadRBXR45Method0 filter","L");
 L1->Draw("same");

}


