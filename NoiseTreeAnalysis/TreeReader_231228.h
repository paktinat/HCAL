//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 20 09:33:43 2015 by ROOT version 5.34/03
// from TTree HcalNoiseTree/Hcal noise tree version 1,2134
// found on file: /dataLOCAL/HCAL/RootFiles/Seema/NoiseTree_Commissionig2014_Cosmics_v3_231228.root
//////////////////////////////////////////////////////////

#ifndef TreeReader_h
#define TreeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Long64_t        RunNumber;
   Long64_t        EventNumber;
   Long64_t        LumiSection;
   Long64_t        Bunch;
   Long64_t        Orbit;
   Long64_t        Time;
   Bool_t          TTrigger[64];
   Bool_t          L1Trigger[128];
   Bool_t          HLTrigger[256];
   Double_t        EBET[2];
   Double_t        EEET[2];
   Double_t        HBET[2];
   Double_t        HEET[2];
   Double_t        HFET[2];
   Double_t        NominalMET[2];
   Double_t        EBSumE;
   Double_t        EESumE;
   Double_t        HBSumE;
   Double_t        HESumE;
   Double_t        HFSumE;
   Double_t        EBSumET;
   Double_t        EESumET;
   Double_t        HBSumET;
   Double_t        HESumET;
   Double_t        HFSumET;
   Int_t           NumberOfGoodTracks;
   Int_t           NumberOfGoodTracks15;
   Int_t           NumberOfGoodTracks30;
   Double_t        TotalPTTracks[2];
   Double_t        SumPTTracks;
   Double_t        SumPTracks;
   Int_t           NumberOfGoodPrimaryVertices;
   Int_t           NumberOfMuonCandidates;
   Int_t           NumberOfCosmicMuonCandidates;
   Int_t           PulseCount;
   Double_t        Charge[5184][10];
   Double_t        Pedestal[5184][10];
   Double_t        Energy[5184];
   Int_t           IEta[5184];
   Int_t           IPhi[5184];
   Int_t           Depth[5184];
   Int_t           Hit_RBXIndex[5184];
   Double_t        RecHitTime[5184];
   UInt_t          FlagWord[5184];
   UInt_t          AuxWord[5184];
   Double_t        RespCorrGain[5184];
   Double_t        fCorr[5184];
   Double_t        SamplesToAdd[5184];
   Double_t        RBXCharge[72][10];
   Double_t        RBXEnergy[72];
   Double_t        RBXCharge15[72][10];
   Double_t        RBXEnergy15[72];
   Int_t           HPDHits;
   Int_t           HPDNoOtherHits;
   Int_t           MaxZeros;
   Double_t        MinE2E10;
   Double_t        MaxE2E10;
   Bool_t          HasBadRBXR45;
   Bool_t          HasBadRBXRechitR45Loose;
   Bool_t          HasBadRBXRechitR45Tight;
   Double_t        LeadingJetEta;
   Double_t        LeadingJetPhi;
   Double_t        LeadingJetPt;
   Double_t        LeadingJetHad;
   Double_t        LeadingJetEM;
   Double_t        FollowingJetEta;
   Double_t        FollowingJetPhi;
   Double_t        FollowingJetPt;
   Double_t        FollowingJetHad;
   Double_t        FollowingJetEM;
   Int_t           JetCount20;
   Int_t           JetCount30;
   Int_t           JetCount50;
   Int_t           JetCount100;
   Double_t        HOMaxEnergyRing0;
   Double_t        HOSecondMaxEnergyRing0;
   Int_t           HOMaxEnergyIDRing0;
   Int_t           HOSecondMaxEnergyIDRing0;
   Int_t           HOHitCount100Ring0;
   Int_t           HOHitCount150Ring0;
   Double_t        HOMaxEnergyRing12;
   Double_t        HOSecondMaxEnergyRing12;
   Int_t           HOMaxEnergyIDRing12;
   Int_t           HOSecondMaxEnergyIDRing12;
   Int_t           HOHitCount100Ring12;
   Int_t           HOHitCount150Ring12;
   Bool_t          OfficialDecision;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_Bunch;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_Time;   //!
   TBranch        *b_TTrigger;   //!
   TBranch        *b_L1Trigger;   //!
   TBranch        *b_HLTrigger;   //!
   TBranch        *b_EBET;   //!
   TBranch        *b_EEET;   //!
   TBranch        *b_HBET;   //!
   TBranch        *b_HEET;   //!
   TBranch        *b_HFET;   //!
   TBranch        *b_NominalMET;   //!
   TBranch        *b_EBSumE;   //!
   TBranch        *b_EESumE;   //!
   TBranch        *b_HBSumE;   //!
   TBranch        *b_HESumE;   //!
   TBranch        *b_HFSumE;   //!
   TBranch        *b_EBSumET;   //!
   TBranch        *b_EESumET;   //!
   TBranch        *b_HBSumET;   //!
   TBranch        *b_HESumET;   //!
   TBranch        *b_HFSumET;   //!
   TBranch        *b_NumberOfGoodTracks;   //!
   TBranch        *b_NumberOfGoodTracks15;   //!
   TBranch        *b_NumberOfGoodTracks30;   //!
   TBranch        *b_TotalPTTracks;   //!
   TBranch        *b_SumPTTracks;   //!
   TBranch        *b_SumPTracks;   //!
   TBranch        *b_NumberOfGoodPrimaryVertices;   //!
   TBranch        *b_NumberOfMuonCandidates;   //!
   TBranch        *b_NumberOfCosmicMuonCandidates;   //!
   TBranch        *b_PulseCount;   //!
   TBranch        *b_Charge;   //!
   TBranch        *b_Pedestal;   //!
   TBranch        *b_Energy;   //!
   TBranch        *b_IEta;   //!
   TBranch        *b_IPhi;   //!
   TBranch        *b_Depth;   //!
   TBranch        *b_Hit_RBXIndex;   //!
   TBranch        *b_RecHitTime;   //!
   TBranch        *b_FlagWord;   //!
   TBranch        *b_AuxWord;   //!
   TBranch        *b_RespCorrGain;   //!
   TBranch        *b_fCorr;   //!
   TBranch        *b_SamplesToAdd;   //!
   TBranch        *b_RBXCharge;   //!
   TBranch        *b_RBXEnergy;   //!
   TBranch        *b_RBXCharge15;   //!
   TBranch        *b_RBXEnergy15;   //!
   TBranch        *b_HPDHits;   //!
   TBranch        *b_HPDNoOtherHits;   //!
   TBranch        *b_MaxZeros;   //!
   TBranch        *b_MinE2E10;   //!
   TBranch        *b_MaxE2E10;   //!
   TBranch        *b_HasBadRBXR45;   //!
   TBranch        *b_HasBadRBXRechitR45Loose;   //!
   TBranch        *b_HasBadRBXRechitR45Tight;   //!
   TBranch        *b_LeadingJetEta;   //!
   TBranch        *b_LeadingJetPhi;   //!
   TBranch        *b_LeadingJetPt;   //!
   TBranch        *b_LeadingJetHad;   //!
   TBranch        *b_LeadingJetEM;   //!
   TBranch        *b_FollowingJetEta;   //!
   TBranch        *b_FollowingJetPhi;   //!
   TBranch        *b_FollowingJetPt;   //!
   TBranch        *b_FollowingJetHad;   //!
   TBranch        *b_FollowingJetEM;   //!
   TBranch        *b_JetCount20;   //!
   TBranch        *b_JetCount30;   //!
   TBranch        *b_JetCount50;   //!
   TBranch        *b_JetCount100;   //!
   TBranch        *b_HOMaxEnergyRing0;   //!
   TBranch        *b_HOSecondMaxEnergyRing0;   //!
   TBranch        *b_HOMaxEnergyIDRing0;   //!
   TBranch        *b_HOSecondMaxEnergyIDRing0;   //!
   TBranch        *b_HOHitCount100Ring0;   //!
   TBranch        *b_HOHitCount150Ring0;   //!
   TBranch        *b_HOMaxEnergyRing12;   //!
   TBranch        *b_HOSecondMaxEnergyRing12;   //!
   TBranch        *b_HOMaxEnergyIDRing12;   //!
   TBranch        *b_HOSecondMaxEnergyIDRing12;   //!
   TBranch        *b_HOHitCount100Ring12;   //!
   TBranch        *b_HOHitCount150Ring12;   //!
   TBranch        *b_OfficialDecision;   //!

   TreeReader(TTree *tree=0);
   virtual ~TreeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(){};
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif
//
//#ifdef TreeReader_cxx
TreeReader::TreeReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/dataLOCAL/HCAL/RootFiles/Seema/NoiseTree_Commissionig2014_Cosmics_v3_231228.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/dataLOCAL/HCAL/RootFiles/Seema/NoiseTree_Commissionig2014_Cosmics_v3_231228.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/dataLOCAL/HCAL/RootFiles/Seema/NoiseTree_Commissionig2014_Cosmics_v3_231228.root:/ExportTree");
      dir->GetObject("HcalNoiseTree",tree);

   }
   Init(tree);
}

TreeReader::~TreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_LumiSection);
   fChain->SetBranchAddress("Bunch", &Bunch, &b_Bunch);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("Time", &Time, &b_Time);
   fChain->SetBranchAddress("TTrigger", TTrigger, &b_TTrigger);
   fChain->SetBranchAddress("L1Trigger", L1Trigger, &b_L1Trigger);
   fChain->SetBranchAddress("HLTrigger", HLTrigger, &b_HLTrigger);
   fChain->SetBranchAddress("EBET", EBET, &b_EBET);
   fChain->SetBranchAddress("EEET", EEET, &b_EEET);
   fChain->SetBranchAddress("HBET", HBET, &b_HBET);
   fChain->SetBranchAddress("HEET", HEET, &b_HEET);
   fChain->SetBranchAddress("HFET", HFET, &b_HFET);
   fChain->SetBranchAddress("NominalMET", NominalMET, &b_NominalMET);
   fChain->SetBranchAddress("EBSumE", &EBSumE, &b_EBSumE);
   fChain->SetBranchAddress("EESumE", &EESumE, &b_EESumE);
   fChain->SetBranchAddress("HBSumE", &HBSumE, &b_HBSumE);
   fChain->SetBranchAddress("HESumE", &HESumE, &b_HESumE);
   fChain->SetBranchAddress("HFSumE", &HFSumE, &b_HFSumE);
   fChain->SetBranchAddress("EBSumET", &EBSumET, &b_EBSumET);
   fChain->SetBranchAddress("EESumET", &EESumET, &b_EESumET);
   fChain->SetBranchAddress("HBSumET", &HBSumET, &b_HBSumET);
   fChain->SetBranchAddress("HESumET", &HESumET, &b_HESumET);
   fChain->SetBranchAddress("HFSumET", &HFSumET, &b_HFSumET);
   fChain->SetBranchAddress("NumberOfGoodTracks", &NumberOfGoodTracks, &b_NumberOfGoodTracks);
   fChain->SetBranchAddress("NumberOfGoodTracks15", &NumberOfGoodTracks15, &b_NumberOfGoodTracks15);
   fChain->SetBranchAddress("NumberOfGoodTracks30", &NumberOfGoodTracks30, &b_NumberOfGoodTracks30);
   fChain->SetBranchAddress("TotalPTTracks", TotalPTTracks, &b_TotalPTTracks);
   fChain->SetBranchAddress("SumPTTracks", &SumPTTracks, &b_SumPTTracks);
   fChain->SetBranchAddress("SumPTracks", &SumPTracks, &b_SumPTracks);
   fChain->SetBranchAddress("NumberOfGoodPrimaryVertices", &NumberOfGoodPrimaryVertices, &b_NumberOfGoodPrimaryVertices);
   fChain->SetBranchAddress("NumberOfMuonCandidates", &NumberOfMuonCandidates, &b_NumberOfMuonCandidates);
   fChain->SetBranchAddress("NumberOfCosmicMuonCandidates", &NumberOfCosmicMuonCandidates, &b_NumberOfCosmicMuonCandidates);
   fChain->SetBranchAddress("PulseCount", &PulseCount, &b_PulseCount);
   fChain->SetBranchAddress("Charge", Charge, &b_Charge);
   fChain->SetBranchAddress("Pedestal", Pedestal, &b_Pedestal);
   fChain->SetBranchAddress("Energy", Energy, &b_Energy);
   fChain->SetBranchAddress("IEta", IEta, &b_IEta);
   fChain->SetBranchAddress("IPhi", IPhi, &b_IPhi);
   fChain->SetBranchAddress("Depth", Depth, &b_Depth);
   fChain->SetBranchAddress("Hit_RBXIndex", Hit_RBXIndex, &b_Hit_RBXIndex);
   fChain->SetBranchAddress("RecHitTime", RecHitTime, &b_RecHitTime);
   fChain->SetBranchAddress("FlagWord", FlagWord, &b_FlagWord);
   fChain->SetBranchAddress("AuxWord", AuxWord, &b_AuxWord);
   fChain->SetBranchAddress("RespCorrGain", RespCorrGain, &b_RespCorrGain);
   fChain->SetBranchAddress("fCorr", fCorr, &b_fCorr);
   fChain->SetBranchAddress("SamplesToAdd", SamplesToAdd, &b_SamplesToAdd);
   fChain->SetBranchAddress("RBXCharge", RBXCharge, &b_RBXCharge);
   fChain->SetBranchAddress("RBXEnergy", RBXEnergy, &b_RBXEnergy);
   fChain->SetBranchAddress("RBXCharge15", RBXCharge15, &b_RBXCharge15);
   fChain->SetBranchAddress("RBXEnergy15", RBXEnergy15, &b_RBXEnergy15);
   fChain->SetBranchAddress("HPDHits", &HPDHits, &b_HPDHits);
   fChain->SetBranchAddress("HPDNoOtherHits", &HPDNoOtherHits, &b_HPDNoOtherHits);
   fChain->SetBranchAddress("MaxZeros", &MaxZeros, &b_MaxZeros);
   fChain->SetBranchAddress("MinE2E10", &MinE2E10, &b_MinE2E10);
   fChain->SetBranchAddress("MaxE2E10", &MaxE2E10, &b_MaxE2E10);
   fChain->SetBranchAddress("HasBadRBXR45", &HasBadRBXR45, &b_HasBadRBXR45);
   fChain->SetBranchAddress("HasBadRBXRechitR45Loose", &HasBadRBXRechitR45Loose, &b_HasBadRBXRechitR45Loose);
   fChain->SetBranchAddress("HasBadRBXRechitR45Tight", &HasBadRBXRechitR45Tight, &b_HasBadRBXRechitR45Tight);
   fChain->SetBranchAddress("LeadingJetEta", &LeadingJetEta, &b_LeadingJetEta);
   fChain->SetBranchAddress("LeadingJetPhi", &LeadingJetPhi, &b_LeadingJetPhi);
   fChain->SetBranchAddress("LeadingJetPt", &LeadingJetPt, &b_LeadingJetPt);
   fChain->SetBranchAddress("LeadingJetHad", &LeadingJetHad, &b_LeadingJetHad);
   fChain->SetBranchAddress("LeadingJetEM", &LeadingJetEM, &b_LeadingJetEM);
   fChain->SetBranchAddress("FollowingJetEta", &FollowingJetEta, &b_FollowingJetEta);
   fChain->SetBranchAddress("FollowingJetPhi", &FollowingJetPhi, &b_FollowingJetPhi);
   fChain->SetBranchAddress("FollowingJetPt", &FollowingJetPt, &b_FollowingJetPt);
   fChain->SetBranchAddress("FollowingJetHad", &FollowingJetHad, &b_FollowingJetHad);
   fChain->SetBranchAddress("FollowingJetEM", &FollowingJetEM, &b_FollowingJetEM);
   fChain->SetBranchAddress("JetCount20", &JetCount20, &b_JetCount20);
   fChain->SetBranchAddress("JetCount30", &JetCount30, &b_JetCount30);
   fChain->SetBranchAddress("JetCount50", &JetCount50, &b_JetCount50);
   fChain->SetBranchAddress("JetCount100", &JetCount100, &b_JetCount100);
   fChain->SetBranchAddress("HOMaxEnergyRing0", &HOMaxEnergyRing0, &b_HOMaxEnergyRing0);
   fChain->SetBranchAddress("HOSecondMaxEnergyRing0", &HOSecondMaxEnergyRing0, &b_HOSecondMaxEnergyRing0);
   fChain->SetBranchAddress("HOMaxEnergyIDRing0", &HOMaxEnergyIDRing0, &b_HOMaxEnergyIDRing0);
   fChain->SetBranchAddress("HOSecondMaxEnergyIDRing0", &HOSecondMaxEnergyIDRing0, &b_HOSecondMaxEnergyIDRing0);
   fChain->SetBranchAddress("HOHitCount100Ring0", &HOHitCount100Ring0, &b_HOHitCount100Ring0);
   fChain->SetBranchAddress("HOHitCount150Ring0", &HOHitCount150Ring0, &b_HOHitCount150Ring0);
   fChain->SetBranchAddress("HOMaxEnergyRing12", &HOMaxEnergyRing12, &b_HOMaxEnergyRing12);
   fChain->SetBranchAddress("HOSecondMaxEnergyRing12", &HOSecondMaxEnergyRing12, &b_HOSecondMaxEnergyRing12);
   fChain->SetBranchAddress("HOMaxEnergyIDRing12", &HOMaxEnergyIDRing12, &b_HOMaxEnergyIDRing12);
   fChain->SetBranchAddress("HOSecondMaxEnergyIDRing12", &HOSecondMaxEnergyIDRing12, &b_HOSecondMaxEnergyIDRing12);
   fChain->SetBranchAddress("HOHitCount100Ring12", &HOHitCount100Ring12, &b_HOHitCount100Ring12);
   fChain->SetBranchAddress("HOHitCount150Ring12", &HOHitCount150Ring12, &b_HOHitCount150Ring12);
   fChain->SetBranchAddress("OfficialDecision", &OfficialDecision, &b_OfficialDecision);
   Notify();
}

Bool_t TreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  if(entry > -100)
   return 1;
  
  return 1;
}
#endif // #ifdef TreeReader_cxx
