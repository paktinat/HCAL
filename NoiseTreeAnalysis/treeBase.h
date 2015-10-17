//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Oct 17 15:53:06 2015 by ROOT version 6.02/05
// from TTree tree/
// found on file: Test.root
//////////////////////////////////////////////////////////

#ifndef treeBase_h
#define treeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class treeBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *EBET;
   vector<double>  *EBSumE;
   vector<double>  *EBSumET;
   vector<double>  *EEET;
   vector<double>  *EESumE;
   vector<double>  *EESumET;
   vector<double>  *HBET;
   vector<double>  *HBSumE;
   vector<double>  *HBSumET;
   vector<double>  *HEET;
   vector<double>  *HESumE;
   vector<double>  *HESumET;
   vector<double>  *HFET;
   vector<double>  *JetEMEB;
   vector<double>  *JetEMEE;
   vector<double>  *JetEMFrac;
   vector<double>  *JetEMHF;
   vector<double>  *JetEta;
   vector<double>  *JetHadFrac;
   vector<double>  *JetHadHB;
   vector<double>  *JetHadHE;
   vector<double>  *JetHadHF;
   vector<double>  *JetPhi;
   vector<double>  *JetPt;
   vector<double>  *NominalMET;
   vector<double>  *CDihitCluE;
   vector<double>  *CDihitIsolEcalE;
   vector<double>  *CDihitIsolHcalE;
   vector<double>  *CDihitIsolTrkE;
   vector<double>  *CDihitTrkFidE;
   vector<double>  *CHPDCluE;
   vector<double>  *CHPDIsolEcalE;
   vector<double>  *CHPDIsolHcalE;
   vector<double>  *CHPDIsolTrkE;
   vector<double>  *CHPDTrkFidE;
   vector<double>  *CMonohitCluE;
   vector<double>  *CMonohitIsolEcalE;
   vector<double>  *CMonohitIsolHcalE;
   vector<double>  *CMonohitIsolTrkE;
   vector<double>  *CMonohitTrkFidE;
   vector<double>  *CRBXCluE;
   vector<double>  *CRBXIsolEcalE;
   vector<double>  *CRBXIsolHcalE;
   vector<double>  *CRBXIsolTrkE;
   vector<double>  *CRBXTrkFidE;
   vector<double>  *HBHERecHitEnergyRaw;
   vector<double>  *IsolatedNoiseSumE;
   vector<double>  *IsolatedNoiseSumEt;
   vector<double>  *MaxE2E10;
   vector<double>  *MinE2E10;
   vector<double>  *NegativeNoiseSumE;
   vector<double>  *NegativeNoiseSumEt;
   vector<double>  *RBXEnergy;
   vector<double>  *RBXEnergy15;
   vector<double>  *SpikeNoiseSumE;
   vector<double>  *SpikeNoiseSumEt;
   vector<vector<double> > *HBHERecHitAuxAllfC;
   vector<vector<double> > *HBHERecHitAuxEnergy;
   vector<vector<double> > *HBHERecHitAuxFC;
   vector<vector<double> > *HBHERecHitAuxGain;
   vector<vector<double> > *HBHERecHitAuxPedFC;
   vector<vector<double> > *HBHERecHitAuxRCGain;
   vector<vector<double> > *RBXCharge;
   vector<vector<double> > *RBXCharge15;
   vector<float>   *HBHERecHitEnergy;
   vector<float>   *HBHERecHitEta;
   vector<float>   *HBHERecHitPhi;
   vector<float>   *HBHERecHitTime;
   vector<int>     *JetN60;
   vector<int>     *JetN90;
   vector<int>     *HBHERecHitAux;
   vector<int>     *HBHERecHitDepth;
   vector<int>     *HBHERecHitFlags;
   vector<int>     *HBHERecHitHPDid;
   vector<int>     *HBHERecHitIEta;
   vector<int>     *HBHERecHitIPhi;
   vector<int>     *HBHERecHitRBXid;
   vector<int>     *CDihitIsIso;
   vector<int>     *CDihitIsTagged;
   vector<int>     *CHPDId;
   vector<int>     *CHPDIsIso;
   vector<int>     *CHPDIsTagged;
   vector<int>     *CHPDNHits;
   vector<int>     *CHPDRBXId;
   vector<int>     *CMonohitIsIso;
   vector<int>     *CMonohitIsTagged;
   vector<int>     *CRBXId;
   vector<int>     *CRBXIsIso;
   vector<int>     *CRBXIsTagged;
   vector<int>     *CRBXNHits;
   vector<int>     *HPDHits;
   vector<int>     *HPDNoOtherHits;
   vector<int>     *HasBadRBXR45;
   vector<int>     *HasBadRBXRechitR45Loose;
   vector<int>     *HasBadRBXRechitR45Tight;
   vector<int>     *IsoNoiseFilterDecision;
   vector<int>     *MaxZeros;
   vector<int>     *NumIsolatedNoiseChannels;
   vector<int>     *NumNegativeNoiseChannels;
   vector<int>     *NumSpikeNoiseChannels;
   vector<int>     *OfficialDecision;
   vector<int>     *OfficialDecisionRun1;
   vector<int>     *OfficialDecisionRun2L;
   vector<int>     *OfficialDecisionRun2T;
   vector<vector<int> > *HBHERecHitAuxADC;
   vector<vector<int> > *HBHERecHitAuxCapID;
   UInt_t          event;
   UInt_t          ls;
   UInt_t          run;

   // List of branches
   TBranch        *b_EBET;   //!
   TBranch        *b_EBSumE;   //!
   TBranch        *b_EBSumET;   //!
   TBranch        *b_EEET;   //!
   TBranch        *b_EESumE;   //!
   TBranch        *b_EESumET;   //!
   TBranch        *b_HBET;   //!
   TBranch        *b_HBSumE;   //!
   TBranch        *b_HBSumET;   //!
   TBranch        *b_HEET;   //!
   TBranch        *b_HESumE;   //!
   TBranch        *b_HESumET;   //!
   TBranch        *b_HFET;   //!
   TBranch        *b_JetEMEB;   //!
   TBranch        *b_JetEMEE;   //!
   TBranch        *b_JetEMFrac;   //!
   TBranch        *b_JetEMHF;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetHadFrac;   //!
   TBranch        *b_JetHadHB;   //!
   TBranch        *b_JetHadHE;   //!
   TBranch        *b_JetHadHF;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_NominalMET;   //!
   TBranch        *b_CDihitCluE;   //!
   TBranch        *b_CDihitIsolEcalE;   //!
   TBranch        *b_CDihitIsolHcalE;   //!
   TBranch        *b_CDihitIsolTrkE;   //!
   TBranch        *b_CDihitTrkFidE;   //!
   TBranch        *b_CHPDCluE;   //!
   TBranch        *b_CHPDIsolEcalE;   //!
   TBranch        *b_CHPDIsolHcalE;   //!
   TBranch        *b_CHPDIsolTrkE;   //!
   TBranch        *b_CHPDTrkFidE;   //!
   TBranch        *b_CMonohitCluE;   //!
   TBranch        *b_CMonohitIsolEcalE;   //!
   TBranch        *b_CMonohitIsolHcalE;   //!
   TBranch        *b_CMonohitIsolTrkE;   //!
   TBranch        *b_CMonohitTrkFidE;   //!
   TBranch        *b_CRBXCluE;   //!
   TBranch        *b_CRBXIsolEcalE;   //!
   TBranch        *b_CRBXIsolHcalE;   //!
   TBranch        *b_CRBXIsolTrkE;   //!
   TBranch        *b_CRBXTrkFidE;   //!
   TBranch        *b_HBHERecHitEnergyRaw;   //!
   TBranch        *b_IsolatedNoiseSumE;   //!
   TBranch        *b_IsolatedNoiseSumEt;   //!
   TBranch        *b_MaxE2E10;   //!
   TBranch        *b_MinE2E10;   //!
   TBranch        *b_NegativeNoiseSumE;   //!
   TBranch        *b_NegativeNoiseSumEt;   //!
   TBranch        *b_RBXEnergy;   //!
   TBranch        *b_RBXEnergy15;   //!
   TBranch        *b_SpikeNoiseSumE;   //!
   TBranch        *b_SpikeNoiseSumEt;   //!
   TBranch        *b_HBHERecHitAuxAllfC;   //!
   TBranch        *b_HBHERecHitAuxEnergy;   //!
   TBranch        *b_HBHERecHitAuxFC;   //!
   TBranch        *b_HBHERecHitAuxGain;   //!
   TBranch        *b_HBHERecHitAuxPedFC;   //!
   TBranch        *b_HBHERecHitAuxRCGain;   //!
   TBranch        *b_RBXCharge;   //!
   TBranch        *b_RBXCharge15;   //!
   TBranch        *b_HBHERecHitEnergy;   //!
   TBranch        *b_HBHERecHitEta;   //!
   TBranch        *b_HBHERecHitPhi;   //!
   TBranch        *b_HBHERecHitTime;   //!
   TBranch        *b_JetN60;   //!
   TBranch        *b_JetN90;   //!
   TBranch        *b_HBHERecHitAux;   //!
   TBranch        *b_HBHERecHitDepth;   //!
   TBranch        *b_HBHERecHitFlags;   //!
   TBranch        *b_HBHERecHitHPDid;   //!
   TBranch        *b_HBHERecHitIEta;   //!
   TBranch        *b_HBHERecHitIPhi;   //!
   TBranch        *b_HBHERecHitRBXid;   //!
   TBranch        *b_CDihitIsIso;   //!
   TBranch        *b_CDihitIsTagged;   //!
   TBranch        *b_CHPDId;   //!
   TBranch        *b_CHPDIsIso;   //!
   TBranch        *b_CHPDIsTagged;   //!
   TBranch        *b_CHPDNHits;   //!
   TBranch        *b_CHPDRBXId;   //!
   TBranch        *b_CMonohitIsIso;   //!
   TBranch        *b_CMonohitIsTagged;   //!
   TBranch        *b_CRBXId;   //!
   TBranch        *b_CRBXIsIso;   //!
   TBranch        *b_CRBXIsTagged;   //!
   TBranch        *b_CRBXNHits;   //!
   TBranch        *b_HPDHits;   //!
   TBranch        *b_HPDNoOtherHits;   //!
   TBranch        *b_HasBadRBXR45;   //!
   TBranch        *b_HasBadRBXRechitR45Loose;   //!
   TBranch        *b_HasBadRBXRechitR45Tight;   //!
   TBranch        *b_IsoNoiseFilterDecision;   //!
   TBranch        *b_MaxZeros;   //!
   TBranch        *b_NumIsolatedNoiseChannels;   //!
   TBranch        *b_NumNegativeNoiseChannels;   //!
   TBranch        *b_NumSpikeNoiseChannels;   //!
   TBranch        *b_OfficialDecision;   //!
   TBranch        *b_OfficialDecisionRun1;   //!
   TBranch        *b_OfficialDecisionRun2L;   //!
   TBranch        *b_OfficialDecisionRun2T;   //!
   TBranch        *b_HBHERecHitAuxADC;   //!
   TBranch        *b_HBHERecHitAuxCapID;   //!
   TBranch        *b_event;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_run;   //!

   treeBase(TTree *tree=0);
   virtual ~treeBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef treeBase_cxx
treeBase::treeBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Test.root:/hcalTupleTree");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

treeBase::~treeBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeBase::LoadTree(Long64_t entry)
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

void treeBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EBET = 0;
   EBSumE = 0;
   EBSumET = 0;
   EEET = 0;
   EESumE = 0;
   EESumET = 0;
   HBET = 0;
   HBSumE = 0;
   HBSumET = 0;
   HEET = 0;
   HESumE = 0;
   HESumET = 0;
   HFET = 0;
   JetEMEB = 0;
   JetEMEE = 0;
   JetEMFrac = 0;
   JetEMHF = 0;
   JetEta = 0;
   JetHadFrac = 0;
   JetHadHB = 0;
   JetHadHE = 0;
   JetHadHF = 0;
   JetPhi = 0;
   JetPt = 0;
   NominalMET = 0;
   CDihitCluE = 0;
   CDihitIsolEcalE = 0;
   CDihitIsolHcalE = 0;
   CDihitIsolTrkE = 0;
   CDihitTrkFidE = 0;
   CHPDCluE = 0;
   CHPDIsolEcalE = 0;
   CHPDIsolHcalE = 0;
   CHPDIsolTrkE = 0;
   CHPDTrkFidE = 0;
   CMonohitCluE = 0;
   CMonohitIsolEcalE = 0;
   CMonohitIsolHcalE = 0;
   CMonohitIsolTrkE = 0;
   CMonohitTrkFidE = 0;
   CRBXCluE = 0;
   CRBXIsolEcalE = 0;
   CRBXIsolHcalE = 0;
   CRBXIsolTrkE = 0;
   CRBXTrkFidE = 0;
   HBHERecHitEnergyRaw = 0;
   IsolatedNoiseSumE = 0;
   IsolatedNoiseSumEt = 0;
   MaxE2E10 = 0;
   MinE2E10 = 0;
   NegativeNoiseSumE = 0;
   NegativeNoiseSumEt = 0;
   RBXEnergy = 0;
   RBXEnergy15 = 0;
   SpikeNoiseSumE = 0;
   SpikeNoiseSumEt = 0;
   HBHERecHitAuxAllfC = 0;
   HBHERecHitAuxEnergy = 0;
   HBHERecHitAuxFC = 0;
   HBHERecHitAuxGain = 0;
   HBHERecHitAuxPedFC = 0;
   HBHERecHitAuxRCGain = 0;
   RBXCharge = 0;
   RBXCharge15 = 0;
   HBHERecHitEnergy = 0;
   HBHERecHitEta = 0;
   HBHERecHitPhi = 0;
   HBHERecHitTime = 0;
   JetN60 = 0;
   JetN90 = 0;
   HBHERecHitAux = 0;
   HBHERecHitDepth = 0;
   HBHERecHitFlags = 0;
   HBHERecHitHPDid = 0;
   HBHERecHitIEta = 0;
   HBHERecHitIPhi = 0;
   HBHERecHitRBXid = 0;
   CDihitIsIso = 0;
   CDihitIsTagged = 0;
   CHPDId = 0;
   CHPDIsIso = 0;
   CHPDIsTagged = 0;
   CHPDNHits = 0;
   CHPDRBXId = 0;
   CMonohitIsIso = 0;
   CMonohitIsTagged = 0;
   CRBXId = 0;
   CRBXIsIso = 0;
   CRBXIsTagged = 0;
   CRBXNHits = 0;
   HPDHits = 0;
   HPDNoOtherHits = 0;
   HasBadRBXR45 = 0;
   HasBadRBXRechitR45Loose = 0;
   HasBadRBXRechitR45Tight = 0;
   IsoNoiseFilterDecision = 0;
   MaxZeros = 0;
   NumIsolatedNoiseChannels = 0;
   NumNegativeNoiseChannels = 0;
   NumSpikeNoiseChannels = 0;
   OfficialDecision = 0;
   OfficialDecisionRun1 = 0;
   OfficialDecisionRun2L = 0;
   OfficialDecisionRun2T = 0;
   HBHERecHitAuxADC = 0;
   HBHERecHitAuxCapID = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EBET", &EBET, &b_EBET);
   fChain->SetBranchAddress("EBSumE", &EBSumE, &b_EBSumE);
   fChain->SetBranchAddress("EBSumET", &EBSumET, &b_EBSumET);
   fChain->SetBranchAddress("EEET", &EEET, &b_EEET);
   fChain->SetBranchAddress("EESumE", &EESumE, &b_EESumE);
   fChain->SetBranchAddress("EESumET", &EESumET, &b_EESumET);
   fChain->SetBranchAddress("HBET", &HBET, &b_HBET);
   fChain->SetBranchAddress("HBSumE", &HBSumE, &b_HBSumE);
   fChain->SetBranchAddress("HBSumET", &HBSumET, &b_HBSumET);
   fChain->SetBranchAddress("HEET", &HEET, &b_HEET);
   fChain->SetBranchAddress("HESumE", &HESumE, &b_HESumE);
   fChain->SetBranchAddress("HESumET", &HESumET, &b_HESumET);
   fChain->SetBranchAddress("HFET", &HFET, &b_HFET);
   fChain->SetBranchAddress("JetEMEB", &JetEMEB, &b_JetEMEB);
   fChain->SetBranchAddress("JetEMEE", &JetEMEE, &b_JetEMEE);
   fChain->SetBranchAddress("JetEMFrac", &JetEMFrac, &b_JetEMFrac);
   fChain->SetBranchAddress("JetEMHF", &JetEMHF, &b_JetEMHF);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetHadFrac", &JetHadFrac, &b_JetHadFrac);
   fChain->SetBranchAddress("JetHadHB", &JetHadHB, &b_JetHadHB);
   fChain->SetBranchAddress("JetHadHE", &JetHadHE, &b_JetHadHE);
   fChain->SetBranchAddress("JetHadHF", &JetHadHF, &b_JetHadHF);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("NominalMET", &NominalMET, &b_NominalMET);
   fChain->SetBranchAddress("CDihitCluE", &CDihitCluE, &b_CDihitCluE);
   fChain->SetBranchAddress("CDihitIsolEcalE", &CDihitIsolEcalE, &b_CDihitIsolEcalE);
   fChain->SetBranchAddress("CDihitIsolHcalE", &CDihitIsolHcalE, &b_CDihitIsolHcalE);
   fChain->SetBranchAddress("CDihitIsolTrkE", &CDihitIsolTrkE, &b_CDihitIsolTrkE);
   fChain->SetBranchAddress("CDihitTrkFidE", &CDihitTrkFidE, &b_CDihitTrkFidE);
   fChain->SetBranchAddress("CHPDCluE", &CHPDCluE, &b_CHPDCluE);
   fChain->SetBranchAddress("CHPDIsolEcalE", &CHPDIsolEcalE, &b_CHPDIsolEcalE);
   fChain->SetBranchAddress("CHPDIsolHcalE", &CHPDIsolHcalE, &b_CHPDIsolHcalE);
   fChain->SetBranchAddress("CHPDIsolTrkE", &CHPDIsolTrkE, &b_CHPDIsolTrkE);
   fChain->SetBranchAddress("CHPDTrkFidE", &CHPDTrkFidE, &b_CHPDTrkFidE);
   fChain->SetBranchAddress("CMonohitCluE", &CMonohitCluE, &b_CMonohitCluE);
   fChain->SetBranchAddress("CMonohitIsolEcalE", &CMonohitIsolEcalE, &b_CMonohitIsolEcalE);
   fChain->SetBranchAddress("CMonohitIsolHcalE", &CMonohitIsolHcalE, &b_CMonohitIsolHcalE);
   fChain->SetBranchAddress("CMonohitIsolTrkE", &CMonohitIsolTrkE, &b_CMonohitIsolTrkE);
   fChain->SetBranchAddress("CMonohitTrkFidE", &CMonohitTrkFidE, &b_CMonohitTrkFidE);
   fChain->SetBranchAddress("CRBXCluE", &CRBXCluE, &b_CRBXCluE);
   fChain->SetBranchAddress("CRBXIsolEcalE", &CRBXIsolEcalE, &b_CRBXIsolEcalE);
   fChain->SetBranchAddress("CRBXIsolHcalE", &CRBXIsolHcalE, &b_CRBXIsolHcalE);
   fChain->SetBranchAddress("CRBXIsolTrkE", &CRBXIsolTrkE, &b_CRBXIsolTrkE);
   fChain->SetBranchAddress("CRBXTrkFidE", &CRBXTrkFidE, &b_CRBXTrkFidE);
   fChain->SetBranchAddress("HBHERecHitEnergyRaw", &HBHERecHitEnergyRaw, &b_HBHERecHitEnergyRaw);
   fChain->SetBranchAddress("IsolatedNoiseSumE", &IsolatedNoiseSumE, &b_IsolatedNoiseSumE);
   fChain->SetBranchAddress("IsolatedNoiseSumEt", &IsolatedNoiseSumEt, &b_IsolatedNoiseSumEt);
   fChain->SetBranchAddress("MaxE2E10", &MaxE2E10, &b_MaxE2E10);
   fChain->SetBranchAddress("MinE2E10", &MinE2E10, &b_MinE2E10);
   fChain->SetBranchAddress("NegativeNoiseSumE", &NegativeNoiseSumE, &b_NegativeNoiseSumE);
   fChain->SetBranchAddress("NegativeNoiseSumEt", &NegativeNoiseSumEt, &b_NegativeNoiseSumEt);
   fChain->SetBranchAddress("RBXEnergy", &RBXEnergy, &b_RBXEnergy);
   fChain->SetBranchAddress("RBXEnergy15", &RBXEnergy15, &b_RBXEnergy15);
   fChain->SetBranchAddress("SpikeNoiseSumE", &SpikeNoiseSumE, &b_SpikeNoiseSumE);
   fChain->SetBranchAddress("SpikeNoiseSumEt", &SpikeNoiseSumEt, &b_SpikeNoiseSumEt);
   fChain->SetBranchAddress("HBHERecHitAuxAllfC", &HBHERecHitAuxAllfC, &b_HBHERecHitAuxAllfC);
   fChain->SetBranchAddress("HBHERecHitAuxEnergy", &HBHERecHitAuxEnergy, &b_HBHERecHitAuxEnergy);
   fChain->SetBranchAddress("HBHERecHitAuxFC", &HBHERecHitAuxFC, &b_HBHERecHitAuxFC);
   fChain->SetBranchAddress("HBHERecHitAuxGain", &HBHERecHitAuxGain, &b_HBHERecHitAuxGain);
   fChain->SetBranchAddress("HBHERecHitAuxPedFC", &HBHERecHitAuxPedFC, &b_HBHERecHitAuxPedFC);
   fChain->SetBranchAddress("HBHERecHitAuxRCGain", &HBHERecHitAuxRCGain, &b_HBHERecHitAuxRCGain);
   fChain->SetBranchAddress("RBXCharge", &RBXCharge, &b_RBXCharge);
   fChain->SetBranchAddress("RBXCharge15", &RBXCharge15, &b_RBXCharge15);
   fChain->SetBranchAddress("HBHERecHitEnergy", &HBHERecHitEnergy, &b_HBHERecHitEnergy);
   fChain->SetBranchAddress("HBHERecHitEta", &HBHERecHitEta, &b_HBHERecHitEta);
   fChain->SetBranchAddress("HBHERecHitPhi", &HBHERecHitPhi, &b_HBHERecHitPhi);
   fChain->SetBranchAddress("HBHERecHitTime", &HBHERecHitTime, &b_HBHERecHitTime);
   fChain->SetBranchAddress("JetN60", &JetN60, &b_JetN60);
   fChain->SetBranchAddress("JetN90", &JetN90, &b_JetN90);
   fChain->SetBranchAddress("HBHERecHitAux", &HBHERecHitAux, &b_HBHERecHitAux);
   fChain->SetBranchAddress("HBHERecHitDepth", &HBHERecHitDepth, &b_HBHERecHitDepth);
   fChain->SetBranchAddress("HBHERecHitFlags", &HBHERecHitFlags, &b_HBHERecHitFlags);
   fChain->SetBranchAddress("HBHERecHitHPDid", &HBHERecHitHPDid, &b_HBHERecHitHPDid);
   fChain->SetBranchAddress("HBHERecHitIEta", &HBHERecHitIEta, &b_HBHERecHitIEta);
   fChain->SetBranchAddress("HBHERecHitIPhi", &HBHERecHitIPhi, &b_HBHERecHitIPhi);
   fChain->SetBranchAddress("HBHERecHitRBXid", &HBHERecHitRBXid, &b_HBHERecHitRBXid);
   fChain->SetBranchAddress("CDihitIsIso", &CDihitIsIso, &b_CDihitIsIso);
   fChain->SetBranchAddress("CDihitIsTagged", &CDihitIsTagged, &b_CDihitIsTagged);
   fChain->SetBranchAddress("CHPDId", &CHPDId, &b_CHPDId);
   fChain->SetBranchAddress("CHPDIsIso", &CHPDIsIso, &b_CHPDIsIso);
   fChain->SetBranchAddress("CHPDIsTagged", &CHPDIsTagged, &b_CHPDIsTagged);
   fChain->SetBranchAddress("CHPDNHits", &CHPDNHits, &b_CHPDNHits);
   fChain->SetBranchAddress("CHPDRBXId", &CHPDRBXId, &b_CHPDRBXId);
   fChain->SetBranchAddress("CMonohitIsIso", &CMonohitIsIso, &b_CMonohitIsIso);
   fChain->SetBranchAddress("CMonohitIsTagged", &CMonohitIsTagged, &b_CMonohitIsTagged);
   fChain->SetBranchAddress("CRBXId", &CRBXId, &b_CRBXId);
   fChain->SetBranchAddress("CRBXIsIso", &CRBXIsIso, &b_CRBXIsIso);
   fChain->SetBranchAddress("CRBXIsTagged", &CRBXIsTagged, &b_CRBXIsTagged);
   fChain->SetBranchAddress("CRBXNHits", &CRBXNHits, &b_CRBXNHits);
   fChain->SetBranchAddress("HPDHits", &HPDHits, &b_HPDHits);
   fChain->SetBranchAddress("HPDNoOtherHits", &HPDNoOtherHits, &b_HPDNoOtherHits);
   fChain->SetBranchAddress("HasBadRBXR45", &HasBadRBXR45, &b_HasBadRBXR45);
   fChain->SetBranchAddress("HasBadRBXRechitR45Loose", &HasBadRBXRechitR45Loose, &b_HasBadRBXRechitR45Loose);
   fChain->SetBranchAddress("HasBadRBXRechitR45Tight", &HasBadRBXRechitR45Tight, &b_HasBadRBXRechitR45Tight);
   fChain->SetBranchAddress("IsoNoiseFilterDecision", &IsoNoiseFilterDecision, &b_IsoNoiseFilterDecision);
   fChain->SetBranchAddress("MaxZeros", &MaxZeros, &b_MaxZeros);
   fChain->SetBranchAddress("NumIsolatedNoiseChannels", &NumIsolatedNoiseChannels, &b_NumIsolatedNoiseChannels);
   fChain->SetBranchAddress("NumNegativeNoiseChannels", &NumNegativeNoiseChannels, &b_NumNegativeNoiseChannels);
   fChain->SetBranchAddress("NumSpikeNoiseChannels", &NumSpikeNoiseChannels, &b_NumSpikeNoiseChannels);
   fChain->SetBranchAddress("OfficialDecision", &OfficialDecision, &b_OfficialDecision);
   fChain->SetBranchAddress("OfficialDecisionRun1", &OfficialDecisionRun1, &b_OfficialDecisionRun1);
   fChain->SetBranchAddress("OfficialDecisionRun2L", &OfficialDecisionRun2L, &b_OfficialDecisionRun2L);
   fChain->SetBranchAddress("OfficialDecisionRun2T", &OfficialDecisionRun2T, &b_OfficialDecisionRun2T);
   fChain->SetBranchAddress("HBHERecHitAuxADC", &HBHERecHitAuxADC, &b_HBHERecHitAuxADC);
   fChain->SetBranchAddress("HBHERecHitAuxCapID", &HBHERecHitAuxCapID, &b_HBHERecHitAuxCapID);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("run", &run, &b_run);
   Notify();
}

Bool_t treeBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeBase_cxx
