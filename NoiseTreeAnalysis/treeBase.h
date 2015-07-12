//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May 30 13:06:39 2015 by ROOT version 5.34/18
// from TTree tree/
// found on file: /dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_2402393Cto3685A126.root
//////////////////////////////////////////////////////////

#ifndef treeBase_h
#define treeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class treeBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   vector<double>  *EBETMethod0;
   vector<double>  *EBSumEMethod0;
   vector<double>  *EBSumETMethod0;
   vector<double>  *EEETMethod0;
   vector<double>  *EESumEMethod0;
   vector<double>  *EESumETMethod0;
   vector<double>  *HBETMethod0;
   vector<double>  *HBSumEMethod0;
   vector<double>  *HBSumETMethod0;
   vector<double>  *HEETMethod0;
   vector<double>  *HESumEMethod0;
   vector<double>  *HESumETMethod0;
   vector<double>  *HFETMethod0;
   vector<double>  *MaxE2E10;
   vector<double>  *MinE2E10;
   vector<double>  *RBXEnergy;
   vector<double>  *RBXEnergy15;
   vector<double>  *MaxE2E10Method0;
   vector<double>  *MinE2E10Method0;
   vector<double>  *RBXEnergy15Method0;
   vector<double>  *RBXEnergyMethod0;
   vector<vector<double> > *RBXCharge;
   vector<vector<double> > *RBXCharge15;
   vector<vector<double> > *RBXCharge15Method0;
   vector<vector<double> > *RBXChargeMethod0;
   vector<float>   *HBHEDigiEta;
   vector<float>   *HBHEDigiPhi;
   vector<float>   *HBHEDigiRecEnergy;
   vector<float>   *HBHEDigiRecTime;
   vector<float>   *HBHERecHitEnergy;
   vector<float>   *HBHERecHitEta;
   vector<float>   *HBHERecHitPhi;
   vector<float>   *HBHERecHitTime;
   vector<float>   *HBHERecHitEnergyMethod0;
   vector<float>   *HBHERecHitEtaMethod0;
   vector<float>   *HBHERecHitPhiMethod0;
   vector<float>   *HBHERecHitTimeMethod0;
   vector<vector<float> > *HBHEDigiAllFC;
   vector<vector<float> > *HBHEDigiEnergy;
   vector<vector<float> > *HBHEDigiFC;
   vector<vector<float> > *HBHEDigiGain;
   vector<vector<float> > *HBHEDigiNomFC;
   vector<vector<float> > *HBHEDigiPedFC;
   vector<vector<float> > *HBHEDigiRCGain;
   vector<int>     *HBHEDigiDepth;
   vector<int>     *HBHEDigiElectronicsID;
   vector<int>     *HBHEDigiFiberIdleOffset;
   vector<int>     *HBHEDigiIEta;
   vector<int>     *HBHEDigiIPhi;
   vector<int>     *HBHEDigiPresamples;
   vector<int>     *HBHEDigiRawID;
   vector<int>     *HBHEDigiSize;
   vector<int>     *HBHEDigiSubdet;
   vector<int>     *HBHERecHitAux;
   vector<int>     *HBHERecHitDepth;
   vector<int>     *HBHERecHitFlags;
   vector<int>     *HBHERecHitIEta;
   vector<int>     *HBHERecHitIPhi;
   vector<int>     *HBHERecHitAuxMethod0;
   vector<int>     *HBHERecHitDepthMethod0;
   vector<int>     *HBHERecHitFlagsMethod0;
   vector<int>     *HBHERecHitIEtaMethod0;
   vector<int>     *HBHERecHitIPhiMethod0;
   vector<int>     *HPDHits;
   vector<int>     *HPDNoOtherHits;
   vector<int>     *HasBadRBXR45;
   vector<int>     *HasBadRBXRechitR45Loose;
   vector<int>     *HasBadRBXRechitR45Tight;
   vector<int>     *MaxZeros;
   vector<int>     *OfficialDecision;
   vector<int>     *HPDHitsMethod0;
   vector<int>     *HPDNoOtherHitsMethod0;
   vector<int>     *HasBadRBXR45Method0;
   vector<int>     *HasBadRBXRechitR45LooseMethod0;
   vector<int>     *HasBadRBXRechitR45TightMethod0;
   vector<int>     *MaxZerosMethod0;
   vector<int>     *OfficialDecisionMethod0;
   vector<vector<int> > *HBHEDigiADC;
   vector<vector<int> > *HBHEDigiCapID;
   vector<vector<int> > *HBHEDigiDV;
   vector<vector<int> > *HBHEDigiER;
   vector<vector<int> > *HBHEDigiFiber;
   vector<vector<int> > *HBHEDigiFiberChan;
   vector<vector<int> > *HBHEDigiLADC;
   vector<vector<int> > *HBHEDigiRaw;
   UInt_t          event;
   UInt_t          ls;
   UInt_t          run;
   vector<unsigned int> *AuxWord;
   vector<unsigned int> *FlagWord;
   vector<unsigned int> *AuxWordMethod0;
   vector<unsigned int> *FlagWordMethod0;

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
   TBranch        *b_EBETMethod0;   //!
   TBranch        *b_EBSumEMethod0;   //!
   TBranch        *b_EBSumETMethod0;   //!
   TBranch        *b_EEETMethod0;   //!
   TBranch        *b_EESumEMethod0;   //!
   TBranch        *b_EESumETMethod0;   //!
   TBranch        *b_HBETMethod0;   //!
   TBranch        *b_HBSumEMethod0;   //!
   TBranch        *b_HBSumETMethod0;   //!
   TBranch        *b_HEETMethod0;   //!
   TBranch        *b_HESumEMethod0;   //!
   TBranch        *b_HESumETMethod0;   //!
   TBranch        *b_HFETMethod0;   //!
   TBranch        *b_MaxE2E10;   //!
   TBranch        *b_MinE2E10;   //!
   TBranch        *b_RBXEnergy;   //!
   TBranch        *b_RBXEnergy15;   //!
   TBranch        *b_MaxE2E10Method0;   //!
   TBranch        *b_MinE2E10Method0;   //!
   TBranch        *b_RBXEnergy15Method0;   //!
   TBranch        *b_RBXEnergyMethod0;   //!
   TBranch        *b_RBXCharge;   //!
   TBranch        *b_RBXCharge15;   //!
   TBranch        *b_RBXCharge15Method0;   //!
   TBranch        *b_RBXChargeMethod0;   //!
   TBranch        *b_HBHEDigiEta;   //!
   TBranch        *b_HBHEDigiPhi;   //!
   TBranch        *b_HBHEDigiRecEnergy;   //!
   TBranch        *b_HBHEDigiRecTime;   //!
   TBranch        *b_HBHERecHitEnergy;   //!
   TBranch        *b_HBHERecHitEta;   //!
   TBranch        *b_HBHERecHitPhi;   //!
   TBranch        *b_HBHERecHitTime;   //!
   TBranch        *b_HBHERecHitEnergyMethod0;   //!
   TBranch        *b_HBHERecHitEtaMethod0;   //!
   TBranch        *b_HBHERecHitPhiMethod0;   //!
   TBranch        *b_HBHERecHitTimeMethod0;   //!
   TBranch        *b_HBHEDigiAllFC;   //!
   TBranch        *b_HBHEDigiEnergy;   //!
   TBranch        *b_HBHEDigiFC;   //!
   TBranch        *b_HBHEDigiGain;   //!
   TBranch        *b_HBHEDigiNomFC;   //!
   TBranch        *b_HBHEDigiPedFC;   //!
   TBranch        *b_HBHEDigiRCGain;   //!
   TBranch        *b_HBHEDigiDepth;   //!
   TBranch        *b_HBHEDigiElectronicsID;   //!
   TBranch        *b_HBHEDigiFiberIdleOffset;   //!
   TBranch        *b_HBHEDigiIEta;   //!
   TBranch        *b_HBHEDigiIPhi;   //!
   TBranch        *b_HBHEDigiPresamples;   //!
   TBranch        *b_HBHEDigiRawID;   //!
   TBranch        *b_HBHEDigiSize;   //!
   TBranch        *b_HBHEDigiSubdet;   //!
   TBranch        *b_HBHERecHitAux;   //!
   TBranch        *b_HBHERecHitDepth;   //!
   TBranch        *b_HBHERecHitFlags;   //!
   TBranch        *b_HBHERecHitIEta;   //!
   TBranch        *b_HBHERecHitIPhi;   //!
   TBranch        *b_HBHERecHitAuxMethod0;   //!
   TBranch        *b_HBHERecHitDepthMethod0;   //!
   TBranch        *b_HBHERecHitFlagsMethod0;   //!
   TBranch        *b_HBHERecHitIEtaMethod0;   //!
   TBranch        *b_HBHERecHitIPhiMethod0;   //!
   TBranch        *b_HPDHits;   //!
   TBranch        *b_HPDNoOtherHits;   //!
   TBranch        *b_HasBadRBXR45;   //!
   TBranch        *b_HasBadRBXRechitR45Loose;   //!
   TBranch        *b_HasBadRBXRechitR45Tight;   //!
   TBranch        *b_MaxZeros;   //!
   TBranch        *b_OfficialDecision;   //!
   TBranch        *b_HPDHitsMethod0;   //!
   TBranch        *b_HPDNoOtherHitsMethod0;   //!
   TBranch        *b_HasBadRBXR45Method0;   //!
   TBranch        *b_HasBadRBXRechitR45LooseMethod0;   //!
   TBranch        *b_HasBadRBXRechitR45TightMethod0;   //!
   TBranch        *b_MaxZerosMethod0;   //!
   TBranch        *b_OfficialDecisionMethod0;   //!
   TBranch        *b_HBHEDigiADC;   //!
   TBranch        *b_HBHEDigiCapID;   //!
   TBranch        *b_HBHEDigiDV;   //!
   TBranch        *b_HBHEDigiER;   //!
   TBranch        *b_HBHEDigiFiber;   //!
   TBranch        *b_HBHEDigiFiberChan;   //!
   TBranch        *b_HBHEDigiLADC;   //!
   TBranch        *b_HBHEDigiRaw;   //!
   TBranch        *b_event;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_run;   //!
   TBranch        *b_AuxWord;   //!
   TBranch        *b_FlagWord;   //!
   TBranch        *b_AuxWordMethod0;   //!
   TBranch        *b_FlagWordMethod0;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_2402393Cto3685A126.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_2402393Cto3685A126.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/dataLOCAL/HCAL/RootFiles/Schenara/data_Commissioning2014_HcalHPDNoise_2402393Cto3685A126.root:/hcalTupleTree");
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
   EBETMethod0 = 0;
   EBSumEMethod0 = 0;
   EBSumETMethod0 = 0;
   EEETMethod0 = 0;
   EESumEMethod0 = 0;
   EESumETMethod0 = 0;
   HBETMethod0 = 0;
   HBSumEMethod0 = 0;
   HBSumETMethod0 = 0;
   HEETMethod0 = 0;
   HESumEMethod0 = 0;
   HESumETMethod0 = 0;
   HFETMethod0 = 0;
   MaxE2E10 = 0;
   MinE2E10 = 0;
   RBXEnergy = 0;
   RBXEnergy15 = 0;
   MaxE2E10Method0 = 0;
   MinE2E10Method0 = 0;
   RBXEnergy15Method0 = 0;
   RBXEnergyMethod0 = 0;
   RBXCharge = 0;
   RBXCharge15 = 0;
   RBXCharge15Method0 = 0;
   RBXChargeMethod0 = 0;
   HBHEDigiEta = 0;
   HBHEDigiPhi = 0;
   HBHEDigiRecEnergy = 0;
   HBHEDigiRecTime = 0;
   HBHERecHitEnergy = 0;
   HBHERecHitEta = 0;
   HBHERecHitPhi = 0;
   HBHERecHitTime = 0;
   HBHERecHitEnergyMethod0 = 0;
   HBHERecHitEtaMethod0 = 0;
   HBHERecHitPhiMethod0 = 0;
   HBHERecHitTimeMethod0 = 0;
   HBHEDigiAllFC = 0;
   HBHEDigiEnergy = 0;
   HBHEDigiFC = 0;
   HBHEDigiGain = 0;
   HBHEDigiNomFC = 0;
   HBHEDigiPedFC = 0;
   HBHEDigiRCGain = 0;
   HBHEDigiDepth = 0;
   HBHEDigiElectronicsID = 0;
   HBHEDigiFiberIdleOffset = 0;
   HBHEDigiIEta = 0;
   HBHEDigiIPhi = 0;
   HBHEDigiPresamples = 0;
   HBHEDigiRawID = 0;
   HBHEDigiSize = 0;
   HBHEDigiSubdet = 0;
   HBHERecHitAux = 0;
   HBHERecHitDepth = 0;
   HBHERecHitFlags = 0;
   HBHERecHitIEta = 0;
   HBHERecHitIPhi = 0;
   HBHERecHitAuxMethod0 = 0;
   HBHERecHitDepthMethod0 = 0;
   HBHERecHitFlagsMethod0 = 0;
   HBHERecHitIEtaMethod0 = 0;
   HBHERecHitIPhiMethod0 = 0;
   HPDHits = 0;
   HPDNoOtherHits = 0;
   HasBadRBXR45 = 0;
   HasBadRBXRechitR45Loose = 0;
   HasBadRBXRechitR45Tight = 0;
   MaxZeros = 0;
   OfficialDecision = 0;
   HPDHitsMethod0 = 0;
   HPDNoOtherHitsMethod0 = 0;
   HasBadRBXR45Method0 = 0;
   HasBadRBXRechitR45LooseMethod0 = 0;
   HasBadRBXRechitR45TightMethod0 = 0;
   MaxZerosMethod0 = 0;
   OfficialDecisionMethod0 = 0;
   HBHEDigiADC = 0;
   HBHEDigiCapID = 0;
   HBHEDigiDV = 0;
   HBHEDigiER = 0;
   HBHEDigiFiber = 0;
   HBHEDigiFiberChan = 0;
   HBHEDigiLADC = 0;
   HBHEDigiRaw = 0;
   AuxWord = 0;
   FlagWord = 0;
   AuxWordMethod0 = 0;
   FlagWordMethod0 = 0;
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
   fChain->SetBranchAddress("EBETMethod0", &EBETMethod0, &b_EBETMethod0);
   fChain->SetBranchAddress("EBSumEMethod0", &EBSumEMethod0, &b_EBSumEMethod0);
   fChain->SetBranchAddress("EBSumETMethod0", &EBSumETMethod0, &b_EBSumETMethod0);
   fChain->SetBranchAddress("EEETMethod0", &EEETMethod0, &b_EEETMethod0);
   fChain->SetBranchAddress("EESumEMethod0", &EESumEMethod0, &b_EESumEMethod0);
   fChain->SetBranchAddress("EESumETMethod0", &EESumETMethod0, &b_EESumETMethod0);
   fChain->SetBranchAddress("HBETMethod0", &HBETMethod0, &b_HBETMethod0);
   fChain->SetBranchAddress("HBSumEMethod0", &HBSumEMethod0, &b_HBSumEMethod0);
   fChain->SetBranchAddress("HBSumETMethod0", &HBSumETMethod0, &b_HBSumETMethod0);
   fChain->SetBranchAddress("HEETMethod0", &HEETMethod0, &b_HEETMethod0);
   fChain->SetBranchAddress("HESumEMethod0", &HESumEMethod0, &b_HESumEMethod0);
   fChain->SetBranchAddress("HESumETMethod0", &HESumETMethod0, &b_HESumETMethod0);
   fChain->SetBranchAddress("HFETMethod0", &HFETMethod0, &b_HFETMethod0);
   fChain->SetBranchAddress("MaxE2E10", &MaxE2E10, &b_MaxE2E10);
   fChain->SetBranchAddress("MinE2E10", &MinE2E10, &b_MinE2E10);
   fChain->SetBranchAddress("RBXEnergy", &RBXEnergy, &b_RBXEnergy);
   fChain->SetBranchAddress("RBXEnergy15", &RBXEnergy15, &b_RBXEnergy15);
   fChain->SetBranchAddress("MaxE2E10Method0", &MaxE2E10Method0, &b_MaxE2E10Method0);
   fChain->SetBranchAddress("MinE2E10Method0", &MinE2E10Method0, &b_MinE2E10Method0);
   fChain->SetBranchAddress("RBXEnergy15Method0", &RBXEnergy15Method0, &b_RBXEnergy15Method0);
   fChain->SetBranchAddress("RBXEnergyMethod0", &RBXEnergyMethod0, &b_RBXEnergyMethod0);
   fChain->SetBranchAddress("RBXCharge", &RBXCharge, &b_RBXCharge);
   fChain->SetBranchAddress("RBXCharge15", &RBXCharge15, &b_RBXCharge15);
   fChain->SetBranchAddress("RBXCharge15Method0", &RBXCharge15Method0, &b_RBXCharge15Method0);
   fChain->SetBranchAddress("RBXChargeMethod0", &RBXChargeMethod0, &b_RBXChargeMethod0);
   fChain->SetBranchAddress("HBHEDigiEta", &HBHEDigiEta, &b_HBHEDigiEta);
   fChain->SetBranchAddress("HBHEDigiPhi", &HBHEDigiPhi, &b_HBHEDigiPhi);
   fChain->SetBranchAddress("HBHEDigiRecEnergy", &HBHEDigiRecEnergy, &b_HBHEDigiRecEnergy);
   fChain->SetBranchAddress("HBHEDigiRecTime", &HBHEDigiRecTime, &b_HBHEDigiRecTime);
   fChain->SetBranchAddress("HBHERecHitEnergy", &HBHERecHitEnergy, &b_HBHERecHitEnergy);
   fChain->SetBranchAddress("HBHERecHitEta", &HBHERecHitEta, &b_HBHERecHitEta);
   fChain->SetBranchAddress("HBHERecHitPhi", &HBHERecHitPhi, &b_HBHERecHitPhi);
   fChain->SetBranchAddress("HBHERecHitTime", &HBHERecHitTime, &b_HBHERecHitTime);
   fChain->SetBranchAddress("HBHERecHitEnergyMethod0", &HBHERecHitEnergyMethod0, &b_HBHERecHitEnergyMethod0);
   fChain->SetBranchAddress("HBHERecHitEtaMethod0", &HBHERecHitEtaMethod0, &b_HBHERecHitEtaMethod0);
   fChain->SetBranchAddress("HBHERecHitPhiMethod0", &HBHERecHitPhiMethod0, &b_HBHERecHitPhiMethod0);
   fChain->SetBranchAddress("HBHERecHitTimeMethod0", &HBHERecHitTimeMethod0, &b_HBHERecHitTimeMethod0);
   fChain->SetBranchAddress("HBHEDigiAllFC", &HBHEDigiAllFC, &b_HBHEDigiAllFC);
   fChain->SetBranchAddress("HBHEDigiEnergy", &HBHEDigiEnergy, &b_HBHEDigiEnergy);
   fChain->SetBranchAddress("HBHEDigiFC", &HBHEDigiFC, &b_HBHEDigiFC);
   fChain->SetBranchAddress("HBHEDigiGain", &HBHEDigiGain, &b_HBHEDigiGain);
   fChain->SetBranchAddress("HBHEDigiNomFC", &HBHEDigiNomFC, &b_HBHEDigiNomFC);
   fChain->SetBranchAddress("HBHEDigiPedFC", &HBHEDigiPedFC, &b_HBHEDigiPedFC);
   fChain->SetBranchAddress("HBHEDigiRCGain", &HBHEDigiRCGain, &b_HBHEDigiRCGain);
   fChain->SetBranchAddress("HBHEDigiDepth", &HBHEDigiDepth, &b_HBHEDigiDepth);
   fChain->SetBranchAddress("HBHEDigiElectronicsID", &HBHEDigiElectronicsID, &b_HBHEDigiElectronicsID);
   fChain->SetBranchAddress("HBHEDigiFiberIdleOffset", &HBHEDigiFiberIdleOffset, &b_HBHEDigiFiberIdleOffset);
   fChain->SetBranchAddress("HBHEDigiIEta", &HBHEDigiIEta, &b_HBHEDigiIEta);
   fChain->SetBranchAddress("HBHEDigiIPhi", &HBHEDigiIPhi, &b_HBHEDigiIPhi);
   fChain->SetBranchAddress("HBHEDigiPresamples", &HBHEDigiPresamples, &b_HBHEDigiPresamples);
   fChain->SetBranchAddress("HBHEDigiRawID", &HBHEDigiRawID, &b_HBHEDigiRawID);
   fChain->SetBranchAddress("HBHEDigiSize", &HBHEDigiSize, &b_HBHEDigiSize);
   fChain->SetBranchAddress("HBHEDigiSubdet", &HBHEDigiSubdet, &b_HBHEDigiSubdet);
   fChain->SetBranchAddress("HBHERecHitAux", &HBHERecHitAux, &b_HBHERecHitAux);
   fChain->SetBranchAddress("HBHERecHitDepth", &HBHERecHitDepth, &b_HBHERecHitDepth);
   fChain->SetBranchAddress("HBHERecHitFlags", &HBHERecHitFlags, &b_HBHERecHitFlags);
   fChain->SetBranchAddress("HBHERecHitIEta", &HBHERecHitIEta, &b_HBHERecHitIEta);
   fChain->SetBranchAddress("HBHERecHitIPhi", &HBHERecHitIPhi, &b_HBHERecHitIPhi);
   fChain->SetBranchAddress("HBHERecHitAuxMethod0", &HBHERecHitAuxMethod0, &b_HBHERecHitAuxMethod0);
   fChain->SetBranchAddress("HBHERecHitDepthMethod0", &HBHERecHitDepthMethod0, &b_HBHERecHitDepthMethod0);
   fChain->SetBranchAddress("HBHERecHitFlagsMethod0", &HBHERecHitFlagsMethod0, &b_HBHERecHitFlagsMethod0);
   fChain->SetBranchAddress("HBHERecHitIEtaMethod0", &HBHERecHitIEtaMethod0, &b_HBHERecHitIEtaMethod0);
   fChain->SetBranchAddress("HBHERecHitIPhiMethod0", &HBHERecHitIPhiMethod0, &b_HBHERecHitIPhiMethod0);
   fChain->SetBranchAddress("HPDHits", &HPDHits, &b_HPDHits);
   fChain->SetBranchAddress("HPDNoOtherHits", &HPDNoOtherHits, &b_HPDNoOtherHits);
   fChain->SetBranchAddress("HasBadRBXR45", &HasBadRBXR45, &b_HasBadRBXR45);
   fChain->SetBranchAddress("HasBadRBXRechitR45Loose", &HasBadRBXRechitR45Loose, &b_HasBadRBXRechitR45Loose);
   fChain->SetBranchAddress("HasBadRBXRechitR45Tight", &HasBadRBXRechitR45Tight, &b_HasBadRBXRechitR45Tight);
   fChain->SetBranchAddress("MaxZeros", &MaxZeros, &b_MaxZeros);
   fChain->SetBranchAddress("OfficialDecision", &OfficialDecision, &b_OfficialDecision);
   fChain->SetBranchAddress("HPDHitsMethod0", &HPDHitsMethod0, &b_HPDHitsMethod0);
   fChain->SetBranchAddress("HPDNoOtherHitsMethod0", &HPDNoOtherHitsMethod0, &b_HPDNoOtherHitsMethod0);
   fChain->SetBranchAddress("HasBadRBXR45Method0", &HasBadRBXR45Method0, &b_HasBadRBXR45Method0);
   fChain->SetBranchAddress("HasBadRBXRechitR45LooseMethod0", &HasBadRBXRechitR45LooseMethod0, &b_HasBadRBXRechitR45LooseMethod0);
   fChain->SetBranchAddress("HasBadRBXRechitR45TightMethod0", &HasBadRBXRechitR45TightMethod0, &b_HasBadRBXRechitR45TightMethod0);
   fChain->SetBranchAddress("MaxZerosMethod0", &MaxZerosMethod0, &b_MaxZerosMethod0);
   fChain->SetBranchAddress("OfficialDecisionMethod0", &OfficialDecisionMethod0, &b_OfficialDecisionMethod0);
   fChain->SetBranchAddress("HBHEDigiADC", &HBHEDigiADC, &b_HBHEDigiADC);
   fChain->SetBranchAddress("HBHEDigiCapID", &HBHEDigiCapID, &b_HBHEDigiCapID);
   fChain->SetBranchAddress("HBHEDigiDV", &HBHEDigiDV, &b_HBHEDigiDV);
   fChain->SetBranchAddress("HBHEDigiER", &HBHEDigiER, &b_HBHEDigiER);
   fChain->SetBranchAddress("HBHEDigiFiber", &HBHEDigiFiber, &b_HBHEDigiFiber);
   fChain->SetBranchAddress("HBHEDigiFiberChan", &HBHEDigiFiberChan, &b_HBHEDigiFiberChan);
   fChain->SetBranchAddress("HBHEDigiLADC", &HBHEDigiLADC, &b_HBHEDigiLADC);
   fChain->SetBranchAddress("HBHEDigiRaw", &HBHEDigiRaw, &b_HBHEDigiRaw);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("AuxWord", &AuxWord, &b_AuxWord);
   fChain->SetBranchAddress("FlagWord", &FlagWord, &b_FlagWord);
   fChain->SetBranchAddress("AuxWordMethod0", &AuxWordMethod0, &b_AuxWordMethod0);
   fChain->SetBranchAddress("FlagWordMethod0", &FlagWordMethod0, &b_FlagWordMethod0);
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
