#define analysisClass_cxx
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
using namespace std;

class analysisClass{

TFile *HBMmapfile;
TFile *HEMmapfile;
TFile *HBPmapfile;
TFile *HEPmapfile;

TTree *HBMmaptree;
TTree *HEMmaptree;
TTree *HBPmaptree;
TTree *HEPmaptree;

int iEta_HBM;
int iPhi_HBM;
int iDepth_HBM;
int RBX_HBM;
int RM_HBM;

int iEta_HEM;
int iPhi_HEM;
int iDepth_HEM;
int RBX_HEM;
int RM_HEM;

int iEta_HBP;
int iPhi_HBP;
int iDepth_HBP;
int RBX_HBP;
int RM_HBP;

int iEta_HEP;
int iPhi_HEP;
int iDepth_HEP;
int RBX_HEP;
int RM_HEP;


void ReadHCALmaps( ){

  //input file
  HBMmapfile=new TFile("HBM_EtaPhiRBXRM.root","READ");
  HEMmapfile=new TFile("HEM_EtaPhiRBXRM.root","READ");
  HBPmapfile=new TFile("HBP_EtaPhiRBXRM.root","READ");
  HEPmapfile=new TFile("HEP_EtaPhiRBXRM.root","READ");
  //
  HBMmaptree=(TTree*)(HBMmapfile->Get("tree"));
  HEMmaptree=(TTree*)(HEMmapfile->Get("tree"));
  HBPmaptree=(TTree*)(HBPmapfile->Get("tree"));
  HEPmaptree=(TTree*)(HEPmapfile->Get("tree"));
  //
  HBMmaptree->SetBranchAddress("iEta",  &iEta_HBM  );
  HBMmaptree->SetBranchAddress("iPhi",  &iPhi_HBM  );
  HBMmaptree->SetBranchAddress("iDepth",&iDepth_HBM);
  HBMmaptree->SetBranchAddress("RBX" ,  &RBX_HBM   );
  HBMmaptree->SetBranchAddress("RM"  ,  &RM_HBM    );
  //
  HEMmaptree->SetBranchAddress("iEta",  &iEta_HEM  );
  HEMmaptree->SetBranchAddress("iPhi",  &iPhi_HEM  );
  HEMmaptree->SetBranchAddress("iDepth",&iDepth_HEM);
  HEMmaptree->SetBranchAddress("RBX" ,  &RBX_HEM   );
  HEMmaptree->SetBranchAddress("RM"  ,  &RM_HEM    );
  //
  HBPmaptree->SetBranchAddress("iEta",  &iEta_HBP  );
  HBPmaptree->SetBranchAddress("iPhi",  &iPhi_HBP  );
  HBPmaptree->SetBranchAddress("iDepth",&iDepth_HBP);
  HBPmaptree->SetBranchAddress("RBX" ,  &RBX_HBP   );
  HBPmaptree->SetBranchAddress("RM"  ,  &RM_HBP    );
  //
  HEPmaptree->SetBranchAddress("iEta",  &iEta_HEP  );
  HEPmaptree->SetBranchAddress("iPhi",  &iPhi_HEP  );
  HEPmaptree->SetBranchAddress("iDepth",&iDepth_HEP);
  HEPmaptree->SetBranchAddress("RBX" ,  &RBX_HEP   );
  HEPmaptree->SetBranchAddress("RM"  ,  &RM_HEP    );
  //
}


int HBM_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBMmaptree->GetEntries(); ientry++){
    HBMmaptree->GetEntry(ientry);
    if( iEta_HBM==abs(ieta) && iPhi_HBM==iphi && iDepth_HBM==idepth ){  output=int((RBX_HBM-int(1))*int(4)+RM_HBM); break; }
  }
  //
  return output;
}
int HBM_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBMmaptree->GetEntries(); ientry++){
    HBMmaptree->GetEntry(ientry);
    if( iEta_HBM==abs(ieta) && iPhi_HBM==iphi && iDepth_HBM==idepth ){  output=int(RBX_HBM); break; }
  }
  //
  return output;
}

int HEM_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEMmaptree->GetEntries(); ientry++){
    HEMmaptree->GetEntry(ientry);
    if( iEta_HEM==abs(ieta) && iPhi_HEM==iphi && iDepth_HEM==idepth ){  output=int((RBX_HEM-int(1))*int(4)+RM_HEM); break; }
  }
  //
  return output;
}
int HEM_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEMmaptree->GetEntries(); ientry++){
    HEMmaptree->GetEntry(ientry);
    if( iEta_HEM==abs(ieta) && iPhi_HEM==iphi && iDepth_HEM==idepth ){  output=int(RBX_HEM); break; }
  }
  //
  return output;
}

int HBP_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBPmaptree->GetEntries(); ientry++){
    HBPmaptree->GetEntry(ientry);
    if( iEta_HBP==ieta && iPhi_HBP==iphi && iDepth_HBP==idepth ){  output=int((RBX_HBP-int(1))*int(4)+RM_HBP); break; }
  }
  //
  return output;
}
int HBP_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBPmaptree->GetEntries(); ientry++){
    HBPmaptree->GetEntry(ientry);
    if( iEta_HBP==ieta && iPhi_HBP==iphi && iDepth_HBP==idepth ){  output=int(RBX_HBP); break; }
  }
  //
  return output;
}


int HEP_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEPmaptree->GetEntries(); ientry++){
    HEPmaptree->GetEntry(ientry);
    if( iEta_HEP==ieta && iPhi_HEP==iphi && iDepth_HEP==idepth ){  output=int((RBX_HEP-int(1))*int(4)+RM_HEP); break; }
  }
  //
  return output;
}
int HEP_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEPmaptree->GetEntries(); ientry++){
    HEPmaptree->GetEntry(ientry);
    if( iEta_HEP==ieta && iPhi_HEP==iphi && iDepth_HEP==idepth ){  output=int(RBX_HEP); break; }
  }
  //
  return output;
}
void CloseHCALmaps( ){

  HBMmapfile->Close();
  HEMmapfile->Close();
  HBPmapfile->Close();
  HEPmapfile->Close();

}

public:

analysisClass(){

HBMmapfile = NULL;
HEMmapfile = NULL;
HBPmapfile = NULL;
HEPmapfile = NULL;

HBMmaptree = NULL;
HEMmaptree = NULL;
HBPmaptree = NULL;
HEPmaptree = NULL;

iEta_HBM = 0;
iPhi_HBM = 0;
iDepth_HBM = 0;
RBX_HBM = 0;
RM_HBM = 0;

iEta_HEM = 0;
iPhi_HEM = 0;
iDepth_HEM = 0;
RBX_HEM = 0;
RM_HEM = 0;

iEta_HBP = 0;
iPhi_HBP = 0;
iDepth_HBP = 0;
RBX_HBP = 0;
RM_HBP = 0;

iEta_HEP = 0;
iPhi_HEP = 0;
iDepth_HEP = 0;
RBX_HEP = 0;
RM_HEP = 0;

ReadHCALmaps();

}

~analysisClass(){
CloseHCALmaps();
/*
delete HBMmapfile;
delete HEMmapfile;
delete HBPmapfile;
delete HEPmapfile;

delete HBMmaptree;
delete HEMmaptree;
delete HBPmaptree;
delete HEPmaptree;
*/
}

int EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=888;
  //
  if(ieta>  0 && ieta< 17 && idepth<3 ) output=HBP_EtaPhitoRBXrm(ieta,iphi,idepth); //HBP
  if(ieta<  0 && ieta>-17 && idepth<3 ) output=HBM_EtaPhitoRBXrm(ieta,iphi,idepth); //HBM
  if(ieta== 16 && idepth==3 ) output=HEP_EtaPhitoRBXrm(ieta,iphi,idepth); //HEP
  if(ieta==-16 && idepth==3 ) output=HEM_EtaPhitoRBXrm(ieta,iphi,idepth); //HEM
  if(ieta> 16 && ieta< 30 ) output=HEP_EtaPhitoRBXrm(ieta,iphi,idepth); //HEP
  if(ieta<-16 && ieta>-30 ) output=HEM_EtaPhitoRBXrm(ieta,iphi,idepth); //HEM
  //
  return output;
}

int EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  // Barrel
  // HBM 1-18
  // HBP 19-36
  //
  // Endcap
  // HEM 37-54
  // HEP 55-72
  //
  int output=888;
  //
  if( ieta<  0  && ieta>-17 && idepth<3  ) output = HBM_EtaPhitoRBX(ieta,iphi,idepth); //HBM
  if( ieta>  0  && ieta< 17 && idepth<3  ) output = HBP_EtaPhitoRBX(ieta,iphi,idepth)+18; //HBP
  if( ieta==-16 &&             idepth==3 ) output = HEM_EtaPhitoRBX(ieta,iphi,idepth)+36; //HEM
  if( ieta<-16  && ieta>-30              ) output = HEM_EtaPhitoRBX(ieta,iphi,idepth)+36; //HEM
  if( ieta== 16 &&             idepth==3 ) output = HEP_EtaPhitoRBX(ieta,iphi,idepth)+54; //HEP
  if( ieta> 16  && ieta< 30              ) output = HEP_EtaPhitoRBX(ieta,iphi,idepth)+54; //HEP
  //
  return output;
}

};
