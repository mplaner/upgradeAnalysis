////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 22 10:15:36 2013 by ROOT version 5.32/03
// from TTree gsf_electron/gsf_electron
// found on file: es_ntuple_electron.root
//////////////////////////////////////////////////////////

#ifndef higgs_gg_h
#define higgs_gg_h


#include <string.h>
#include <TROOT.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TTree.h>
#include <TProfile.h>
#include <TChain.h>
#include <TFile.h>
#include <cmath>
#include <TCanvas.h>
#include <TStyle.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
//make two trees, one of w/ es and one w/o es

class higgs_gg {
 public :
  TTree          *fChain;   //!pointer to the analyzed  TTree or TChain
  Int_t          fCurrent; //!current Tree number in  TChain
  
  // Declaration of leaf types
  Float_t         ecalenergy[200];
  Float_t         photonenergy[200];
  Float_t         pt[200];
  Double_t        esenergy[200];
  Double_t        rawenergy[200];
  Double_t        eta[200];
  Double_t        phi[200];
  Float_t         e3[200];
  Float_t         e5[200];
  Float_t         r9[200];
  Float_t         sieie[200];
  Float_t         ecalIso[200];
  Float_t         hcalIso[200];
  Float_t         trkIso[200];
  Int_t           p_size;

  
  
  Double_t        GEN_eta[5];
  Double_t        GEN_phi[5];
  Float_t         GEN_pt[5];
  Int_t           GEN_id[5];
  Int_t           NGEN;
  Int_t           run_number;
  Int_t           lumi_number;
  Int_t           event_number;
  
  // List of branches
  TBranch        *b_esenergy;
  TBranch        *b_rawenergy;
  TBranch        *b_photonenergy;
  TBranch        *b_ecalenergy;   //!
  TBranch        *b_pt;
  TBranch        *b_5x5_energy;
  TBranch        *b_3x3_energy;
  TBranch        *b_r9;
  TBranch        *b_sieie;
  TBranch        *b_ecalIso;
  TBranch        *b_hcalIso;
  TBranch        *b_trkIso;
  TBranch        *b_p_size;
  TBranch        *b_g_eta;
  TBranch        *b_g_phi;
  TBranch        *b_g_pt;
  TBranch        *b_g_size;
  TBranch        *b_eta;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_pdg_id;   //!
  TBranch        *b_run_number;
  TBranch        *b_lumi_number;
  TBranch        *b_event_number;
  
  
  higgs_gg(TTree *tree);
  virtual ~higgs_gg();
  virtual Int_t     GetEntry(Long64_t entry);
  virtual Long64_t  LoadTree(Long64_t entry);
  virtual void      Init(TTree *tree); //particle is particle id
  virtual int       GetEtaBin(double eta);
  };

higgs_gg::higgs_gg(TTree *tree) : fChain(0)  
{
  // used to generate this class and read the Tree.
  //for some reason, must open file first
  //std::cout << "in constructor " << std::endl;
  //  if (tree == 0) 
  if (tree == 0) 
    {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntuple_PU50_age_0.root");
      if (!f || !f->IsOpen()) 
	{
	  f = new TFile("ntuple_PU50_age_0.root");
	}
      TDirectory * dir = (TDirectory*)f->Get("ntuple_PU50_age_0.root:/demo");
      dir->GetObject("photon",tree);
      std::cout <<tree->GetEntries() << std::endl;
    }
  std::cout << tree << std::endl;
  Init(tree);
  
  std::cout << "initialized tree " << std::endl;
}

higgs_gg::~higgs_gg()
{
  std::cout << "         fchain " << fChain << std::endl;
  if (!fChain) return;
  fChain->GetCurrentFile()->ls();
  delete fChain;
  //fChain->GetCurrentFile()->Close("R");
}

Int_t higgs_gg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t higgs_gg::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void higgs_gg::Init(TTree *tree)
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

   fChain->SetBranchAddress("pEsenergy", &esenergy, &b_esenergy);
   fChain->SetBranchAddress("pRawenergy", &rawenergy, &b_rawenergy);
   fChain->SetBranchAddress("pSC_energy", &ecalenergy, &b_ecalenergy);
   fChain->SetBranchAddress("pPhoton_energy", &photonenergy, &b_photonenergy);
   fChain->SetBranchAddress("pPt", &pt, &b_pt);
   fChain->SetBranchAddress("p5x5_energy", &e5, &b_5x5_energy);
   fChain->SetBranchAddress("p3x3_energy", &e3, &b_3x3_energy);
   fChain->SetBranchAddress("pR9", &r9, &b_r9);
   fChain->SetBranchAddress("pSigmaIetaIeta", &sieie, &b_sieie);
   fChain->SetBranchAddress("p_ecalRecHitSumEtConeDR03", &ecalIso, &b_ecalIso);
   fChain->SetBranchAddress("p_hcalTowerSumEtConeDR03", &hcalIso, &b_hcalIso);
   fChain->SetBranchAddress("p_trkSumPtSolidConeDR03", &trkIso, &b_trkIso);
   fChain->SetBranchAddress("pEta", &eta, &b_eta);
   fChain->SetBranchAddress("pPhi", &phi, &b_phi);
   fChain->SetBranchAddress("psize", &p_size, &b_p_size);
      
   
   fChain->SetBranchAddress("gEta", &GEN_eta, &b_g_eta);
   fChain->SetBranchAddress("gPhi", &GEN_phi, &b_g_phi);
   
   fChain->SetBranchAddress("gsize", &NGEN, &b_g_size);
   fChain->SetBranchAddress("pdgid", &GEN_id, &b_pdg_id);
   fChain->SetBranchAddress("gPt", &GEN_pt, &b_g_pt);
      
   std::cout << "branches loaded " << std::endl;
}
#endif
