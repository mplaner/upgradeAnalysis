#define zAlexey_cxx
#include "photon.h"
#include "params.h"
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void Loop();

int main()
{
  Loop();
}

void Loop()
{
  TH1F* drHist;
  TH1F *etaHist;
  TH1F *ratio_ptHist[nEta];
  vector <float> ratio_ptVect[nEta];
  TTree *b;

  higgs_gg * h = new higgs_gg(b);
  
  int PARENT   =0;
  int CHILD    =0;
  if(type==HGG)
    {
      int  PARENT  = 25;
      int  CHILD  = 22;
    }
  char hname[40];
  for(int i=0;i<nEta;i++)
    {
      sprintf(hname,"ratio_pt_#eta_%1.2f",etaVal[i]);
      ratio_ptHist[i] = new TH1F(hname,hname,100,0,3);
    }
  if (h->fChain == 0) return;

   Long64_t nentries = h->fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   const int nPhotons=2; //hardcoded for two h->gg photons for now
   double g_pt[nPhotons], p_energy[nPhotons];//recophotons matched to genPhotons  
   double p_eta[nPhotons],g_eta[nPhotons];
   double p_phi[nPhotons], g_phi[nPhotons];
   double p_r9[nPhotons];
   double p_deltaR[nPhotons];
   double temp_deltaR[nPhotons];
   int subDet[nPhotons];
   
   
   int nGEB=0;
   int nGEE=0;
   std::cout << "number of entries: " << nentries << std::endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
     {
       Long64_t ientry = h->LoadTree(jentry);
       if (ientry < 0) break;
       nb = h->fChain->GetEntry(jentry);   nbytes += nb;
       /***************initialize electron values********************/
       for(int i=0; i<nPhotons; i++) 
	 {
	   subDet[i]=0;
	   p_energy[i]=0;
	   g_eta[i]=-99;
	   g_phi[i]=-99;
	   p_eta[i]=99;
	   p_phi[i]=99;
	   g_pt[i]=0;
	   p_r9[i]=0;
	   p_deltaR[i]=drCut;
	   temp_deltaR[i]=0;
	 }
       int nGen=0;
       
       for(int i=0; i< h->NGEN; i++)
	 {
	   /****Gets GEN electrons from electron*****/
	   //etaHist->Fill(h->GEN_phi[i]);
	   //std::cout << "h->GEN_id[i] " << h->GEN_id[i] << " gen pt " << h->GEN_pt[i] << std::endl;	   
	   if(std::abs(h->GEN_id[i])==CHILD)
		 if(h->GEN_pt[i]>g_ptCut)
		   {
		     nGen++;
		     if(g_pt[0])
		       {
			 g_pt[1] = h->GEN_pt[i];
			 g_eta[1] = h->GEN_eta[i];
			 g_phi[1] = h->GEN_phi[i];
		       }
		     else
		       {
			 g_pt[0] = h->GEN_pt[i];
			 g_eta[0] = h->GEN_eta[i];
			 g_phi[0] = h->GEN_phi[i];
		       }
		     if(nGen>2)
		       std::cout << "too many decay particles: " << nGen << std::endl;
		   }
	 }
       //std::cout << "jentry: " << jentry << " nGen " << nGen << std::endl;
	      

       /******checks position of GEN electrons in detector*******/
       for(int i=0; i<nPhotons; i++)
	 {
	   /***sets pt=0 if in gap****/
	   if(std::abs(g_eta[i])>EBMAX&& std::abs(g_eta[i])<EEMIN)
	     g_pt[i]=0;
	   /***set pt=0 if outside detector***/
	   if(std::abs(g_eta[i])>EEMAX)
	     g_pt[i]=0;
	   //	   if(std::abs(g_eta[i])<EBMAX)
	   // g_pt[i]=0;
	   if(g_pt[i]>0)
	     {
	       if(std::abs(g_eta[i])>EEMIN)
		 nGEE++;
	       if(std::abs(g_eta[i])<EBMAX)
		 nGEB++;
	     }
	 }

       /******matching Algorithm**********/

       for(int l=0;l<nPhotons;l++)
	 for(int i=0; i<h->p_size; i++)
	   {
	     if(!g_pt[l]) //if good gen phoctron
	       continue;
	     if(h->ecalenergy[i]/cosh(h->eta[i])<p_ptCut) //apply ptCut
	       continue;
	     if(GetEtaBin(h->eta[i])==-1) //checks photon is in detector
	       continue;
	     //	     std::cout << "nEntry : " << jentry << " nGen " << l << std::endl;
	     /*****fix dR later *********/	
	     temp_deltaR[l] = dR(g_eta[l],g_phi[l],h->eta[i],h->phi[i]);
	     //std::cout << "temp dr: " << temp_deltaR[l] << " dr cut: " << p_deltaR[l] << std::endl;
	     //drHist->Fill(temp_deltaR[l]);
	     if(p_deltaR[l]> temp_deltaR[l]) //if better dR
	       {
		   //     std::cout << jentry << " ngen : " << l << " delta r: " << temp_deltaR[l] << std::endl;
		   p_deltaR[l]    = temp_deltaR[l];
		   p_energy[l]    = h->ecalenergy[i];
		   p_eta[l]       = h->eta[i];
		   p_phi[l]       = h->phi[i];
		   //p_r9[l]        = h->PHO_r9[i];
		   if(std::abs(p_eta[l])>EEMIN&&std::abs(p_eta[l])<EEMAX)
		     subDet[l]=EE;
		   if(std::abs(p_eta[l])>EBMIN&&std::abs(p_eta[l])<EBMAX)
		     subDet[l]=EB;
		 }
	   }
       for(int l=0; l<nPhotons; l++)
	 {
	   /****check that p_energy is filled*****/
	   if(p_energy[l]) 
	     {
	       ratio_ptHist[GetEtaBin(p_eta[l])]->Fill(p_energy[l]/g_pt[l]/cosh(g_eta[l]));
	       ratio_ptVect[GetEtaBin(p_eta[l])].push_back(p_energy[l]/g_pt[l]/cosh(g_eta[l]));
	     }
	 }
     }
   TFile * temp = new TFile("testout.root","RECREATE");
   ofstream fileA;
   fileA.open("plots.txt");
   if(fileA.is_open())
     std::cout << "writing to file" << std::endl;
   fileA << "starting energy resolution" << std::endl;
   fileA << "etaVal " << " min " << " max " << " effSigma " << std::endl;
   for(int i=0;i<nEta;i++)
     {
       float min=0,max=0, effSigma=0;
       effSigma = GetEffSigma(.68,ratio_ptVect[i], min, max);
       fileA << etaVal[i] << " " << min << " " << max << " " << effSigma << std::endl;
       ratio_ptHist[i]->Write();
     }
   std::cout << "nGEB: " << nGEB << " nGEE: " << nGEE << std::endl;
}
