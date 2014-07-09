#define zAlexey_cxx
#include "photon.h"
#include "params.h"
#include "formats.h"
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
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
  
  TH1F * r9Hist[2];
//  TH1F * r9Hist= new TH1F("0/fb","barrel Zpeak",20,0,1);
//  TH1F * r9Hist1= new TH1F("fb","barrel Zpeak",20,0,1);
  TH1F* drHist;
  TH1F *etaHist;
  TH1F *ratio_ptHist[nEta];
  vector <float> ratio_ptVect[nEta];
  TTree *b=0;

  higgs_gg * h = new higgs_gg(b);
  
  int PARENT   =0;
  int CHILD    =0;
  if(type==HGG)
    {
      int  PARENT  = 25;
      int  CHILD  = 22;
    }
  char hname[40];
  etaHist = new TH1F("pEta","pEta",100,-3.2,3.2);
  drHist = new TH1F("dr","dr",10000,0,1);
  r9Hist[0]=new TH1F("hname","hname",20,0,1.2);
  r9Hist[1]=new TH1F("hname1","hname1",20,0,1.2);
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
   double p_3x3[nPhotons], p_5x5[nPhotons];
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
           p_3x3[i]=-99;
           p_5x5[i]=-99;
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
	     if(h->ecalenergy[i]/cosh(h->eta[i])<p_ptCut) //apply ptCut
	       continue;
             r9Hist[0]->Fill(h->e3[i]/h->e5[i]);
	     if(!g_pt[l]) //if good gen phoctron
	       continue;
	     if(GetEtaBin(h->eta[i])==-1) //checks photon is in detector
	       continue;
	     temp_deltaR[l] = dR(g_eta[l],g_phi[l],h->eta[i],h->phi[i]);
	     drHist->Fill(temp_deltaR[l]);
	     if(p_deltaR[l]> temp_deltaR[l]) //if better dR
	       {
		 p_deltaR[l]    = temp_deltaR[l];
		 p_energy[l]    = h->ecalenergy[i];
		 p_eta[l]       = h->eta[i];
		 p_phi[l]       = h->phi[i];
                 p_3x3[l]       =h->e3[i];
                 p_5x5[l]       =h->e5[i];		
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
	       etaHist->Fill(p_eta[l]);
	       ratio_ptHist[GetEtaBin(p_eta[l])]->Fill(p_energy[l]/g_pt[l]/cosh(g_eta[l]));
	       ratio_ptVect[GetEtaBin(p_eta[l])].push_back(p_energy[l]/g_pt[l]/cosh(g_eta[l]));
	       r9Hist[1]->Fill(p_3x3[l]/p_5x5[l]);
	     }
	 }
     }
   TFile * temp = new TFile("testout.root","RECREATE");
   std::ofstream fileA;
   fileA.open("plots.txt");
   if(fileA.is_open())
     std::cout << "writing to file" << std::endl;
   fileA << "starting energy resolution" << std::endl;
   fileA << "etaVal " << " min " << " max " << " effSigma " << std::endl;
   
   //   string legName[] = {"<PU>: 50, int lumi: 0fb^{-1}"};
   string legName[nEta];
   int colors[nEta] = {1,2,4,5,6,8};
   for(int i=0;i<nEta;i++)
     {
       float min=0,max=0, effSigma=0;
       effSigma = GetEffSigma(.68,ratio_ptVect[i], min, max);
       fileA << etaVal[i] << " " << min << " " << max << " " << effSigma << std::endl;
       ratio_ptHist[i]->Write();
       char tempName[40];
       sprintf(tempName,"eta %1.2f",etaVal[i]);
       legName[i] = tempName;
     }
   string legendname[2]={"RECON_PHOTON", "GEN_MATCH"};
   /**format hist example: formatHisto(uselogY,"name for plot", array of strings for legend names, "xtitle",xmin,xmax,"ytitle",array of ints for colors, pointer to histogram, number of histograms to be plotted on same canvas*/
   formatHisto(0,"<PU>: 50, int lumi: 0fb^{-1}", legName, "Ratio (reco pt)/(gen pt)",.5,1.3,"PDF", colors,&ratio_ptHist[0],nEta);
   legName[0]="#eta distribution";
   formatHisto(0,"<PU>: 50, int lumi: 0fb^{-1}", legName, "#eta of matched reco photons",-3.2,3.2,"PDF", colors,&etaHist,1);
   formatHisto(0,"<PU>: 50, R9 DISTRIBUTION",legendname," #R9 of photons",0,1.2,"PDF", colors,&r9Hist[0],2);
 //  r9Hist->Write();
   //formatHisto(0, "legendary", legName);
    
//   r9Hist1->Draw();
 //  r9Hist->Draw();
   temp->Close();
   std::cout << "nGEB: " << nGEB << " nGEE: " << nGEE << std::endl;
}
