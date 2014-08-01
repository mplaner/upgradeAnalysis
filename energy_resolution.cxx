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
  
  TH1F* drHist;
  TH1F *etaHist;
  TH1F *ratio_ptHist[nEta];
  TH1F *energy_Hist[nEta];
  TH1F *pt_Hist[nEta];
  vector <float> ratio_ptVect[nEta];
  vector <float> e3x3_ratioVect[nEta];
  vector <float> e5x5_ratioVect[nEta];
  TTree *b=0;

  higgs_gg * h = new higgs_gg(b);
  
  int CONV_CHECK =0;
  int PARENT   =0;
  int CHILD    =0;
  if(Type==HGG)
    {
      CONV_CHECK=0;
      PARENT  = 25;
      CHILD  = 22;
    }
  if(Type==GJETS)
    {
      CONV_CHECK = -1;
    }
  char hname[40];
  etaHist = new TH1F("pEta","pEta",100,-3.2,3.2);
  drHist = new TH1F("dr","dr",10000,0,1);
  
  for(int i=0;i<nEta;i++)
    {
      sprintf(hname,"ratio_pt_#eta_%1.2f",etaVal[i]);
      ratio_ptHist[i] = new TH1F(hname,hname,100,0,3);
      sprintf(hname,"energy_#eta_%1.2f",etaVal[i]);
      energy_Hist[i] = new TH1F(hname,hname,100,0,300);
      sprintf(hname,"pt_#eta_%1.2f",etaVal[i]);
      pt_Hist[i] = new TH1F(hname,hname,100,0,300);
    }
  if (h->fChain == 0) return;

  Long64_t nentries = h->fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  const int nPhotons=200; 
   double g_pt[nPhotons], p_energy[nPhotons];//recophotons matched to genPhotons  
   double p_eta[nPhotons],g_eta[nPhotons];
   double p_phi[nPhotons], g_phi[nPhotons];
   double p_3x3[nPhotons], p_5x5[nPhotons];
   double p_r9[nPhotons];
   double p_hovere[nPhotons];
   double p_sigmaIetaIeta[nPhotons];
   double p_deltaR[nPhotons];
   double temp_deltaR[nPhotons];
   int subDet[nPhotons];
   int nGEB=0;
   int nGEE=0;
   std::cout << "number of entries: " << nentries << std::endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
     //for (Long64_t jentry=0; jentry<1000;jentry++) 
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
           p_sigmaIetaIeta[i]=-99;
           p_hovere[i]=-99;
	   temp_deltaR[i]=0;
	 }
       int nGen=0;
       
       for(int i=0; i< h->NGEN; i++)
	 {
	   //etaHist->Fill(h->GEN_phi[i]);
	   //std::cout << "h->GEN_id[i] " << h->GEN_id[i] << " gen pt " << h->GEN_pt[i] << std::endl;	   
	   //if(h->GEN_fromH[i]) //check if its a higgs decay particle
	   if(Type==HGG)
	     {
	       if(!h->GEN_fromH[i]) //check if it's a higgs decay
		 continue;
	       //if(!h->isConv[i]==1) //if notconverted
	       if(!h->isConv[i]==0) // if converted 
		 continue;  //continue
	       if(!(std::abs(h->GEN_id[i])==CHILD)) //check if it's a photon (not really necessary...)
		 continue;
	     }
	   if(Type==GJETS)
	     {
	       if(h->GEN_status[i]==-2)
		 continue;
	     }
	   
	   //if(( h->GEN_pt[i]*cosh(h->GEN_eta[i]) )<100)
	   //if(( h->GEN_pt[i]*cosh(h->GEN_eta[i]) )>70)
	     if(h->GEN_pt[i]>g_ptCut)
	     {
	       g_eta[nGen] = h->GEN_eta[i];
	       //    std::cout << "eta: " << g_eta[nGen] << std::endl;
	       if(fabs(g_eta[nGen])>EBMAX)
		 if(fabs(g_eta[nGen])<EEMIN) //in crack
		   continue;
	       if(fabs(g_eta[nGen])>EEMAX) //outside detector
		 continue;
	       if(g_eta[nGen]< -EBMAX)
		 {
		   //std::cout << "ee- : " << g_eta[nGen] << std::endl;
		   continue;
		 }
	       g_pt[nGen] = h->GEN_pt[i];
	       g_phi[nGen] = h->GEN_phi[i];
	       nGen++;
	     }
	   if(Type==HGG&&nGen>2)
	     std::cout << "too many decay particles: " << nGen << std::endl;
	 }
       /******checks position of GEN electrons in detector*******/
       for(int i=0; i<nGen; i++)
	 {
	   if(GetEtaBin(g_eta[i])==-1) //checks photon is in detector
	     g_pt[i]=0;
	   if(g_pt[i]>0)
	     {
	       if((g_eta[i])>EEMIN) //for now EE+ only
		 nGEE++;
	       if(fabs(g_eta[i])<EBMAX) //for now EE+ only
		 nGEB++;
	     }
	 }
       /******Reco PT cut selection*****************/
       /******matching Algorithm**********/
       for(int l=0;l<nGen;l++)
	 for(int i=0; i<h->p_size; i++)
	   {
	     if(h->e5[i]/cosh(h->eta[i])<p_ptCut) //apply ptCut
	       continue;
	     
	     if(!g_pt[l]) //if good gen photon
	       continue;
	     //if(GetEtaBin(h->eta[i])==-1) //checks photon is in detector
	     // continue;
	     temp_deltaR[l] = dR(g_eta[l],g_phi[l],h->eta[i],h->phi[i]);
	     drHist->Fill(temp_deltaR[l]);
	     if(p_deltaR[l]> temp_deltaR[l]) //if better dR
	       {
		 p_deltaR[l]    = temp_deltaR[l];
		 p_energy[l]    = h->ecalenergy[i];
		 //p_energy[l]    = h->e5[i];  //5x5 has corrections and es added
		 //p_energy[l]    = h->e3[i] + h->esenergy[i]; //has no crack correction/es added
		 p_eta[l]       = h->eta[i];
		 p_phi[l]       = h->phi[i];
		 p_3x3[l]       =h->e3[i]+ h->esenergy[i];
		 p_5x5[l]       =h->e5[i];
		 p_hovere[l]    =h->hovere[i];
		 p_sigmaIetaIeta[l]=h->sieie[i];	
		 //p_r9[l]        = h->PHO_r9[i];
		 if(std::abs(p_eta[l])>EEMIN&&std::abs(p_eta[l])<EEMAX)
		   subDet[l]=EE;
		 if(std::abs(p_eta[l])>EBMIN&&std::abs(p_eta[l])<EBMAX)
		   subDet[l]=EB;
	       }
	   }
	 for(int l=0; l<nGen; l++)
	   {
	     /****check that p_energy is filled*****/
	     if(p_energy[l]) 
	       {
		 etaHist->Fill(p_eta[l]);
		 energy_Hist[GetEtaBin(g_eta[l])]->Fill(p_5x5[l]);
		 pt_Hist[GetEtaBin(g_eta[l])]->Fill(p_5x5[l]/cosh(p_eta[l]));
		 ratio_ptHist[GetEtaBin(g_eta[l])]->Fill(p_5x5[l]/g_pt[l]/cosh(g_eta[l]));
		 ratio_ptVect[GetEtaBin(g_eta[l])].push_back(p_energy[l]/g_pt[l]/cosh(g_eta[l]));
		 e3x3_ratioVect[GetEtaBin(g_eta[l])].push_back(p_3x3[l]/g_pt[l]/cosh(g_eta[l]));
		 e5x5_ratioVect[GetEtaBin(g_eta[l])].push_back(p_5x5[l]/g_pt[l]/cosh(g_eta[l]));
		 if(subDet[l]==EB)
		   {
		     //EB stuff
		   }
		 if(subDet[l]==EE)
		   {
		     //EE stuff
		   }	
	       }
	   }
}
   TFile * temp = new TFile("testout1.root","RECREATE");
   std::ofstream fileA;
   std::ofstream fileB;
   std::ofstream fileC;
   fileA.open("full_range_sc_error.txt", std::ofstream::app);
   fileB.open("full_range_3x3_error.txt",std::ofstream::app);
   fileC.open("full_range_5x5_error.txt",std::ofstream::app);
   if(fileA.is_open())
     std::cout << "writing to file" << std::endl;
   //   fileA << "starting energy resolution" << std::endl;
   //fileA << "etaVal " << " min " << " max " << " effSigma " << std::endl;
   const char * scenarioName[4] = {"<PU>: 70, int lumi: 0fb^{-1}","<PU>: 140, int lumi: 1000fb^{-1}","<PU>: 140, int lumi: 3000fb^{-1}", "Run 1 MC"};
   fileA << scenarioName[FILETYPE-1] << std::endl;   
   fileB << scenarioName[FILETYPE-1] << std::endl;   
   fileC << scenarioName[FILETYPE-1] << std::endl;   
   
   string legName[nEta];
int colors[nEta] = {1,2,4,5,6,8};
   for(int i=0;i<nEta;i++)
     {
       float min=0,max=0, effSigma=0, eHigh, eLow;
       effSigma = GetEffSigma(.68, ratio_ptVect[i], min, max, eHigh, eLow);
       fileA << etaVal[i] << " " << min << " " << max << " " << effSigma << " " << eLow << " " << eHigh << std::endl;
       effSigma = GetEffSigma(.68, e3x3_ratioVect[i], min, max, eHigh, eLow);
       fileB << etaVal[i] << " " << min << " " << max << " " << effSigma << " " << eLow << " " << eHigh << std::endl;
       effSigma = GetEffSigma(.68, e5x5_ratioVect[i], min, max, eHigh, eLow);
       fileC << etaVal[i] << " " << min << " " << max << " " << effSigma << " " << eLow << " " << eHigh << std::endl;
       ratio_ptHist[i]->Write();
       char tempName[40];
       sprintf(tempName,"eta %1.2f",etaVal[i]);
       legName[i] = tempName;
     }
   //string legendname[2]={"RECO_PHOTON", "GEN_MATCH"};
   /**format hist example: formatHisto(uselogY,"name for plot", array of strings for legend names, "xtitle",xmin,xmax,"ytitle",array of ints for colors, pointer to histogram, number of histograms to be plotted on same canvas*/
   formatHisto(0, scenarioName[FILETYPE-1], legName, "#eta of matched reco photons",-3.2,3.2,"PDF", colors,&etaHist,1, "etaplot");
   formatHisto(0, scenarioName[FILETYPE-1], legName, "Ratio (reco pt)/(gen pt)",.5,1.3,"PDF", colors,&ratio_ptHist[0],nEta,"ptRatio");
   formatHisto(0, scenarioName[FILETYPE-1], legName, "reco pt",0,300,"PDF", colors,&pt_Hist[0],nEta,"pt"); 
   formatHisto(0, scenarioName[FILETYPE-1], legName, "reco energy",0,300,"PDF", colors,&energy_Hist[0],nEta,"energy");

   legName[0]="#eta distribution";
   //formatHisto(0, "legendary", legName);
   temp->Close();
   std::cout << "nGEB: " << nGEB << " nGEE: " << nGEE << std::endl;
}
