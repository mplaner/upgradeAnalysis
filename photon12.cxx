#define zAlexey_cxx
#include "photon9.h"
#include "params.h"
#include "formats10.h"
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>

int photon12()
{ 
  TFile *file1=TFile::Open("hgg_70_age0_ntuple.root");
  TTree *t1 = (TTree*)file1->Get("ntuple/photon");
  TFile *file2=TFile::Open("gjet_70_age0_ntuple.root");
  TTree *t2 = (TTree*)file2->Get("ntuple/photon");
  TFile *file3=TFile::Open("hgg_140_age1000_ntuple.root");
  TTree *t3 = (TTree*)file3->Get("ntuple/photon");
  TFile *file4=TFile::Open("gjet_140_age1000_ntuple.root");
  TTree *t4 = (TTree*)file4->Get("ntuple/photon");
  TFile *file5=TFile::Open("hgg_140_age3000_ntuple.root");
  TTree *t5 = (TTree*)file5->Get("ntuple/photon");
  TFile *file6=TFile::Open("gjet_140_age3000_ntuple.root");
  TTree *t6 = (TTree*)file6->Get("ntuple/photon");
 // int HIK=0;
//  string filename[3]={"<PU>: 50, int lumi: 0fb^{-1}","<PU>: 100, int lumi: 0fb^{-1}","<PU>: 140, int lumi: 0fb^{-1}"};
  //Loop(t1);
 char *legname="<PU>70 age0";
 char *dirctoryname="PLOTS/%s";   
//  onegroup(t1,t2,legname,dirctoryname);
//  legname="<PU>140 age1000";
//  dirctoryname="PLOTS140/%s";
//   onegroup(t3,t4,legname,dirctoryname);
//  legname="<PU>140 age3000";
//  dirctoryname="PLOTS140_3000/%s";
//  onegroup(t5,t6,legname,dirctoryname);
  threegroup2(t1,t2,t3,t4,t5,t6,dirctoryname);
//  threegroup(t1,t2,t3,t4,t5,t6,dirctoryname);
// threegroup1(t1,t2,t3,t4,t5,t6,char *dirctoryname);
}

void Loop(int filenum,TTree *tree,TH1F *r9EEHist1,TH1F *r9EBHist1,TH1F *hovereEBHist1,TH1F *hovereEEHist1,TH1F * sigmIetIetaEBHist1,TH1F * sigmIetIetaEEHist1)
{
    
  TH1F* drHist;
  TH1F *etaHist;
  TH1F *ratio_ptHist[nEta];
  vector <float> ratio_ptVect[nEta];
  higgs_gg * h = new higgs_gg(tree);
  
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
  if (h->fChain == 0) return;

   Long64_t nentries = h->fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   if(filenum==0)
      const int nPhotons=30; //hardcoded for two h->gg photons for now
   if(filenum==1)
      {
        const int nPhotons=30; 
        CHILD  = 111;
      }
   double g_pt[nPhotons], p_energy[nPhotons];//recophotons matched to genPhotons  
   double p_eta[nPhotons],g_eta[nPhotons];
   double p_phi[nPhotons], g_phi[nPhotons];
   double p_3x3[nPhotons], p_5x5[nPhotons];
   double p_r9[nPhotons];
   double p_hovere[nPhotons];
   double p_sigmaIetaIeta[nPhotons];
   double p_deltaR[nPhotons];
   double temp_deltaR[nPhotons];
 //  float r9[nPhotons];
   int subDet[nPhotons];
   std::cout << "number of entries: " << nentries << std::endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
     {
      
      int nGEB=0;
      int nGEE=0;
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
	   p_r9[i]=99;
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
	   /****Gets GEN electrons from electron*****/
	   if(std::abs(h->GEN_id[i])==CHILD)
		 if(h->GEN_pt[i]>g_ptCut)
		   {
                      
			 g_pt[i] = h->GEN_pt[i];
			 g_eta[i] = h->GEN_eta[i];
			 g_phi[i] = h->GEN_phi[i];
		     nGen++;
		   }
	 }
	      
       /******checks position of GEN electrons in detector*******/
       for(int i=0; i<nPhotons; i++)
	 {
	   /***sets pt=0 if in gap****/
	   if(std::fabs(g_eta[i])>EBMAX&& std::fabs(g_eta[i])<EEMIN)
	     g_pt[i]=0;
	   /***set pt=0 if outside detector***/
	   if(std::fabs(g_eta[i])>EEMAX)
	     g_pt[i]=0;
	   //	   if(std::abs(g_eta[i])<EBMAX)
	   // g_pt[i]=0;
	   if((g_eta[i])<-EBMAX)
	     g_pt[i]=0;
	   if(g_pt[i]>0)
	     {
	       if(std::fabs(g_eta[i])>EEMIN)
                     nGEE++;
	       if(std::fabs(g_eta[i])<EBMAX)
 		    nGEB++;
	     }
	 }
       /******Reco PT cut selection*****************/
       /******matching Algorithm**********/

       for(int l=0;l<nPhotons;l++)
        {
          double temp_ptdiff=30;
   
	 for(int i=0; i<h->p_size; i++)
	   {
	     if(h->pt[i]<p_ptCut) //apply ptCut
	       continue;
	     if(!g_pt[l]) //if good gen phoctron
	       continue;
             if(filenum==1)
              {
	           if(h->GEN_ePt[l]>10)
			continue; 
               //    if(h->GEN_pPt[l]>10)
               //         continue;     
              }
		temp_deltaR[l] = dR(g_eta[l],g_phi[l],h->eta[i],h->phi[i]);
	     drHist->Fill(temp_deltaR[l]);
	     if(p_deltaR[l]> temp_deltaR[l]) //if better dR
	       {
                if(std::fabs(h->pt[i]-g_pt[l])<temp_ptdiff)
                 {
                  temp_ptdiff   =std::fabs(h->pt[i]-g_pt[l]);
	          // p_deltaR[l]    = temp_deltaR[l];
		  p_energy[l]    = h->ecalenergy[i];
		  p_eta[l]       = h->eta[i];
		  p_phi[l]       = h->phi[i];
                  p_3x3[l]       =h->e3[i];
                  p_5x5[l]       =h->e5[i];	
                
                  p_hovere[l]    =h->hovere[i];
                  p_sigmaIetaIeta[l]=h->sieie[i];	
                  p_r9[l]        = h->r9[i];
		 if(std::fabs(p_eta[l])>EEMIN&&std::fabs(p_eta[l])<EEMAX)
		   subDet[l]=EE;
		 if(std::fabs(p_eta[l])>EBMIN&&std::fabs(p_eta[l])<EBMAX)
		   subDet[l]=EB;
	        }
              }
	    }
        }
        for(int l=0; l<nPhotons; l++)
 	  {
	   /****check that p_energy is filled*****/
	   if(p_energy[l]) 
	     {
	       etaHist->Fill(p_eta[l]);
	       if(subDet[l]==EB)
                 {
	          // r9EBHist1->Fill(p_3x3[l]/p_5x5[l]);
	           r9EBHist1->Fill(p_r9[l]);
	           hovereEBHist1->Fill(p_hovere[l]);
                   sigmIetIetaEBHist1->Fill(p_sigmaIetaIeta[l]);
 	         }
	       if(subDet[l]==EE)
               {
	        // r9EEHist1->Fill(p_3x3[l]/p_5x5[l]);
	         r9EEHist1->Fill(p_r9[l]);
	         hovereEEHist1->Fill(p_hovere[l]);
                 sigmIetIetaEEHist1->Fill(p_sigmaIetaIeta[l]);
	       }	
	     }
	 }
        
       /******Reco PT cut selection*****************/
    }
   /**format hist example: formatHisto(uselogY,"name for plot", array of strings for legend names, "xtitle",xmin,xmax,"ytitle",array of ints for colors, pointer to histogram, number of histograms to be plotted on same canvas*/
}
//ok
void integrateHis(TH1F *hist1,TH1F *hist2)
 {
   Float_t sum = 0;
   int size=hist1->GetSize();
//   cout<<"the size of all bin"<<size<<endl;
   for (Int_t i=1;i<=size-2;i++)
   {
       sum +=hist1->GetBinContent(i);
       hist2->SetBinContent(i,sum);
   }
 }
void inverseintegrateHis(TH1F *hist1,TH1F *hist2)
 {
   Float_t sum = 0;
   int size=hist1->GetSize();
   for (Int_t i=1;i<=size-2;i++)
   {
       sum +=hist1->GetBinContent(i);
       hist2->SetBinContent(i,1-sum);
   }
 }

void normalize (TH1F *hist1)
{
 Double_t integral = hist1->TH1F::Integral();
 if (integral > 0)
 hist1->Scale(1/integral);
}

void efficencycal(TH1F *reco_hist,TH1F *gen_hist,TH1F *efficenHist)
{
 int size=reco_hist->GetSize();
// float retio=0;
// cout<<size<<endl;
 for(int i=1;i<=size-2;i++)
  {
    float genvalue=gen_hist->GetBinContent(i);
    if(genvalue!=0)
     {
       efficenHist->SetBinContent(i,(reco_hist->GetBinContent(i)/genvalue));
     }

  }
 
}
void hispasscontent(TH1F *hist1,TH1F *hist2)
{
int size=hist1->GetSize();
for(int i=1;i<size-2;i++)
    hist2->SetBinContent(i,hist1->GetBinContent(i));

}
void inverseefficencycal(TH1F *reco_hist,TH1F *gen_hist,TH1F *efficenHist)
{
 int size=reco_hist->GetSize();
// float retio=0;
// cout<<size<<endl;
 for(int i=1;i<=size-2;i++)
  {
    float genvalue=gen_hist->GetBinContent(i);
    if(genvalue!=0)
     {
       efficenHist->SetBinContent(i,(1-reco_hist->GetBinContent(i)/genvalue));
     }

  }
 
}



//TGraph *roc[3];
//void ROCCOVE(TH1F *hist1,TH1F *hist2,TGraph *roc,int binnum)
void ROCCOVE(TH1F *hist1,TH1F *hist2,TH1F *hist3,TH1F *hist4,TH1F *hist5,TH1F *hist6,TGraph *roc1,TGraph *roc2,TGraph *roc3,int binnum)
    {
      int k=0;
   //   int l=i*2;
   //   double x[binnum],y[binnum];
      //roc = new TGraph(binnum,x,y);
      for(int j=0;j<binnum;j++)
	{
	  {
	    roc1->SetPoint(k,hist1->GetBinContent(j),hist2->GetBinContent(j));
	    roc2->SetPoint(k,hist3->GetBinContent(j),hist4->GetBinContent(j));
	    roc3->SetPoint(k,hist5->GetBinContent(j),hist6->GetBinContent(j));
	    k++;
	  }
	}
    }


//void onegroup(TTree *tree1,TTree *tree2)
void onegroup(TTree *tree1,TTree *tree2,char *legname,char *dirctoryname)
{

  TH1F * r9EBHist[2];
  TH1F * r9EBsleHist[2];
  TH1F * r9EEHist[2];
  TH1F * r9EEsleHist[2];
  TH1F * hovereEBHist[2];
  TH1F * hovereEEHist[2];
  TH1F * hovereEBsleHist[2];
  TH1F * hovereEEsleHist[2];
  TH1F * sigmIetIetaEBHist[2];
  TH1F * sigmIetIetaEEHist[2];
  TH1F * sigmIetIetaEBsleHist[2];
  TH1F * sigmIetIetaEEsleHist[2];
  r9EBHist[0]=new TH1F("hname","hname",100,0,1);
  r9EBHist[1]=new TH1F("hname1","hname1",100,0,1);
  r9EBsleHist[0]=new TH1F("hname11","hname_",100,0,1);
  r9EBsleHist[1]=new TH1F("hname111","hname_",100,0,1);
  r9EEHist[0]=new TH1F("hname1111","hname_",100,0,1);
  r9EEsleHist[0]=new TH1F("hname11111","hname_",100,0,1);
  r9EEHist[1]=new TH1F("hname111111","hname_1",100,0,1);
  r9EEsleHist[1]=new TH1F("hname1111111","hname_",100,0,1);
  hovereEBHist[0]=new TH1F("hname2","hname2",200,-1,1);
  hovereEBHist[1]=new TH1F("hname22","hname3",200,-1,1);
  hovereEEHist[0]=new TH1F("hname222","hname_2",200,-1,1);
  hovereEEHist[1]=new TH1F("hname2222","hname_3",200,-1,1);
  hovereEBsleHist[0]=new TH1F("hname22222","hname_3",200,-1,1);
  hovereEBsleHist[1]=new TH1F("hname222222","hname_3",200,-1,1);
  hovereEEsleHist[0]=new TH1F("hname2222222","hname_3",200,-1,1);
  hovereEEsleHist[1]=new TH1F("hname22222222","hname_3",200,-1,1);
  sigmIetIetaEBHist[0]=new TH1F("hname3","hname4",100,0,0.1);
  sigmIetIetaEBHist[1]=new TH1F("hname33","hname5",100,0,0.1);
  sigmIetIetaEEHist[0]=new TH1F("hname_333","hname_4",100,0,0.1);
  sigmIetIetaEEHist[1]=new TH1F("hname_3333","hname_5",100,0,0.1);
  sigmIetIetaEBsleHist[0]=new TH1F("hname33333","hname4",100,0,0.1);
  sigmIetIetaEBsleHist[1]=new TH1F("hname333333","hname4",100,0,0.1);
  sigmIetIetaEEsleHist[0]=new TH1F("hname3333333","hname4",100,0,0.1);
  sigmIetIetaEEsleHist[1]=new TH1F("hname33333333","hname4",100,0,0.1);
  string legendname[2]={"RECON_PHOTON EFFICIENCY", "JET_BACKGROUND REJECTION"};
  int colors[nEta] = {1,2,4,5,6,8};
  Loop(0,tree1,r9EEHist[0],r9EBHist[0],hovereEBHist[0],hovereEEHist[0],sigmIetIetaEBHist[0],sigmIetIetaEEHist[0]);
  Loop(1,tree2,r9EEHist[1],r9EBHist[1],hovereEBHist[1],hovereEEHist[1],sigmIetIetaEBHist[1],sigmIetIetaEEHist[1]);
  normalize(r9EBHist[0]);
  normalize(r9EBHist[1]);
  normalize(r9EEHist[0]);
  normalize(r9EEHist[1]);
  normalize(hovereEBHist[0]);
  normalize(hovereEBHist[1]);
  normalize(hovereEEHist[0]);
  normalize(hovereEEHist[1]);
  normalize(sigmIetIetaEBHist[0]);
  normalize(sigmIetIetaEBHist[1]);
  normalize(sigmIetIetaEEHist[0]);
  normalize(sigmIetIetaEEHist[1]);
//  hovereEBHist[0]->Draw();
//  hovereEBHist[1]->Draw("SAME");;
  inverseintegrateHis(r9EBHist[0],r9EBsleHist[0]); 
  inverseintegrateHis(r9EEHist[0],r9EEsleHist[0]); 
  integrateHis(hovereEBHist[0],hovereEBsleHist[0]); 
  integrateHis(hovereEEHist[0],hovereEEsleHist[0]);
  integrateHis(sigmIetIetaEBHist[0],sigmIetIetaEBsleHist[0]); 
  integrateHis(sigmIetIetaEEHist[0],sigmIetIetaEEsleHist[0]); 
//  integrateHis(r9EEHist[1],r9EEsleHist[1]); 
  integrateHis(r9EBHist[1],r9EBsleHist[1]); 
  integrateHis(r9EEHist[1],r9EEsleHist[1]);
  inverseintegrateHis(hovereEBHist[1],hovereEBsleHist[1]); 
  inverseintegrateHis(hovereEEHist[1],hovereEEsleHist[1]); 
  inverseintegrateHis(sigmIetIetaEBHist[1],sigmIetIetaEBsleHist[1]); 
  inverseintegrateHis(sigmIetIetaEEHist[1],sigmIetIetaEEsleHist[1]); 
  //char *dirc="PLOTS/%s";
   formatHisto(0,legname,legendname," #R9 of photons",0,1,"efficiency/rejection", colors,&r9EBsleHist[0],2,"r9EBsleplot2.png",dirctoryname);
   formatHisto(0,legname,legendname," #R9 of photons",0,1,"efficiency/rejection", colors,&r9EEsleHist[0],2,"r9EEsleplot2.png",dirctoryname);
   formatHisto(0,legname,legendname," #HoverE of photons",-0.1,0.2,"efficiency/rejection", colors,&hovereEBsleHist[0],2,"hovereEBsleplot2.png",dirctoryname);
   formatHisto(0,legname,legendname," #HoverE of photons",-0.1,0.4,"efficiency/rejection", colors,&hovereEEsleHist[0],2,"hovereEEsleplot2.png",dirctoryname);
   formatHisto(0,legname,legendname," #SIGMAIetaIeta of photons",0,1,"efficiency/rejection", colors,&sigmIetIetaEBsleHist[0],2,"sigmaIetaIetaEB2.png",dirctoryname);
   formatHisto(0,legname,legendname," #SIGMAIetaIeta of photons",0,1,"efficiency/rejection", colors,&sigmIetIetaEEsleHist[0],2,"sigmaIetaIetaEE2.png",dirctoryname);
//  hovereEBsleHist[0]->Draw();
delete[] r9EBHist; delete[] r9EBsleHist;delete[] r9EEHist;delete[] r9EEsleHist;delete[] hovereEBHist;
delete[] hovereEEHist; delete[] hovereEBsleHist;delete[] hovereEEsleHist;delete[]  sigmIetIetaEBHist;
delete[] sigmIetIetaEEHist; delete[] sigmIetIetaEBsleHist; delete[] sigmIetIetaEEsleHist;
//r9EEsleHist[0]->Draw();
//  r9EEsleHist[1]->SetLineColor(6);
 // r9EEsleHist[1]->Draw();
  //r9EEHist[1]->Draw();
 // delete  
}
void threegroup(TTree *tree1,TTree *tree2,TTree *tree3,TTree *tree4,TTree *tree5,TTree *tree6,char *dirctoryname)
{
 float r9EB[3]={0.65,0.66,0.66},r9EE[3]={0.75,0.73,0.5},HoverE_EB[3]={0.12,0.16,0.17},HoverE_EE[3]={0.20,0.32,0.335},SigmaIetaIetaEB[3]={0.011,0.014,0.018},SigmaIetaIetaEE[3]={0.033,0.052,0.061};
// float r9EB[3]={0,0,0},r9EE[3]={0,0,0},HoverE_EB[3]={1,1,1},HoverE_EE[3]={1,1,1},SigmaIetaIetaEB[3]={1,1,1},SigmaIetaIetaEE[3]={1,1,1};
 
 TH1F *eta_RECO[6];
 TH1F *eta_GEN[6];
 TH1F *eta_efficenHist[6];
 TH1F *pt_RECOEB[6];
 TH1F *pt_GENEB[6];
 TH1F *pt_efficencyEB[6];
 TH1F *pt_RECOEE[6];
 TH1F *pt_GENEE[6];
 TH1F *pt_efficencyEE[6];
 eta_RECO[0]=new TH1F("etareco0","etareco",100,-EEMAX,EEMAX);
 eta_RECO[1]=new TH1F("etareco1","etareco",100,-EEMAX,EEMAX);
 eta_GEN[0]= new TH1F("etagen0","etagen",100,-EEMAX,EEMAX);
 eta_GEN[1]= new TH1F("etagen1","etagen",100,-EEMAX,EEMAX);
 eta_efficenHist[0]=new TH1F("etaeffi","etaeffi",100,-EEMAX,EEMAX);
 eta_efficenHist[1]=new TH1F("etaeffi1","etaeffi",100,-EEMAX,EEMAX);
 pt_RECOEE[0] = new TH1F("ptreco0EE","pereco",100,0,200);
 pt_RECOEE[1] = new TH1F("ptreco1EE","pereco",100,0,200);
 pt_GENEE[0]  = new TH1F("ptgen0EE","ptgen",100,0,200);
 pt_GENEE[1]  = new TH1F("ptgen1EE","ptgen",100,0,200);
 pt_efficencyEE[0]  = new TH1F("pteffi0EE","pteffi",100,0,200);
 pt_efficencyEE[1]  = new TH1F("pteffi1EE","pteffi",100,0,200);
 pt_RECOEB[0] = new TH1F("ptreco0EB","pereco",100,0,200);
 pt_RECOEB[1] = new TH1F("ptreco1EB","pereco",100,0,200);
 pt_GENEB[0]  = new TH1F("ptgen0EB","ptgen",100,0,200);
 pt_GENEB[1]  = new TH1F("ptgen1EB","ptgen",100,0,200);
 pt_efficencyEB[0]  = new TH1F("pteffi0EB","pteffi",100,0,200);
 pt_efficencyEB[1]  = new TH1F("pteffi1EB","pteffi",100,0,200);
 string legendname[6]={"PU70 age0 reco efficiency","PU70 age0 jet REJECTION","PU140 age1000 reco efficiency","PU140 age1000 jet REJECTION","PU140 age3000 reco efficiency","PU140 age3000 jet REJECTION"};
 int colors[nEta] = {1,2,4,5,6,8};
 Loop1(0,tree1,eta_RECO[0],eta_GEN[0],pt_RECOEB[0],pt_GENEB[0],pt_RECOEE[0],pt_GENEE[0],r9EB[0],r9EE[0],HoverE_EB[0],HoverE_EE[0],SigmaIetaIetaEB[0],SigmaIetaIetaEE[0]);
 Loop1(1,tree2,eta_RECO[1],eta_GEN[1],pt_RECOEB[1],pt_GENEB[1],pt_RECOEE[1],pt_GENEE[1],r9EB[0],r9EE[0],HoverE_EB[0],HoverE_EE[0],SigmaIetaIetaEB[0],SigmaIetaIetaEE[0]);
 efficencycal(eta_RECO[0],eta_GEN[0],eta_efficenHist[0]);
 efficencycal(pt_RECOEB[0],pt_GENEB[0],pt_efficencyEB[0]);
 efficencycal(pt_RECOEE[0],pt_GENEE[0],pt_efficencyEE[0]);
 //eta_RECO[0]->Draw();
// eta_GEN[0]->Draw("SAME");
 inverseefficencycal(eta_RECO[1],eta_GEN[1],eta_efficenHist[1]);
 inverseefficencycal(pt_RECOEB[1],pt_GENEB[1],pt_efficencyEB[1]);
 inverseefficencycal(pt_RECOEE[1],pt_GENEE[1],pt_efficencyEE[1]);
delete  eta_RECO[0]; delete eta_GEN[0];delete pt_RECOEB[0];
delete pt_GENEB[0];
delete pt_RECOEE[0];
delete pt_GENEE[0];
//cout<<"where is it0"<<endl;
delete  eta_RECO[1]; delete eta_GEN[1];delete pt_RECOEB[1];
//delete pt_GENEB[1];
delete pt_RECOEE[1];
//delete pt_GENEE[1];

//cout<<"where is it1"<<endl;
 eta_RECO[2]=new TH1F("etareco2","etareco",100,-EEMAX,EEMAX);
 eta_RECO[3]=new TH1F("etareco3","etareco",100,-EEMAX,EEMAX);

 eta_GEN[2]= new TH1F("etagen2","etagen",100,-EEMAX,EEMAX);
 eta_GEN[3]= new TH1F("etagen3","etagen",100,-EEMAX,EEMAX);

 eta_efficenHist[2]=new TH1F("etaeffi2","etaeffi",100,-EEMAX,EEMAX);
 eta_efficenHist[3]=new TH1F("etaeffi3","etaeffi",100,-EEMAX,EEMAX);

 pt_RECOEE[2] = new TH1F("ptreco2EE","pereco",100,0,200);
 pt_RECOEE[3] = new TH1F("ptreco3EE","pereco",100,0,200);

 pt_GENEE[2]  = new TH1F("ptgen2EE","ptgen",100,0,200);
 pt_GENEE[3]  = new TH1F("ptgen3EE","ptgen",100,0,200);

 pt_efficencyEE[2]  = new TH1F("pteffi2EE","pteffi",100,0,200);
 pt_efficencyEE[3]  = new TH1F("pteffi3EE","pteffi",100,0,200);

 pt_RECOEB[2] = new TH1F("ptreco2EB","pereco",100,0,200);
 pt_RECOEB[3] = new TH1F("ptreco3EB","pereco",100,0,200);


 pt_GENEB[2]  = new TH1F("ptgen2EB","ptgen",100,0,200);
 pt_GENEB[3]  = new TH1F("ptgen3EB","ptgen",100,0,200);
 

 pt_efficencyEB[2]  = new TH1F("pteffi2EB","pteffi",100,0,200);
 pt_efficencyEB[3]  = new TH1F("pteffi3EB","pteffi",100,0,200);
 
 Loop1(0,tree3,eta_RECO[2],eta_GEN[2],pt_RECOEB[2],pt_GENEB[2],pt_RECOEE[2],pt_GENEE[2],r9EB[1],r9EE[1],HoverE_EB[1],HoverE_EE[1],SigmaIetaIetaEB[1],SigmaIetaIetaEE[1]);

 Loop1(1,tree4,eta_RECO[3],eta_GEN[3],pt_RECOEB[3],pt_GENEB[3],pt_RECOEE[3],pt_GENEE[3],r9EB[1],r9EE[1],HoverE_EB[1],HoverE_EE[1],SigmaIetaIetaEB[1],SigmaIetaIetaEE[1]);

 efficencycal(eta_RECO[2],eta_GEN[2],eta_efficenHist[2]);
 efficencycal(pt_RECOEB[2],pt_GENEB[2],pt_efficencyEB[2]);
 efficencycal(pt_RECOEE[2],pt_GENEE[2],pt_efficencyEE[2]);
 inverseefficencycal(eta_RECO[3],eta_GEN[3],eta_efficenHist[3]);
 inverseefficencycal(pt_RECOEB[3],pt_GENEB[3],pt_efficencyEB[3]);
 inverseefficencycal(pt_RECOEE[3],pt_GENEE[3],pt_efficencyEE[3]);
delete  eta_RECO[2]; delete eta_GEN[2];delete pt_RECOEB[2];delete pt_GENEB[2];delete pt_RECOEE[2];delete pt_GENEE[2];
delete  eta_RECO[3]; delete eta_GEN[3];delete pt_RECOEB[3];
//delete pt_GENEB[3];
delete pt_RECOEE[3];
//delete pt_GENEE[3];


 eta_RECO[4]=new TH1F("etareco4","etareco",100,-EEMAX,EEMAX);
 eta_RECO[5]=new TH1F("etareco5","etareco",100,-EEMAX,EEMAX);

 eta_GEN[4]= new TH1F("etagen4","etagen",100,-EEMAX,EEMAX);
 eta_GEN[5]= new TH1F("etagen5","etagen",100,-EEMAX,EEMAX);

 eta_efficenHist[4]=new TH1F("etaeffi4","etaeffi",100,-EEMAX,EEMAX);
 eta_efficenHist[5]=new TH1F("etaeffi5","etaeffi",100,-EEMAX,EEMAX);

 pt_RECOEE[4] = new TH1F("ptreco4EE","pereco",100,0,200);
 pt_RECOEE[5] = new TH1F("ptreco5EE","pereco",100,0,200);

 pt_GENEE[4]  = new TH1F("ptgen4EE","ptgen",100,0,200);
 pt_GENEE[5]  = new TH1F("ptgen5EE","ptgen",100,0,200);

 pt_efficencyEE[4]  = new TH1F("pteffi4EE","pteffi",100,0,200);
 pt_efficencyEE[5]  = new TH1F("pteffi5EE","pteffi",100,0,200);

 pt_RECOEB[4] = new TH1F("ptreco4EB","pereco",100,0,200);
 pt_RECOEB[5] = new TH1F("ptreco5EB","pereco",100,0,200);


 pt_GENEB[4]  = new TH1F("ptgen4EB","ptgen",100,0,200);
 pt_GENEB[5]  = new TH1F("ptgen5EB","ptgen",100,0,200);
 

 pt_efficencyEB[4]  = new TH1F("pteffi4EB","pteffi",100,0,200);
 pt_efficencyEB[5]  = new TH1F("pteffi5EB","pteffi",100,0,200);
 
 Loop1(0,tree5,eta_RECO[4],eta_GEN[4],pt_RECOEB[4],pt_GENEB[4],pt_RECOEE[4],pt_GENEE[4],r9EB[2],r9EE[2],HoverE_EB[2],HoverE_EE[2],SigmaIetaIetaEB[2],SigmaIetaIetaEE[2]);

 Loop1(1,tree6,eta_RECO[5],eta_GEN[5],pt_RECOEB[5],pt_GENEB[5],pt_RECOEE[5],pt_GENEE[5],r9EB[2],r9EE[2],HoverE_EB[2],HoverE_EE[2],SigmaIetaIetaEB[2],SigmaIetaIetaEE[2]);

 efficencycal(eta_RECO[4],eta_GEN[4],eta_efficenHist[4]);
 efficencycal(pt_RECOEB[4],pt_GENEB[4],pt_efficencyEB[4]);
 efficencycal(pt_RECOEE[4],pt_GENEE[4],pt_efficencyEE[4]);
 inverseefficencycal(eta_RECO[5],eta_GEN[5],eta_efficenHist[5]);
 inverseefficencycal(pt_RECOEB[5],pt_GENEB[5],pt_efficencyEB[5]);
 inverseefficencycal(pt_RECOEE[5],pt_GENEE[5],pt_efficencyEE[5]);
delete  eta_RECO[4]; delete eta_GEN[4];delete pt_RECOEB[4];delete pt_GENEB[4];delete pt_RECOEE[4];delete pt_GENEE[4];
delete  eta_RECO[5]; delete eta_GEN[5];delete pt_RECOEB[5];
//delete pt_GENEB[5];
delete pt_RECOEE[5];
//delete pt_GENEE[5];
 TH1F *GEN_PTEB[3];
GEN_PTEB[0]=new TH1F("genptEB","genpt",100,0,200);
GEN_PTEB[1]=new TH1F("genptEB1","genpt",100,0,200);
GEN_PTEB[2]=new TH1F("genptEB2","genpt",100,0,200);
hispasscontent(pt_GENEB[1],GEN_PTEB[0]);
hispasscontent(pt_GENEB[3],GEN_PTEB[1]);
hispasscontent(pt_GENEB[5],GEN_PTEB[2]);
/*  int binnum1=100;
  double x[100],y[100];
  roc[0] = new TGraph(binnum1,x,y);
  roc[1] = new TGraph(binnum1,x,y);
  roc[2] = new TGraph(binnum1,x,y);*/
  string legendname1[3]={"PU70 aging0","PU140 aging 1000","PU140 aging 3000"}; 
 formatHisto(0,"jet Gen PT EB distribution",legendname1," pt of photons",0,200,"pdf", colors,&GEN_PTEB[0],3,"JETGENEB_ptcut10.png",dirctoryname);
delete GEN_PTEB[0];delete GEN_PTEB[1];delete GEN_PTEB[2];delete pt_GENEB[1];delete pt_GENEB[3];delete pt_GENEB[5];
 TH1F *GEN_PTEE[3];
 GEN_PTEE[0]=new TH1F("genptEE","genpt",100,0,200);
 GEN_PTEE[1]=new TH1F("genptEE1","genpt",100,0,200);
 GEN_PTEE[2]=new TH1F("genptEE2","genpt",100,0,200);
 hispasscontent(pt_GENEE[1],GEN_PTEE[0]);
 hispasscontent(pt_GENEE[3],GEN_PTEE[1]);
 hispasscontent(pt_GENEE[5],GEN_PTEE[2]);
 formatHisto(0,"jet Gen PT EE distribution",legendname1," pt of photons",0,200,"pdf", colors,&GEN_PTEE[0],3,"JETGENEE_ptcut10.png",dirctoryname);
 delete GEN_PTEE[0];delete GEN_PTEE[1];delete GEN_PTEE[2];delete pt_GENEE[1];delete pt_GENEE[3];delete pt_GENEE[5];
//void ROCCOVE(TH1F *hist1,TH1F *hist2,TGraph *roc,int binnum)
//void ROCCOVE(TH1F *hist1,TH1F *hist2,TH1F *hist3,TH1F *hist4,TH1F *hist5,TH1F *hist6,TGraph *roc,int binnum)
//void ROCCOVE(TH1F *hist1,TH1F *hist2,TH1F *hist3,TH1F *hist4,TH1F *hist5,TH1F *hist6,TGraph *roc1,TGraph *roc2,TGraph *roc3,int binnum)
//ROCCOVE(eta_efficenHist[0],eta_efficenHist[1],eta_efficenHist[2],eta_efficenHist[3],eta_efficenHist[4],eta_efficenHist[5],roc[0],roc[1],roc[2],binnum1);
//************************************
/*  TMultiGraph *mg = new TMultiGraph();
  mg->Add(roc[0],"AC");
  mg->Add(roc[1],"AC");
  mg->Add(roc[2],"AC");
  mg->SetTitle("title;signal efficiency;background rejection");

  roc[0]->SetTitle("8TeV rho corrected");
  roc[1]->SetTitle("PU40bx25 rho corrected");
  roc[2]->SetTitle("PU40bx25 rho timecut=15ns");
  
  roc[0]->SetLineColor(1);
  roc[1]->SetLineColor(2);
  roc[2]->SetLineColor(3);

//  roc[0]->SetFillColor(0);
//  roc[1]->SetFillColor(0);
//  roc[2]->SetFillColor(0);
  
  roc[0]->SetLineWidth(3);
  roc[1]->SetLineWidth(3);
  roc[2]->SetLineWidth(3);
  


    mg->Draw();*/
   //   mg->Draw("AC");
  //
  
// formatHisto1(0," ETA Roc curve of reco signals and jets",legendname1,"identified jets rejection rate",0.5,1,"identified signals reconstruction rate ", colors,&roc[0],3,"etaRoc.png",dirctoryname);
 formatHisto(0,"ETA distribution of signal and jets",legendname," eta of photons",0,EEMAX,"efficiency/rejection", colors,&eta_efficenHist[0],6,"Etaefficency.png",dirctoryname);
 formatHisto(0,"EB PT distribution of signal and jets",legendname," pt of photons",0,200,"efficiency/rejection", colors,&pt_efficencyEB[0],6,"PtefficencyEB.png",dirctoryname);
 formatHisto(0,"EE PT distribution of signal and jets",legendname," pt of photons",0,200,"efficiency/rejection", colors,&pt_efficencyEE[0],6,"PtefficencyEE.png",dirctoryname);
// h1.Divide(h2)
//  eta_RECO[0]->Divide(eta_GEN[0]); 
//  eta_GEN[1]->Draw(); 
// eta_RECO[1]->Draw(); 
// eta_efficenHist[1]->Draw();
}//ok3
void Loop1(int filenum,TTree *tree,TH1F * eta_recoHist1,TH1F * eta_genHist1,TH1F *pt_recoEBHist1,TH1F *pt_genEBHist1,TH1F *pt_recoEEHist1,TH1F *pt_genEEHist1,float r9EB,float r9EE,float HoverE_EB,float HoverE_EE,float SigmaIetaIetaEB,float SigmaIetaIetaEE)
//void Loop1(int filenum,TTree *tree,TH1F * eta_recoHist1,TH1F * eta_genHist1,float r9EB,float r9EE,float HoverE_EB,float HoverE_EE,float SigmaIetaIetaEB,float SigmaIetaIetaEE)
{
    
  TH1F* drHist;
//  TH1F *etaHist;
  TH1F *ratio_ptHist[nEta];
  vector <float> ratio_ptVect[nEta];
  higgs_gg * h = new higgs_gg(tree);
  
  int PARENT   =0;
  int CHILD    =0;
  if(type==HGG)
    {
      int  PARENT  = 25;
      int  CHILD  = 22;
    }
  char hname[40];
//  etaHist = new TH1F("pEta","pEta",100,-3.2,3.2);
  drHist = new TH1F("dr","dr",10000,0,1);
  if (h->fChain == 0) return;

   Long64_t nentries = h->fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   if(filenum==0)
      const int nPhotons=100; //hardcoded for two h->gg photons for now
   if(filenum==1)
      {
        const int nPhotons=100; 
        CHILD  = 111;
      }
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
   float p_pt[nPhotons];
   std::cout << "number of entries: " << nentries << std::endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
     {
      
   int nGEB=0;
   int nGEE=0;

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
           p_pt[i]=0;
	 }
       int nGen=0;
       for(int i=0; i< h->NGEN; i++)
	 {
	   /****Gets GEN electrons from electron*****/
	   if(std::abs(h->GEN_id[i])==CHILD)
		 if(h->GEN_pt[i]>g_ptCut)
		   {
                      
			 g_pt[i] = h->GEN_pt[i];
			 g_eta[i] = h->GEN_eta[i];
			 g_phi[i] = h->GEN_phi[i];
		     nGen++;
		   }
	 }
 //       cout<<"nGen is "<<nGen<<"child is "<<CHILD<<endl;	      

       /******checks position of GEN electrons in detector*******/
       for(int i=0; i<nPhotons; i++)
	 {
	   /***sets pt=0 if in gap****/
	   if(std::fabs(g_eta[i])>EBMAX&& std::fabs(g_eta[i])<EEMIN)
	     g_pt[i]=0;
	   /***set pt=0 if outside detector***/
	   if(std::fabs(g_eta[i])>EEMAX)
	     g_pt[i]=0;
	   if((g_eta[i])<-EBMAX)
	     g_pt[i]=0;
	   if(g_pt[i]>0)
	     {
               eta_genHist1->Fill(std::fabs(g_eta[i]));

	       if(std::fabs(g_eta[i])>EEMIN)
                 {  
                     nGEE++;
                     pt_genEEHist1->Fill(g_pt[i]);
                 }
	       if(std::fabs(g_eta[i])<EBMAX)
 		 {
                    nGEB++;
                    pt_genEBHist1->Fill(g_pt[i]);

                 }
	     }
	 }
     //  cout<<"nGEE and nGEB are "<<nGEE+nGEB<<endl;
       /******Reco PT cut selection*****************/
       /******matching Algorithm**********/

       for(int l=0;l<nPhotons;l++)
        {
	  if(g_pt[l]==0) //if good gen phoctron
	      continue;
          float temptdiff=100;
	 for(int i=0; i<h->p_size; i++)
	   {
	     if(h->ecalenergy[i]/cosh(h->eta[i])<p_ptCut) //apply ptCut
	     { 
//	 	cout<<p_ptCut<<endl;
                 continue;
 
 	     }
	     if(!g_pt[l]) //if good gen phoctron
	       continue;
/**************************************photon in detector?*****************************************/
	     //if(GetEtaBin(h->eta[i])==-1) //checks photon is in detector
	     //  continue;
             if(filenum==1)
              {
	           if(h->GEN_ePt[l]>10)
			continue; 
                 //  if(h->GEN_pPt[l]>10)
                 //       continue;     
              }
	     temp_deltaR[l] = dR(g_eta[l],g_phi[l],h->eta[i],h->phi[i]);
	    // drHist->Fill(temp_deltaR[l]);
	     if(p_deltaR[l]> temp_deltaR[l]) //if better dR
	       {
	//	 p_deltaR[l]    = temp_deltaR[l];
                if(std::fabs(g_pt[l]-p_pt[i])<temptdiff)
                {
                 temptdiff=std::fabs(g_pt[l]-p_pt[i]);
		 p_energy[l]    = h->ecalenergy[i];
		 p_eta[l]       = g_eta[l];//to gen eta
                 p_pt[l]        = g_pt[l];//gen pt

  //   cout<<"pt"<<p_pt[l]<<endl;
		 p_phi[l]       = h->phi[i];
                 p_3x3[l]       =h->e3[i];
                 p_5x5[l]       =h->e5[i];	
                 p_hovere[l]    =h->hovere[i];
                 p_sigmaIetaIeta[l]=h->sieie[i];
                 p_r9[l]=h->r9[i];
          //       if(std::fabs(p_eta[l])>2.5)
	//		cout<<"larger "<<p_eta[l]<<endl;
                  
          //       cout<<"nGen"<<nGen<<"eta"<<p_eta[l]<<endl;	
		 if(std::fabs(p_eta[l])>EEMIN&&std::fabs(p_eta[l])<EEMAX)
		   subDet[l]=EE;
		 if(std::fabs(p_eta[l])>EBMIN&&std::fabs(p_eta[l])<EBMAX)
		   subDet[l]=EB;
                }
	       }
	    }
        }
        for(int l=0; l<nPhotons; l++)
 	  {
	   /****check that p_energy is filled*****/
	   if(p_energy[l]) 
	     {
	    //   if(std::fabs(p_eta[l])<3)
             //    {
             //      eta_recoHist1->Fill(p_eta[l]);//HERE________________________
              //   }
                  //  cout<<"ha"<<endl;
	     //  etaHist->Fill(p_eta[l]);
	       if(subDet[l]==EB)
                 {
                   if(p_r9[l]=r9EB&&p_hovere[l]<=HoverE_EB&&p_sigmaIetaIeta[l]<=SigmaIetaIetaEB)
                     {
			 eta_recoHist1->Fill(std::fabs(p_eta[l]));
                         pt_recoEBHist1->Fill(p_pt[l]);
                     }
	       //    r9EBHist1->Fill(p_3x3[l]/p_5x5[l]);
	       //    hovereEBHist1->Fill(p_hovere[l]);
               //    sigmIetIetaEBHist1->Fill(p_sigmaIetaIeta[l]);
 	         }
	       if(subDet[l]==EE)
               {
                   if(p_3x3[l]/p_5x5[l]>=r9EE&&p_hovere[l]<=HoverE_EE&&p_sigmaIetaIeta[l]<=SigmaIetaIetaEE)
			{                   

                             eta_recoHist1->Fill(std::fabs(p_eta[l]));
			     pt_recoEEHist1->Fill(std::fabs(p_pt[l]));	
                         }
	       //  r9EEHist1->Fill(p_3x3[l]/p_5x5[l]);
	       //  hovereEEHist1->Fill(p_hovere[l]);
               //  sigmIetIetaEEHist1->Fill(p_sigmaIetaIeta[l]);
	       }	
	     }
	 }
        
       /******Reco PT cut selection*****************/
    }
}
//
/*void signaltune(TTree *tree1,TTree *tree2,TTree *tree2)
{

 TH1F *eta_RECO[3];
 TH1F *eta_GEN[3];
 TH1F *eta_efficenHist[3];
 eta_RECO[0]=new TH1F("etaroco","eta",0,2.5,50);
 eta_RECO[1]=new TH1F("etaroco1","eta",0,2.5,50);
 eta_RECO[2]=new TH1F("etaroco2","eta",0,2.5,50);
 eta_GEN[0]=new

 Loop2(0,tree1,eta_recoHist1,eta_genHist1,r9EB,r9EE,HoverE_EB,HoverE_EE,SigmaIetaIetaEB,SigmaIetaIetaEE)


}*/
void Loop2(int filenum,TTree *tree,float r9EB1,float r9EB2,float r9EE1,float r9EE2,float HoverE_EB1,float HoverE_EB2,float HoverE_EE1,float HoverE_EE2,float SigmaIetaIetaEB1,float SigmaIetaIetaEB2,float SigmaIetaIetaEE1,float SigmaIetaIetaEE2,char *outputfilename)
{
// float r9EBg1=0,r9EEg1=0,HoverE_EBg1=0,HoverE_EEg1=0,SigmaIetaIetaEBg1=0,SigmaIetaIetaEEg1=0;  
 float intervalr9EB=(r9EB1-r9EB2)/5.0, intervalr9EE=(r9EE1-r9EE2)/5.0,intervalHoverE_EB=(HoverE_EB1-HoverE_EB2)/5.0,intervalHoverE_EE=(HoverE_EE1-HoverE_EE2)/5,intervalSigmaIetaIetaEB=(SigmaIetaIetaEB1-SigmaIetaIetaEB2)/5.0,intervalSigmaIetaIetaEE=(SigmaIetaIetaEE1-SigmaIetaIetaEE2)/5.0;
 // cout<<"intervalr9EB"<<intervalr9EB<<endl;
 // cout<<"intervalr9EB"<<intervalr9EB<<endl;
// float r9EBg3=0,r9EEg3=0,HoverE_EBg3=0,HoverE_EEg3=0,SigmaIetaIetaEBg3=0,SigmaIetaIetaEEg3=0;
// float r9EB=0,r9EE=0,HoverE_EB=0,HoverE_EE=0,SigmaIetaIetaEB=0,SigmaIetaIetaEE=0;   
		 
 std::ofstream fileA;  
 fileA.open(outputfilename, std::ofstream::app);
   if(fileA.is_open())                 
     {
       cout<<outputfilename<<"is opened"<<endl;
     // if(first_file)
            fileA <<outputfilename<<"is following"<<endl;
       fileA <<"r9EB"<<"  "<<"HoverE_EB"<<"  "<<"SigmaIetaIetaEB"<<"  "<<"EBrecobin1effi"<< "  " <<" EBrecobin2effi" << "         "<<"r9EE"<<"  "<<"HoverE_EE"<<"  "<<"SigmaIetaIetaEE"<<"  "<<"EErecobin1effi"<<"  "<<"EErecobin2effi"<<endl;
    }
 
  higgs_gg * h = new higgs_gg(tree);
  
 for(int j1=0;j1<=5;j1++)
 {
    float r9EB=0,r9EE=0;
    r9EB=r9EB2+j1*intervalr9EB;
    r9EE=r9EE2+j1*intervalr9EE;
     for(int j2=0;j2<=5;j2++)
      {
          float HoverE_EB=0,HoverE_EE=0;
          HoverE_EB=HoverE_EB2+j2*intervalHoverE_EB;
          HoverE_EE=HoverE_EE2+j2*intervalHoverE_EE;
         
          for(int j3=0;j3<=5;j3++)
          {
            float SigmaIetaIetaEB=0,SigmaIetaIetaEE=0;           
            float EBrecobin1=0,EBrecobin2=0,EErecobin1=0,EErecobin2=0,EBgenbin1=0,EBgenbin2=0,EEgenbin1=0,EEgenbin2=0;
            SigmaIetaIetaEB=SigmaIetaIetaEB2+j3*intervalSigmaIetaIetaEB;
            SigmaIetaIetaEE=SigmaIetaIetaEE2+j3*intervalSigmaIetaIetaEE;
            cout<<"SigmaIetaIetaEE"<<SigmaIetaIetaEE<<endl;
  int PARENT   =0;
  int CHILD    =0;
  if(type==HGG)
    {
      int  PARENT  = 25;
      int  CHILD  = 22;
    }
  char hname[40];
  if (h->fChain == 0) return;

   Long64_t nentries = h->fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   if(filenum==0)
      const int nPhotons=100; //hardcoded for two h->gg photons for now
   if(filenum==1)
      {
        const int nPhotons=100; 
        CHILD  = 111;
      }
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
   float p_pt[nPhotons];
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
           p_sigmaIetaIeta[i]=-99;
           p_hovere[i]=-99;
	   temp_deltaR[i]=0;
           p_pt[i]=0;
	 }
       int nGen=0;
          //  cout<<"NGEN"<<h->NGEN<<endl;
       for(int i=0; i< h->NGEN; i++)
	 {
	   /****Gets GEN electrons from electron*****/
	   if(std::abs(h->GEN_id[i])==CHILD)
		 if(h->GEN_pt[i]>g_ptCut)
		   {
                      
			 g_pt[i] = h->GEN_pt[i];
			 g_eta[i] = h->GEN_eta[i];
			 g_phi[i] = h->GEN_phi[i];
		     nGen++;
		   }
	 }
 //       cout<<"nGen is "<<nGen<<"child is "<<CHILD<<endl;	      

       /******checks position of GEN electrons in detector*******/
       for(int i=0; i<nPhotons; i++)
	 {
	   /***sets pt=0 if in gap****/
	   if(std::fabs(g_eta[i])>EBMAX&& std::fabs(g_eta[i])<EEMIN)
	     g_pt[i]=0;
	   /***set pt=0 if outside detector***/
	   if(std::fabs(g_eta[i])>EEMAX)
	     g_pt[i]=0;
	   if((g_eta[i])<-EBMAX)
	     g_pt[i]=0;
	   if(g_pt[i]>0)
	     {
		if(std::fabs(g_eta[i])>0&&std::fabs(g_eta[i])<=0.5)
                     EBgenbin1++;  
		if(std::fabs(g_eta[i])>0.5&&std::fabs(g_eta[i])<=1.44)
                     EBgenbin2++;
		if(std::fabs(g_eta[i])>1.567&&std::fabs(g_eta[i])<=2)
                     EEgenbin1++;
		if(std::fabs(g_eta[i])>2&&std::fabs(g_eta[i])<=2.5)
                     EEgenbin2++;
	     }
        }
	 
       /******Reco PT cut selection*****************/
       /******matching Algorithm**********/

       for(int l=0;l<nPhotons;l++)
        {
	  if(g_pt[l]==0) //if good gen phoctron
	      continue;
          float temptdiff=100;
	 for(int i=0; i<h->p_size; i++)
	   {
	     if(h->ecalenergy[i]/cosh(h->eta[i])<p_ptCut) //apply ptCut
	     { 
//	 	cout<<p_ptCut<<endl;
                 continue;
 
 	     }
             if(h->TrkIso[i]>5)
                 continue;
	     if(!g_pt[l]) //if good gen phoctron
	       continue;
/**************************************photon in detector?*****************************************/
	     //if(GetEtaBin(h->eta[i])==-1) //checks photon is in detector
	     //  continue;
             if(filenum==1)
              {
	           if(h->GEN_ePt[l]>10)
			continue; 
                 //  if(h->GEN_pPt[l]>10)
                 //       continue;     
              }
	     temp_deltaR[l] = dR(g_eta[l],g_phi[l],h->eta[i],h->phi[i]);
	     if(p_deltaR[l]> temp_deltaR[l]) //if better dR
	       {
                if(std::fabs(g_pt[l]-p_pt[i])<temptdiff)
                 {
	//	 p_deltaR[l]    = temp_deltaR[l];
                 temptdiff=std::fabs(g_pt[l]-p_pt[i]);
		 p_energy[l]    = h->ecalenergy[i];
		 p_eta[l]       = g_eta[l];//to gen eta
                 p_pt[l]        = g_pt[l];//gen pt
		 p_phi[l]       = h->phi[i];
                 p_3x3[l]       =h->e3[i];
                 p_5x5[l]       =h->e5[i];	
                 p_hovere[l]    =h->hovere[i];
                 p_sigmaIetaIeta[l]=h->sieie[i];
                 p_r9[l]=h->r9[i];
		 if(std::fabs(p_eta[l])>EEMIN&&std::fabs(p_eta[l])<EEMAX)
		   subDet[l]=EE;
		 if(std::fabs(p_eta[l])>EBMIN&&std::fabs(p_eta[l])<EBMAX)
		   subDet[l]=EB;
	       }
             }
	   }
        }
        for(int l=0; l<nPhotons; l++)
 	  {
	   /****check that p_energy is filled*****/
	   if(p_energy[l]) 
	     {
	       if(subDet[l]==EB)
                 {
                   if(p_r9[l]>=r9EB&&p_hovere[l]<=HoverE_EB&&p_sigmaIetaIeta[l]<=SigmaIetaIetaEB)
                     {
		       if(std::fabs(p_eta[l])>0&&std::fabs(p_eta[l])<=0.5)
                           EBrecobin1++;
		       if(std::fabs(p_eta[l])>0.5&&std::fabs(p_eta[l])<=1.44)
                           EBrecobin2++;
                     }
 	         }
	       if(subDet[l]==EE)
               {
                   if(p_r9[l]>=r9EE&&p_hovere[l]<=HoverE_EE&&p_sigmaIetaIeta[l]<=SigmaIetaIetaEE)
			{                   
		       if(std::fabs(p_eta[l])>1.57&&std::fabs(p_eta[l])<=2)
                           EErecobin1++;
		       if(std::fabs(p_eta[l])>2&&std::fabs(p_eta[l])<=2.5)
                           EErecobin2++;

                         }
	       }	
	     }
	 }
        
       /******Reco PT cut selection*****************/
    }
            //float EBrecobin1=0,EBrecobin2=0,EEreocbin1=0,EErecobin2=0,EBgenbin1=0,EBgenbin2=0,EEgenbin1=0,EEgenbin2=0;
    cout<<"the EB1 effficiency"<<EBrecobin1/EBgenbin1<<endl;
    cout<<"the EB2 effficiency"<<EBrecobin2/EBgenbin2<<endl;
    cout<<"the EE1 effficiency"<<EErecobin1/EEgenbin1<<endl;
    cout<<"the EE2 effficiency"<<EErecobin2/EEgenbin2<<endl;
                 //  if(p_r9[l]>=r9EB&&p_hovere[l]<=HoverE_EB&&p_sigmaIetaIeta[l]<=SigmaIetaIetaEB)
                //   if(p_r9[l]>=r9EE&&p_hovere[l]<=HoverE_EE&&p_sigmaIetaIeta[l]<=SigmaIetaIetaEE)
    fileA <<r9EB<<"    "<<HoverE_EB<<"    "<<SigmaIetaIetaEB<<"    "<<EBrecobin1/EBgenbin1<< "    " << EBrecobin2/EBgenbin2 << "                             "<<r9EE<<"    "<<HoverE_EE<<"    "<<SigmaIetaIetaEE<<"    "<<EErecobin1/EEgenbin1<<"    "<<EErecobin2/EEgenbin2<<endl;
    //  else if(not_last_file)
  }
  }
  }
   fileA.close();
}
void threegroup2(TTree *tree1,TTree *tree2,TTree *tree3,TTree *tree4,TTree *tree5,TTree *tree6,char *dirctoryname)
{
 float r9EB1=0.8,r9EB2=0.4,r9EE1=0.85,r9EE2=0.5,HoverE_EB1=0.14,HoverE_EB2=0.025,HoverE_EE1=0.2,HoverE_EE2=0.08,SigmaIetaIetaEB1=0.013,SigmaIetaIetaEB2=0.008,SigmaIetaIetaEE1=0.04,SigmaIetIetaEE2=0.025;
 char *outputfilename="pile70signal.txt";
 Loop2(0,tree1,r9EB1,r9EB2,r9EE1,r9EE2,HoverE_EB1,HoverE_EB2,HoverE_EE1,HoverE_EE2,SigmaIetaIetaEB1,SigmaIetaIetaEB2,SigmaIetaIetaEE1,SigmaIetIetaEE2,outputfilename);


/* outputfilename="pile70background.txt";
 Loop2(1,tree2,r9EB1,r9EB2,r9EE1,r9EE2,HoverE_EB1,HoverE_EB2,HoverE_EE1,HoverE_EE2,SigmaIetaIetaEB1,SigmaIetaIetaEB2,SigmaIetaIetaEE1,SigmaIetIetaEE2,outputfilename);

 outputfilename="pile140aging1000signal.txt";
 r9EB1=0.8,r9EB2=0.4,r9EE1=0.85,r9EE2=0.4,HoverE_EB1=0.18,HoverE_EB2=0.07,HoverE_EE1=0.32,HoverE_EE2=0.15,SigmaIetaIetaEB1=0.013,SigmaIetaIetaEB2=0.01,SigmaIetaIetaEE1=0.031,SigmaIetIetaEE2=0.05;
 Loop2(0,tree3,r9EB1,r9EB2,r9EE1,r9EE2,HoverE_EB1,HoverE_EB2,HoverE_EE1,HoverE_EE2,SigmaIetaIetaEB1,SigmaIetaIetaEB2,SigmaIetaIetaEE1,SigmaIetIetaEE2,outputfilename);
 outputfilename="pile140aging1000background.txt";
 Loop2(1,tree4,r9EB1,r9EB2,r9EE1,r9EE2,HoverE_EB1,HoverE_EB2,HoverE_EE1,HoverE_EE2,SigmaIetaIetaEB1,SigmaIetaIetaEB2,SigmaIetaIetaEE1,SigmaIetIetaEE2,outputfilename);


 r9EB1=0.77,r9EB2=0.4,r9EE1=0.8,r9EE2=0.5,HoverE_EB1=0.16,HoverE_EB2=0.06,HoverE_EE1=0.3,HoverE_EE2=0.15,SigmaIetaIetaEB1=0.018,SigmaIetaIetaEB2=0.013,SigmaIetaIetaEE1=0.042,SigmaIetIetaEE2=0.06;
 outputfilename="pile140aging3000signal.txt";
 Loop2(0,tree5,r9EB1,r9EB2,r9EE1,r9EE2,HoverE_EB1,HoverE_EB2,HoverE_EE1,HoverE_EE2,SigmaIetaIetaEB1,SigmaIetaIetaEB2,SigmaIetaIetaEE1,SigmaIetIetaEE2,outputfilename);
 outputfilename="pile140aging3000background.txt";
 Loop2(1,tree6,r9EB1,r9EB2,r9EE1,r9EE2,HoverE_EB1,HoverE_EB2,HoverE_EE1,HoverE_EE2,SigmaIetaIetaEB1,SigmaIetaIetaEB2,SigmaIetaIetaEE1,SigmaIetIetaEE2,outputfilename);*/
}

