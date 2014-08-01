#include "TH1F.h"
#include <vector>
#include "root.cc"
#include "params.h"
#include "formats.h"
#include <fstream>
#include <iostream>
using namespace std;



//.L combine_histos.C
//setTDRStyle1(1,1,1,1); 
//TCanvas c1("c1","c1",700,600); 
TFile *f1,*f2,*f3,*f4, *f5, *f6;
//setTDRStyle1(1,1,1,1); 
//int RUNTYPE=1; //1 =qcd 0=signal

//combine_histos("Eiso","EB");
//char name[40], subDet[40],tHisto[50]; sprintf(name,"Eiso");sprintf(subDet,"EB");
void load_file(int comp_type)
{
  if(comp_type==1)
    {
      f1 = new TFile("output.root.PU70age0");
      f2 = new TFile("output.root.GJETS.PU70age0");
      f3 = new TFile("output.root.PU140age1k");
      f4 = new TFile("output.root.GJETS.PU140age1k");
      f5 = new TFile("output.root.PU140age3k");
      f6 = new TFile("output.root.GJETS.PU140age3k");
    }
  else if(comp_type==2)
    {
      
      f1 = new TFile("output.root.GJETS.PU70age0");
      f2 = new TFile("output.root.GJETS.PU70age0");
      
      f3 = new TFile("output.root.GJETS.PU140age1k");
      f4 = new TFile("output.root.GJETS.PU140age1k");
      f5 = new TFile("output.root.GJETS.PU140age3k");
      f6 = new TFile("output.root.GJETS.PU140age3k");
    }
  else
    {
      f1 = new TFile("output.root.CON.PU70age0");
      f2 = new TFile("output.root.UCON.PU70age0");
      f3 = new TFile("output.root.CON.PU140age1k");
      f4 = new TFile("output.root.UCON.PU140age1k");
      f5 = new TFile("output.root.CON.PU140age3k");
      f6 = new TFile("output.root.UCON.PU140age3k");
    }
  if(f1)
    std::cout << "loaded file 1" << std::endl;
  if(f2)
    std::cout << "loaded file 2" << std::endl;
  if(f3)
    std::cout << "loaded file 3" << std::endl;
  if(f4)
    std::cout << "loaded file 4" << std::endl;
  if(f5)
    std::cout << "loaded file 5" << std::endl;
  if(f6)
    std::cout << "loaded file 6" << std::endl;
}

int combine_histos(char* name, int comp_type)
{
  
  char nameB[70];
  if(comp_type==2)
    sprintf(nameB,"Unweighted_%s",name);
  else
    sprintf(nameB,"%s",name);
  load_file(comp_type);
  TH1F *HISTO[6];
  string AlegNames[6] = {"Photons_PU70age0","Jets_PU70age0","Photons_PU140age1k","Jets_PU140age1k","Photons_PU140age3k","Jets_PU140age3k"};
  string BlegNames[6] = {"Converted_PU70age0","Unconverted_PU70age0","Converted_PU140age1k","Unconverted_PU140age1k","Converted_PU140age3k","Unconverted_PU140age3k"};
  string ClegNames[6] = {"unweighted_PU70age0","Weighted_PU70age0","Unweighted_PU140age1k","Weighted_PU140age1k","Unweighted_PU140age3k","Weighted_PU140age3k"};
  std::cout << "legName[0] : " << AlegNames[0] << "  " << AlegNames[1] << std::endl;
  std::cout << "nameA : " << name << "  " << nameB << std::endl;
  
  HISTO[0] = (TH1F*)f1->Get(nameB);
  HISTO[0]->SetTitle(AlegNames[0].c_str());
  HISTO[0]->SetName(AlegNames[0].c_str()); 
  
  HISTO[1] = (TH1F*)f2->Get(name);
  HISTO[1]->SetTitle(AlegNames[1].c_str());
  HISTO[1]->SetName(AlegNames[1].c_str());

  HISTO[2] = (TH1F*)f3->Get(nameB);
  HISTO[2]->SetTitle(AlegNames[2].c_str());
  HISTO[2]->SetName(AlegNames[2].c_str());

  HISTO[3] = (TH1F*)f4->Get(name);
  HISTO[3]->SetTitle(AlegNames[3].c_str());
  HISTO[3]->SetName(AlegNames[3].c_str());
			
  HISTO[4] = (TH1F*)f5->Get(nameB);
  HISTO[4]->SetTitle(AlegNames[4].c_str());
  HISTO[4]->SetName(AlegNames[4].c_str());
  
  HISTO[5] = (TH1F*)f6->Get(name);
  HISTO[5]->SetTitle(AlegNames[5].c_str());
  HISTO[5]->SetName(AlegNames[5].c_str());
  
  //void formatHisto(int logY, const char* legendTitle, string * legendNames, const char* xTitle, float xMin, float xMax, const char* yTitle, int *colors, TH1F** hist, const int nPoints,const char* filename)
  if(comp_type==1)
    formatHisto(0,name, AlegNames, name,0,0,"PDF",etaColor,&HISTO[0],6,name);
  else if(comp_type==2)
    formatHisto(0,name, ClegNames, name,0,0,"PDF",etaColor,&HISTO[0],6,name);
  else
    formatHisto(0,name, BlegNames, name,0,0,"PDF",etaColor,&HISTO[0],6,name);
}
int combine_efficieny(char* name, int comp_type, int low_to_high)
{

 load_file(comp_type);
  TH1F *HISTO[6];
  string AlegNames[6] = {"Photons_PU70age0","Jets_PU70age0","Photons_PU140age1k","Jets_PU140age1k","Photons_PU140age3k","Jets_PU140age3k"};
  string BlegNames[6] = {"Converted_PU70age0","Unconverted_PU70age0","Converted_PU140age1k","Unconverted_PU140age1k","Converted_PU140age3k","Unconverted_PU140age3k"};
  std::cout << "legName[0] : " << AlegNames[0] << "  " << AlegNames[1] << std::endl;
  
  HISTO[0] = (TH1F*)f1->Get(name);
  HISTO[0]->SetTitle(AlegNames[0].c_str());
  HISTO[0]->SetName(AlegNames[0].c_str()); 
  
  HISTO[1] = (TH1F*)f2->Get(name);
  HISTO[1]->SetTitle(AlegNames[1].c_str());
  HISTO[1]->SetName(AlegNames[1].c_str());

  HISTO[2] = (TH1F*)f3->Get(name);
  HISTO[2]->SetTitle(AlegNames[2].c_str());
  HISTO[2]->SetName(AlegNames[2].c_str());

  HISTO[3] = (TH1F*)f4->Get(name);
  HISTO[3]->SetTitle(AlegNames[3].c_str());
  HISTO[3]->SetName(AlegNames[3].c_str());
			
  HISTO[4] = (TH1F*)f5->Get(name);
  HISTO[4]->SetTitle(AlegNames[4].c_str());
  HISTO[4]->SetName(AlegNames[4].c_str());
  
  HISTO[5] = (TH1F*)f6->Get(name);
  HISTO[5]->SetTitle(AlegNames[5].c_str());
  HISTO[5]->SetName(AlegNames[5].c_str());
  
  TH1F *eff[6];
  for(int i=0;i<6;i++) //calculate efficiencies                                                                         
    {
      char title[50];
      sprintf(title,"eff_%s",HISTO[i]->GetTitle());
      int nBin = HISTO[i]->GetNbinsX();
      double xmin = HISTO[i]->GetXaxis()->GetXmin();
      double xmax = HISTO[i]->GetXaxis()->GetXmax();
      eff[i]= new TH1F(title,title,nBin,xmin,xmax);
      //input[i]->Draw();                                                                                               
      float integral = HISTO[i]->Integral(0,nBin);
      if(low_to_high)
	for(int j=0;j<nBin;j++)
	  {
	    eff[i]->SetBinContent(j,HISTO[i]->GetBinContent(j));
	    if(j>0)
	      eff[i]->AddBinContent(j,eff[i]->GetBinContent(j-1));
	  }
      else
	for(int j=nBin-1;j>0;j--)
	  {
	    eff[i]->SetBinContent(j,HISTO[i]->GetBinContent(j));
	    if(j>0)
	      eff[i]->AddBinContent(j,eff[i]->GetBinContent(j+1));
	  }
	
      eff[i]->Scale(1.0/integral);
    }
  if(comp_type==1)
   formatEfficiency(0,name, AlegNames, name,0,0,"PDF",etaColor,&eff[0],6,name);
  else
   formatEfficiency(0,name, BlegNames, name,0,0,"PDF",etaColor,&eff[0],6,name);
  /*
    TGraph *roc[3];
    for(int i=0;i<3;i++)
    {
      int k=0;
      int l=i*2;
      double x[800],y[800];
      roc[i] = new TGraph(MAX,x,y);
      for(int j=0;j<MAX;j++)
        {
          {
            roc[i]->SetPoint(k,eff[l]->GetBinContent(j),1-eff[l+1]->GetBinContent(j));
            k++;
          }
        }
    }
  c3.cd();
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(roc[0],"AC");
  mg->Add(roc[1],"AC");
  mg->Add(roc[2],"AC");
  mg->SetTitle("title;signal efficiency;background rejection");

  roc[0]->SetTitle("8TeV rho corrected");
  roc[1]->SetTitle("PU40bx25 rho corrected");
  roc[2]->SetTitle("PU40bx25 rho timecut=15ns");

  roc[0]->SetLineColor(1);
  roc[1]->SetLineColor(2);
  roc[2]->SetLineColor(kViolet-3);

  roc[0]->SetFillColor(0);
  roc[1]->SetFillColor(0);
  roc[2]->SetFillColor(0);

  roc[0]->SetLineWidth(3);
  roc[1]->SetLineWidth(3);
  roc[2]->SetLineWidth(3);
  mg->Draw("AC");
  c3.BuildLegend(.45,.60,.95,.92);
  */
}


int make_weights(char* name, char* fname)
{
  load_file(1); //1 for signal vs background
  TH1F *HISTO[6];
  string legNames[6] = {"Photons_PU70age0","Jets_PU70age0","Photons_PU140age1k","Jets_PU140age1k","Photons_PU140age3k","Jets_PU140age3k"};
  std::cout << "legName[0] : " << legNames[0] << "  " << legNames[1] << std::endl;
  
  HISTO[0] = (TH1F*)f1->Get(name);
  HISTO[0]->SetTitle(legNames[0].c_str());
  HISTO[0]->SetName(legNames[0].c_str()); 
  
  HISTO[1] = (TH1F*)f2->Get(name);
  HISTO[1]->SetTitle(legNames[1].c_str());
  HISTO[1]->SetName(legNames[1].c_str());

  HISTO[2] = (TH1F*)f3->Get(name);
  HISTO[2]->SetTitle(legNames[2].c_str());
  HISTO[2]->SetName(legNames[2].c_str());

  HISTO[3] = (TH1F*)f4->Get(name);
  HISTO[3]->SetTitle(legNames[3].c_str());
  HISTO[3]->SetName(legNames[3].c_str());
			
  HISTO[4] = (TH1F*)f5->Get(name);
  HISTO[4]->SetTitle(legNames[4].c_str());
  HISTO[4]->SetName(legNames[4].c_str());
  
  HISTO[5] = (TH1F*)f6->Get(name);
  HISTO[5]->SetTitle(legNames[5].c_str());
  HISTO[5]->SetName(legNames[5].c_str());
  
  std::ofstream fileA;
  fileA.open(fname);
  if(fileA.is_open())
    std::cout << "writing weigths to file" << std::endl;

  for(int i=0;i<6;i+=2) //calculate weights                                                                         
    {
      int nBin = HISTO[i]->GetNbinsX();
      int nBinB = HISTO[i+1]->GetNbinsX();
      double xmin = HISTO[i]->GetXaxis()->GetXmin();
      double xmax = HISTO[i]->GetXaxis()->GetXmax();
      float integral = HISTO[i]->Integral(0,nBin);
      float integralB = HISTO[i+1]->Integral(0,nBinB);
      double weight =1;
      std::cout << legNames[i] << std::endl;
      std::cout  << "{";
      for(int j=0;j<nBin;j++)
	{
	  /*****weight is (#signal)/(#background) in a bin******/
	  weight = HISTO[i]->GetBinContent(j)/HISTO[i+1]->GetBinContent(j)*integral/integralB;
	  if(HISTO[i+1]->GetBinContent(j)==0)
	    weight =1; 
	  //std::cout << HISTO[i]->GetBinContent(j) << "/" << HISTO[i+1]->GetBinContent(j) << " == ";
	  std::cout << weight << ",";
	}
      std::cout << "};" << std::endl;
    }  
  fileA.close();
  
  //void formatHisto(int logY, const char* legendTitle, string * legendNames, const char* xTitle, float xMin, float xMax, const char* yTitle, int *colors, TH1F** hist, const int nPoints,const char* filename)
  formatHisto(0,name, legNames, name,0,0,"PDF",etaColor,&HISTO[0],6,name);
}
