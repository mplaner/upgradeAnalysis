#include "root.cc"
#include "params.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <vector>

void formatHisto(int logY, const char* legendTitle, string * legendNames, const char* xTitle, float xMin, float xMax, const char* yTitle, int *colors, TH1F** hist, const int nPoints, const char* filename);
void formatHisto(int logY, const char* legendTitle, string * legendNames, const char* xTitle, float xMin, float xMax, const char* yTitle, int *colors, TH1F** hist, const int nPoints,const char* filename)
{
  setTDRStyle1(1,1,0,logY);//gridX/Y, logX/Y
  gStyle->TStyle::SetOptStat(0);
  gStyle->TStyle::SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
  TLegend* leg = new TLegend(0.20,0.65,0.45,0.90,legendTitle);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03679144);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  //leg->SetFillColor(0);
  
  for(int j=0; j<nPoints;j++)
    {
      hist[j]->SetLineColor(colors[j]);
      hist[j]->SetLineWidth(3);
      hist[j]->SetMarkerStyle(20);
      hist[j]->SetMarkerSize(0.9);
      hist[j]->SetMarkerColor(colors[j]);
      
      if(j==0)
	{
	  hist[j]->GetXaxis()->SetNdivisions(505);
	  hist[j]->GetYaxis()->SetNdivisions(505);
	  hist[j]->GetXaxis()->SetTitle(xTitle);
	  hist[j]->GetYaxis()->SetTitle(yTitle);
	  hist[j]->GetYaxis()->SetTitleOffset(1.3);
	  if(xMin!=xMax)
	    hist[j]->GetXaxis()->SetRangeUser(xMin,xMax); //make this automatic
	  /*set y range automatically too*/
	  hist[j]->DrawNormalized();
	}
      else
	hist[j]->DrawNormalized("same");
      TString ll=legendNames[j];
      leg->AddEntry(hist[j],ll,"pl");
    }
  
  leg->Draw("same");
  
  c1->cd();
  TLatex T1;
  T1.SetTextSize(.04014599);
  T1.SetNDC();
  T1.DrawLatex(0.45,0.95,"CMS Simulation Preliminary"); //follows x-axis and y-axis values for position...
  /*
  T1.SetTextSize(.035);
  T1.DrawLatex(0.2,0.5,"Aging model sys uncert.: ~30%");
  T1.DrawLatex(0.2,0.45,"Run1 reconstruction");
  T1.DrawLatex(0.2,0.40,"Prompt photons");
  */
  char  ss1[50];
  sprintf(ss1,"PLOTS/%s.png",filename);
  std::cout << ss1 << std::endl;
  c1->SaveAs(ss1);
  sprintf(ss1,"PLOTS/C/%s.C",filename);
  c1->SaveAs(ss1);
  //c1->Clear();
}

//"Energy resolution, #sigma_{eff}(E)/E"

/*used to format any graph which is produced with nEta points on the x-axis */

void formatEtaGraph(int logY, const char* legendTitle, string * legendNames,const char * yTitle, int * colors, std::vector< std::vector<float> > &data, std::vector< std::vector<float> > &highErrors, std::vector< std::vector<float> > &lowErrors, const int nPoints, const char* fileName );
void formatEtaGraph(int logY, const char* legendTitle, string * legendNames,const char * yTitle, int * colors, std::vector< std::vector<float> > &data, std::vector< std::vector<float> > &highErrors, std::vector< std::vector<float> > &lowErrors, const int nPoints, const char* fileName )
{

  setTDRStyle1(1,1,0,logY);//gridX/Y, logX/Y
  gStyle->TStyle::SetOptStat(0);
  gStyle->TStyle::SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
  TLegend* leg = new TLegend(0.22,0.65,0.55,0.90,legendTitle);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03679144);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  /*
  leg->SetLineColor(kBlue);
  leg->SetLineColorAlpha(kBlue,0.0);
  */
  TGraphAsymmErrors *gr_res[nPoints];
  for (int j=0; j<nPoints; j++)
    {
      float yHigh[nEta];
      float yLow[nEta];
      float y[nEta];
      for (int i=0; i<nEta; i++)
	{
	  yHigh[i]=highErrors[j][i];
	  yLow[i]=lowErrors[j][i];
	  //yHigh[i]=.01;
	  //yLow[i]=.01;
	  y[i]=data[j][i];
	    std::cout << "data: " << y[i] << std::endl;
	}                               //replace with nEta
      gr_res[j] = new TGraphAsymmErrors(6,etaVal,y,etaWidthLow,etaWidthHigh,yHigh,yHigh);
      gr_res[j]->GetY();
      gr_res[j]->SetLineColor(colors[j]);
      gr_res[j]->SetLineWidth(2);
      gr_res[j]->SetMarkerStyle(20);
      gr_res[j]->SetMarkerSize(0.9);
      gr_res[j]->SetMarkerColor(colors[j]);
      if (j==0)
	{
	  gr_res[j]->SetMinimum(0.0);
	  gr_res[j]->SetMaximum(0.25);
	  gr_res[j]->GetXaxis()->SetNdivisions(505);
	  gr_res[j]->GetXaxis()->SetRangeUser(0,2.5);
	  gr_res[j]->GetYaxis()->SetNdivisions(505);
	  gr_res[j]->GetXaxis()->SetTitle("#eta");
	  gr_res[j]->GetYaxis()->SetTitle(yTitle);
	  gr_res[j]->Draw("AP");
	}
      else
	{
	  gr_res[j]->Draw("P");
	}
      TString ll=legendNames[j];
      leg->AddEntry(gr_res[j],ll,"pl");
    }
  leg->Draw("same");
  TLatex T1;
  T1.SetNDC();
  T1.SetTextSize(.04014599);
  T1.DrawLatex(.45,0.95,"CMS Simulation Preliminary"); //follows x-axis and y-axis values for position...
  T1.SetTextSize(.035);
  T1.DrawLatex(0.2,0.5,"Aging model sys uncert.: ~30%");
  T1.DrawLatex(0.2,0.45,"Run1 reconstruction");
  T1.DrawLatex(0.2,0.40,"Prompt photons, pT>30GeV");
  
  char ss1[50];
  sprintf(ss1,"PLOTS/%s.png",fileName);
  c1->SaveAs(ss1);
  sprintf(ss1,"PLOTS/C/%s.C",fileName);
  c1->SaveAs(ss1);
  //c1->Clear();
}

void formatEfficiency(int logY, const char* legendTitle, string * legendNames, const char* xTitle, float xMin, float xMax, const char* yTitle, int *colors, TH1F** hist, const int nPoints, const char* filename);
void formatEfficiency(int logY, const char* legendTitle, string * legendNames, const char* xTitle, float xMin, float xMax, const char* yTitle, int *colors, TH1F** hist, const int nPoints,const char* filename)
{
  setTDRStyle1(1,1,0,logY);//gridX/Y, logX/Y
  gStyle->TStyle::SetOptStat(0);
  gStyle->TStyle::SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
  TLegend* leg = new TLegend(0.20,0.65,0.45,0.90,legendTitle);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03679144);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  //leg->SetFillColor(0);
  
  for(int j=0; j<nPoints;j++)
    {
      hist[j]->SetLineColor(colors[j/2]);
      hist[j]->SetLineWidth(3);
      if(j%2)
	hist[j]->SetLineStyle(2);
      hist[j]->SetMarkerStyle(20);
      hist[j]->SetMarkerSize(0.9);
      hist[j]->SetMarkerColor(colors[j/2]);
      if(j==0)
	{
	  hist[j]->GetXaxis()->SetNdivisions(505);
	  hist[j]->GetYaxis()->SetNdivisions(505);
	  hist[j]->GetXaxis()->SetTitle(xTitle);
	  hist[j]->GetYaxis()->SetTitle("Relative ID cut efficiency");
	  hist[j]->GetYaxis()->SetTitleOffset(1.3);
	  if(xMin!=xMax)
	    hist[j]->GetXaxis()->SetRangeUser(xMin,xMax); //make this automatic
	  /*set y range automatically too*/
	  hist[j]->Draw();
	}
      else
	hist[j]->Draw("same");
      
      TString ll=legendNames[j];
      leg->AddEntry(hist[j],ll,"pl");
    }
  
  leg->Draw("same");
  
  c1->cd();
  TLatex T1;
  T1.SetTextSize(.04014599);
  T1.SetNDC();
  T1.DrawLatex(0.45,0.95,"CMS Simulation Preliminary"); //follows x-axis and y-axis values for position...
  /*
  T1.SetTextSize(.035);
  T1.DrawLatex(0.2,0.5,"Aging model sys uncert.: ~30%");
  T1.DrawLatex(0.2,0.45,"Run1 reconstruction");
  T1.DrawLatex(0.2,0.40,"Prompt photons");
  */
  char  ss1[50];
  sprintf(ss1,"PLOTS/%s.png",filename);
  std::cout << ss1 << std::endl;
  c1->SaveAs(ss1);
  sprintf(ss1,"PLOTS/C/%s.C",filename);
  c1->SaveAs(ss1);
  //c1->Clear();
}


void formatROCGraph(int logY, const char* legendTitle, string * legendNames, int * colors, std::vector <std::vector<float> >  &xVal, std::vector <std::vector<float> > &yVal, const char* fileName );
void formatROCGraph(int logY, const char* legendTitle, string * legendNames, int * colors, std::vector <std::vector<float> >  &xVal, std::vector <std::vector<float> > &yVal, const char* fileName )
//void formatROCGraph(int logY, const char* legendTitle, string * legendNames, int * colors, std::vector<float>  &xVal, std::vector<float> &yVal, const char* fileName );
//void formatROCGraph(int logY, const char* legendTitle, string * legendNames, int * colors, TGraph &graph[3], const char* fileName )
{

  setTDRStyle1(1,1,0,logY);//gridX/Y, logX/Y
  gStyle->TStyle::SetOptStat(0);
  gStyle->TStyle::SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
  TLegend* leg = new TLegend(0.22,0.20,0.55,0.45,legendTitle);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03679144);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  /*
  leg->SetLineColor(kBlue);
  leg->SetLineColorAlpha(kBlue,0.0);
  */
  
  TGraph *gr_res[3];
  for(int i=0;i<3;i++)
    {
      std::cout << xVal[i].size() << " y " << xVal[i][xVal[i].size()-1] << std::endl;
      gr_res[i] = new TGraph(xVal[1].size()-1,&xVal[i][0],&yVal[i][0]);
      gr_res[i]->GetY();
      gr_res[i]->SetLineColor(colors[i]);
      gr_res[i]->SetLineWidth(2);
      gr_res[i]->SetMarkerStyle(20);
      gr_res[i]->SetMarkerSize(0.9);
      gr_res[i]->SetMarkerColor(colors[i]);
      if(i==0)
	{
	  gr_res[i]->SetMinimum(0.0);
	  gr_res[i]->SetMaximum(1.0);
	  gr_res[i]->GetXaxis()->SetNdivisions(505);
	  gr_res[i]->GetXaxis()->SetRangeUser(0,1.0);
	  gr_res[i]->GetYaxis()->SetRangeUser(0.0,1.0);
	  gr_res[i]->GetYaxis()->SetNdivisions(505);
	  gr_res[i]->GetXaxis()->SetTitle("Signal Efficiency");
	  gr_res[i]->GetYaxis()->SetTitle("background rejection");
	  gr_res[i]->Draw("AP");
	}
      else
	 gr_res[i]->Draw("P");
      TString ll=legendNames[i];
      leg->AddEntry(gr_res[i],ll,"pl");
      leg->Draw("same");
    }
  /*
  TLatex T1;
  T1.SetNDC();
  T1.SetTextSize(.04014599);
  T1.DrawLatex(.45,0.95,"CMS Simulation Preliminary"); //follows x-axis and y-axis values for position...
  T1.SetTextSize(.035);
  T1.DrawLatex(0.2,0.3,"Aging model sys uncert.: ~30%");
  T1.DrawLatex(0.2,0.25,"Run1 reconstruction");
  T1.DrawLatex(0.2,0.20,"Prompt photons, pT>30GeV");
  */
  char ss1[50];
  sprintf(ss1,"PLOTS/%s.png",fileName);
  c1->SaveAs(ss1);
  sprintf(ss1,"PLOTS/C/%s.C",fileName);
  c1->SaveAs(ss1);
  //c1->Clear();
}
