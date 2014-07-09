#include "root.cc"
#include "params.h"

void formatHisto(int logY, char* legendTitle, string * legendNames, char* xTitle, float xMin, float xMax, char* yTitle, int *colors, TH1F** hist, const int nPoints);
void formatHisto(int logY, char* legendTitle, string * legendNames, char* xTitle, float xMin, float xMax,char* yTitle, int *colors, TH1F** hist, const int nPoints)
{
  setTDRStyle1(1,1,0,logY);//gridX/Y, logX/Y
  gStyle->TStyle::SetOptStat(0);
  gStyle->TStyle::SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
  TLegend* leg = new TLegend(0.20,0.65,0.45,0.90,legendTitle);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03679144);
  leg->SetFillColor(0);
  
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
  T1.DrawLatex(0.8,0.43,"CMS Simulation Preliminary"); //follows x-axis and y-axis values for position...
  T1.SetTextSize(.035);
  T1.DrawLatex(0.5,0.18,"Aging model sys uncert.: ~30%");
  T1.DrawLatex(0.5,0.15,"Run1 reconstruction");
  T1.DrawLatex(0.5,0.12,"Prompt photons");
  
  string ss1="PLOTS/histoTest.png";
  c1->SaveAs(ss1.c_str());
  ss1="PLOTS/C/histoTest.C";
  c1->SaveAs(ss1.c_str());
  c1->Clear();
}

//"Energy resolution, #sigma_{eff}(E)/E"
/*used to format any graph which is produced with nEta points on the x-axis */
void formatEtaGraph(int logY, char* legendTitle, string * legendNames, char * yTitle, int * colors, float *data, const int nPoints );
void formatEtaGraph(int logY, char* legendTitle, string * legendNames, char * yTitle, int * colors, float *data, const int nPoints )
{

  setTDRStyle1(1,1,0,logY);//gridX/Y, logX/Y
  gStyle->TStyle::SetOptStat(0);
  gStyle->TStyle::SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
  TLegend* leg = new TLegend(0.25,0.65,0.55,0.90,legendTitle);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03679144);
  leg->SetFillColor(0);

  TGraph *gr_res[nPoints];
  for (int j=0; j<nPoints; j++)
    {
      float y[nEta];
      for (int i=0; i<nEta; i++)
	{
	  y[i]=data[i][j];
	  //  std::cout << "data: " << y[i] << std::endl;
	}
      gr_res[j] = new TGraph(nEta,etaVal,y);
      gr_res[j]->GetY();
      gr_res[j]->SetLineColor(colors[j]);
      gr_res[j]->SetLineWidth(3);
      gr_res[j]->SetMarkerStyle(20);
      gr_res[j]->SetMarkerSize(0.9);
      gr_res[j]->SetMarkerColor(colors[j]);
      if (j==0)
	{
	  gr_res[j]->SetMinimum(0.0);
	  gr_res[j]->SetMaximum(0.25);
	  gr_res[j]->GetXaxis()->SetNdivisions(505);
	  gr_res[j]->GetYaxis()->SetNdivisions(505);
	  gr_res[j]->GetXaxis()->SetTitle("#eta");
	  gr_res[j]->GetYaxis()->SetTitle(yTitle);
	  gr_res[j]->Draw("ALP");
	}
      else
	{
	  gr_res[j]->Draw("LP");
	}
      TString ll=legendNames[j];
      leg->AddEntry(gr_res[j],ll,"pl");
    }
  leg->Draw("same");
  TLatex * latex2 = new TLatex(1.0,0.25,"CMS Simulation Preliminary");
  latex2->SetTextSize(0.04014599);
  latex2->Draw("same");
  TLatex * latex2 = new TLatex(0.2,0.08,"Aging model sys uncert.: ~30%");
  latex2->SetTextSize(0.035);
  latex2->Draw("same");
  TLatex * latex2 = new TLatex(0.2,0.07,"Run1 reconstruction");
  latex2->SetTextSize(0.035);
  latex2->Draw("same");
  TLatex * latex2 = new TLatex(0.2,0.06,"Prompt photons");
  latex2->SetTextSize(0.035);
  latex2->Draw("same");
  string ss1="PLOTS/graph_test.png";
  c1->SaveAs(ss1.c_str());
  ss1="PLOTS/C/graph_test.C";
  c1->SaveAs(ss1.c_str());
  c1->Clear();
}
