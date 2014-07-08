#include "root.cc"
#include "params.h"
using namespace std;

int main()
{
  const int plot_resolution_michael=1;
  
  if (plot_resolution_michael)
    {
      ifstream myfile;
      myfile.open("plots.txt");
      if(myfile.is_open())
	std::cout << "plotting file opened" << std::endl;
      
      int not_present=1;
      while(!myfile.eof()) 
	{
	  string line;
	  std::getline(myfile,line);
	  if(!line.compare("starting energy resolution"))
	    {
	      std::getline(myfile,line);
	      std::cout << line << std::endl;
	      not_present=0;
	      break;
	    }
	}
      if(not_present)
	return(0);
      const int nPoints=1;
      std::string points[nPoints];
      //for (int i=0; i<nPoints; i++) 
      //	myfile >> points[i];
      float data[nEta][nPoints];
      std::string temp;
      for (int i=0; i<nEta; i++) 
	{
	  myfile >> temp;
	  std::cout << temp << std::endl;
	  for (int j=0; j<nPoints; j++) 
	    {
	      myfile >> temp;
	      myfile >> temp;
	      myfile >> data[i][j]; 
	      std::cout << data[i][j] << std::endl;
	    }
	}
      myfile.close();
      // CHECK:
      for (int i=0; i<nEta; i++)
	{
	  cout << etaVal[i];
	  for (int j=0; j<nPoints; j++) 
	    {
	      cout << "\t" << data[i][j];
	    }
	  cout << endl;
	}
      int logOY=0;
      setTDRStyle1(1,1,0,logOY);//gridX/Y, logX/Y
      gStyle->TStyle::SetOptStat(0);
      gStyle->TStyle::SetOptFit(0);
      TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
      TLegend* leg = new TLegend(0.25,0.65,0.55,0.90,"Integrated lumi  nPU");//"");
      leg->SetTextFont(62);
      leg->SetTextSize(0.03679144);
      leg->SetFillColor(0);
      string legendNames[nPoints] = {"test"};
      int color[nPoints]={1};
      TGraph *gr_res[nPoints];// = new TGraph(n,xx,yy4);
      for (int j=0; j<nPoints; j++)
	{
	  float y[nEta];
	  for (int i=0; i<nEta; i++) 
	    {
	      y[i]=data[i][j];
	    }
	  gr_res[j] = new TGraph(nEta,etaVal,y);
	  gr_res[j]->GetY();
	  gr_res[j]->SetLineColor(color[j]);
	  gr_res[j]->SetLineWidth(3);
	  gr_res[j]->SetMarkerStyle(20);
	  gr_res[j]->SetMarkerSize(0.9);
	  gr_res[j]->SetMarkerColor(color[j]);
	  
	  if (j==0)
	    {
	      gr_res[j]->SetMinimum(0.0);
	      gr_res[j]->SetMaximum(0.25);
	      gr_res[j]->GetXaxis()->SetNdivisions(505);
	      gr_res[j]->GetYaxis()->SetNdivisions(505);
	      gr_res[j]->GetXaxis()->SetTitle("#eta");
	      gr_res[j]->GetYaxis()->SetTitle("Energy resolution, #sigma_{eff}(E)/E");
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
      {
	TLatex * latex2 = new TLatex(1.0,0.25,"CMS Simulation Preliminary");
	latex2->SetTextSize(0.04014599);
	latex2->Draw("same");
      }
      {
	TLatex * latex2 = new TLatex(0.2,0.08,"Aging model sys uncert.: ~30%");
	latex2->SetTextSize(0.035);
	latex2->Draw("same");
      }
      {
	TLatex * latex2 = new TLatex(0.2,0.07,"Run1 reconstruction");
	latex2->SetTextSize(0.035);
	latex2->Draw("same");
      }
      {
	TLatex * latex2 = new TLatex(0.2,0.06,"Prompt photons");
	latex2->SetTextSize(0.035);
	latex2->Draw("same");
      }
      string ss1="PLOTS/gr_photon_energy_resolution_vs_eta_for_michael.png";
      c1->SaveAs(ss1.c_str());
      ss1="PLOTS/C/gr_photon_energy_resolution_vs_eta_for_michael.C";
      c1->SaveAs(ss1.c_str());
      c1->Clear();
    }
}

