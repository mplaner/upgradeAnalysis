#include "root.cc"
#include "params.h"
#include "formats.h"
using namespace std;

int PlotFromText(const char * filename, const int CLUSTERTYPE)
{
  const int plot_resolution_michael=1;
  //const int CLUSTERTYPE=2;
  if (plot_resolution_michael)
    {
      ifstream myfile;
      if(CLUSTERTYPE==2)
	myfile.open("full_range_sc_error.txt");
      else if(CLUSTERTYPE==1)
	myfile.open("full_range_3x3_error.txt");
      else
	myfile.open("full_range_5x5_error.txt");
      
      if(myfile.is_open())
	std::cout << "plotting file opened" << std::endl;
      else
	{
	  std::cout << "cannot open txt file" << std::endl;
	  return();
	}
      
      int not_present=1;
      /*
      while(!myfile.eof()) 
	{
	  string line;
	  std::getline(myfile,line);
	  if(!line.compare("starting energy resolution"))
	    {
	      //std::getline(myfile,line);
	      std::cout << line << std::endl;
	      not_present=0;
	      break;
	    }
	}
      if(not_present)
	return(0);
      */
      const int nPoints=4;
      std::string points[nPoints];
      //for (int i=0; i<nPoints; i++) 
      //	myfile >> points[i];
      std::vector < std::vector<float> > data;
      std::vector < std::vector<float> > lError;
      std::vector < std::vector<float> > hError;
      //float data[nEta][nPoints];
      float data_val;
      float high_err;
      float low_err;
      std::string temp;
      for (int j=0; j<nPoints; j++) 
	{
	  std::vector < float> row;
	  std::vector < float> hrow;
	  std::vector < float> lrow;
	  string line;
	  std::getline(myfile, line);
	  std::cout << "scenario: " << line << std::endl;
	  for (int i=0; i<nEta; i++) 
	    {
	      myfile >> temp;
	      myfile >> temp;
	      myfile >> temp;
	      //myfile >> data[i][j]; 
	      myfile >> data_val;
	      myfile >> low_err;
	      myfile >> high_err;

	      std::cout << data_val << " ";
	      hrow.push_back(high_err);
	      lrow.push_back(low_err);
	      row.push_back(data_val);
	      //  std::cout << data[i][j] << std::endl;
	    }
	  std::cout << std::endl;
	  std::getline(myfile, line);
	  data.push_back(row);
	  lError.push_back(lrow);
	  hError.push_back(hrow);
	}
      myfile.close();
      // CHECK:
      for (int i=0; i<nEta; i++)
	{
	  cout << etaVal[i];
	  for (int j=0; j<nPoints; j++) 
	    {
	      cout << "\t" << data[j][i];
	    }
	  cout << endl;
	}
      
      string legName[4] = {"Run1 MC","<PU>: 70, int lumi: 0fb^{-1}","<PU>: 140, int lumi: 1000fb^{-1}","<PU>: 140, int lumi: 3000fb^{-1}"};
      int colors[4] = {1,2,4,8};
      if(CLUSTERTYPE==2)
	formatEtaGraph(0,"energy (SC) resolution", legName, "Energy resolution, #sigma_{eff}(E)/E" ,colors,data,lError, lError, nPoints,filename);
      else if(CLUSTERTYPE==1)
	formatEtaGraph(0,"energy (3x3) resolution", legName, "Energy resolution, #sigma_{eff}(E)/E" ,colors,data,lError, lError, nPoints,filename);
      else
	formatEtaGraph(0,"energy (5x5) resolution", legName, "Energy resolution, #sigma_{eff}(E)/E" ,colors,data,lError, lError, nPoints,filename);
      return(0);
    }
}



int PlotClusterTypes(const char * filename, const int RUNTYPE)
{
  ifstream myfile;
  if(RUNTYPE==1)
    myfile.open("run_1_all_clusters.txt");
  else
    myfile.open("PU140_age3k_all_clusters.txt");
      
      if(myfile.is_open())
	std::cout << "plotting file opened" << std::endl;
      else
	{
	  std::cout << "cannot open txt file" << std::endl;
	  return();
	}
      
      int not_present=1;
      while(!myfile.eof()) 
	{
	  string line;
	  std::getline(myfile,line);
	  if(!line.compare("starting energy resolution"))
	    {
	      //std::getline(myfile,line);
	      std::cout << line << std::endl;
	      not_present=0;
	      break;
	    }
	}
      if(not_present)
	return(0);
      const int nPoints=3;
      std::string points[nPoints];
      //for (int i=0; i<nPoints; i++) 
      //	myfile >> points[i];
      std::vector < std::vector<float> > data;
      std::vector < std::vector<float> > lError;
      std::vector < std::vector<float> > hError;
      //float data[nEta][nPoints];
      float data_val;
      float high_err;
      float low_err;
      std::string temp;
      for (int j=0; j<nPoints; j++) 
	{
	  std::vector < float> row;
	  std::vector < float> hrow;
	  std::vector < float> lrow;
	  string line;
	  std::getline(myfile, line);
	  std::cout << "scenario: " << line << std::endl;
	  for (int i=0; i<nEta; i++) 
	    {
	      myfile >> temp;
	      myfile >> temp;
	      myfile >> temp;
	      //myfile >> data[i][j]; 
	      myfile >> data_val;
	      myfile >> low_err;
	      myfile >> high_err;

	      std::cout << data_val << " ";
	      hrow.push_back(high_err);
	      lrow.push_back(low_err);
	      row.push_back(data_val);
	      //  std::cout << data[i][j] << std::endl;
	    }
	  std::cout << std::endl;
	  std::getline(myfile, line);
	  data.push_back(row);
	  lError.push_back(lrow);
	  hError.push_back(hrow);
	}
      myfile.close();
      // CHECK:
      for (int i=0; i<nEta; i++)
	{
	  cout << etaVal[i];
	  for (int j=0; j<nPoints; j++) 
	    {
	      cout << "\t" << data[j][i];
	    }
	  cout << endl;
	}
      
      string legName[4] = {"e3x3","e5x5","SC"};
      int colors[4] = {1,2,4,8};
      if(RUNTYPE==1)
	formatEtaGraph(0,"energy resolution Run 1", legName, "Energy resolution, #sigma_{eff}(E)/E" ,colors,data,lError, lError, nPoints,filename);
      else
	formatEtaGraph(0,"energy resolution <PU>: 140, int lumi: 3ab^{-1}", legName, "Energy resolution, #sigma_{eff}(E)/E" ,colors,data,lError, lError, nPoints,filename);
      return(0);
}
