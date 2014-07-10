#include "root.cc"
#include "params.h"
#include "formats.h"
using namespace std;

int PlotFromText(const char * filename)
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
      std::vector < std::vector<float> > data;
      //float data[nEta][nPoints];
      float data_val;
      std::string temp;
      for (int i=0; i<nEta; i++) 
	{
	  std::vector < float> row;
	  myfile >> temp;
	  //std::cout << temp << std::endl;
	  for (int j=0; j<nPoints; j++) 
	    {
	      myfile >> temp;
	      myfile >> temp;
	      //myfile >> data[i][j]; 
	      myfile >> data_val; 
	      row.push_back(data_val);
	      //  std::cout << data[i][j] << std::endl;
	    }
	  data.push_back(row);
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
      
      string legName[3] = {"<PU>: 50, int lumi: 0fb^{-1}","<PU>: 140, int lumi: 1000fb^{-1}","<PU>: 140, int lumi: 3000fb^{-1}"};
      int colors[3] = {1,2,4};
      formatEtaGraph(0,"energy (SC) resolution", legName, "Energy resolution, #sigma_{eff}(E)/E" ,colors,data,nPoints,filename);
      return(0);
    }
}
