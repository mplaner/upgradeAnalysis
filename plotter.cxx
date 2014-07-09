#include "root.cc"
#include "params.h"
#include "formats.h"
using namespace std;

//int main()
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
	  //std::cout << temp << std::endl;
	  for (int j=0; j<nPoints; j++) 
	    {
	      myfile >> temp;
	      myfile >> temp;
	      myfile >> data[i][j]; 
	      //  std::cout << data[i][j] << std::endl;
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
      
      string legName[1] = {"<PU>: 50, int lumi: 0fb^{-1}"};
      int colors[1] = {1};
      formatEtaGraph(0,"energy resolution", legName, "Energy resolution, #sigma_{eff}(E)/E" ,colors,&data[0][0],nPoints);
      return(0);
    }
}
