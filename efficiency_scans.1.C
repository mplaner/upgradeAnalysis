//#include "root.cc"
#include <iostream>
#include <fstream>
#include "params.h"
#include "formats.h"
#include <algorithm>
#include <string>


using namespace std;


bool mSort(const vector<float>& vec1, const vector<float>& vec2);
bool mSort(const vector<float>& vec1, const vector<float>& vec2)
{
  if(vec1[3]<vec2[3])
    return true;
  else if(vec1[3]>vec2[3])
    return false;
  else
    return vec1[4]<vec2[4];
  //return vec1[3]<vec2[3];
};

int PlotFromText(int etaBin, const char * filename)
{
  std::ifstream files[6];
  
  files[0].open("pile70signal.txt");
  files[1].open("pile70background.txt");
  files[2].open("pile140aging1000signal.txt");
  files[3].open("pile140aging1000background.txt");
  files[4].open("pile140aging3000signal.txt");
  files[5].open("pile140aging3000background.txt");
  
  for(int i =0; i<6;i++)
    if(files[i].is_open())
      std::cout << "plotting file " << i << " opened" << std::endl;
    else
      {
	std::cout << "cannot open txt file #" << i << std::endl;
	return(0);
      }
  std::vector< std::vector <float> >  xVAL;
  std::vector <std::vector <float> >  yVAL;
  for(int i=0;i<6;i+=2)
    {
      std::vector <float> xval;
      std::vector <float> yval;
      //while(!files[i].eof()) 
      //{
      string line;
      std::getline(files[i],line);
      std::cout << line << std::endl;
      std::getline(files[i+1],line);
      std::cout << line << std::endl;
      std::getline(files[i],line);
      std::cout << line << std::endl;
      std::getline(files[i+1],line);
      std::cout << line << std::endl;
      
      std::vector < std::vector<float> > dataEB1;
      std::vector < std::vector<float> > dataEE1;
      std::vector < std::vector<float> > dataEB2;
      std::vector < std::vector<float> > dataEE2;
      while(!files[i].eof()||!files[i+1].eof())
	{
	  float data_val;
	  std::string temp;
	  std::vector < float> rowEB1;
	  std::vector < float> rowEB2;
	  std::vector < float> rowEE1;
	  std::vector < float> rowEE2;
	  /*********EB ppart***********/
	  files[i] >> data_val;
	  files[i+1] >> data_val;
	  if(files[i].eof()||files[i+1].eof()) //deals with newline at the end of the file...
	    break;
	  //std::cout << "ID #1 : " << data_val;
	  rowEB1.push_back(data_val); //adds first ID paramater
	  rowEB2.push_back(data_val); //adds first ID paramater
	  files[i] >> data_val;
	  files[i+1] >> data_val;
	  rowEB1.push_back(data_val); //adds second ID parameter
	  rowEB2.push_back(data_val); //adds first ID paramater
	  
	  files[i] >> data_val;
	  files[i+1] >> data_val;
	  rowEB1.push_back(data_val); //adds third ID parameter
	  rowEB2.push_back(data_val); //adds first ID paramater
	  files[i] >> data_val;
	  data_val = floorf(data_val*100)/100.0;
	  if(etaBin==0)
	    xval.push_back(data_val);
	  rowEB1.push_back(data_val); //add first efficiency
	  files[i+1] >> data_val;
	  if(etaBin==0)
	    yval.push_back(1-data_val);
	  rowEB1.push_back(data_val); //add first background rejection
	  
	  files[i] >> data_val;
	  data_val = floorf(data_val*100)/100.0;
	  if(etaBin==1)
	    xval.push_back(data_val);
	  rowEB2.push_back(data_val); //add 2nd efficiency
	  files[i+1] >> data_val;
	  if(etaBin==1)
	    yval.push_back(1-data_val);
	  rowEB2.push_back(data_val); //add 2nd background rejection
	  
	  /******EE part*********/
	  files[i] >> data_val;
	  files[i+1] >> data_val;
	  rowEE1.push_back(data_val); //adds first ID paramater
	  rowEE2.push_back(data_val); //adds first ID paramater
	  files[i] >> data_val;
	  files[i+1] >> data_val;
	  rowEE1.push_back(data_val); //adds second ID parameter
	  rowEE2.push_back(data_val); //adds first ID paramater
	  files[i] >> data_val;
	  files[i+1] >> data_val;
	  rowEE1.push_back(data_val); //adds third ID parameter
	  rowEE2.push_back(data_val); //adds first ID paramater
	  
	  files[i] >> data_val;
	  data_val = floorf(data_val*100)/100;
	  if(etaBin==2)
	    xval.push_back(data_val);
	  rowEE1.push_back(data_val); //add 3rd efficiency
	  files[i+1] >> data_val;
	  if(etaBin==2)
	    yval.push_back(1-data_val);
	  rowEE1.push_back(data_val); //add 3rd background rejection
	  
	  files[i] >> data_val;
	  data_val = floorf(data_val*100)/100.0;
	  if(etaBin==3)
	    xval.push_back(data_val);
	  rowEE2.push_back(data_val); //add 4th efficiency
	  files[i+1] >> data_val;
	  if(etaBin==3)
	    yval.push_back(1-data_val);
	  rowEE2.push_back(data_val); //add 4th background rejection
	    
	  
	  dataEB1.push_back(rowEB1); //adds line to vector for sorting
	  dataEB2.push_back(rowEB2);
	  dataEE1.push_back(rowEE1); //adds line to vector for sorting
	  dataEE2.push_back(rowEE2);
	  
	}
      std::cout << "starting to sort" << std::endl;
      std::sort(dataEB1.begin(),dataEB1.end(),mSort);
      std::sort(dataEB2.begin(),dataEB2.end(),mSort);
      std::sort(dataEE1.begin(),dataEE1.end(),mSort);
      std::sort(dataEE2.begin(),dataEE2.end(),mSort);
      //std::cout << "finished file # " << i << std::endl;
      //std::cout.width(5);
      std::cout.precision(3);
      std::ofstream ofile;
      if(i==0)
	ofile.open("PU70age0_rocs.txt", std::ofstream::out);
      else 
	if(i==2)
	  ofile.open("PU140age1k_rocs.txt", std::ofstream::out);
      else
	ofile.open("PU140age3k_rocs.txt", std::ofstream::out);
      for(int k=0; k< dataEB1.size(); k++)  //change to file out once it's clear everything's working
	{
	  for(int j =0; j<5;j++)
	    {
	      //std::cout << dataEB1[k][j] << " ";
	      ofile << dataEB1[k][j] << " ";
	    }
	  
	  for(int j =0; j<5;j++)
	    {
	      //std::cout << dataEB2[k][j] << " ";
	      ofile << dataEB2[k][j] << " ";
	    }
	  
	  for(int j =0; j<5;j++)
	    {
	      //std::cout << dataEE1[k][j] << " ";
	      ofile << dataEE1[k][j] << " ";
	    }
	  
	  for(int j =0; j<5;j++)
	    {
	      //std::cout << dataEE2[k][j] << " ";
	      ofile << dataEE2[k][j] << " ";
	    }
	  //std::cout  << std::endl;
	  ofile  << std::endl;
	}
      std::cout << "x1 " << xval[0] << " y1 " << yval[0] << std::endl;
      xVAL.push_back(xval);
      yVAL.push_back(yval);
      ofile.close();
      files[i].close();
      files[i+1].close();
    }
  //finished running over all 6 files
  string legName[3] = {"<PU>: 70, int lumi: 0fb^{-1}","<PU>: 140, int lumi: 1000fb^{-1}","<PU>: 140, int lumi: 3000fb^{-1}"};
  int colors[3] = {1,2,4};
  
  if(etaBin==0)
    formatROCGraph(0,"|#eta|<0.5", legName,colors,xVAL,yVAL,filename);
  if(etaBin==1)
    formatROCGraph(0,"0.5<|#eta|<.1.44", legName,colors,xVAL,yVAL,filename);
  if(etaBin==2)
    formatROCGraph(0,"1.566<|#eta|<2.0", legName,colors,xVAL,yVAL,filename);
  if(etaBin==3)
    formatROCGraph(0,"2.0<#eta<2.5", legName,colors,xVAL,yVAL,filename);

  return(0);
}

