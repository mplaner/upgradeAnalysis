#ifndef params
#define params
#include <math.h>
#include <vector>
#include <iostream>

const double EEMIN= 1.566,EBMAX=1.444; //cuts for gap between EE and EB 
const double EEMAX= 3.0, EBMIN=0; //edge of detector (can be used to limit eta range for comparison plots)
const double drCut= .15; //cut on deltaR 
const double g_ptCut = 30; //cut on GEN level pt 
const double g_eCut  = 0;
const double p_ptCut = 10; //cut on RECO level pt
const double p_r9Cut =0;   //cut on RECO level r9  
const double p_sieieCut =1; //cut on RECO level sIetaIeta
const double p_hoverE   =1; //cut on RECO level H/E     
const int HGG    =  1;
const int GJETS  =  3;
const int EE     =  2;
const int EB     =  1;
const int FILETYPE   =  1;

const int Type = GJETS;
const int nEta = 2;
const float etaMin[nEta] = {0.000, 1.566};
const float etaMax[nEta] = {1.440, 2.500};
const float etaVal[nEta] = {0.000, 2.033};
const float etaWidthHigh[nEta] = {1.44, 0.467};
const float etaWidthLow[nEta] =  {0.00, 0.467};
//const int  etaColor[nEta] = {1,2};
const int  etaColor[7] = {1,2,4,6,8,9,41};
/*
const float etaMin[nEta] = {0.000, 0.500, 1.000, 1.566, 1.870, 2.180, 2.50};
const float etaMax[nEta] = {0.500, 1.000, 1.444, 1.870, 2.180, 2.500, 3.00};
const float etaVal[nEta] = {0.00,  0.75,  1.222, 1.718, 2.025, 2.340, 2.75};
const float etaWidthHigh[nEta] = {0.5, 0.250, 0.222, 0.152, 0.155, 0.160, 0.250};
const float etaWidthLow[nEta] =  {0.0, 0.250, 0.222, 0.152, 0.155, 0.160, 0.250};
const int  etaColor[nEta] = {1,2,4,6,8,9,41};
*/

int GetEtaBin(double ETA);
int GetEtaBin(double ETA)
{
  for(int i=0; i<nEta;i++)
    {
      if(fabs(ETA)>etaMin[i]&&fabs(ETA)<etaMax[i])
	return(i);
    }
  return(-1);
  
}
float dR(double eta1, double phi1, double eta2, double phi2);
float dR(double eta1, double phi1, double eta2, double phi2)
{
  double pi = 3.141596;
  //check if phi complete rotation of phi is needed
  if (fabs(phi1-phi2) > pi)
    //calculate dr after rotating phi by 2pi
    return sqrt((eta1-eta2)*(eta1-eta2) + pow((2.0*pi - fabs(phi1-phi2)),2.0));
  //calculate dr
  return sqrt((eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
}

float GetEffSigmaErrHigh(int max, std::vector<float> v, float iMin, float xWidth);
float GetEffSigmaErrHigh(int max, std::vector<float> v, float iMin, float xWidth)
{
  float percent=.05;
  int total = v.size();
  int variation = (int)(percent*total);
  float high_error =0;
  if(iMin+max+variation< total)
    high_error = v[iMin + max + variation] - v[iMin + variation];
  else
    high_error = v[total-1] - v[iMin + variation];
  high_error = (high_error - xWidth)/2;
  return(xWidth/sqrt((double) variation));
  //return(high_error);
}

float GetEffSigmaErrLow(int max, std::vector<float> v, float iMin, float xWidth);
float GetEffSigmaErrLow(int max, std::vector<float> v, float iMin, float xWidth)
{
  float percent=.05;
  int total = v.size();
  int variation = (int)(percent*total);
  float low_error  =0;
  low_error  = v[iMin + max - variation] - v[iMin - variation];
  low_error = (low_error - xWidth)/2;
  return(low_error);
}
int GetPtWeightBin(float pt);
int GetPtWeightBin(float pt)
{
  int bin = pt +1;
  if(bin>0 && bin<100)
    return(bin);
  return(0); //bin of 1 gives a weight of 1, so safety catch
}

int GetEtaWeightBin(double eta);
int GetEtaWeightBin(double eta)
{
  int bin = 100.0*(eta+3.2)/6.4;
  if(bin>0&& bin<100)
    return(bin);
  return(0); //bin of 1 gives a weight of 1, so safety catch
}

float GetEffSigma(float interval, std::vector<float> v, float &xmin, float &xmax, float &highError, float &lowError);
float GetEffSigma(float interval, std::vector<float> v, float &xmin, float &xmax,  float &highError, float &lowError)
{
  //sort vector low to high
  std::sort(v.begin(),v.end());
  //get number of points
  int total = v.size();
  //get number of points in interval
  int max = (int)(interval*total);
  //initialize width to first point
  float width = v[max]-v[0];
  int imin = 0;
  //search for smallest xinterval containing 'interval' percent of points
  for(int i=1; i<(total-max)-1;i++)
    {
      //if new width is larger, continue
      if(v[i+max]-v[i]> width)
	continue;
      //if new width is smaller set xmin and width
      xmin = v[i];
      imin = i;
      width = v[i+max]-xmin;
    }
  //set xmax after the minimization is finished to save time
  xmax = xmin+width;
  //float widthError = 0;
  highError = GetEffSigmaErrHigh(max, v, imin, width); 
  lowError = GetEffSigmaErrLow(max, v, imin, width);
  //std::cout << "width: " << width/2 << std::endl;
  //return effective sigma
  return(width/2);  
}



#endif
