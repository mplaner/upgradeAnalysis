const double EEMIN= 1.566,EBMAX=1.444; //cuts for gap between EE and EB 
const double EEMAX= 2.5, EBMIN=0; //edge of detector (can be used to limit eta range for comparison plots)
const double drCut= .15; //cut on deltaR (optimized for EB pu140 cut =.21) 
const double g_ptCut = 30; //cut on GEN level pt 
const double p_ptCut = 20; //cut on RECO level pt
const double p_r9Cut =0;   //cut on RECO level r9  
const double p_sieieCut =1; //cut on RECO level sIetaIeta
const double p_hoverE   =1; //cut on RECO level H/E     
const int HGG    =  1;
const int GJETS  =  3;
const int EE     =  2;
const int EB     =  1;

const int type = HGG;
const int nEta = 6;
const double etaMin[nEta] = {0.000, 0.500, 1.000, 1.566, 1.870, 2.180};
const double etaMax[nEta] = {0.500, 1.000, 1.444, 1.870, 2.180, 2.500};
const double etaVal[nEta] = {0.00,  0.75,  1.222, 1.718, 2.025, 2.340};
const int  etaColor[nEta] = {1,2,4,6,8,9};

int GetEtaBin(double ETA);
int GetEtaBin(double ETA)
{
  for(int i=0; i<nEta;i++)
    {
      if(ETA>etaMin[i]&&ETA<etaMax[i])
	return(i);
    }
  return(-1);
  
}
float dR(double eta1, double phi1, double eta2, double phi2);
float dR(double eta1, double phi1, double eta2, double phi2)
{
  double pi = 3.141596;
  //check if phi complete rotation of phi is needed
  if (std::fabs(phi1-phi2) > pi)
    //calculate dr after rotating phi by 2pi
    return sqrt((eta1-eta2)*(eta1-eta2) + pow((2.0*pi - std::fabs(phi1-phi2)),2.0));
  //calculate dr
  return sqrt((eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
}

float GetEffSigma(float interval, vector<float> v, float &xmin, float &xmax);
float GetEffSigma(float interval, vector<float> v, float &xmin, float &xmax)
{
  //sort vector low to high
  std::sort(v.begin(),v.end());
  //get number of points
  int total = v.size();
  //get number of points in interval
  int max = (int)(interval*total);
  //initialize width to first point
  float width = v[max]-v[0];
  //search for smallest xinterval containing 'interval' percent of points
  for(int i=1; i<(total-max);i++)
    {
      //if new width is larger, continue
      if(v[i+max]-v[i]> width)
	continue;
      //if new width is smaller set xmin and width
      xmin = v[i];
      width = v[i+max]-xmin;
    }
  //set xmax after the minimization is finished to save time
  xmax = xmin+width;
  //return effective sigma
  return(width/2);  
}
