#include <vector>
#include <algorithm>
#include <iostream>


double sigmaEff(std::vector<double> v, float threshold,float& xmin, float& xmax)
{
  
  std::sort(v.begin(),v.end());
  
  //  threshold=threshold/2.0;
  
  int total = v.size();
  int max = (int)(threshold * total);

  std::cout << "vector size " << v.size() << std::endl;
  
  std::vector<double>  start;
  std::vector<double>  stop;
  std::vector<double>  width;
    
  unsigned i = 0 ;
  while (i != v.size()-1){
    
    unsigned j = i+max;
    if(j>=v.size())
      break;
    if ( j != v.size()-1)
      {
	start.push_back(v[i]);
	stop.push_back(v[j]);
	width.push_back(v[j] - v[i]);
      }
    ++i;
    //    std::cout << "first " << v[i] << " width " << v[j]-v[i] <<std::endl;
  }
  
  double minwidth = *min_element(width.begin(), width.end());
  
  unsigned pos = min_element(width.begin(),width.end()) - width.begin();    
  
  xmin = start[pos];
  xmax = stop [pos];
  return minwidth;

}
