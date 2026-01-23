#include <iostream>
#include <cmath>

int main() {
  const int nBinY = 14;
  //const int nBinY = 20;
  float binsY[nBinY+1];
  float firstBinYWidth = 0.001;
  float binYWidthRatio = 1.27;

  float v = 0; 
  std::cout<<v<<", ";
  for (int i=1; i<=nBinY; i++) {
    v = v + firstBinYWidth*pow(binYWidthRatio,i-1);
    std::cout<<v<<", ";    
  }
  std::cout<<std::endl;
}
