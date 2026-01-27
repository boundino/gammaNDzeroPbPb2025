#ifndef __BINS__
#define __BINS__

namespace bins {
  const std::vector<float> ybins = { -2., -1., 0., 1., 2. },
    yinclbins = { -1., 1. },
    ptbins = { 2., 3., 4., 5. };

#ifdef __BINS_EQ__
  const int ny = 4; const float miny = -2, maxy = 2;
  const int npt = 1; const float minpt = 2, maxpt = 5;
#endif

#ifdef __BINS_DCA__
  const std::vector<float> ip3dbins = { 0, 0.001, 0.00227, 0.0038829, 0.00593128, 0.00853273, 0.0118366, 0.0160324, 0.0213612, 0.0281287, 0.0367235, 0.0476388, 0.0615013, 0.0791066, 0.101465 };
  // const std::vector<float> ip3dsigbins = {  };
#endif

}
#endif
