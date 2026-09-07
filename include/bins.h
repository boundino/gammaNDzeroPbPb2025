#ifndef __BINS__
#define __BINS__

namespace bins {

#ifdef __BINS_PTY_PLACEHOLDER__
  std::vector<double> ybins;
  std::vector<double> ptbins;
#endif

#ifdef __BINS_PTY_BDT__
  const std::vector<double> ybins = { -2., -1., 0., 1., 2. };
  const std::vector<double> ptbins = { 2., 5. };
#endif

#ifdef __BINS_PTY_EFF__
  const int ny = 48; const float miny = -2.4, maxy = 2.4;
  const int npt = 30; const float minpt = 2, maxpt = 5;
#endif

#ifdef __BINS_PTY_DATAMCCOMP__
  std::vector<double> ybins = { -2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2. }; // can be redefined in .cc
  const int npt = 1; const float minpt = 2, maxpt = 5;
#endif

#ifdef __BINS_MULT__
  const int nmult = 50; const float minmult = 0, maxmult = 50;
#endif

#ifdef __BINS_MASS__
  const int nmass = 80; const float minmass = 1.66, maxmass = 2.06;
#endif

  // used in comparison/savehist.cc -> to update
#ifdef __BINS_PTY_EQ__
  const int ny = 4; const float miny = -2, maxy = 2;
  const int npt = 1; const float minpt = 2, maxpt = 5;
#endif
  
}

#endif
