#pragma once

namespace global {
  struct variable
  {
    std::string varname = "";
    std::string var = "";
    std::string vartex = "";
    float varmin = 0;
    float varmax = 1;
    int nbin = 1;
    int logy = 0;
    int isbranch = 1;
    std::vector<double> bins = {};
    std::string note = "";
  };
}

