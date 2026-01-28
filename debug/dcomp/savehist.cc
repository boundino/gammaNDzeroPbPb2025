#include <TH3F.h>

#include "config.h"
#include "xjjanauti.h"
#include "xjjcuti.h"

#define __BINS_EQ__
#include "../include/bins.h"

#include "Variable.h"

int macro(std::string configfile, std::string var) {
  std::cout<<std::endl;
  xjjc::config conf(configfile);

  auto inputfiles = conf.get_vec("Inputs");
  auto cut = conf.get("Cut");
  auto name = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);

  auto idx = i_variables(var, globals::variables);
  if (idx < 0) {
    std::cout<<"error: "<<var<<" is not in variable list."<<std::endl;
    return 2;
  }
  const auto v = globals::variables[idx];

  auto* tr = xjjana::get_tree_multifiles(inputfiles, "Tree");
  if (!tr) { std::cout<<"error: bad input file "<<conf.get("Inputs")<<"."<<std::endl; return 2; }

  auto* outf = xjjroot::newfile("rootfiles/" + var + "_" + name + ".root");
  auto* h3 = new TH3F("h3", Form(";y;%s;p_{T} (GeV)", v.latex.c_str()),
                      bins::ny, bins::miny, bins::maxy,
                      v.nbin, v.minbin, v.maxbin,
                      bins::npt, bins::minpt, bins::maxpt);

  std::cout<<std::endl<<"\e[2mFilling\e[0m "<<name<<std::endl
           <<"    \e[2mCut\e[0m "<<cut<<std::endl
           <<"    \e[2mVar\e[0m "<<v.var<<std::endl;

  tr->Project(h3->GetName(), Form("Dpt:%s:Dy", v.var.c_str()), cut.c_str());
  xjjroot::writehist(h3);

  outf->Close();

  return 0;
  
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
