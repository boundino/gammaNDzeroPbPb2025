#include <TH3F.h>

#include "config.h"
#include "xjjanauti.h"
#include "xjjcuti.h"

#define __BINS_EQ__
#include "../include/bins.h"

#include "Variable.h"

int macro(std::string var, std::string configfile) {
  xjjc::config conf(configfile); std::cout<<std::endl;

  auto inputfiles = conf.get_vec("Inputs");
  auto cut = conf.get("Cut");
  auto outputname = conf.has("Output") ? conf.get("Output") : xjjc::str_gettag_from_file(configfile);

  auto idx = i_variables(var, globals::variables);
  if (idx < 0) {
    std::cout<<"error: "<<var<<" is not in variable list."<<std::endl;
    return 2;
  }
  const auto v = globals::variables[idx];

  auto* outf = xjjroot::newfile("rootfiles/" + var + "_" + outputname + ".root");

  for (const auto& i : inputfiles) {
    auto ipar = xjjc::str_divide_trim(i, ":");
    if (ipar.size() < 5) { std::cout<<"warning: skip bad input "<<i<<std::endl; continue; }
    auto* inf = TFile::Open(ipar[0].c_str());
    if (!inf) { std::cout<<"warning: skip bad file "<<ipar[0]<<"."<<std::endl; continue; }
    auto* tr = (TTree*)inf->Get("Tree");
    if (!tr) { std::cout<<"warning: skip bad file "<<ipar[0]<<"."<<std::endl; continue; }
    auto ip_name = ipar[1],
      add_cut = ipar[4], ip_cut = cut + " && " + add_cut;
    auto* h3 = new TH3F(Form("h3_%s", ip_name.c_str()), Form(";y;%s;p_{T} (GeV)", v.latex.c_str()),
                        bins::ny, bins::miny, bins::maxy,
                        v.nbin, v.minbin, v.maxbin,
                        bins::npt, bins::minpt, bins::maxpt);

    std::cout<<"\e[2mFilling\e[0m "<<ipar[0]<<std::endl
             <<"    \e[2mCut\e[0m "<<ip_cut<<std::endl
             <<"    \e[2mVar\e[0m "<<v.var<<std::endl;
    // tr->Project(h3->GetName(), Form("Dpt:%s:Dy", v.var.c_str()), cut.c_str(), "", 1.e7);
    tr->Project(h3->GetName(), Form("Dpt:%s:Dy", v.var.c_str()), ip_cut.c_str());
    outf->cd();
    xjjroot::writehist(h3);
  }

  outf->Close();

  return 0;
  
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
