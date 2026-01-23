#include <TH3F.h>

#include "config.h"
#include "xjjanauti.h"
#include "xjjcuti.h"

#define __BINS_DCA__
#include "../include/bins.h"

int macro(std::string configfile) {
  xjjc::config conf(configfile); std::cout<<std::endl;

  auto inputfiles = conf.get_vec("Inputs");
  auto inputfiles_prompt = conf.get_vec("Inputs_prompt");
  auto inputfiles_nonprompt = conf.get_vec("Inputs_nonprompt");
  auto selevt = conf.get("Sel_event");
  auto cut = conf.get("Cut");
  auto outputname = conf.has("Output") ? conf.get("Output") : xjjc::str_gettag_from_file(configfile);
  auto label_1 = conf.has("Label_1") ? conf.get("Label_1") : "";

  auto* tdata = xjjana::get_tree_multifiles(inputfiles, "Tree");
  auto* tmc_prompt = xjjana::get_tree_multifiles(inputfiles_prompt, "Tree");
  auto* tmc_nonprompt = xjjana::get_tree_multifiles(inputfiles_nonprompt, "Tree");

  outputname = "rootfiles/" + outputname + "_project.root";
  auto* outf = xjjroot::newfile(outputname);

#define PROJECT(l, tr)                                                  \
  std::cout<<std::endl<<"\e[2m" << #l << ":\e[0m "<<cut##l<<std::endl;  \
  auto* h3_##l = new TH3F("h3_" #l, "y;DCA (cm);p_{T} (GeV)", yinclbins.size()-1, yinclbins.data(), ip3dbins.size()-1, ip3dbins.data(), ptbins.size()-1, ptbins.data()); \
  tr->Project(h3_##l->GetName(), "Dpt:DsvpvDistance*sin(Dalpha):Dy", cut##l.c_str(), "", 1.e7); \
  xjjroot::writehist(h3_##l);
  
  auto cutmc_prompt = cut + " && Dgen==23333";
  PROJECT(mc_prompt, tmc_prompt);
  auto cutmc_nonprompt = cut + " && Dgen==23333";
  PROJECT(mc_nonprompt, tmc_nonprompt);

  // auto cutdata = cut + " && fabs(Dmass-1.8649)<0.025";
  // auto* hdata = new TH3F("hdata", "y;DCA (cm);p_{T} (GeV)", ip3dbins.size()-1, ip3dbins.data(), yinclbins.size()-1, yinclbins.data(), ptbins.size()-1, ptbins.data());
  // std::cout<<" cutdata: "<<cutdata<<std::endl;
  // t->Project(hdata->GetName(), "Dpt:DsvpvDistance*TMath::Sin(Dalpha):Dy", cutdata.c_str());
  // auto* hdata_sideband = new TH3F("hdata_sideband", "y;DCA (cm);p_{T} (GeV)", ip3dbins.size()-1, ip3dbins.data(), yinclbins.size()-1, yinclbins.data(), ptbins.size()-1, ptbins.data());
  // std::cout<<" cutdata_sideband: "<<cutdata_sideband<<std::endl;
  // t->Project(hdata_sideband->GetName(), "Dpt:DsvpvDistance*TMath::Sin(Dalpha):Dy", cutdata_sideband.c_str());

  outf->Close();

  return 0;
  
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}
