#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"

#define __BINS_PTY_EQ__
#define __BINS_MASS__
#include "../../include/bins.h"

#include "mvaroot.h"
#include "save_calc.h"

int macro(std::string configfile_train, std::string configfile_cut) {
  std::cout<<std::endl;
  xjjc::config conf_t(configfile_train), conf_c(configfile_cut);

  // train
  auto dir = conf_t.has("Name") ? conf_t.get("Name") : xjjc::str_tag_from_file(configfile_train);
  auto* rmva = new mytmva::mvaroot(conf_t.get("Rootoutput"));
  if (!rmva->good()) {
    __XJJLOG << "!! bad input " << rmva->info("rootfile") << std::endl;
    return 2;
  }
  // rmva->print();
  auto cuts = rmva->info("cuts"), cutb = rmva->info("cutb");
  auto inputSname = conf_t.get("Inputs"), inputBname = conf_t.get("Inputb");    
  // cut
  auto cut = conf_c.get("Cut");
  auto name = conf_c.has("Name") ? conf_c.get("Name") : xjjc::str_tag_from_file(configfile_cut);
  std::cout<<"\e[1m"<<dir<<" / "<<name<<"\e[0m"<<std::endl;

  auto* trS = xjjana::chain_files({ inputSname }, "Tree");
  if (!trS) { __XJJLOG<<"!! bad input file "<<inputSname<<std::endl; return 2; }
  save::mask_branch(trS);
  auto* trB = xjjana::chain_files({ inputBname }, "Tree");
  if (!trB) { __XJJLOG<<"!! bad input file "<<inputBname<<std::endl; return 2; }
  save::mask_branch(trB);

  auto* outf = xjjroot::newfile("rootfiles/" + dir + "/" + name + ".root");
  std::map<std::string, TH3F*> h3;
  auto project = [&](TChain* tr, std::string key, std::string icut) {
    h3[key] = new TH3F(Form("h3_cut%s", key.c_str()), ";y;m_{K#pi} (GeV);Pass",
                       bins::ny, bins::miny, bins::maxy,
                       bins::nmass, bins::minmass, bins::maxmass,
                       2, 0, 2);
    __XJJLOG << ">> "<<h3[key]->GetName()<<" [ "<<cut<<" ] \e[2m"<<icut<<"\e[0m"<<std::endl;
    tr->Project(h3[key]->GetName(), Form("%s:Dmass:Dy", cut.c_str()), icut.c_str());
    xjjroot::writehist(h3[key]);
    const auto eff = calc::eff_cut(h3[key]);
    __XJJLOG << ">> eff: [ "<<eff<<" ] rej: [ "<<1-eff<<" ]" << std::endl;;
  };
  project(trS, "_S", cuts);
  project(trB, "_B", cutb);

  auto* t = new TTree("info", "");
  t->Branch("cuts", &cuts);
  t->Branch("cutb", &cutb);
  t->Branch("cut", &cut);
  t->Branch("inputSname", &inputSname);
  t->Branch("inputBname", &inputBname);
  auto rootfile = rmva->info("rootfile");
  t->Branch("rootfile", &rootfile);
  t->Branch("config_train", &configfile_train);
  t->Branch("config_cut", &configfile_cut);
  t->Fill();
  t->Write();

  outf->mkdir("dataset")->cd();
  for (auto& h : rmva->hrocs("MVA_.*_rejBvsS"))
    xjjroot::writehist(h);
  outf->cd();
  
  outf->Close();
  
  return 0;
  
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
