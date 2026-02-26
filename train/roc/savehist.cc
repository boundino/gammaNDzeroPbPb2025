#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"

#define __BINS_PTY_EQ__
#define __BINS_MASS__
#include "../../include/bins.h"

#include "mvaroot.h"
#include "save_calc.h"

int macro(std::string config_train, std::string config_sig) {
  std::cout<<std::endl;
  xjjc::config conf_t(config_train), conf_s(config_sig);

  // train
  auto* rmva = new mytmva::mvaroot(conf_t.get("Rootoutput"));
  if (!rmva->good()) {
    __XJJLOG << "!! bad input " << rmva->info("rootfile") << std::endl;
    return 2;
  }
  // rmva->print();
  auto dir = conf_t.has("Name") ? conf_t.get("Name") : xjjc::str_tag_from_file(config_train);
  auto cuts = rmva->info("cuts"), cutb = rmva->info("cutb");

  // significance
  std::map<std::string, TChain*> trs;
  auto prepare_tree = [&trs, &conf_s](const std::string &key, const std::string &Input) {
    trs[key] = xjjana::chain_files(conf_s.get_vec(Input), "Tree");
    if (!trs[key]) { __XJJLOG<<"!! bad input file "<<conf_s.get(Input)<<std::endl; return std::string(""); }
    save::mask_branch(trs[key]);
    return conf_s.get(Input);
  };
  auto method = conf_s.get("Method");
  auto name = conf_s.has("Name") ? conf_s.get("Name") : xjjc::str_tag_from_file(config_sig);
  std::cout<<"\e[1m"<<dir<<" / "<<name<<" [ "<<method<<" ]\e[0m"<<std::endl;
  if (!rmva->hroc("MVA_" + method + "_S_high")) {
    __XJJLOG << "!! no hist "<<("MVA_" + method + "_S_high")<<std::endl; return 2;
  }
  auto inputSname = prepare_tree("_S", "Inputs"); if (inputSname.empty()) return 2;
  auto inputBname = prepare_tree("_B", "Inputb"); if (inputBname.empty()) return 2;

  auto* outf = xjjroot::newfile("rootfiles/" + dir + "/" + name + ".root");

  const auto h_dump = rmva->hroc("MVA_" + method + "_S_high");
  std::map<std::string, TH3F*> h3_meths;
  auto make_h3 = [&h3_meths, &method, &trs, &h_dump](const std::string &key, const std::string &icut) { // can add href
    h3_meths[key+"_high"] = new TH3F(Form("h3_%s%s_high", method.c_str(), key.c_str()), Form(";y;m_{K#pi} (GeV);%s", method.c_str()),
                                     bins::ny, bins::miny, bins::maxy,
                                     bins::nmass, bins::minmass, bins::maxmass,
                                     h_dump->GetXaxis()->GetNbins(), h_dump->GetXaxis()->GetXmin(), h_dump->GetXaxis()->GetXmax());
    trs[key]->Project(h3_meths[key+"_high"]->GetName(), "Dmva:Dmass:Dy", icut.c_str());
    xjjroot::writehist(h3_meths[key+"_high"]);
  };
  make_h3("_S", cuts);
  make_h3("_B", cutb);

  outf->mkdir("datasets")->cd();
  rmva->write_hrocs(".*_" + method + "_.*");
  outf->cd();

  auto* t = new TTree("info", "");
  t->Branch("cuts", &cuts);
  t->Branch("cutb", &cutb);
  t->Branch("inputSname", &inputSname);
  t->Branch("inputBname", &inputBname);
  t->Branch("method", &method);
  std::string rootfile = rmva->info("rootfile");
  t->Branch("rootfile", &rootfile);
  t->Branch("config_train", &config_train);
  t->Branch("config_sig", &config_sig);
  t->Fill();
  t->Write();

  outf->Close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
