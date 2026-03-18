#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"

#define __BINS_PTY_EQ__
#define __BINS_MASS__
#include "../../include/bins.h"

#include "mvaroot.h"
#include "save_calc.h"

int macro(std::string config_train, std::string config_cut) {
  std::cout<<std::endl;
  xjjc::config conf_t(config_train), conf_c(config_cut);

  auto dir = conf_t.has("Name") ? conf_t.get("Name") : xjjc::str_tag_from_file(config_train);
  auto name = conf_c.has("Name") ? conf_c.get("Name") : xjjc::str_tag_from_file(config_cut);
  std::cout<<"\e[1m"<<dir<<" / "<<name<<"\e[0m"<<std::endl;

  // train
  std::vector<mytmva::mvaroot*> rmva;
  for (const auto& rfname : conf_t.get_vec("Rootoutputs")) {
    __XJJLOG << ">> " << rfname << std::endl;
    auto* rr = new mytmva::mvaroot(rfname);
    if (!rr->good()) {
      __XJJLOG << "!! bad mva input " << rfname << ", abort." << std::endl;
      return 2;
    }
    if (!( rr->has_info<float>("ptmin") && rr->has_info<float>("ptmax") && rr->has_info<float>("ymin") && rr->has_info<float>("ymax") )) {
      __XJJLOG << "?? " << rfname << " doesn't have [ptmin, ptmax, ymin, ymax], set to inclusive" << std::endl;
      rr->add_info<float>("ptmin", 2); rr->add_info<float>("ptmax", 5); rr->add_info<float>("ymin", -2); rr->add_info<float>("ymax", 2);
    }
    // rr->print();
    rmva.push_back(rr);
  }
  if (rmva.empty()) {
    __XJJLOG << "!! no good mva root files, abort." << std::endl;
    return 2;
  }

  // cut
  std::map<std::string, std::string> cut;
  for (const auto& rr : rmva) {
    cut["S"] = (cut.find("S")==cut.end() ? rr->info("cutS") : save::common_cut(cut.at("S"), rr->info("cutS")));
    cut["B"] = (cut.find("B")==cut.end() ? rr->info("cutB") : save::common_cut(cut.at("B"), rr->info("cutB")));
    auto rrcutMass = (rr->has_info("cutData") ? rr->info("cutData") : xjjc::str_eraseall(rr->info("cutB"), " && TMath::Abs(Dmass-1.8648) > 0.05 && TMath::Abs(Dmass-1.8648) < 0.12"));
    cut["Mass"] = (cut.find("Mass")==cut.end() ? rrcutMass : save::common_cut(cut.at("Mass"), rrcutMass));
  }
  cut["Swap"] = xjjc::str_replaceall(cut.at("S"), "23333", "23344"); // !!
  cut["Onep"] = conf_c.get("Cut");
  xjjc::print_tab(cut, -1);

  std::map<std::string, TChain*> trs; std::map<std::string, std::string> inputstr;
  auto prepare_tree = [&trs, &conf_c, &inputstr](const std::string &key, const std::string &Input) {
    trs[key] = xjjana::chain_files(conf_c.get_vec(Input), "Tree");
    if (!trs[key]) { __XJJLOG<<"!! bad input file "<<conf_c.get(Input)<<std::endl; return 2; }
    save::mask_branch(trs[key]);
    inputstr[key] = conf_c.get(Input); //
    return 0;
  };
  if (prepare_tree("S", "InputsS")) return 2;
  if (prepare_tree("B", "InputsB")) return 2;
  if (prepare_tree("Mass", conf_c.has("InputsMass") ? "InputsMass" : "InputsB")) return 2;

  auto* outf = xjjroot::newfile("rootfiles/" + dir + "/" + name + ".root");
  
  std::map<std::string, TH3F*> h3;
  auto project = [&h3, &cut](TChain* tr, std::string key, std::string icut) {
    h3[key] = new TH3F(Form("h3_cut_%s", key.c_str()), ";y;m_{K#pi} (GeV);Pass",
                       bins::ny, bins::miny, bins::maxy,
                       bins::nmass, bins::minmass, bins::maxmass,
                       2, 0, 2);
    __XJJLOG << ">> "<<h3[key]->GetName()<<" [ "<<cut.at("Onep")<<" ] \e[2m"<<icut<<"\e[0m"<<std::endl;
    tr->Project(h3[key]->GetName(), Form("%s:Dmass:Dy", cut.at("Onep").c_str()), icut.c_str());
    xjjroot::writehist(h3[key]);
    const auto eff = calc::eff_cut(h3[key]);
    __XJJLOG << ">> eff: [ "<<eff<<" ] rej: [ "<<1-eff<<" ]" << std::endl;;
  };
  project(trs.at("S"), "S", cut.at("S"));
  project(trs.at("B"), "B", cut.at("B"));
  project(trs.at("S"), "Swap", cut.at("Swap"));
  project(trs.at("Mass"), "Mass", cut.at("Mass"));

  for (auto& rr : rmva) {
    auto* dr = outf->mkdir(Form("pt-%s-%s_y-%s-%s",
                                xjjc::number_to_string(rr->info<float>("ptmin")).c_str(),
                                xjjc::number_to_string(rr->info<float>("ptmax")).c_str(),
                                xjjc::number_to_string(rr->info<float>("ymin")).c_str(),
                                xjjc::number_to_string(rr->info<float>("ymax")).c_str()
                                ));
    dr->mkdir("dataset")->cd();
    rr->write_hrocs();
    auto* t = new TTree("tmvainfo", "");
    rr->write_info(t);
    t->Fill();
    t->Write();
    outf->cd();
  }
  
  auto* t = new TTree("info", "");
  auto branchrval = [&t](const std::string& bname, const std::string& content) {
    auto* tr = new std::string(content);
    t->Branch(bname.c_str(), tr);
  };
  for (auto& [k, c] : cut)
    t->Branch(Form("cut%s", k.c_str()), &c);
  for (auto& [k, c] : inputstr)
    t->Branch(Form("input%sname", k.c_str()), &c);
  branchrval("label", conf_t.get("Label"));
  branchrval("label2", conf_t.get("Label2"));
  branchrval("labelcut", conf_c.get("Label"));
  t->Branch("config_train", &config_train);
  t->Branch("config_cut", &config_cut);
  t->Fill();
  t->Write();

  outf->mkdir("dataset")->cd();
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
