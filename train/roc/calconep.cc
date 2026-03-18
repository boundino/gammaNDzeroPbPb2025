#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

#include "mvaroot.h"
#include "save_calc.h"

int macro(std::string inputname) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  if (!inf) {
    __XJJLOG << "!! bad input file, abort." << std::endl; return 2;
  }
  std::vector<mytmva::mvaroot*> rmva;
  TIter next(inf->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    auto* obj = key->ReadObj();
    if (!obj->InheritsFrom(TDirectory::Class())) continue;
    if (!std::regex_match(obj->GetName(), std::regex("pt-[0-9]+-[0-9]+_y-[M0-9]+-[M0-9]+"))) continue;
    __XJJLOG << ">> " << obj->GetName() << " >> " << static_cast<TDirectory*>(obj)->GetPath() << std::endl;
    auto* r = new mytmva::mvaroot(static_cast<TDirectory*>(obj));
    if (!r->good()) continue;
    r->add_info("dir", obj->GetName());
    rmva.push_back(r);
    // r->print();
    // r->print_info();
  }
  if (rmva.empty()) {
    __XJJLOG << "!! no good mva root files, abort." << std::endl;
    return 2;
  }
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  xjjc::print_tab(info, -1);

  std::map<std::string, TH3F*> h3s;
  std::map<std::string, std::vector<TH1D*>> h1ys;
  for (const std::string &t : { "S", "B", "Swap", "Mass" }) {
    h3s[t] = xjjana::getobj<TH3F>(inf, "h3_cut_"+t); //
    h3s.at(t)->Sumw2();

    for (int i=0; i<h3s.at(t)->GetXaxis()->GetNbins(); i++) {
      std::pair<int, int> jcut = { 0, -1 }; // doesn't matter if the cut is passed or not
      if (t=="Mass") jcut = { 2, 2 };
      auto* h1_mass = h3s.at(t)->ProjectionY(Form("h1_mass_%s__y-%d", t.c_str(), i),
                                             i+1, i+1,
                                             jcut.first, jcut.second,
                                             "e");
      h1ys["mass_"+t].push_back(h1_mass);
    }
  }
  xjjroot::print_tab(h3s, 0);
  xjjroot::print_tab(h1ys, 0);

  // roc for the one point
  __XJJLOG << "-- reference point" << std::endl;
  std::map<std::string, std::vector<TGraph*>> grys;
  for (int i=0; i<h3s.at("S")->GetXaxis()->GetNbins(); i++) {
    auto effS = calc::eff_cut(h3s.at("S"), i+1, i+1), effB = calc::eff_cut(h3s.at("B"), i+1, i+1);
    auto* grcut = xjjroot::drawpoint(effS, 1-effB, kBlack, 47, 4, "goff");
    grcut->SetName(Form("gr_roc_cut__y-%d", i));
    xjjroot::print_gr(grcut);
    grys["Onep"].push_back(grcut);
  }

  auto* outf = xjjroot::newfile(xjjc::str_getdir(inputname) + "/calc_" + xjjc::str_erasestar(inputname, "*/"));
  for (const auto& [_, h] : h3s) xjjroot::writehist(h);
  for (const auto& [_, hh] : h1ys)
    for (const auto& h : hh) xjjroot::writehist(h);
  for (const auto& [_, hh] : grys)
    for (const auto& h : hh) xjjroot::writehist(h);
  for (auto& rr : rmva) {
    auto* dr = outf->mkdir(rr->info("dir").c_str());
    rr->rm_info("dir");
    dr->mkdir("dataset")->cd();
    rr->write_hrocs(".+_rejBvsS");
    auto* t = new TTree("tmvainfo", "");
    rr->write_info(t);
    t->Fill();
    t->Write();
    outf->cd();
  }
  auto* t = new TTree("info", "");
  for (auto& [key, content] : info) {
    t->Branch(key.c_str(), &content);
  }
  t->Fill();
  t->Write();

  outf->Close();
  return 0;  

}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

