#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

#include "mvaroot.h"

int ngroup(const mytmva::mvaroot* m);

const float eff_anchor = 0.3;

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
  const auto method = info.at("method");

  auto* h_dump_rejBvsS = static_cast<TH1D*>(rmva.front()->hrocs("MVA_.*_rejBvsS").front()->Clone("h_dump_rejBvsS")); h_dump_rejBvsS->Reset();

  // h3
  std::map<std::string, TH3F*> h3s;
  TH3F* h3_dump = nullptr;
  for (const std::string &t : { "S", "B", "Swap", "Mass" }) {
    h3s[t] = xjjana::getobj_regexp<TH3F>(inf, ".*_"+t+"_high").front(); //
    h3s.at(t)->Sumw2();
    if (!h3_dump) { h3_dump = (TH3F*)h3s.at(t)->Clone("h3_dump"); h3_dump->Reset(); }
  }
  xjjroot::print_tab(h3s, 0);

  std::map<std::string, std::vector<TH1D*>> h1ys;
  // training binning
  for (const std::string &t : { "S_trainbin", "B_trainbin" }) {
    // h3_S_trainbin__y-2_high
    for (auto& h3 : xjjana::getobj_regexp<TH3F>(inf, ".*_"+t+"__y-[0-9]*_high")) {
      int iy_train = atoi(xjjc::str_eraseall(xjjc::str_erasestar(h3->GetName(), "*__y-"), "_high").c_str());
      h3->Sumw2();
      // loop ny
      for (int i=0; i<h3->GetXaxis()->GetNbins(); i++) {
        if (iy_train != i) continue;
        auto* h1_high = h3->ProjectionZ(Form("h1_%s_high__y-%d", t.c_str(), i), // Z : mva
                                        i+1, i+1, // X : y
                                        0, -1, "e"); // Y : mass
        h1ys[t+"_high"].push_back(h1_high); //

        auto ngroup_rebin = ngroup(rmva[i]);
        auto* h1_rebin = (TH1D*)h1_high->Rebin(ngroup_rebin, xjjc::str_eraseall(h1_high->GetName(), "_high").c_str());
        h1ys[t].push_back(h1_rebin); //
      }
    }
  }
  
  // efficiency
  for (const auto& [t, h3] : h3s) {
    if (t=="S" || t=="B") {
      // loop ny
      for (int i=0; i<h3->GetXaxis()->GetNbins(); i++) {
        auto* h1_high = h3->ProjectionZ(Form("h1_%s_high__y-%d", t.c_str(), i), // Z : mva
                                        i+1, i+1, // X : y
                                        0, -1, "e"); // Y : mass
        h1ys[t+"_high"].push_back(h1_high); //

        // eff
        if (t=="S" || t=="B") {    
          auto* h1_eff = (TH1D*)h1_high->Clone(xjjc::str_replaceall(h1_high->GetName(), t+"_high", "eff"+t).c_str());
          float sum = 0;
          for (int j=h1_high->GetXaxis()->GetNbins()+1; j>=0; j--) {
            sum += h1_high->GetBinContent(j);
            h1_eff->SetBinContent(j, sum);
            h1_eff->SetBinError(j, 0);
          }
          h1_eff->Scale(1./sum);
          if (h1_eff->GetBinContent(0) != h1_eff->GetBinContent(1))
            __XJJLOG << "!! underflow: " << h1_eff->GetName() << " [ "<<h1_eff->GetBinContent(0)<<" vs "<<h1_eff->GetBinContent(1)<<" ]" << std::endl;
          if (h1_eff->GetBinContent(h1_eff->GetXaxis()->GetNbins()+1) != h1_eff->GetBinContent(h1_eff->GetXaxis()->GetNbins()))
            __XJJLOG << "!! overflow: " << h1_eff->GetName() << " [ "<<h1_eff->GetBinContent(h1_eff->GetXaxis()->GetNbins()+1)<<" vs "<<h1_eff->GetBinContent(h1_eff->GetXaxis()->GetNbins())<<" ]" << std::endl;
          h1ys["eff"+t].push_back(h1_eff); //
        }
      } // for (int i=0; i<h3->GetXaxis()->GetNbins(); i++) {
    }

    if (t=="S" || t=="Swap") {
      for (int i=0; i<h3->GetXaxis()->GetNbins(); i++) {
        auto* h1_mass = h3->ProjectionY(Form("h1_massTemplate_%s__y-%d", t.c_str(), i), i+1, i+1, 0, -1, "e");
        h1ys["mass_"+t].push_back(h1_mass);
      }
    }
  } // for (const auto& [t, h3] : h3s) {

  std::map<std::string, TH1D*> h1s;
  h1s["mvaAnchor"] = new TH1D("h1_mvaAnchor", Form(";MVA;jMVA at efficiency ~ %.2f", eff_anchor), h3_dump->GetXaxis()->GetNbins(), h3_dump->GetXaxis()->GetXmin(), h3_dump->GetXaxis()->GetXmax());
  std::map<std::string, std::vector<TGraph*>> grys;
  // loop ny
  for (int i=0; i<h3_dump->GetXaxis()->GetNbins(); i++) {
    const auto *h1_effS = h1ys["effS"].at(i), *h1_effB = h1ys["effB"].at(i); // eff vs MVA
    auto *h1_rejBvsS = (TH1D*)h_dump_rejBvsS->Clone(Form("h1_rejBvsS__y-%d", i));
    auto *h1_effvsMVA = (TH1D*)h_dump_rejBvsS->Clone(Form("h1_effvsMVA__y-%d", i)); h1_effvsMVA->GetYaxis()->SetTitle(method.c_str());
    auto *gr_rejBvsS_high = new TGraph(h1_effS->GetXaxis()->GetNbins()); gr_rejBvsS_high->SetName(Form("gr_rejBvsS_high__y-%d", i));
    int eff_bin_now = 1, jmva_now = -1; float eff_dev_now = 1.e+5;
    int jmva_anchor = 1; float eff_dev_anchor = 1.e+5; // !!
    for (int j=h1_effS->GetXaxis()->GetNbins(); j>=1; j--) {
      const float effS = h1_effS->GetBinContent(j), effB = h1_effB->GetBinContent(j);
      gr_rejBvsS_high->SetPoint(j-1, effS, 1-effB);
      // anchor point
      if (fabs(effS-eff_anchor) <= eff_dev_anchor) {
        jmva_anchor = j;
        eff_dev_anchor = fabs(effS-eff_anchor);
      }
      // 100 effS bin
      const float eff_target_now = h1_rejBvsS->GetBinCenter(eff_bin_now);
      if (fabs(effS-eff_target_now) > eff_dev_now) {
        h1_rejBvsS->SetBinContent(eff_bin_now, 1-effB);
        h1_effvsMVA->SetBinContent(eff_bin_now, h1_effS->GetBinCenter(jmva_now));
        // move to next target eff
        eff_bin_now++;
      }
      
      eff_dev_now = fabs(effS - h1_rejBvsS->GetBinCenter(eff_bin_now));
      jmva_now = j;
    }
    __XJJLOG << ">> jMVA for eff = "<<eff_anchor<<": "<<jmva_anchor<<" ("<<h1_effS->GetBinCenter(jmva_anchor)<<")"<<std::endl;
    h1s["mvaAnchor"]->SetBinContent(i+1, h1_effS->GetBinCenter(jmva_anchor));
    auto* hanchor = h3s["Mass"]->ProjectionY(Form("h1_massAnchor__y-%d", i), i+1, i+1, jmva_anchor, h3s["Mass"]->GetZaxis()->GetNbins(), "e");
    
    h1ys["rejBvsS"].push_back(h1_rejBvsS);
    h1ys["effvsMVA"].push_back(h1_effvsMVA);
    h1ys["massAnchor"].push_back(hanchor);
    grys["rejBvsS_high"].push_back(gr_rejBvsS_high);
  }
  
  xjjroot::print_tab(h1s, 0);
  xjjroot::print_tab(h1ys, 0);
  xjjroot::print_tab(grys, 0);

  auto* outf = xjjroot::newfile(xjjc::str_getdir(inputname) + "/calc_" + xjjc::str_erasestar(inputname, "*/"));
  for (const auto& [_, h] : h3s) xjjroot::writehist(h);
  for (const auto& [_, h] : h1s) xjjroot::writehist(h);
  for (const auto& [_, hh] : h1ys)
    for (const auto& h : hh) xjjroot::writehist(h);
  for (const auto& [_, hh] : grys)
    for (const auto& h : hh) xjjroot::writehist(h);
  for (auto& rr : rmva) {
    auto* dr = outf->mkdir(rr->info("dir").c_str());
    rr->rm_info("dir");
    dr->mkdir("dataset")->cd();
    rr->write_hrocs();
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

int ngroup(const mytmva::mvaroot* m) {
  const auto meth = m->methods().front();
  const auto *h_high = m->hroc("MVA_"+meth+"_S_high"),
    *h_rebin = m->hroc("MVA_"+meth+"_S");
  if (h_high && h_rebin) {
    // __XJJLOG << ">> h_rebin variable array: " << h_rebin->GetXaxis()->GetXbins()->GetSize() << std::endl;
    if (h_high->GetXaxis()->GetNbins() % h_rebin->GetXaxis()->GetNbins() != 0)
      __XJJLOG << "!! bad rebin ngroup" << std::endl;
    return h_high->GetXaxis()->GetNbins() / h_rebin->GetXaxis()->GetNbins();
  }
  return 0;
}

