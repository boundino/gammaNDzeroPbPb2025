#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

#include "mvaroot.h"
#include "save_calc.h"

TLegend* stylehists(std::vector<TH1D*>&);
int macro(std::string inputname) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  const auto& conf = xjjana::getstr_regexp((TTree*)inf->Get("info"), ".*");
  xjjc::print_tab(conf, -1);
  std::map<std::string, TH3F*> h3;
  h3["_S"] = xjjana::getobj<TH3F>(inf, "h3_cut_S");
  h3["_B"] = xjjana::getobj<TH3F>(inf, "h3_cut_B");
  auto effS = calc::eff_cut(h3["_S"]), effB = calc::eff_cut(h3["_B"]);
  auto* grcut = xjjroot::drawpoint(effS, 1-effB, kBlack, 47, 4, "goff");
  grcut->SetName("grcut");
  xjjroot::print_gr(grcut);
  auto hrejBvsS = xjjana::getobj_regexp<TH1D>(inf->GetDirectory("dataset"), "MVA_.*_rejBvsS", "TH1");
  
  auto* leg = stylehists(hrejBvsS);
  auto* hempty = new TH2F("hempty", ";Signal efficiency;Background rejection", 10, 0, 1, 10, 0, 1);
  xjjroot::sethempty(hempty);

  xjjroot::setgstyle(1);  
  auto* pdf = new xjjroot::mypdf("figpdfs/" + xjjc::str_erasestar(xjjc::str_divide_once(inputname, "/").back(), ".*") + ".pdf");
  pdf->prepare();
  hempty->Draw("axis");
  for (const auto& h : hrejBvsS) {
    h->Draw("c same");
  }
  grcut->Draw("psame");
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, "2025 PbPb (5.36 TeV)");
  pdf->write();

  // zoom
  hempty = new TH2F("hempty_zoom", ";Signal efficiency;Background rejection", 10, 0, 0.6, 10, 0.9, 1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");
  for (const auto& h : hrejBvsS) {
    h->Draw("c same");
  }
  grcut->Draw("psame");
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, "2025 PbPb (5.36 TeV)");
  pdf->write();

  pdf->draw_cover( {
      "#bf{Signal cut} " + conf.at("cuts"),
      "#bf{Background cut} " + conf.at("cutb"),
      "#bf{Single working point} " + conf.at("cut"),
      "#bf{Signal sample} " + conf.at("inputSname"),
      "#bf{Background sample} " + conf.at("inputBname"),
    }, 0.03);
  
  pdf->close();
  
  return 0;
  
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

TLegend* stylehists(std::vector<TH1D*>& hists) {
  auto* leg = new TLegend(0.25, 0.5-0.04*hists.size(), 0.5, 0.5);
  xjjroot::setleg(leg, 0.038);
  for (std::size_t i=0; i<hists.size(); i++) {
    auto& h = hists.at(i);
    xjjroot::setthgrstyle(h, xjjroot::colorlist_middle[i], 21, 1, xjjroot::colorlist_middle[i], 1, 4);
    auto name = xjjc::str_eraseall(h->GetName(), "MVA_");
    name = xjjc::str_erasestar(name, "_*");
    leg->AddEntry(h, name.c_str(), "l");
  }
  return leg;
}

