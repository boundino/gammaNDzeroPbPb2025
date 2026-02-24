#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

// #define __BINS_PTY_EQ__
// #define __BINS_MASS__
// #include "../../include/bins.h"

TLegend* stylehists(std::vector<TH1D*>&);

int macro(std::string inputname) {
  std::cout<<std::endl;
  auto* inf = TFile::Open(inputname.c_str());
  auto* tr = (TTree*)inf->Get("info");
  const auto& conf = xjjana::getstr_regexp(tr, ".*");
  for (const auto& [t, val] : conf) std::cout<<t<<" \e[2m"<<val<<"\e[0m"<<std::endl;
  auto hrocs = xjjana::getobj_regexp<TH1D>(inf, ".*", "TH1");
  auto* grcut = xjjana::getobj<TGraph>(inf, "grcut");
  xjjroot::setmarkerstyle(grcut, kBlack, -1, 4);
  xjjroot::print_gr(grcut);

  auto* leg = stylehists(hrocs);
  auto* hempty = new TH2F("hempty", ";Signal efficiency;Background rejection", 10, 0, 1, 10, 0, 1);
  xjjroot::sethempty(hempty);

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figpdfs/" + conf.at("dir") + "/" + conf.at("name") + ".pdf");
  pdf->prepare();
  hempty->Draw("axis");
  for (const auto& h : hrocs) {
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

