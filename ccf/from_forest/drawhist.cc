#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "config.h"
#include "xjjmypdf.h"

#include "clusComp.h"

int macro(std::string configfile) {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto outputname = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);
  auto* inf = TFile::Open(Form("rootfiles/%s.root", outputname.c_str()));
  auto* h_qual_nhit = xjjana::getobj<TH2F>(inf, "h_qual_nhit");
  xjjroot::sethempty(h_qual_nhit, 0, -0.2);
  auto* h_qual_nhit_pass = xjjana::getobj<TH2F>(inf, "h_qual_nhit_pass");
  xjjroot::sethempty(h_qual_nhit_pass, 0, -0.2);

#define DRAW_COMMON() \
  xjjroot::drawCMS(xjjroot::CMS::internal, "2025 PbPb (5.36 TeV)");
  
  xjjroot::setgstyle(1, 2, 1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  pdf->prepare();
  h_qual_nhit->Draw("colz");
  globals::drawline(h_qual_nhit);
  // xjjana::drawhoutline(h_qual_nhit_pass, kRed, 1, 2);
  DRAW_COMMON();
  pdf->write();

  pdf->prepare();
  h_qual_nhit_pass->Draw("colz");
  globals::drawline(h_qual_nhit_pass);
  DRAW_COMMON();
  pdf->write();

  pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

