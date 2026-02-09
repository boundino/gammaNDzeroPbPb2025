#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

#include "clusComp.h"

int macro(std::string configfile) {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto outputname = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);
  auto* inf = TFile::Open(Form("rootfiles/%s.root", outputname.c_str()));
  auto labels = conf.get_vec("Label_2");
  auto ismc = conf.has("isMC") ? conf.get<bool>("isMC") : false;

#define GET_AND_DRAW(q)                                                 \
  auto* h_qual_nhit##q = xjjana::getobj<TH2F>(inf, "h_qual_nhit" #q);   \
  xjjroot::sethempty(h_qual_nhit##q, 0, -0.3);                          \
  h_qual_nhit##q->Draw("colz");                                         \
  globals::drawline(h_qual_nhit##q);                                    \
  xjjroot::drawCMS(ismc ? xjjroot::CMS::simulation : xjjroot::CMS::internal, conf.get("Label_1")); \
  xjjroot::drawtexgroup(0.83, 0.85, labels, 0.04, 33, 42, 1.15);

  // xjjana::drawhoutline(h_qual_nhit_pass, kRed, 1, 2);
  xjjroot::setgstyle(1, 2, 1);
  xjjroot::adjust_margin(1, 3, 1, 0.7);
  TGaxis::SetMaxDigits(4);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  pdf->prepare();
  GET_AND_DRAW();
  pdf->write();

  labels.push_back("CCF applied manually");
  pdf->prepare();
  GET_AND_DRAW(_pass);
  pdf->write();

  labels.pop_back();
  labels.push_back("CCF in forest applied");
  pdf->prepare();
  GET_AND_DRAW(_filter);
  pdf->write();

  pdf->draw_cover({
      "#bf{" + outputname + "}",
      "#bf{Inputs} " + conf.get("Inputs"),
      "#bf{Base cut} " + conf.get("Cut"),
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

