#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

void print_ratio(std::vector<TH1F*> hists) {
  if (hists.empty()) return;
  auto denom = hists.front()->Integral(0, hists.front()->GetXaxis()->GetNbins()+1);
  if (!denom) return;
  for (int i=0; i<hists.size(); i++) {
    if (!i) continue;
    auto r = hists[i]->Integral(0, hists[i]->GetXaxis()->GetNbins()+1) / denom;
    xjjroot::drawtex(0.55, 0.85-0.043*4-0.1-0.042*(i-1), Form("Ratio = %.2f", r), 0.04, 13, 42, hists[i]->GetLineColor());
  }
}

int macro(std::string configfile) {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto name = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);
  // std::cout<<"\e[1m[ "<<name<<" ]\e[0m"<<std::endl;
  auto* inf = TFile::Open(Form("rootfiles/%s.root", name.c_str()));
  auto labels = conf.get_vec("Label_2");
  const auto& isgammaN = conf.get<int>("IsgammaN");
  const std::string& v_gamma = isgammaN ? "Plus" : "Minus",
    v_N = isgammaN ? "Minus" : "Plus";
  const auto& cut = conf.has("Cut") ? conf.get("Cut") : "(1)";
  const auto& tPFptmin = conf.get("PFptmin");

#define GET_HIST(l)                                                     \
  auto* hHFEMax##l = xjjana::getobj<TH1F>(inf, "hHFEMax" #l);           \
  xjjroot::sethempty(hHFEMax##l, 0, 0);                                 \
  if (xjjc::str_contains(hHFEMax##l->GetName(), "_forest")) {           \
    xjjroot::setthgrstyle(hHFEMax##l, kGray, 21, 1, kGray, 1, 3, kGray, 0.5); \
  } else if (xjjc::str_contains(hHFEMax##l->GetName(), "_eta5")) {      \
    xjjroot::setthgrstyle(hHFEMax##l, xjjroot::mycolor_satmiddle["azure"], 21, 1, xjjroot::mycolor_satmiddle["azure"], 1, 3, xjjroot::mycolor_satmiddle["azure"], 1., 3054); \
  } else {                                                              \
    xjjroot::setthgrstyle(hHFEMax##l, xjjroot::mycolor_satmiddle["red"], 21, 1, xjjroot::mycolor_satmiddle["red"], 1, 3, xjjroot::mycolor_satmiddle["red"], 1., 3045); \
  }

  GET_HIST(_gamma);
  GET_HIST(_gamma_forest);
  GET_HIST(_gamma_eta5);
  GET_HIST(_N);
  GET_HIST(_N_forest);
  GET_HIST(_N_eta5);
  GET_HIST(_gamma_gap);
  GET_HIST(_gamma_forest_gap);
  GET_HIST(_gamma_eta5_gap);
  GET_HIST(_gamma_gaploose);
  GET_HIST(_gamma_forest_gaploose);
  GET_HIST(_gamma_eta5_gaploose);

  auto* leg = new TLegend(0.55, 0.85-0.043*4, 0.8, 0.85);
  leg->SetHeader("Leading tower in PF", "L");
  xjjroot::setleg(leg, 0.04);
  leg->AddEntry(hHFEMax_gamma_forest, "p_{T} > 0, |#eta| < 6", "f");
  leg->AddEntry(hHFEMax_gamma, Form("p_{T} > %s, |#eta| < 5.2", tPFptmin.c_str()), "l");
  leg->AddEntry(hHFEMax_gamma_eta5, Form("p_{T} > %s, |#eta| < 5", tPFptmin.c_str()), "l");

#define DRAW_GLOBAL()                                                   \
  xjjroot::drawCMS(xjjroot::CMS::internal, conf.get("Label_1"));        \
  leg->Draw();

  xjjroot::setgstyle(1);
  TGaxis::SetMaxDigits(4);
  TGaxis::SetExponentOffset(-0.1, 0., "y");
  auto* pdf = new xjjroot::mypdf("figspdf/" + name + ".pdf");
  pdf->prepare();
  hHFEMax_gamma_forest->Draw("hist");
  hHFEMax_gamma_eta5->Draw("hist same");
  hHFEMax_gamma->Draw("hist same");
  print_ratio({ hHFEMax_gamma_forest, hHFEMax_gamma, hHFEMax_gamma_eta5 });
  DRAW_GLOBAL();
  pdf->write();

  pdf->prepare();
  hHFEMax_gamma_forest_gap->Draw("hist");
  hHFEMax_gamma_eta5_gap->Draw("hist same");
  hHFEMax_gamma_gap->Draw("hist same");  
  print_ratio({ hHFEMax_gamma_forest_gap, hHFEMax_gamma_gap, hHFEMax_gamma_eta5_gap });
  DRAW_GLOBAL();
  pdf->write();

  pdf->prepare();
  hHFEMax_gamma_forest_gaploose->Draw("hist");
  hHFEMax_gamma_eta5_gaploose->Draw("hist same");
  hHFEMax_gamma_gaploose->Draw("hist same");  
  print_ratio({ hHFEMax_gamma_forest_gaploose, hHFEMax_gamma_gaploose, hHFEMax_gamma_eta5_gaploose });
  DRAW_GLOBAL();
  pdf->write();

  pdf->prepare();
  hHFEMax_N_forest->Draw("hist");
  hHFEMax_N_eta5->Draw("hist same");
  hHFEMax_N->Draw("hist same");
  DRAW_GLOBAL();
  pdf->write();

  pdf->draw_cover({
      "#bf{" + name + "}",
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

// bool filter(const double& quality, const double& nPxlHits) {
//   bool accept = true;
//   // construct polynomial cut on cluster vertex quality vs. npixelhits
//   double polyCut = 0;
//   for (unsigned int i = 0; i < globals::clusterPars_.size(); i++) {
//     polyCut += globals::clusterPars_[i] * std::pow((double)nPxlHits, (int)i);
//   }
//   if (nPxlHits < globals::nhitsTrunc_)
//     polyCut = 0;  // don't use cut below nhitsTrunc_ pixel hits
//   if (polyCut > globals::clusterTrunc_ && globals::clusterTrunc_ > 0)
//     polyCut = globals::clusterTrunc_;  // no cut above clusterTrunc_

//   if (quality < polyCut)
//     accept = false;

//   // return with final filter decision
//   return accept;
// }

