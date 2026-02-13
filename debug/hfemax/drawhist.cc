#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

void print_ratio(std::vector<TH1F*> hists) {
  if (hists.empty()) return;
  auto denom = hists.front()->Integral(0, hists.front()->GetXaxis()->GetNbins()+1);
  if (!denom) return;
  xjjroot::drawtex(0.52, 0.5, "Ratio to", 0.04, 13);
  xjjroot::drawtex(0.66, 0.5, "forest event tree", 0.04, 13, 42, hists.front()->GetLineColor());
  for (int i=0; i<hists.size(); i++) {
    if (!i) continue;
    auto r = hists[i]->Integral(0, hists[i]->GetXaxis()->GetNbins()+1) / denom;
    xjjroot::drawtex(0.52, 0.5-0.042*i, Form("%.2f", r), 0.04, 13, 42, hists[i]->GetLineColor());
  }
}

struct Hist{
  TH1F* h_gamma;
  TH1F* h_gamma_gap;
  TH1F* h_gamma_gaploose;
  TH1F* h_N;
};

xjjroot::thgrstyle style_by_var(const std::string& v) {
  if (xjjc::str_contains(v, "_forest")) {
    return xjjroot::thgrstyle({ kGray+1, 21, 1, kGray+1, 1, 3, kGray+1, 0.2, 1001 });
  } else if (xjjc::str_contains(v, "_eta5")) {
    return xjjroot::thgrstyle({ xjjroot::mycolor_satmiddle["red"], 21, 1, xjjroot::mycolor_satmiddle["red"], 1, 3, 0, 0, 0 });
  } else if (xjjc::str_contains(v, "_pt0p1")) {
    return xjjroot::thgrstyle({ xjjroot::mycolor_middle["olive"], 21, 1, xjjroot::mycolor_middle["olive"], 1, 3, 0, 0, 0 });
  }
  return xjjroot::thgrstyle({ xjjroot::mycolor_satmiddle["azure"], 21, 1, xjjroot::mycolor_satmiddle["azure"], 1, 3, xjjroot::mycolor_satmiddle["azure"], 0.5, 3354 });
}

std::string leg_by_var(const std::string& v, float ptmin_sample, float etamax_sample) {
  if (xjjc::str_contains(v, "_forest")) {
    return "p_{T} > 0, |#eta| < 6";
  } else if (xjjc::str_contains(v, "_eta5")) {
    return Form("p_{T} > %s, |#eta| < %s", xjjc::number_remove_zero(ptmin_sample).c_str(), (etamax_sample<5?xjjc::number_remove_zero(etamax_sample).c_str():"5"));
  } else if (xjjc::str_contains(v, "_pt0p1")) {
    return Form("p_{T} > %s, |#eta| < %s", (ptmin_sample>0.1?xjjc::number_remove_zero(ptmin_sample).c_str():"0.1"), xjjc::number_remove_zero(etamax_sample).c_str());
  }
  return Form("p_{T} > %s, |#eta| < %s", xjjc::number_remove_zero(ptmin_sample).c_str(), xjjc::number_remove_zero(etamax_sample).c_str());
}

void style_hist(Hist& h, xjjroot::thgrstyle s) {
  xjjroot::setthgrstyle(h.h_gamma, s);
  xjjroot::setthgrstyle(h.h_gamma_gap, s);
  xjjroot::setthgrstyle(h.h_gamma_gaploose, s);
  xjjroot::setthgrstyle(h.h_N, s);
}

int macro(std::string configfile, std::vector<std::string> vars) {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto name = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);
  // std::cout<<"\e[1m[ "<<name<<" ]\e[0m"<<std::endl;
  auto labels = conf.get_vec("Label_2");
  const auto& isgammaN = conf.get<int>("IsgammaN");
  const std::string& v_gamma = isgammaN ? "Plus" : "Minus",
    v_N = isgammaN ? "Minus" : "Plus";
  const auto& cut = conf.has("Cut") ? conf.get("Cut") : "(1)";
  const auto &PFptmin = conf.get<float>("PFptmin_sample"),
    &PFetamax = conf.get<float>("PFetamax_sample");

  std::map<std::string, Hist> h1s;
  for (const auto& v : vars) {
    auto* inf = TFile::Open(Form("rootfiles/%s%s.root", name.c_str(), v.c_str()));
    if (!inf) continue;
    Hist hs = {
      .h_gamma = xjjana::getobj<TH1F>(inf, "hHFEMax_gamma"),
      .h_gamma_gap = xjjana::getobj<TH1F>(inf, "hHFEMax_gamma_gap"),
      .h_gamma_gaploose = xjjana::getobj<TH1F>(inf, "hHFEMax_gamma_gaploose"),
      .h_N = xjjana::getobj<TH1F>(inf, "hHFEMax_N")
    };
    style_hist(hs, style_by_var(v));
    h1s[v] = hs; //
  }

  std::erase_if(vars, [&h1s](const auto& v) {
    return h1s.find(v) == h1s.end();
  });
  xjjc::print_vec_h(vars);

  auto* leg = new TLegend(0.52, 0.85-0.043*(h1s.size()+1), 0.8, 0.85);
  leg->SetHeader("Leading tower in PF", "L");
  xjjroot::setleg(leg, 0.04);
  for (const auto& v : vars) {
    leg->AddEntry(h1s.at(v).h_gamma, leg_by_var(v, PFptmin, PFetamax).c_str(), (h1s.at(v).h_gamma->GetFillStyle()?"f":"l"));
  }

  xjjroot::setgstyle(1);
  TGaxis::SetMaxDigits(4);
  TGaxis::SetExponentOffset(-0.08, 0., "y");
  auto* pdf = new xjjroot::mypdf("figspdf/" + name + ".pdf");

  TH1F* hempty = nullptr;
  std::vector<TH1F*> vhs;

#define DRAW_HIST(l)                                                    \
  vhs.clear();                                                          \
  for (const auto& v : vars) {                                          \
    vhs.push_back(h1s.at(v).h##l);                                      \
  }                                                                     \
  xjjana::sethsmax(vhs, 1.2);                                           \
  hempty = (TH1F*)h1s.at("_forest").h##l->Clone();                      \
  xjjroot::sethempty(hempty, 0);                                        \
  hempty->GetXaxis()->SetTitle(xjjc::str_erasestar(hempty->GetXaxis()->GetTitle(), "_*").c_str()); \
  pdf->prepare();                                                       \
  hempty->Draw("axis");                                                 \
  for (const auto& v : vars) {                                          \
    h1s.at(v).h##l->Draw("hist same");                                  \
  }                                                                     \
  xjjroot::drawCMS(xjjroot::CMS::internal, conf.get("Label_1"));        \
  leg->Draw();

  DRAW_HIST(_gamma);
  pdf->write();
  DRAW_HIST(_gamma_gap);
  print_ratio(vhs);
  pdf->write();
  DRAW_HIST(_gamma_gaploose);
  print_ratio(vhs);
  pdf->write();
  DRAW_HIST(_N);
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
  if (argc == 2)
    return macro(argv[1], {"_forest", "", "_eta5"});
  if (argc == 3)
    return macro(argv[1], xjjc::str_divide_trim(argv[2], ","));

  return 1;
}

