#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

#include "mvaroot.h"
#include "../../include/dfitter.h"

template<class T> void style_hist(T*);
TGraph* coarse_bin(TH1D* h, int step = 100, const std::string &name = "");
int macro(std::string inputname) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  if (!inf) {
    __XJJLOG << "!! bad input file, abort." << std::endl; return 2;
  }
  auto rmva = mytmva::read_rmvas(inf, "pt-[0-9]+-[0-9]+_y-[M0-9]+-[M0-9]+");
  if (rmva.empty()) {
    __XJJLOG << "!! no good mva root files, abort." << std::endl;
    return 2;
  }
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  xjjc::print_tab(info, -1);
  const auto method = info.at("method"), mprefix = "MVA_" + method;
  
  auto* h3_dump = static_cast<TH3F*>(xjjana::getobj_regexp<TH3F>(inf, "h3_.+").front()->Clone("h3_dump")); h3_dump->Reset();
  for (int i=0; i<rmva.size(); i++) {
    for (auto& h : rmva.at(i)->hrocs(".+")) {
      style_hist(h);
      if (rmva.size() > 1 && xjjc::str_contains(h->GetName(), "rej"))
        xjjroot::setthgrstyle(h, xjjroot::colorlist_middle.at(i), -1, -1, xjjroot::colorlist_middle.at(i), -1, -1);
    }
  }

  // application
  std::map<std::string, std::vector<TH1D*>> h1ys;
  for (const std::string &t : { "S_trainbin", "B_trainbin", "effS", "effB", "effvsMVA", "rejBvsS", "massAnchor", "massTemplate_S", "massTemplate_Swap" }) {
    h1ys[t] = xjjana::getobj_regexp<TH1D>(inf, "h1_" + t + "__y-[0-9]*");
    for (auto &h : h1ys[t]) { style_hist(h); }
  } 
  xjjroot::print_tab(h1ys, 0);
  std::map<std::string, TH1D*> h1s;
  for (const std::string &t : { "mvaAnchor" }) {
    h1s[t] = xjjana::getobj<TH1D>(inf, "h1_" + t);
    style_hist(h1s[t]);
  }
  xjjroot::print_tab(h1s, 0);
  std::map<std::string, std::vector<TGraph*>> grys;

  auto* df = new xjjroot::dfitter("S3");

  auto label_y = [&h3_dump](int i) {
    const auto ymin = h3_dump->GetXaxis()->GetBinLowEdge(i+1), ymax = h3_dump->GetXaxis()->GetBinUpEdge(i+1);
    return xjjc::number_range_string(ymin, ymax, "y", -1.e1);
  };
  auto leg_y = [&label_y](TLegend* leg, const std::vector<TH1D*>& h1legs, const std::string& opt) {
    xjjroot::setleg(leg, 0.038);
    for (int i=0; i<h1legs.size(); i++) {
      leg->AddEntry(h1legs.at(i), label_y(i).c_str(), opt.c_str());
    }
  };
  auto* leg = new TLegend(0.25, 0.5-0.04*h1ys.at("rejBvsS").size(), 0.5, 0.5);
  // leg->AddEntry(rmva->hroc(".+_rejBvsS"), "Train", "l");
  leg_y(leg, h1ys.at("rejBvsS"), "l");

  auto* legr = new TLegend(0.25, 0.5-0.04*rmva.size(), 0.5, 0.5);
  xjjroot::setleg(legr, 0.038);
  for (const auto& rr : rmva) {
    legr->AddEntry(rr->hroc(".+_rejBvsS"), xjjc::number_range_string(rr->info<float>("ymin"), rr->info<float>("ymax"), "y", (float)-1.e1).c_str(), "l");
  }

  auto* hempty = new TH2F("hempty", "", 10, 0, 1, 10, 0, 1);
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf(xjjc::str_replaceall(xjjc::str_replaceall(inputname, "rootfiles/", "figspdf/"), ".root", ".pdf"));

  if (h1ys.at("S_trainbin").size() != rmva.size()) {
    __XJJLOG << "!! application trainbin hists have different binning number with mva size" << std::endl;
  } else {
    for (int i=0; i<h1ys.at("S_trainbin").size(); i++) {
      h1ys.at("S_trainbin")[i]->Scale(1./h1ys.at("S_trainbin")[i]->Integral(), "width");
      h1ys.at("B_trainbin")[i]->Scale(1./h1ys.at("B_trainbin")[i]->Integral(), "width");
      auto hmvas = std::vector<std::pair<std::string, TH1D*>>({
          { "Train", rmva[i]->hroc(mprefix + "_S") },
          { "Application", h1ys.at("S_trainbin")[i] },
          { "Train", rmva[i]->hroc(mprefix + "_B") },
          { "Application", h1ys.at("B_trainbin")[i] },
        });
      auto hymax = xjjana::sethsmax(hmvas);
      delete hempty; hempty = new TH2F("hempty", Form(";%s;Distribution", method.c_str()), 10, rmva[i]->hroc(mprefix + "_S")->GetXaxis()->GetXmin(), rmva[i]->hroc(mprefix + "_S")->GetXaxis()->GetXmax(), 10, 0, hymax*1.3);
      xjjroot::sethempty(hempty);
      auto* legm = new TLegend(0.5, 0.85-0.04*2, 0.85, 0.85);
      legm->SetNColumns(2);
      xjjroot::setleg(legm, 0.038);
      pdf->prepare();
      hempty->Draw("axis");
      for (const auto& [t, h] : hmvas) {
        h->Draw(Form("%s same", xjjc::str_contains(t, "rain") ? "hist" : "pe"));
        legm->AddEntry(h, t.c_str(), xjjc::str_contains(t, "rain") ? "f" : "p");
      }
      legm->Draw();
      xjjroot::drawtexgroup(0.47, 0.85, { "Signal", "Background" }, 0.038, 33, 42, 1.15, 1, 0, { h1ys.at("S_trainbin")[i]->GetLineColor(), h1ys.at("B_trainbin")[i]->GetLineColor() });
      xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
      xjjroot::drawtex(0.85, 0.75, xjjc::number_range_string(rmva[i]->info<float>("ymin"), rmva[i]->info<float>("ymax"), "y", (float)-1.e3).c_str(), 0.038, 33);
      pdf->write();
    }
  }
  
  delete hempty; hempty = new TH2F("hempty", ";Signal Efficiency;Background Rejection", 10, 0, 1, 10, 0, 1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");
  for (int i=0; i<h1ys.at("rejBvsS").size(); i++) {
    h1ys.at("rejBvsS").at(i)->Draw("c same");
  }
  auto tex_apply = xjjroot::drawtex(0.25, 0.51+0.04, "Apply curves", 0.038, 11);
  auto tex_label2 = xjjroot::drawtex(0.25, 0.51, info.at("label2").c_str(), 0.038, 11);
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  pdf->write();

  pdf->prepare();
  hempty->Draw("axis");
  for (auto& rr : rmva)
    rr->hroc(".+_rejBvsS")->Draw("c same");
  auto tex_train = xjjroot::drawtex(0.25, 0.51+0.04, "Training curves", 0.038, 11);
  tex_label2->Draw();
  legr->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  pdf->write();

  delete hempty; hempty = new TH2F("hempty", ";Signal Efficiency;Background Rejection", 10, 0, 0.6, 10, 0.9, 1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");
  for (const auto& h : h1ys.at("rejBvsS"))
    h->Draw("c same");
  tex_apply->Draw();
  tex_label2->Draw();
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  pdf->write();

  pdf->prepare();
  hempty->Draw("axis");
  for (auto& rr : rmva)
    rr->hroc(".+_rejBvsS")->Draw("c same");
  tex_train->Draw();
  tex_label2->Draw();
  legr->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  pdf->write();

  // significance
  const auto scaleto25 = std::atof(info.at("scaleto25").c_str());
  __XJJLOG << ">> scale to 25 = "<<scaleto25<<std::endl;
  for (int i=0; i<h1ys.at("massAnchor").size(); i++) {
    auto *h1_effS = h1ys["effS"].at(i), *h1_effB = h1ys["effB"].at(i);

    // Anchor point
    const auto jmva_anchor = h1_effS->FindBin(h1s["mvaAnchor"]->GetBinContent(i+1));
    const auto& h = h1ys.at("massAnchor").at(i),
      hmc = h1ys.at("massTemplate_S").at(i), hmcswap = h1ys.at("massTemplate_Swap").at(i);

    pdf->prepare();
    df->fit(h, hmc, hmcswap);
    xjjroot::drawtexgroup(0.25, 0.86, { label_y(i), Form("%s > %.2f", method.c_str(), h1_effS->GetBinLowEdge(jmva_anchor)) }, 0.035, 13);
    df->draw_leg();
    df->draw_result(0.25, 0.86-2*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
    pdf->write();

    if (h1_effS->GetBinContent(jmva_anchor) <= 0 || h1_effB->GetBinContent(jmva_anchor) <= 0) {
      __XJJLOG << "!! anchor "<<i<<"; jmva = "<<jmva_anchor<<"; effS = "<<h1_effS->GetBinContent(jmva_anchor)<<"; effB = "<<h1_effB->GetBinContent(jmva_anchor)<<std::endl;
      return 3;
    }
    auto S0 = df->S() / h1_effS->GetBinContent(jmva_anchor) * scaleto25,
      B0 = df->B() / h1_effB->GetBinContent(jmva_anchor) * scaleto25;

    auto* hsig25 = (TH1D*)h1_effS->Clone(Form("h1_sig25__y-%d", i)); hsig25->Reset();
    auto* hsig = (TH1D*)h1_effS->Clone(Form("h1_sig__y-%d", i)); hsig->Reset();
    hsig25->SetTitle(Form(";%s;Expected S /#scale[0.5]{ }#sqrt{S+B} in 2025", method.c_str()));
    hsig->SetTitle(Form(";%s;Expected S /#scale[0.5]{ }#sqrt{S+B} in Sample", method.c_str()));
    for (int j=0; j<hsig25->GetXaxis()->GetNbins(); j++) {
      if ((h1_effS->GetBinContent(j+1)>0 || h1_effB->GetBinContent(j+1)>0)) {
        hsig25->SetBinContent(j+1, S0*h1_effS->GetBinContent(j+1) / TMath::Sqrt(S0*h1_effS->GetBinContent(j+1) + B0*h1_effB->GetBinContent(j+1)));
        hsig->SetBinContent(j+1, S0/scaleto25*h1_effS->GetBinContent(j+1) / TMath::Sqrt(S0/scaleto25*h1_effS->GetBinContent(j+1) + B0/scaleto25*h1_effB->GetBinContent(j+1)));
      }
    }

    hsig25->SetMinimum(0);
    xjjroot::sethempty(hsig25);
    style_hist(hsig25);
    h1ys["sig25"].push_back(hsig25); //
    grys["sig25"].push_back(coarse_bin(hsig25, 100)); //
    hsig->SetMinimum(0);
    xjjroot::sethempty(hsig);
    style_hist(hsig);
    h1ys["sig"].push_back(hsig); //
    grys["sig"].push_back(coarse_bin(hsig, 100)); //

    pdf->prepare();
    hmc->Draw("pe");
    hmcswap->Draw("pe same");
    df->draw_fmc();
    xjjroot::drawtex(0.25, 0.86, label_y(i).c_str(), 0.035, 13);
    df->draw_params(0.25, 0.86-1*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::simulation, "P#scale[0.9]{YTHIA}8#scale[0.5]{ }#gammaN (5.02 TeV)");
    pdf->write();
  }
  xjjroot::print_tab(grys, 0);

  auto* legp = new TLegend(0.25, 0.80-0.04*h1ys.at("sig25").size(), 0.5, 0.8);
  leg_y(legp, h1ys.at("sig25"), "p");

  xjjana::sethsmax(h1ys.at("sig25"), 1.3);
  pdf->prepare();
  h1ys.at("sig25").front()->Draw("axis");
  std::map<std::string, std::vector<TF1*>> tfys;
  h1s["xopt"] = (TH1D*)h1s.at("mvaAnchor")->Clone("h1_xopt"); h1s.at("xopt")->Reset();
  for (int i=0; i<h1ys.at("sig25").size(); i++) {
    auto* &h = h1ys.at("sig25")[i];
    // h->Sumw2();
    auto jmax = h->GetMaximumBin();
    auto xmax = h->GetBinCenter(jmax), frange = (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin()) / 25.;
    __XJJLOG << ">> " << i << " -> xmax = " << xmax << ", frange = " << frange << std::endl;
    auto* f1 = new TF1(Form("f_sig25__y-%d", i), "gaus(0)", xmax-frange, xmax+frange);
    xjjroot::setthgrstyle(f1, -1, -1, -1, h->GetLineColor(), h->GetLineStyle(), h->GetLineWidth());
    f1->SetParameters(1, xmax, frange);
    grys.at("sig25")[i]->Draw("p same");
    grys.at("sig25")[i]->Fit(f1, "R q");
    tfys["sig25"].push_back(f1);
    auto xopt = f1->GetParameter(1);
    h1s.at("xopt")->SetBinContent(i+1, xopt);
    h1s.at("xopt")->SetBinError(i+1, f1->GetParError(1));
    __XJJLOG << "   ^ xopt = " << Form("%.4f", xopt) << std::endl;
    xjjroot::drawline(xopt, 0, xopt, f1->Eval(xopt), f1->GetLineColor(), 2, f1->GetLineWidth(), 0.5);
  }
  legp->Draw();
  xjjroot::drawtex(0.25, 0.81, info.at("label2").c_str(), 0.038, 11);
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  pdf->write();

  for (auto& f1 : tfys.at("sig25")) {
    auto* f1_2 = (TF1*)f1->Clone(xjjc::str_replaceall(f1->GetName(), "sig25", "sig").c_str());
    f1_2->SetParameters(f1->GetParameter(0)/TMath::Sqrt(scaleto25), f1->GetParameter(1), f1->GetParameter(2));
    tfys["sig"].push_back(f1_2);
  }
  xjjana::sethsmax(h1ys.at("sig"), 1.3);
  pdf->prepare();
  h1ys.at("sig").front()->Draw("axis");
  for (int i=0; i<h1ys.at("sig").size(); i++) {
    grys.at("sig")[i]->Draw("p same");
    auto* f1 =  tfys.at("sig")[i];
    f1->Draw("same"); 
    auto xopt = h1s.at("xopt")->GetBinContent(i+1);
    xjjroot::drawline(xopt, 0, xopt, f1->Eval(xopt), f1->GetLineColor(), 2, f1->GetLineWidth(), 0.5);
  }
  legp->Draw();
  xjjroot::drawtex(0.25, 0.81, info.at("label2").c_str(), 0.038, 11);
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  pdf->write();

  delete hempty; hempty = new TH2F("hempty", Form(";%s;Signal Efficiency", method.c_str()), 10, h1s["S_high"]->GetXaxis()->GetXmin(), h1s["S_high"]->GetXaxis()->GetXmax(), 10, 0, 1.1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");  
  for (const auto& h : h1ys.at("effS")) {
    auto* gr = coarse_bin(h);
    gr->Draw("c same");
    // rmva->hroc(".+_effS")->Draw("c same");
  }
  tex_label2->Draw();
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  pdf->write();

  // delete hempty; hempty = new TH2F("hempty", Form(";%s;Background Efficiency", method.c_str()), 10, h1s["S_high"]->GetXaxis()->GetXmin(), h1s["S_high"]->GetXaxis()->GetXmax(), 10, 0, 1.1);
  // xjjroot::sethempty(hempty);
  // pdf->prepare();
  // hempty->Draw("axis");  
  // for (const auto& h : h1ys.at("effB"))
  //   h->Draw("c same");
  // // rmva->hroc(".+_effB")->Draw("c same");
  // xjjroot::movetex_n_draw(tex_label2, 0.55, 0.81, 11);
  // xjjroot::moveleg_n_draw(leg, 0.55, 0.8);
  // xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
  // pdf->write();

  pdf->draw_cover( {
      "#bf{Mass files} " + info.at("inputMassname"),
      "#bf{Scale to 25 stats} " + info.at("scaleto25"),
      "#bf{Mass cut} " + info.at("cutMass"),
      "#bf{Background cut} " + info.at("cutB"),
      "#bf{Background sample} " + info.at("inputBname"),
      "#bf{Signal cut} " + info.at("cutS"),
      "#bf{Signal sample} " + info.at("inputSname"),
    }, 0.025);
  
  pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

template<class T>
void style_hist(T* h) {
  std::string t = xjjc::str_eraseall(h->GetName(), "MVA_");
  if (xjjc::str_contains(t, "mass")) {
    if constexpr (std::is_same_v<T, TH1D>) {
      xjjroot::dfitter::set_hist(h);
      return;
    }
  }
  if (std::regex_match(t, std::regex(".+__y-[0-9]+"))) {
    int i = std::atoi(xjjc::str_erasestar(t, "*__y-").c_str());
    if (xjjc::str_contains(t, "rejBvsS") || xjjc::str_contains(t, "eff") || xjjc::str_contains(t, "sig")) {
      xjjroot::setthgrstyle(h, xjjroot::colorlist_middle.at(i), 53, 1.2, xjjroot::colorlist_middle.at(i), 1, 4); return;
    }
  }
  if (xjjc::str_contains(t, "_S") && !xjjc::str_contains(t, "_B")) {
    xjjroot::setthgrstyle(h, xjjroot::mycolor_satmiddle.at("blue"), 20, 1.3, xjjroot::mycolor_satmiddle.at("blue"), 1, 1, xjjroot::mycolor_satmiddle.at("blue"), 0.5, 1001); return;
  } else if (xjjc::str_contains(t, "_B") && !xjjc::str_contains(t, "_S")) {
    xjjroot::setthgrstyle(h, xjjroot::mycolor_satmiddle2.at("red"), 20, 1.3, xjjroot::mycolor_satmiddle2.at("red"), 1, 1, xjjroot::mycolor_satmiddle2.at("red"), 1, 3354); return;
  }
  if (xjjc::str_contains(h->GetName(), "MVA_")) {
    xjjroot::setlinestyle(h, kBlack, 4, 4); return;
  }
}

TGraph* coarse_bin(TH1D* h, int step, const std::string& name) {
  auto *gr = new TGraph();
  if (name.empty()) 
    gr->SetName(xjjc::str_replaceall(h->GetName(), "h1_", "gr_").c_str());
  else
    gr->SetName(name.c_str());

  for (int j=0; j<h->GetXaxis()->GetNbins(); j++) {
    if (j%step == 0 && h->GetBinContent(j+1) > 0) {
      gr->SetPoint(gr->GetN(), h->GetBinCenter(j+1), h->GetBinContent(j+1));
    }
  }
  xjjroot::setthgrstyle(gr, h->GetMarkerColor(), h->GetMarkerStyle(), h->GetMarkerSize(),
                        h->GetLineColor(), h->GetLineStyle(), h->GetLineWidth());
  return gr;
}
