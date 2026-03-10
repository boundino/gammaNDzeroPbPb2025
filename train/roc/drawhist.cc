#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

#include "mvaroot.h"
#include "../../include/dfitter.h"

void style_hist(TH1D*);

int macro(std::string inputname) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  // train
  auto* rmva = new mytmva::mvaroot(inf);
  rmva->print_info();
  auto* h3_dump = (TH3F*)xjjana::getobj_regexp<TH3F>(inf, "h3_.+").front()->Clone("h3_dump"); h3_dump->Reset();
  const auto method = rmva->info("method"), mprefix = "MVA_" + method;
  for (auto& h : rmva->hrocs(".+"))
    style_hist(h);
  // application
  std::map<std::string, std::vector<TH1D*>> h1ys;
  for (const std::string &t : { "S", "B", "effS", "effB", "effvsMVA", "rejBvsS", "massAnchor", "massTemplate_S", "massTemplate_Swap" }) {
    h1ys[t] = xjjana::getobj_regexp<TH1D>(inf, "h1_" + t + "__y-[0-9]*");
    for (auto &h : h1ys[t]) { style_hist(h); }
  } 
  xjjroot::print_tab(h1ys, 0);
  std::map<std::string, TH1D*> h1s;
  for (const std::string &t : { "mvaAnchor", "S", "B", "S_high", "B_high" }) {
    h1s[t] = xjjana::getobj<TH1D>(inf, "h1_" + t);
    style_hist(h1s[t]);
  }
  xjjroot::print_tab(h1s, 0);

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
  auto* leg = new TLegend(0.25, 0.5-0.04*(h1ys.at("rejBvsS").size()+1), 0.5, 0.5);
  leg->AddEntry(rmva->hroc(".+_rejBvsS"), "Train", "l");
  leg_y(leg, h1ys.at("rejBvsS"), "l");

  auto* hempty = new TH2F("hempty", "", 10, 0, 1, 10, 0, 1);
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf(xjjc::str_replaceall(xjjc::str_replaceall(inputname, "rootfiles/", "figspdf/"), ".root", ".pdf"));

  h1s.at("S")->Scale(1./h1s.at("S")->Integral(), "width");
  h1s.at("B")->Scale(1./h1s.at("B")->Integral(), "width");
  auto hmvas = std::vector<std::pair<std::string, TH1D*>>({
      { "Train", rmva->hroc(mprefix + "_S") },
      { "Application", h1s.at("S") },
      { "Train", rmva->hroc(mprefix + "_B") },
      { "Application", h1s.at("B") },
    });
  auto hymax = xjjana::sethsmax(hmvas);
  delete hempty; hempty = new TH2F("hempty", Form(";%s;Distribution", method.c_str()), 10, rmva->hroc(mprefix + "_S")->GetXaxis()->GetXmin(), rmva->hroc(mprefix + "_S")->GetXaxis()->GetXmax(), 10, 0, hymax*1.3);
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
  xjjroot::drawtexgroup(0.47, 0.85, { "Signal", "Background" }, 0.038, 33, 42, 1.15, 1, 0, { h1s.at("S")->GetLineColor(), h1s.at("B")->GetLineColor() });
  xjjroot::drawCMS(xjjroot::CMS::internal, rmva->info("label"));
  pdf->write();  
  
  delete hempty; hempty = new TH2F("hempty", ";Signal Efficiency;Background Rejection", 10, 0, 1, 10, 0, 1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");
  for (int i=0; i<h1ys.at("rejBvsS").size(); i++) {
    h1ys.at("rejBvsS").at(i)->Draw("c same");
  }
  rmva->hroc(".+_rejBvsS")->Draw("c same");
  auto tex_label2 = xjjroot::drawtex(0.25, 0.51, rmva->info("label2").c_str(), 0.038, 11);
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, rmva->info("label"));
  pdf->write();

  delete hempty; hempty = new TH2F("hempty", ";Signal Efficiency;Background Rejection", 10, 0, 0.6, 10, 0.9, 1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");
  for (const auto& h : h1ys.at("rejBvsS"))
    h->Draw("c same");
  rmva->hroc(".+_rejBvsS")->Draw("c same");
  tex_label2->Draw();
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, rmva->info("label"));
  pdf->write();

  // significance
  const auto scaleto25 = std::atof(rmva->info("scaleto25").c_str());
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
    xjjroot::drawCMS(xjjroot::CMS::internal, rmva->info("label"));
    pdf->write();

    if (h1_effS->GetBinContent(jmva_anchor) <= 0 || h1_effB->GetBinContent(jmva_anchor) <= 0) {
      __XJJLOG << "!! anchor "<<i<<"; jmva = "<<jmva_anchor<<"; effS = "<<h1_effS->GetBinContent(jmva_anchor)<<"; effB = "<<h1_effB->GetBinContent(jmva_anchor)<<std::endl;
      return 3;
    }
    auto S0 = df->S() / h1_effS->GetBinContent(jmva_anchor) * scaleto25,
      B0 = df->B() / h1_effB->GetBinContent(jmva_anchor) * scaleto25;

    auto* hsig = (TH1D*)h1_effS->Clone(Form("h1_sig25__y-%d", i)); hsig->Reset();
    hsig->SetTitle(Form(";%s;Expected S /#scale[0.5]{ }#sqrt{S+B} in 2025", method.c_str()));
    for (int j=0; j<hsig->GetXaxis()->GetNbins(); j++) {
      if (j%100 == 0 && (h1_effS->GetBinContent(j+1)>0 || h1_effB->GetBinContent(j+1)>0))
        hsig->SetBinContent(j+1, S0*h1_effS->GetBinContent(j+1) / TMath::Sqrt(S0*h1_effS->GetBinContent(j+1) + B0*h1_effB->GetBinContent(j+1)));
    }
    hsig->SetMinimum(0);
    xjjroot::sethempty(hsig);
    style_hist(hsig);
    h1ys["sig25"].push_back(hsig); //

    pdf->prepare();
    hmc->Draw("pe");
    hmcswap->Draw("pe same");
    df->draw_fmc();
    xjjroot::drawtex(0.25, 0.86, label_y(i).c_str(), 0.035, 13);
    df->draw_params(0.25, 0.86-1*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::simulation, "P#scale[0.9]{YTHIA}8#scale[0.5]{ }#gammaN (5.02 TeV)");
    pdf->write();
  }

  auto* legp = new TLegend(0.25, 0.80-0.04*h1ys.at("sig25").size(), 0.5, 0.8);
  leg_y(legp, h1ys.at("sig25"), "p");
  auto ymax = xjjana::sethsmax(h1ys.at("sig25"), 1.3);
  pdf->prepare();
  h1ys.at("sig25").front()->Draw("axis");
  for (const auto& h : h1ys.at("sig25")) {
    h->Draw("p same");
  }
  legp->Draw();
  xjjroot::drawtex(0.25, 0.81, rmva->info("label2").c_str(), 0.038, 11);
  xjjroot::drawCMS(xjjroot::CMS::internal, rmva->info("label"));
  pdf->write();
  
  delete hempty; hempty = new TH2F("hempty", Form(";%s;Signal Efficiency", method.c_str()), 10, h1s["S_high"]->GetXaxis()->GetXmin(), h1s["S_high"]->GetXaxis()->GetXmax(), 10, 0, 1.1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");  
  for (const auto& h : h1ys.at("effS"))
    h->Draw("c same");
  rmva->hroc(".+_effS")->Draw("c same");
  tex_label2->Draw();
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, rmva->info("label"));
  pdf->write();

  delete hempty; hempty = new TH2F("hempty", Form(";%s;Background Efficiency", method.c_str()), 10, h1s["S_high"]->GetXaxis()->GetXmin(), h1s["S_high"]->GetXaxis()->GetXmax(), 10, 0, 1.1);
  xjjroot::sethempty(hempty);
  pdf->prepare();
  hempty->Draw("axis");  
  for (const auto& h : h1ys.at("effB"))
    h->Draw("c same");
  rmva->hroc(".+_effB")->Draw("c same");
  xjjroot::movetex_n_draw(tex_label2, 0.55, 0.81, 11);
  xjjroot::moveleg_n_draw(leg, 0.55, 0.8);
  xjjroot::drawCMS(xjjroot::CMS::internal, rmva->info("label"));
  pdf->write();

  pdf->draw_cover( {
      "#bf{Application files} " + rmva->info("inputMname"),
      "#bf{Scale to 25 stats} " + rmva->info("scaleto25"),
      "#bf{Background cut} " + rmva->info("cutb"),
      "#bf{Background sample} " + rmva->info("inputBname"),
      "#bf{Signal cut} " + rmva->info("cuts"),
      "#bf{Signal sample} " + rmva->info("inputSname"),
      "#bf{Training output} " + rmva->info("rootfile"),
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

void style_hist(TH1D* h) {
  std::string t = xjjc::str_eraseall(h->GetName(), "MVA_");
  if (xjjc::str_contains(t, "mass")) {
    xjjroot::dfitter::set_hist(h);
  } else if (std::regex_match(t, std::regex(".+__y-[0-9]+")) &&
             (xjjc::str_contains(t, "rejBvsS") || xjjc::str_contains(t, "eff") || xjjc::str_contains(t, "sig"))) {
    int i = std::atoi(xjjc::str_erasestar(t, "*__y-").c_str());
    xjjroot::setthgrstyle(h, xjjroot::colorlist_middle.at(i), 53, 1.2, xjjroot::colorlist_middle.at(i), 1, 4);
  } else if (xjjc::str_contains(t, "_S") && !xjjc::str_contains(t, "_B")) {
    xjjroot::setthgrstyle(h, xjjroot::mycolor_satmiddle.at("blue"), 20, 1.3, xjjroot::mycolor_satmiddle.at("blue"), 1, 1, xjjroot::mycolor_satmiddle.at("blue"), 0.5, 1001);
  } else if (xjjc::str_contains(t, "_B") && !xjjc::str_contains(t, "_S")) {
    xjjroot::setthgrstyle(h, xjjroot::mycolor_satmiddle2.at("red"), 20, 1.3, xjjroot::mycolor_satmiddle2.at("red"), 1, 1, xjjroot::mycolor_satmiddle2.at("red"), 1, 3354);
  } else if (xjjc::str_contains(h->GetName(), "MVA_")) {
    xjjroot::setlinestyle(h, kBlack, 4, 4);
  }
}
