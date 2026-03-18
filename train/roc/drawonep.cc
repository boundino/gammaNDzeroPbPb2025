#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

#include "mvaroot.h"
// #include "save_calc.h"
#include "../../include/dfitter.h"

void style_hist(TH1D*);
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

  std::vector<std::string> methods;
  std::map<std::string, Color_t> colors;
  for (int i=0; i<rmva.size(); i++) {
    for (auto& h : rmva.at(i)->hrocs()) {
      auto method = xjjc::str_eraseall(xjjc::str_eraseall(h->GetName(), "MVA_"), "_rejBvsS");
      if (std::ranges::find(methods, method) == methods.end()) {
        methods.push_back(method);
        colors[method] = xjjroot::colorlist_middle.at(methods.size()-1);
      }
      xjjroot::setthgrstyle(h, colors.at(method), 21, 1, colors.at(method), 1, 4);
    }
  }
  
  std::map<std::string, std::vector<TH1D*>> h1ys;
  for (const std::string &t : { "mass_Mass", "mass_S", "mass_Swap" }) {
    h1ys[t] = xjjana::getobj_regexp<TH1D>(inf, "h1_" + t + "__y-[0-9]*");
    for (auto &h : h1ys[t]) { style_hist(h); }
  }
  xjjroot::print_tab(h1ys, 0);

  std::map<std::string, std::vector<TGraph*>> grys;
  for (const std::string &t : { "roc_cut" }) {
    grys[t] = xjjana::getobj_regexp<TGraph>(inf, "gr_" + t + "__y-[0-9]*");
  }
  xjjroot::print_tab(grys, 0);

  auto* h3_dump = static_cast<TH3F*>(xjjana::getobj_regexp<TH3F>(inf, "h3_.+").front()->Clone("h3_dump")); h3_dump->Reset();
  auto label_y = [&h3_dump](int i) {
    const auto ymin = h3_dump->GetXaxis()->GetBinLowEdge(i+1), ymax = h3_dump->GetXaxis()->GetBinUpEdge(i+1);
    return xjjc::number_range_string(ymin, ymax, "y", -1.e1);
  };
  
  auto* df = new xjjroot::dfitter("S3");
  auto* hempty = new TH2F("hempty", "", 10, 0, 1, 10, 0, 1);
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf(xjjc::str_replaceall(xjjc::str_replaceall(inputname, "rootfiles/", "figspdf/"), ".root", ".pdf"));

  for (int i=0; i<h1ys.at("mass_Mass").size(); i++) {
    const auto& h = h1ys.at("mass_Mass").at(i),
      hmc = h1ys.at("mass_S").at(i), hmcswap = h1ys.at("mass_Swap").at(i);
    pdf->prepare();
    df->fit(h, hmc, hmcswap);
    xjjroot::drawtexgroup(0.25, 0.86, { label_y(i), "#bf{" + info.at("labelcut") + "}" }, 0.035, 13);
    df->draw_leg();
    df->draw_result(0.25, 0.86-2*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
    pdf->write();

    pdf->prepare();
    hmc->Draw("pe");
    hmcswap->Draw("pe same");
    df->draw_fmc();
    xjjroot::drawtex(0.25, 0.86, label_y(i).c_str(), 0.035, 13);
    df->draw_params(0.25, 0.86-1*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::simulation, "P#scale[0.9]{YTHIA}8#scale[0.5]{ }#gammaN (5.36 TeV)");
    pdf->write();
  }

  auto* leg = new TLegend(0.25, 0.5-0.04*(methods.size()+1), 0.5, 0.5);
  xjjroot::setleg(leg, 0.038);
  for (const auto& [mm, cc] : colors) {
    xjjroot::addentrybystyle(leg, mm, "l", xjjroot::thgrstyle{ .lcolor = cc, .lstyle = 1, .lwidth = 4 });
  }
  xjjroot::addentrybystyle(leg, info.at("labelcut"), "p", xjjroot::thgrstyle{ .mcolor = grys.at("roc_cut").front()->GetMarkerColor(), .mstyle = grys.at("roc_cut").front()->GetMarkerStyle(), .msize = grys.at("roc_cut").front()->GetMarkerSize() });
  auto tex = xjjroot::drawtex(0.25, 0.51, info.at("label2").c_str(), 0.038, 11);

  delete hempty; hempty = new TH2F("hempty", ";Signal Efficiency;Background Rejection", 10, 0, 1, 10, 0, 1);
  xjjroot::sethempty(hempty);
  for (int i=0; i<rmva.size(); i++) {
    pdf->prepare();
    hempty->Draw("axis");
    for (auto& h : rmva.at(i)->hrocs()) {
      h->Draw("c same");
    }
    grys.at("roc_cut").at(i)->Draw("p same");
    leg->Draw();
    xjjroot::drawtexgroup(0.25, 0.51, { label_y(i) }, 0.038, 11);
    xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
    pdf->write();
  } //    
  
  delete hempty; hempty = new TH2F("hempty", ";Signal Efficiency;Background Rejection", 10, 0, 0.8, 10, 0.85, 1);
  xjjroot::sethempty(hempty);
  for (int i=0; i<rmva.size(); i++) {
    pdf->prepare();
    hempty->Draw("axis");
    for (auto& h : rmva.at(i)->hrocs()) {
      h->Draw("c same");
    }
    grys.at("roc_cut").at(i)->Draw("p same");
    leg->Draw();
    xjjroot::drawtexgroup(0.25, 0.51, { label_y(i) }, 0.038, 11);
    xjjroot::drawCMS(xjjroot::CMS::internal, info.at("label"));
    pdf->write();
  } //
  
  pdf->draw_cover( {
      "#bf{" + info.at("labelcut") + "}",
      "#bf{Working point} " + info.at("cutOnep"),
      "#bf{Base cut} " + info.at("cutMass"),
      "#bf{Sample} " + info.at("inputMassname"),
    }, 0.03);

  pdf->draw_cover( {
      "#bf{Signal sample} " + info.at("inputSname"),
      "#bf{Signal cut} " + info.at("cutS"),
      "#bf{Background sample} " + info.at("inputBname"),
      "#bf{Background cut} " + info.at("cutB"),
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
