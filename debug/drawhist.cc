#include <TH3F.h>
#include <TH1F.h>

#include "config.h"
#include "xjjanauti.h"
#include "xjjcuti.h"
#include "xjjmypdf.h"

struct Input {
  TH3F* h3;
  std::vector<std::vector<TH1F*>> h1var, h1var_norm, h1var_lumi;
  std::string latex;
  float lumi;
  std::vector<std::pair<float, float>> xaxis, zaxis;
};

int macro(std::string configfile) {
  xjjc::config conf(configfile); std::cout<<std::endl;
  // auto var = conf.get("Variable");
  auto inputfiles = conf.get_vec("Inputs");
  auto label_1 = conf.has("Label_1") ? conf.get("Label_1") : "";
  auto outputname = conf.has("Output") ? conf.get("Output") : xjjc::str_gettag_from_file(configfile);
  auto inputname = "rootfiles/" + outputname + ".root";

  auto* inf = TFile::Open(inputname.c_str());

  std::vector<int> colors = { xjjroot::mycolor_satmiddle2["red"], xjjroot::mycolor_satmiddle2["blue"] };
  std::vector<Input> ips;
  for (const auto& input : inputfiles) {
    std::cout<<input<<std::endl;
    auto ipar = xjjc::str_divide_trim(input, ":");
    if (ipar.size() < 5) { std::cout<<"warning: skip bad input "<<input<<std::endl; continue; }
    auto ip_name = ipar[1];
    auto h3 = xjjana::getobj<TH3F>(inf, Form("h3_%s", ip_name.c_str()));
    if (!h3) { std::cout<<"warning: no hist "<<Form("h3_%s", ip_name.c_str())<<". skip..."<<std::endl; continue; }
    auto nx = h3->GetXaxis()->GetNbins(), nz = h3->GetZaxis()->GetNbins();

    Input ip = {
      .h3 = h3,
      .h1var = xjjc::array2d<TH1F*>(nx, nz),
      .h1var_norm = xjjc::array2d<TH1F*>(nx, nz),
      .h1var_lumi = xjjc::array2d<TH1F*>(nx, nz),
      .latex = ipar[2],
      .lumi = (float)std::atof(ipar[3].c_str()),
    };

    ip.h3->Sumw2();
    for (int i = 0; i < nx; i++) {
      ip.xaxis.push_back( { ip.h3->GetXaxis()->GetBinLowEdge(i+1), ip.h3->GetXaxis()->GetBinUpEdge(i+1) } );
      for (int j = 0; j < nz; j++) {
        if (!i) ip.zaxis.push_back( { ip.h3->GetZaxis()->GetBinLowEdge(j+1), ip.h3->GetZaxis()->GetBinUpEdge(j+1) } );
        ip.h1var[i][j] = (TH1F*)ip.h3->ProjectionY(Form("h1var_%s_%d_%d", ip.h3->GetName(), i, j),
                                                   i+1, i+1,
                                                   j+1, j+1);
        ip.h1var_norm[i][j] = (TH1F*)ip.h1var[i][j]->Clone(Form("%s_norm", ip.h1var[i][j]->GetName()));
        ip.h1var_norm[i][j]->Scale(1./ip.h1var_norm[i][j]->Integral(), "width");
        ip.h1var_lumi[i][j] = (TH1F*)ip.h1var[i][j]->Clone(Form("%s_lumi", ip.h1var[i][j]->GetName()));
        ip.h1var_lumi[i][j]->Scale(1./ip.lumi, "width");

        xjjroot::setthgrstyle(ip.h1var[i][j], colors[ips.size()], 20, 1., colors[ips.size()], 1, 1);
        xjjroot::setthgrstyle(ip.h1var_norm[i][j], colors[ips.size()], 20, 1., colors[ips.size()], 1, 1);
        xjjroot::setthgrstyle(ip.h1var_lumi[i][j], colors[ips.size()], 20, 1., colors[ips.size()], 1, 1);
      }
    }
    
    ips.push_back(ip);
  }

  auto* leg = new TLegend(0.60, 0.72-0.045*ips.size(), 0.87, 0.72);
  xjjroot::setleg(leg, 0.042);
  for (const auto& ip : ips) {
    leg->AddEntry(ip.h1var.front().front(), ip.latex.c_str(), "pl");
  }
  
  auto nx = ips.front().xaxis.size(), nz = ips.front().zaxis.size();
  std::cout<<"nx "<<nx<<" ; nz "<<nz<<std::endl;
  auto ip_dump = ips.front();

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf", "c", 600*nx, 600*nz);

#define DRAWCOMP(l, t)                                                  \
  pdf->prepare();                                                       \
  pdf->getc()->Divide(nx, nz);                                          \
  for (int i=0; i<nx; i++) {                                            \
    for (int j=0; j<nz; j++) {                                          \
      pdf->getc()->cd(i+1+j*nz);                                        \
      std::vector<TH1F*> h1s;                                           \
      double ymax = 0;                                                  \
      for (const auto& ip : ips) ymax = std::max(ymax, ip.h1##l[i][j]->GetMaximum()); \
      auto* hempty##l = new TH2F(Form("hempty" #l  "_%d_%d", i, j), Form(";%s;Entries %s", ip_dump.h3->GetXaxis()->GetTitle(), #t), \
                                 10, ip_dump.h1##l[i][j]->GetXaxis()->GetXmin(), ip_dump.h1##l[i][j]->GetXaxis()->GetXmax(), \
                                 10, 0, ymax*1.25);                     \
      xjjroot::sethempty(hempty##l);                                    \
      hempty##l->Draw("axis");                                          \
                                                                        \
      for (auto& ip : ips)                                              \
        ip.h1##l[i][j]->Draw("pe same");                                \
                                                                        \
      xjjroot::drawtexgroup(0.88, 0.80, {                               \
          xjjc::number_range_string(ip_dump.zaxis[j].first, ip_dump.zaxis[j].second, "p_{T}", (float)0, (float)1.e+3, "GeV"), \
          xjjc::number_range_string(ip_dump.xaxis[i].first, ip_dump.xaxis[i].second, "y", (float)-1.e+3, (float)1.e+3), \
        }, 0.042, 33, 42, 1.25);                                        \
      leg->Draw();                                                      \
      xjjroot::drawCMS(xjjroot::CMS::internal, label_1);                \
    }                                                                   \
  }                                                                     \
  pdf->write();                                                         \


  DRAWCOMP(var,);
  DRAWCOMP(var_lumi, / lumi);
  DRAWCOMP(var_norm, Self-normalized);
  
  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}
