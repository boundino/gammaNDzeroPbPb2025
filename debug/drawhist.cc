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

namespace globals {
  std::vector<std::pair<std::string, std::vector<std::vector<float>>>> cuts = {
    { "Dalpha", {
        { 0.2, 0.4, 0.4, 0.2 },
      }},
    { "Dchi2cl", {
        { 0.1, 0.1, 0.1, 0.1 },
      }},
    { "Ddtheta", {
        { 0.3, 0.5, 0.5, 0.3 },
      }},
    { "dls", {
        { 2.5, 2.5, 2.5, 2.5 },
      }},
    { "Dtrk1Pt", {
        { 1, 1, 1, 1 },
      }},
    { "Dtrk2Pt", {
        { 1, 1, 1, 1 },
      }},
    // { "", {
    //     { , , , },
    //   }},
  };
}

int idx_varcut(std::string var,
               const std::vector<std::pair<std::string, std::vector<std::vector<float>>>>& cuts = globals::cuts) {
  for (int i=0; i<cuts.size(); i++) {
    if (cuts[i].first == var) {
      return i;
    }
  }
  return -1;
}

int macro(std::string configfiles, std::string var) {

  std::vector<int> colors = {
    xjjroot::mycolor_satmiddle2["blue"],
    xjjroot::mycolor_satmiddle["green"],
    xjjroot::mycolor_satmiddle2["red"],
  };
  auto idx_cut = idx_varcut(var, globals::cuts);

  std::string outputname, label_1;
  std::vector<std::string> comments;
  std::vector<Input> ips;
  for (const auto& c : xjjc::str_divide_trim(configfiles, ",")) {
    std::cout<<std::endl<<c<<std::endl;
    auto* conf = new xjjc::config(c);
    auto name = conf->has("Name") ? conf->get("Name") : xjjc::str_gettag_from_file(c);
    if (label_1.empty()) label_1 = conf->get("Label_1");
    auto* h3 = xjjana::getobj<TH3F>("rootfiles/" + var + "_" + name + ".root::h3");
    if (!h3) { std::cout<<"warning: no hist h3 in "<<c<<". skip..."<<std::endl; continue; }
    h3->SetName(Form("h3_%s", name.c_str()));
    auto nx = h3->GetXaxis()->GetNbins(), nz = h3->GetZaxis()->GetNbins();

    xjjc::vec_append(comments, { "#bf{"+conf->get("Latex")+"}", conf->get("Inputs"), conf->get("Cut"), "Lumi = " + conf->get("Lumi") });
    
    Input ip = {
      .h3 = h3,
      .h1var = xjjc::array2d<TH1F*>(nx, nz),
      .h1var_norm = xjjc::array2d<TH1F*>(nx, nz),
      .h1var_lumi = xjjc::array2d<TH1F*>(nx, nz),
      .latex = conf->get("Latex"),
      .lumi = (float)std::atof(conf->get("Lumi").c_str()),
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

        xjjroot::setthgrstyle(ip.h1var[i][j], colors[ips.size()], 20, 1.2, colors[ips.size()], 1, 2);
        xjjroot::setthgrstyle(ip.h1var_norm[i][j], colors[ips.size()], 20, 1.2, colors[ips.size()], 1, 2);
        xjjroot::setthgrstyle(ip.h1var_lumi[i][j], colors[ips.size()], 20, 1.2, colors[ips.size()], 1, 2);
      }
    }

    if (outputname.empty()) outputname += (var + "__" + name);
    else outputname += ("-" + name);

    ips.push_back(ip);
  }

  if (ips.empty()) { std::cout<<"error: no valid inputs."<<std::endl; return 2; }
  
  auto* leg = new TLegend(0.55, 0.72-0.045*ips.size(), 0.89, 0.72);
  xjjroot::setleg(leg, 0.042);
  auto* leg_wlumi = new TLegend(0.55, 0.72-0.045*ips.size(), 0.89, 0.72);
  xjjroot::setleg(leg_wlumi, 0.042);
  for (const auto& ip : ips) {
    leg->AddEntry(ip.h1var.front().front(), ip.latex.c_str(), "pl");
    leg_wlumi->AddEntry(ip.h1var.front().front(), Form("%s#scale[0.3]{ }#scale[0.7]{#it{%s/#mub}}", ip.latex.c_str(), xjjc::number_remove_zero(ip.lumi).c_str()), "pl");
  }
  
  auto ip_dump = ips.front();
  auto nx = ip_dump.xaxis.size(), nz = ip_dump.zaxis.size();
  std::cout<<std::endl;
  xjjc::print_tab({
      { "y (nx)", Form("%d", nx) },
      { "p_T (nz)", Form("%d", nz) },
    });

  xjjroot::setgstyle(1);
  gStyle->SetLineScalePS(1);

  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf", "c", 600*nx, 600*nz);
  std::vector<TH2F*> trash;

#define DRAWCOMP(l, t, le, order)                                       \
  pdf->prepare();                                                       \
  pdf->getc()->Divide(nx, nz);                                          \
  for (int i=0; i<nx; i++) {                                            \
    for (int j=0; j<nz; j++) {                                          \
      pdf->getc()->cd(order + j*nz);                                    \
      std::vector<TH1F*> h1s;                                           \
      double ymax = 0;                                                  \
      for (const auto& ip : ips) ymax = std::max(ymax, ip.h1##l[i][j]->GetMaximum()); \
      auto* hempty = new TH2F(Form("hempty_%d_%d", i, j), Form(";%s;Entries%s", ip_dump.h3->GetYaxis()->GetTitle(), t), \
                              10, ip_dump.h1##l[i][j]->GetXaxis()->GetXmin(), ip_dump.h1##l[i][j]->GetXaxis()->GetXmax(), \
                              10, 0, ymax*1.25);                        \
      xjjroot::sethempty(hempty);                                       \
      hempty->Draw("axis");                                             \
      if (idx_cut >= 0) {                                               \
        xjjroot::drawline_vertical(globals::cuts[idx_cut].second[j][i], \
                                   hempty, kBlack, 2, 2, 1.);           \
      }                                                                 \
      for (auto& ip : ips)                                              \
        ip.h1##l[i][j]->Draw("ple same");                               \
                                                                        \
      xjjroot::drawtexgroup(0.88, 0.80, {                               \
          xjjc::number_range_string(ip_dump.zaxis[j].first, ip_dump.zaxis[j].second, "p_{T}", (float)0, (float)1.e+3, "GeV"), \
          xjjc::number_range_string(ip_dump.xaxis[i].first, ip_dump.xaxis[i].second, "y", (float)-1.e+3, (float)1.e+3), \
        }, 0.042, 33, 42, 1.25);                                        \
      le->Draw();                                                       \
      xjjroot::drawCMS(xjjroot::CMS::internal, label_1);                \
                                                                        \
      trash.push_back(hempty);                                          \
    }                                                                   \
  }                                                                     \
  pdf->write();                                                         \
  for (auto& h : trash) delete h;                                       \
  trash.clear();


  pdf->draw_cover(comments, 0.03);
  
  DRAWCOMP(var, "", leg_wlumi, i+1);
  DRAWCOMP(var_lumi, " / lumi", leg, i+1);
  DRAWCOMP(var_norm, " (self-normalized)", leg, i+1);

  pdf->draw_cover({ "#bf{Reverse rapidity bins}" }, 0.05);
  
  DRAWCOMP(var, "", leg, nx-i);
  DRAWCOMP(var_lumi, " / lumi", leg, nx-i);
  DRAWCOMP(var_norm, " (self-normalized)", leg, nx-i);

  pdf->close();

  return 0;
}
       
int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
