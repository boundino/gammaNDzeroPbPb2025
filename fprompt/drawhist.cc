#include <TH3F.h>

#include "config.h"
#include "xjjanauti.h"
#include "xjjcuti.h"
#include "xjjmypdf.h"

int macro(std::string configfile) {
  xjjc::config conf(configfile); std::cout<<std::endl;
  auto label_1 = conf.has("Label_1") ? conf.get("Label_1") : "";
  auto outputname = conf.has("Output") ? conf.get("Output") : xjjc::str_gettag_from_file(configfile);
  auto inputname = "rootfiles/" + outputname + "_project.root";
  auto* inf = TFile::Open(inputname.c_str());

#define READHIST(l)                                                     \
  auto* h3_##l = xjjana::getobj<TH3F>(inf, "h3_" #l);                   \
  if (!h3_##l) return 2;                                                \
  h3_##l->Sumw2();                                                      \
  auto nx_##l = h3_##l->GetXaxis()->GetNbins(), nz_##l = h3_##l->GetZaxis()->GetNbins(); \
  auto h1dca_##l = xjjc::array2d<TH1F*>(nx_##l, nz_##l);                \
  for (int i=0; i<nx_##l; i++) {                                        \
    for (int j=0; j<nz_##l; j++) {                                      \
      h1dca_##l[i][j] = (TH1F*)h3_##l->ProjectionY(Form("h1dca_" #l "_%d_%d", i, j), \
                                                   i+1, i+1,            \
                                                   j+1, j+1);           \
      h1dca_##l[i][j]->Scale(1./h1dca_##l[i][j]->Integral(), "width");  \
    }                                                                   \
  }                                                                     \

  READHIST(mc_prompt);
  for (auto& hh : h1dca_mc_prompt) {
    for (auto& h : hh)
      xjjroot::setthgrstyle(h, xjjroot::mycolor_satmiddle2["red"], 20, 2., xjjroot::mycolor_satmiddle2["red"], 1, 2);
  }
  READHIST(mc_nonprompt);
  for (auto& hh : h1dca_mc_nonprompt) {
    for (auto& h : hh)
      xjjroot::setthgrstyle(h, xjjroot::mycolor_satmiddle2["blue"], 20, 2., xjjroot::mycolor_satmiddle2["blue"], 1, 2);
  }

  auto* leg = new TLegend(0.60, 0.72-0.045*2, 0.87, 0.72);
  xjjroot::setleg(leg, 0.042);
  leg->AddEntry(h1dca_mc_prompt.front().front(), "Prompt", "pl");
  leg->AddEntry(h1dca_mc_nonprompt.front().front(), "Nonprompt", "pl");
  leg->Draw();
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  pdf->getc()->SetLogy();

  for (int i=0; i<nx_mc_prompt; i++) {
    for (int j=0; j<nz_mc_prompt; j++) {
      float xmin = h3_mc_prompt->GetXaxis()->GetBinLowEdge(i+1), xmax = h3_mc_prompt->GetXaxis()->GetBinUpEdge(i+1),
        zmin = h3_mc_prompt->GetZaxis()->GetBinLowEdge(j+1), zmax = h3_mc_prompt->GetZaxis()->GetBinUpEdge(j+1);
      auto* hempty = new TH2F(Form("hempty_%d_%d", i, j), ";DCA (cm); Self normalized / cm",
                              10, h1dca_mc_prompt[i][j]->GetXaxis()->GetXmin(), h1dca_mc_prompt[i][j]->GetXaxis()->GetXmax(),
                              10, h1dca_mc_nonprompt[i][j]->GetMinimum()/5., h1dca_mc_prompt[i][j]->GetMaximum()*2);
      xjjroot::sethempty(hempty);

      pdf->prepare();
      hempty->Draw("axis");
      h1dca_mc_nonprompt[i][j]->Draw("pe same");
      h1dca_mc_prompt[i][j]->Draw("pe same");

      xjjroot::drawtexgroup(0.88, 0.80, {
          xjjc::number_range_string(zmin, zmax, "p_{T}", (float)0, (float)1.e+3, "GeV"),
          xjjc::number_range_string(xmin, xmax, "y", (float)-1.e+3, (float)1.e+3),
        }, 0.042, 33, 42, 1.25);
      leg->Draw();
      xjjroot::drawCMS(xjjroot::CMS::simulation, label_1);
      pdf->write();
    }
  }
  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}
