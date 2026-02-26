// https://root.cern.ch/doc/v610/TMVAGui_8cxx_source.html
// https://root.cern.ch/doc/v608/efficiencies_8cxx_source.html
#include <iostream>
#include <string>
#include <algorithm>

#include <TKey.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TH2F.h>
#include "TMVA/tmvaglob.h"

#include "TMVAClassification.h"
#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "xjjanauti.h"

namespace mytmva
{
  void guieffvar(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage);
  void effvar(std::vector<std::string> outfname, std::string method);
}

void mytmva::guieffvar(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage)
{
  TMVA::TMVAGlob::SetTMVAStyle();

  std::vector<std::string> outfname;
  while (xjjc::str_contains(stage, ",")) {
    outfname.push_back(mytmva::mkname(outputname, mymethod, stage, ptmin, ptmax));
    stage.erase(stage.rfind(','));
  }
  std::string outputstr(xjjc::str_tag_from_file(outfname.front()));
  std::sort(outfname.begin(), outfname.end(), xjjc::sortbydes<std::string>);

  auto methods = xjjc::str_divide_trim(mymethod, ",");
  for (auto& mm : methods) {
    mytmva::effvar(outfname, mm); 
    gSystem->Exec(Form("mkdir -p dataset/plots/%s", outputstr.c_str()));
    gSystem->Exec(Form("mv rejBvsSvar.png dataset/plots/%s/rejBvsSvar_%s.png", outputstr.c_str(), mm.c_str()));
  }
}

void mytmva::effvar(std::vector<std::string> outfname, std::string method)
{
  Float_t x1 = 0, x2 = 1, y1 = 0.2, y2 = 1.0;

  // frame
  TCanvas* c = new TCanvas( "c", "", 200, 0, 650, 500 );
  c->SetGrid();
  c->SetTicks();
  // c->SetFrameFillColor(10);

  TString xtit = "Signal efficiency";
  TString ytit = "Background rejection";
  TString ftit = ytit + " versus " + xtit;
  TString hNameRef = "rejBvsS";
  if(gROOT->FindObject("frame")!=0) gROOT->FindObject("frame")->Delete();
  c->cd();
  TH2F* frame = new TH2F( "frame", ftit, 500, x1, x2, 500, y1, y2 );
  frame->GetXaxis()->SetTitle( xtit );
  frame->GetYaxis()->SetTitle( ytit );
  TMVA::TMVAGlob::SetFrameStyle( frame, 1.0 );
  frame->Draw();
 
  // legend
  Float_t x0L = 0.15, dxL = 0.457-x0L, dyH = 0.19, y0H = 1 - 0.899 + dyH + 0.07; // dyH = 0.22
  TLegend *legend = new TLegend( x0L, y0H-dyH, x0L+dxL, y0H );
  legend->SetTextSize( 0.038 );
  legend->SetHeader( "MVA Variables:" );
  legend->SetMargin( 0.4 );

  Int_t color = 1;
  std::vector<TH1D*> hists; 
  std::string lastlegtex("");
  TString dataset("");
  for(auto& ff : outfname)
    {
      dataset = "";
      // set up dataset
      std::cout << "==> " << __FUNCTION__ << ": Attaching file:" << ff << std::endl; 
      TFile* file = TFile::Open( ff.c_str() );
      if(!file)
        { std::cout << "==> Abort " << __FUNCTION__ << ", please verify filename." << std::endl; continue; }
      if(file->GetListOfKeys()->GetEntries()<=0)
        { std::cout << "==> Abort " << __FUNCTION__ << ", please verify if dataset exist." << std::endl; continue; }
      // --->>
      if( (dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()==1))
        { dataset = ((TKey*)file->GetListOfKeys()->At(0))->GetName(); }
      // <<---
      else { return; }

      TTree* info = (TTree*)file->Get(dataset+"/tmvainfo");
      std::string* varinfo = 0;
      info->SetBranchAddress("var", &varinfo);
      info->GetEntry(0);

      std::string varlegtex = xjjc::str_replaceall(*varinfo, ";", "+= ");
      
      std::string histname(dataset + "/Method_" + method + "/" + method + "/" + "MVA_" + method + "_" + hNameRef);
      TH1D* h = (TH1D*)file->Get(histname.c_str());
      if(!h) { std::cout << "==> " << __FUNCTION__ << ": No hist "<< histname << std::endl; continue; }
      TString hname = h->GetName();
      h->SetLineWidth(3);
      h->SetLineColor(color);
      color++; if(color == 5 || color == 10 || color == 11) color++;
      h->Draw("csame");
      hists.push_back(xjjana::copyobject(h, ff));
      legend->AddEntry(hists[hists.size()-1], xjjc::str_replaceall(varlegtex, lastlegtex, "").c_str(), "l");
      lastlegtex = varlegtex;
    }
  int nmva = hists.size();
  if(!nmva) { std::cout << "==> Abort "<< __FUNCTION__ << ", no valid files." << std::endl; delete c; return; }

  dyH *= (Float_t(TMath::Min(10,nmva) - 3.0)/4.0);
  legend->SetY2( y0H + dyH);

  frame->Draw("sameaxis");
  legend->Draw("same");
  xjjroot::drawtex(1 - c->GetRightMargin(), 1 - c->GetTopMargin()*0.6, method.c_str(), frame->GetTitleSize(), 31, 62, TMVA::TMVAGlob::getTitleBox());
 
  c->Update();
  c->SaveAs("rejBvsSvar.png");

  for(auto& h : hists) { delete h; }
  delete frame;
  delete c;
}

int main(int argc, char* argv[])
{
  if(argc==4)
    {
      for(int i=0; i<mytmva::nptbins; i++)
        {
          mytmva::guieffvar(argv[1], mytmva::ptbins[i], mytmva::ptbins[i+1], argv[2], argv[3]);
        }
      return 0; 
    }
  return 1;
}
