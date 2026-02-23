#ifndef _MVAEFFS_H_
#define _MVAEFFS_H_

#include "TMVA/mvaeffs.h"
#include <TFormula.h>
#include <TString.h>
#include <TList.h>
#include <TDirectory.h>
#include <TKey.h>
#include <iostream>
#include <iomanip>

namespace mytmva
{
  class mvaeffs
  {
  public:
    mvaeffs(TString ds, TString fin, Float_t ns=0, Float_t nb=0, TString formula="S/sqrt(S+B)", Bool_t useTMVAStyle=kTRUE);
    bool isvalid() { return fvalid; }

  private:
    Float_t fNSignal;
    Float_t fNBackground;  
    TString fFormula;
    TString dataset;
    TList* fInfoList;
    bool fvalid;
    bool fsignificance;
 
    Int_t maxLenTitle;

    TString GetFormula();
    void ReadHistograms(TFile* file);
    void UpdateSignificanceHists();
    void DrawHistograms();

  };
}

mytmva::mvaeffs::mvaeffs(TString ds, TString fin, Float_t ns, Float_t nb, TString formula, Bool_t useTMVAStyle) :
  fNSignal(ns),
  fNBackground(nb),
  fFormula(formula),
  dataset(ds),
  fInfoList(0),
  fsignificance(true),
  fvalid(true),
  maxLenTitle(0)
{
  TMVA::TMVAGlob::Initialize( useTMVAStyle );

  TFile* inf = TFile::Open(fin.Data());
  if(!inf) { std::cout << "==> Abort " << __FUNCTION__ << ": invalid input file." << std::endl; fvalid = false; return; }

  if(ns<=0 || nb<=0) { fsignificance = false; }

  ReadHistograms(inf);
  UpdateSignificanceHists();
  DrawHistograms();
}

void mytmva::mvaeffs::ReadHistograms(TFile* file) 
{
  if(fInfoList)
    { 
      TIter next(fInfoList);
      TMVA::MethodInfo *info(0);
      while ( (info = (TMVA::MethodInfo*)next()) ) {
        delete info;
      }
      delete fInfoList;
      fInfoList = 0;
    }
  fInfoList = new TList;
 
  // search for the right histograms in full list of keys
  TIter next(file->GetDirectory(dataset.Data())->GetListOfKeys());
  TKey *key(0);   
  while( (key = (TKey*)next()) ) {
 
    if (!TString(key->GetName()).BeginsWith("Method_")) continue;
    if( ! gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory") ) continue;
 
    std::cout << "--- Found directory: " << ((TDirectory*)key->ReadObj())->GetName() << std::endl;
 
    TDirectory* mDir = (TDirectory*)key->ReadObj();
 
    TIter keyIt(mDir->GetListOfKeys());
    TKey *titkey;
    while((titkey = (TKey*)keyIt())) {
      if( ! gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory") ) continue;
         
      TMVA::MethodInfo* info = new TMVA::MethodInfo();
      TDirectory* titDir = (TDirectory *)titkey->ReadObj();
 
      TMVA::TMVAGlob::GetMethodName(info->methodName,key);
      TMVA::TMVAGlob::GetMethodTitle(info->methodTitle,titDir);        
      if (info->methodTitle.Length() > maxLenTitle) maxLenTitle = info->methodTitle.Length();
      TString hname = "MVA_" + info->methodTitle;
         
      std::cout << "--- Classifier: " << info->methodTitle << std::endl;
         
      info->sig = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
      info->bgd = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));
      info->origSigE = dynamic_cast<TH1*>(titDir->Get( hname + "_effS" ));
      info->origBgdE = dynamic_cast<TH1*>(titDir->Get( hname + "_effB" ));      
      if (info->origSigE==0 || info->origBgdE==0) { delete info; continue; }
 
      info->SetResultHists();
      fInfoList->Add(info);
    }
  }
  return;
}


void mytmva::mvaeffs::DrawHistograms() 
{
  // counter variables
  Int_t countCanvas = 0;
 
  // define Canvas layout here!
  const Int_t width = 600;   // size of canvas
  Int_t signifColor = TColor::GetColor( "#00aa00" );
 
  TIter next(fInfoList);
  TMVA::MethodInfo* info(0);
  while ( (info = (TMVA::MethodInfo*)next()) ) 
    { 
      // create new canvas
      TCanvas* c = new TCanvas( Form("canvas%d", countCanvas+1), 
                                Form("Cut efficiencies for %s classifier",info->methodTitle.Data()), 
                                countCanvas*50+200, countCanvas*20, width, Int_t(width*0.78) ); 
      info->canvas = c;
      c->cd();
 
      // draw grid
      c->SetGrid(1);
      c->SetTickx(0);
      c->SetTicky(0);

      // TMVA::TMVAGlob::SetTMVAStyle(); 
      TStyle *TMVAStyle = gROOT->GetStyle("Plain"); // our style is based on Plain
      TMVAStyle->SetLineStyleString( 5, "[32 22]" );
      TMVAStyle->SetLineStyleString( 6, "[12 22]" );
      // TMVAStyle->SetLineStyleString( 7, "[22 10 7 10]" );
      c->SetTopMargin(.2);
      c->SetFrameFillColor(10);
      c->SetFrameLineWidth(2.);

      // TH1* effcanvas = fsignificance?info->effpurS:info->sigE;
      TH1* effcanvas = info->effpurS;
       
      // and the signal purity and quality
      effcanvas->SetTitle("Cut efficiencies and optimal cut value");
      if (info->methodTitle.Contains("Cuts")) {
        effcanvas->GetXaxis()->SetTitle( "Signal Efficiency" );
      }
      else {
        effcanvas->GetXaxis()->SetTitle( TString("Cut value applied on ") + info->methodTitle + " output" );
      }
      effcanvas->GetYaxis()->SetTitle( "Efficiency (Purity)" );

      TMVA::TMVAGlob::SetFrameStyle( effcanvas );
 
      c->SetTicks(0,0);
      c->SetRightMargin ( 2.0 );
 
      effcanvas->SetMaximum(1.1);
      effcanvas->Draw("histl");

      //
      if(fsignificance) info->effpurS->Draw("samehistl");
      if(fsignificance) info->purS->Draw("samehistl");      
 
      // overlay signal and background histograms
      info->sigE->Draw("samehistl");
      info->bgdE->Draw("samehistl");
 
      info->sSig->SetLineColor( signifColor );
      if(fsignificance) info->sSig->Draw("samehistl");
 
      // redraw axes
      // effcanvas->Draw( "sameaxis" );

      // Draw legend               
      TLegend *legend1= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin(), 
                                     c->GetLeftMargin() + 0.4, 1 - c->GetTopMargin() + 0.12 );
      // legend1->SetFillStyle( 1 );
      legend1->AddEntry(info->sigE,"Signal efficiency","L");
      legend1->AddEntry(info->bgdE,"Background efficiency","L");
      legend1->SetBorderSize(0); // 1
      legend1->SetMargin( 0.3 );
      legend1->Draw("same");
 
      TLegend *legend2= new TLegend( c->GetLeftMargin() + 0.4, 1 - c->GetTopMargin(), 
                                     1 - c->GetRightMargin(), 1 - c->GetTopMargin() + 0.12 );
      // legend2->SetFillStyle( 1 );
      legend2->AddEntry(info->purS,"Signal purity","L");
      legend2->AddEntry(info->effpurS,"Signal efficiency*purity","L");
      legend2->AddEntry(info->sSig,fFormula.Data(),"L");
      legend2->SetBorderSize(0); // 1
      legend2->SetMargin( 0.3 );
      if(fsignificance) legend2->Draw("same");
          
      // line to indicate maximum efficiency
      TLine* effline = new TLine( info->sSig->GetXaxis()->GetXmin(), 1, info->sSig->GetXaxis()->GetXmax(), 1 );
      effline->SetLineWidth( 3 );
      effline->SetLineColor( 1 );
      effline->SetLineStyle( 2 );
      effline->Draw();
 
      TLatex tl;
      tl.SetNDC();
      tl.SetTextSize( 0.033 );

      // print comments
      if(fsignificance)
        {
          Int_t maxbin = info->sSig->GetMaximumBin();
          info->line1 = tl.DrawLatex( 0.15, 0.23, Form("For %1.0f signal and %1.0f background", fNSignal, fNBackground));
          tl.DrawLatex( 0.15, 0.19, "events the maximum "+fFormula+" is");
 
          if (info->maxSignificanceErr > 0) {
            info->line2 = tl.DrawLatex( 0.15, 0.15, Form("%5.2f +- %4.2f when cutting at %5.2f", 
                                                         info->maxSignificance, 
                                                         info->maxSignificanceErr, 
                                                         info->sSig->GetXaxis()->GetBinCenter(maxbin)) );
          }
          else {
            info->line2 = tl.DrawLatex( 0.15, 0.15, Form("%4.2f when cutting at %5.2f", 
                                                         info->maxSignificance, 
                                                         info->sSig->GetXaxis()->GetBinCenter(maxbin)) );
          }

          // Draw second axes
          info->rightAxis = new TGaxis(effcanvas->GetXaxis()->GetXmax(), c->GetUymin(),
                                       effcanvas->GetXaxis()->GetXmax(), c->GetUymax(),0,1.1*info->maxSignificance,510,"+L");
          info->rightAxis->SetLineColor ( signifColor );
          info->rightAxis->SetLabelColor( signifColor );
          info->rightAxis->SetTitleColor( signifColor );
 
          info->rightAxis->SetTitleSize( info->sSig->GetXaxis()->GetTitleSize() );
          info->rightAxis->SetTitle( "Significance" );
          info->rightAxis->Draw();
        } 

      if (info->methodTitle.Contains("Cuts")){
        tl.DrawLatex( 0.13, 0.77, "Method Cuts provides a bundle of cut selections, each tuned to a");
        tl.DrawLatex(0.13, 0.74, "different signal efficiency. Shown is the purity for each cut selection.");
      }

      // save canvas to file
      c->Update();
 
      // switche
      const Bool_t Save_Images = kTRUE;
 
      if (Save_Images) {
        TMVA::TMVAGlob::imgconv( c, Form("%s/plots/mvaeffs_%s%s",dataset.Data(), info->methodTitle.Data(), fsignificance?"":"_nosig") ); 
      }
      countCanvas++;
      delete c;
    }
}

void mytmva::mvaeffs::UpdateSignificanceHists() 
{
  TFormula f("sigf", GetFormula());
  TIter next(fInfoList);
  TMVA::MethodInfo* info(0);
  TString cname = "Classifier";
  if (cname.Length() >  maxLenTitle)  maxLenTitle = cname.Length();
  TString str = Form( "%*s   (  #signal, #backgr.)  Optimal-cut  %s      NSig      NBkg   EffSig   EffBkg", 
                      maxLenTitle, cname.Data(), fFormula.Data() );
  std::cout << "--- " << std::setfill('=') << std::setw(str.Length()) << "" << std::setfill(' ') << std::endl;
  std::cout << "--- " << str << std::endl;
  std::cout << "--- " << std::setfill('-') << std::setw(str.Length()) << "" << std::setfill(' ') << std::endl;
  Double_t maxSig    = -1;
  Double_t maxSigErr = -1;
  while ((info = (TMVA::MethodInfo*)next())) {
    for (Int_t i=1; i<=info->origSigE->GetNbinsX(); i++) {
      Float_t eS = info->origSigE->GetBinContent( i );
      Float_t S = eS * fNSignal;
      Float_t B = info->origBgdE->GetBinContent( i ) * fNBackground;
      info->purS->SetBinContent( i, (S+B==0)?0:S/(S+B) );
      
      Double_t sig = (S+B==0)?0:f.Eval(S,B);
      if (sig > maxSig) {
        maxSig    = sig;
        if (fFormula== "S/sqrt(B)" && (S>0 && B>0)) {
          maxSigErr = sig * sqrt( 1./S + 1./(2.*B));
        }
      }
      info->sSig->SetBinContent( i, sig );
      info->effpurS->SetBinContent( i, eS*info->purS->GetBinContent( i ) );
    }
       
    info->maxSignificance    = info->sSig->GetMaximum();
    info->maxSignificanceErr = (maxSigErr > 0) ? maxSigErr : 0;
    info->sSig->Scale(1/info->maxSignificance);
 
    // PrintResults( info );
  }
  std::cout << "--- " << std::setfill('-') << std::setw(str.Length()) << "" << std::setfill(' ') << std::endl << std::endl;
}

TString mytmva::mvaeffs::GetFormula() 
{
  TString f = fFormula;
  f.ReplaceAll("S","x");
  f.ReplaceAll("B","y");
  return f;
}

#endif
