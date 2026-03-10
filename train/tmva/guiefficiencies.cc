// https://root.cern.ch/doc/v610/TMVAGui_8cxx_source.html
// https://root.cern.ch/doc/v608/efficiencies_8cxx_source.html
#include <iostream>
#include <string>

#include <TKey.h>
#include <TList.h>
#include "TMVA/efficiencies.h"
#include "TMVA/mvas.h"
#include "TMVA/correlations.h"

#include "mvaeffs.h"
#include "TMVAClassification.h"
#include "xjjcuti.h"

namespace mytmva
{
  void guiefficiencies(std::string outputname, float ptmin, float ptmax, float ymin, float ymax, std::string mymethod, std::string stage = "0,1,2,3");
  void efficiencies(std::string outfname);
}

void mytmva::guiefficiencies(std::string outputname, float ptmin, float ptmax, float ymin, float ymax, std::string mymethod, std::string stage/* = "0,1,2,3"*/) {
  std::string mvafname = mytmva::mkname(outputname, mymethod, stage, ptmin, ptmax, ymin, ymax);
  __XJJLOG << "++ [ " << ptmin << ", " << ptmax << ", " << ymin << ", " << ymax << "] " << mvafname << std::endl;
  mytmva::efficiencies(mvafname);
}

void mytmva::efficiencies(std::string mvafname) {
  // Set up dataset
  TString dataset("");
  auto* file = TFile::Open(mvafname.c_str());
  if (!file) {
    __XJJLOG << "!! Abort, please verify filename." << std::endl; return; }
  if (file->GetListOfKeys()->GetEntries() <= 0) {
    __XJJLOG << "!! Abort, please verify if dataset exist." << std::endl; return; }
  if ((dataset == "" || dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries() == 1)) {
    dataset = ((TKey*)file->GetListOfKeys()->At(0))->GetName();
  } else if ((dataset == "" || dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries() > 1)) {
    __XJJLOG << "!! Abort, more than 1 dataset." << std::endl;
    for (int i=0; i<file->GetListOfKeys()->GetEntries(); i++) {
      dataset = ((TKey*)file->GetListOfKeys()->At(i))->GetName();
      __XJJLOG << "    " << dataset << std::endl;
    }
    return;
  }
  else { return; }

  TMVA::efficiencies(dataset.Data(), mvafname.c_str(), 2); // 1,2,3
  TMVA::mvas(dataset.Data(), mvafname.c_str(), TMVA::kCompareType);
  TMVA::correlations(dataset.Data(), mvafname.c_str());
  mytmva::mvaeffs(dataset.Data(), mvafname.c_str());
  mytmva::mvaeffs(dataset.Data(), mvafname.c_str(), 4.e+2, 1.e+5);

  const auto outputstr = xjjc::str_tag_from_file(mvafname);
  // gSystem->Exec(Form("rm %s/plots/*.png", dataset.Data()));
  gSystem->Exec(Form("mkdir -p %s/plots/%s", dataset.Data(), outputstr.c_str()));
  gSystem->Exec(Form("mv %s/plots/*.png %s/plots/%s/", dataset.Data(), dataset.Data(), outputstr.c_str()));
}

int main(int argc, char* argv[]) {
  if (argc==4) {
    for (int i=0; i<mytmva::ptbins.size()-1; i++) {
      for (int j=0; j<mytmva::ybins.size()-1; j++) {
        mytmva::guiefficiencies(argv[1], mytmva::ptbins[i], mytmva::ptbins[i+1], mytmva::ybins[j], mytmva::ybins[j+1], argv[2], argv[3]);
      }
    }
    return 0; 
  }
  return 1;
}
