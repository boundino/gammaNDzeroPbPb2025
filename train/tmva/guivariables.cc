// https://root.cern.ch/doc/v610/TMVAGui_8cxx_source.html
// https://root.cern.ch/doc/v608/variables_8cxx_source.html
#include <iostream>
#include <string>

#include <TKey.h>
#include <TList.h>
#include <TObjString.h>
#include "TMVA/variables.h"

#include "TMVAClassification.h"
#include "xjjcuti.h"

namespace mytmva
{
  static TList* TMVAGui_keyContent;
  void guivariables(std::string outputname, float ptmin, float ptmax, float ymin, float ymax, std::string mymethod, std::string stage = "0,1,2,3");
  void variables(const std::string&);
  TList* GetKeyList(const std::string& pattern);
}

void mytmva::guivariables(std::string outputname, float ptmin, float ptmax, float ymin, float ymax, std::string mymethod, std::string stage/* = "0,1,2,3"*/) {
  std::string mvafname = mytmva::mkname(outputname, mymethod, stage, ptmin, ptmax, ymin, ymax);
  __XJJLOG << "++ [ " << ptmin << ", " << ptmax << ", " << ymin << ", " << ymax << "] " << mvafname << std::endl;
  mytmva::variables(mvafname);
}

void mytmva::variables(const std::string& mvafname) {  
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

  TMVAGui_keyContent = (TList*)file->GetDirectory(dataset.Data())->GetListOfKeys()->Clone(); //
  
  TList* keylist = GetKeyList("InputVariables");
  TListIter it( keylist );
  TObjString* str = 0;
  while ((str = (TObjString*)it())) {
    auto tmp = str->GetString();
    std::string title = Form("Input variables '%s'-transformed (training sample)", 
                             tmp.ReplaceAll("InputVariables_","").Data() );
    if (tmp.Contains("Id")) title = "Input variables (training sample)";
    // void TMVA::variables (TString dataset, TString fin="TMVA.root", TString dirName="InputVariables_Id", TString title="TMVA Input Variables", Bool_t isRegression=kFALSE, Bool_t useTMVAStyle=kTRUE)
    TMVA::variables(dataset.Data(), mvafname.c_str(), str->GetString().Data(), title.c_str());
  }

  const auto outputstr = xjjc::str_tag_from_file(mvafname);
  // gSystem->Exec(Form("rm %s/plots/*.png", dataset.Data()));
  gSystem->Exec(Form("mkdir -p %s/plots/%s", dataset.Data(), outputstr.c_str()));
  gSystem->Exec(Form("mv %s/plots/*.png %s/plots/%s/", dataset.Data(), dataset.Data(), outputstr.c_str()));
}

TList* mytmva::GetKeyList(const std::string& pattern) {
  TList* list = new TList();
  TIter next( TMVAGui_keyContent );
  TKey* key(0);
  while ((key = (TKey*)next())) {
    if (xjjc::str_contains(key->GetName(), pattern)) { list->Add( new TObjString(key->GetName()) ); }
  }
  return list;
}

int main(int argc, char* argv[]) {
  if (argc==4) {
    for (int i=0; i<mytmva::ptbins.size()-1; i++) {
      for (int j=0; j<mytmva::ybins.size()-1; j++) {
        mytmva::guivariables(argv[1], mytmva::ptbins[i], mytmva::ptbins[i+1], mytmva::ybins[j], mytmva::ybins[j+1], argv[2], argv[3]); 
      }
    }
    return 0; 
  }
  return 1;
}
