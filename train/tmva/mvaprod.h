#ifndef __MVAPROD_H_
#define __MVAPROD_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "TMVAClassification.h"
#include "xjjcuti.h"

namespace mytmva
{
  std::string titlecolor = "\e[34m", nocolor = "\e[0m", contentcolor = "\e[34m", errorcolor = "\e[31m";

  class mvaprod
  {
  public:
    mvaprod(std::string weight);
    void getnote(const std::string& trainrootfile = "");
    void writenote(TTree* info, const std::string& suffix);
    float evalmva(mytmva::varval* values, int j);
    std::string method() { return method_; }
    bool valid() { return valid_; }
    
  private:
    bool valid_;
    std::string xmlname_;
    TMVA::Reader* reader_;
    std::string method_;
    std::vector<std::string> varnames_, specnames_;
    std::map<std::string, float> varval_, specval_;
    std::string cuts_, cutb_, varinfo_, varnote_;
  };
}

mytmva::mvaprod::mvaprod(std::string weight) :
  valid_(true),
  xmlname_(weight),
  reader_(new TMVA::Reader("!Color:!Silent")) {

  // read weight file
  void *doc = TMVA::gTools().xmlengine().ParseFile(xmlname_.c_str(), TMVA::gTools().xmlenginebuffersize());
  if (!doc) { std::cout<<__FUNCTION__<<": failed to parse weight xml file: "<<xmlname_<<std::endl; valid_ = false; }
  if (valid_) {
    void* rootnode = TMVA::gTools().xmlengine().DocGetRootElement(doc); // node "MethodSetup"
    // method
    std::string fullmethodname("");
    TMVA::gTools().ReadAttr(rootnode, "Method", fullmethodname);
    method_ = fullmethodname; //
    method_.erase(0, fullmethodname.find("::") + 2);
    std::cout<<std::left<<std::setw(10)<<method_<<" // "<<fullmethodname<<mytmva::nocolor<<std::endl;
    // variable
    void* variables = TMVA::gTools().GetChild(rootnode, "Variables");
    UInt_t NVar = 0;
    TMVA::gTools().ReadAttr(variables, "NVar", NVar);
    void* var = TMVA::gTools().GetChild(variables, "Variable");
    for (unsigned int k=0; k<NVar; k++) {
      std::string varlabel("");
      TMVA::gTools().ReadAttr(var, "Label", varlabel);
      varnames_.push_back(varlabel); //
      var = TMVA::gTools().GetNextChild(var);
    }
    // spectator
    void* spectators = TMVA::gTools().GetChild(rootnode, "Spectators");
    UInt_t NSpec = 0;
    TMVA::gTools().ReadAttr(spectators, "NSpec", NSpec);
    void* spec = TMVA::gTools().GetChild(spectators, "Spectator");
    for (unsigned int k=0; k<NSpec; k++) {
      std::string speclabel("");
      TMVA::gTools().ReadAttr(spec, "Label", speclabel);
      specnames_.push_back(speclabel); //
      spec = TMVA::gTools().GetNextChild(spec);
    }
  }
  valid_ = (varnames_.size() > 0 && !method_.empty());

  if (valid_) {
    // reader
    std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Add variable:"<<mytmva::nocolor<<std::endl;
    for (const auto& v : varnames_) {
      varval_[v] = 0; //
      std::cout<<std::left<<std::setw(10)<<v<<" // "<<mytmva::findvar(v)->var.c_str()<<std::endl;
      reader_->AddVariable(mytmva::findvar(v)->var.c_str(), &(varval_.at(v)));
      varnote_ += (";" + v);
    }
    std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Add spectator:"<<mytmva::nocolor<<std::endl;
    for (const auto& v : specnames_) {
      specval_[v] = 0; //
      std::cout<<std::left<<std::setw(10)<<v<<" // "<<mytmva::findvar(v)->var.c_str()<<std::endl;
      reader_->AddSpectator(mytmva::findvar(v)->var.c_str(), &(specval_.at(v)));
    }
    //
    std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Book method:"<<mytmva::nocolor<<std::endl;
    std::string methodtag(method_ + " method");
    reader_->BookMVA(methodtag.c_str(), xmlname_.c_str()); // ~
  }
}

void mytmva::mvaprod::getnote(const std::string& trainrootfile) {
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Training root file:"<<mytmva::nocolor<<std::endl;
  std::cout<<trainrootfile<<std::endl;
  auto* inf = TFile::Open(trainrootfile.c_str());
  if (inf) {
    auto* rinfo = (TTree*)inf->Get("dataset/tmvainfo");
    if (rinfo) { 
      TString *cuts = nullptr, *cutb = nullptr; std::string *varinfo = nullptr;
      rinfo->SetBranchAddress("cuts", &cuts);
      rinfo->SetBranchAddress("cutb", &cutb);
      rinfo->SetBranchAddress("var", &varinfo);
      std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": mva info:"<<mytmva::nocolor<<std::endl;
      rinfo->Show(0); std::cout<<std::endl;
      rinfo->GetEntry(0);
      cuts_ = *cuts;
      cutb_ = *cutb;
      varinfo_ = *varinfo;
      inf->Close();
    }
  }
}

void mytmva::mvaprod::writenote(TTree* info, const std::string& suffix) {
  std::cout<<__FUNCTION__<<std::endl;
  info->Branch(Form("cuts%s", suffix.c_str()), &(cuts_));
  info->Branch(Form("cutb%s", suffix.c_str()), &(cutb_));
  info->Branch(Form("varinfo%s", suffix.c_str()), &(varinfo_));
  info->Branch(Form("varnote%s", suffix.c_str()), &(varnote_));
}

float mytmva::mvaprod::evalmva(mytmva::varval* values, int j) {
  int badinput = 0;
  for (auto& [v, val] : varval_) {
    val = values->getval(v, j);
    if (std::isnan(val)) badinput++;
  }
  for (auto& [v, val] : specval_) {
    val = values->getval(v, j);
  }
  // xjjc::print_tab(varval_, 0);
  float result = (badinput ? -999 : reader_->EvaluateMVA(Form("%s method", method_.c_str())));
  return result;
}

#endif
