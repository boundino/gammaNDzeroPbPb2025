#ifndef __MVAPROD_H_
#define __MVAPROD_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
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
    mvaprod(std::string weight, bool verbose = true);
    void getnote(const std::string& trainrootfile = "");
    void writenote(TTree* info, const std::string& suffix);
    float evalmva(mytmva::varval* values, int j);
    bool valid() { return valid_; }
    // print info
    std::string method() const { return method_; }
    std::string xmlname() const { return xmlname_; }
    std::string briefing() const;
    
  private:
    bool valid_, verbose_;
    std::string xmlname_;
    TMVA::Reader* reader_;
    std::string method_;
    std::vector<std::string> varnames_, specnames_;
    std::map<std::string, float> varval_, specval_;
    std::string cuts_, cutb_, varinfo_, varnote_;
  };

  xjjc::array2D<std::map<std::string, mvaprod*>> dir_to_weights(const std::string& prefix, std::vector<float> ybins, bool verbose = false);
}

mytmva::mvaprod::mvaprod(std::string weight, bool verbose) :
  valid_(true),
  verbose_(verbose),
  xmlname_(weight) {

  reader_ = new TMVA::Reader(Form("%s", verbose_ ? "!Silent" : "!Color:Silent"));
  
  // read weight file
  void *doc = TMVA::gTools().xmlengine().ParseFile(xmlname_.c_str(), TMVA::gTools().xmlenginebuffersize());
  if (!doc) { __XJJLOG << "!! failed to parse weight xml file: " << xmlname_ << std::endl; valid_ = false; }
  if (valid_) {
    void* rootnode = TMVA::gTools().xmlengine().DocGetRootElement(doc); // node "MethodSetup"
    // method
    std::string fullmethodname("");
    TMVA::gTools().ReadAttr(rootnode, "Method", fullmethodname);
    method_ = fullmethodname; //
    method_.erase(0, fullmethodname.find("::") + 2);
    if (verbose_) std::cout<<std::left<<std::setw(10)<<method_<<" // "<<fullmethodname<<mytmva::nocolor<<std::endl;
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
    if (verbose_) __XJJLOG << mytmva::titlecolor << "++ Add variable:" << mytmva::nocolor << std::endl;
    for (const auto& v : varnames_) {
      varval_[v] = 0; //
      if (verbose_) __XJJLOG << ">> " << std::left << std::setw(10) << v << " // " << mytmva::findvar(v)->var.c_str() << std::endl;
      reader_->AddVariable(mytmva::findvar(v)->var.c_str(), &(varval_.at(v)));
      varnote_ += (";" + v);
    }
    if (verbose_) __XJJLOG<<mytmva::titlecolor<<"++ Add spectator:"<<mytmva::nocolor<<std::endl;
    for (const auto& v : specnames_) {
      specval_[v] = 0; //
      if (verbose_) __XJJLOG<<std::left<<std::setw(10)<<v<<" // "<<mytmva::findvar(v)->var.c_str()<<std::endl;
      reader_->AddSpectator(mytmva::findvar(v)->var.c_str(), &(specval_.at(v)));
    }
    //
    if (verbose_) __XJJLOG<<mytmva::titlecolor<<"++ Book method:"<<mytmva::nocolor<<std::endl;
    std::string methodtag(method_ + " method");
    reader_->BookMVA(methodtag.c_str(), xmlname_.c_str()); // ~
  }
}

void mytmva::mvaprod::getnote(const std::string& trainrootfile) {
  if (verbose_) __XJJLOG<<mytmva::titlecolor<<"++ Training root file:"<<mytmva::nocolor<<std::endl<<"    "<<trainrootfile<<std::endl;
  auto* inf = TFile::Open(trainrootfile.c_str());
  if (inf) {
    auto* rinfo = (TTree*)inf->Get("dataset/tmvainfo");
    if (rinfo) {
      TString *cuts = nullptr, *cutb = nullptr; std::string *varinfo = nullptr;
      rinfo->SetBranchAddress("cuts", &cuts);
      rinfo->SetBranchAddress("cutb", &cutb);
      rinfo->SetBranchAddress("var", &varinfo);
      __XJJLOG<<mytmva::titlecolor<< "++ mva info:"<<mytmva::nocolor<<std::endl;
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

xjjc::array2D<std::map<std::string, mytmva::mvaprod*>> mytmva::dir_to_weights(const std::string& prefix, std::vector<float> ybins, bool verbose) {
  xjjc::print_vec_h(ybins, 0);
  if (ybins.empty()) {
    __XJJLOG << "!! bad ybins, mva deactivated." << std::endl;
    return {};
  }
  auto prods = xjjc::array2d<std::map<std::string, mvaprod*>>(mytmva::ptbins.size()-1, ybins.size()-1);
  for (int i=0; i<mytmva::ptbins.size()-1; i++) {
    for (int j=0; j<ybins.size()-1; j++) {
      mvaprod* pd = nullptr;
      auto& iprod = prods.at(i).at(j);
      if (xjjc::str_contains(prefix, ".weights.xml")) { // single weight for all kinematic bins
        pd = new mvaprod(prefix, verbose);
        if (!pd->valid()) { __XJJLOG << "!! has a bad weight for " << prefix << ", mva deactivated." << std::endl;
          pd = nullptr;
        } else if (iprod.find(pd->method()) != iprod.end()) { __XJJLOG << "!! repeated method " << pd->method() << std::endl;
          continue;
        }
        iprod[pd->method()] = pd; //
      } else { // different weights 
        auto weightdir = prefix + "_" + mytmva::mkname_pt(mytmva::ptbins.at(i), mytmva::ptbins.at(i+1)) + "_" + mytmva::mkname_y(ybins.at(j), ybins.at(j+1));
        __XJJLOG << "++ found weight files in " << weightdir << std::endl;
        for (const auto& entry : fs::directory_iterator(weightdir)) {
          std::string entrypath(entry.path());
          if (!xjjc::str_contains(entrypath, ".weights.xml")) continue;
          if (verbose) __XJJLOG<<"   >> "<<entrypath<<std::endl;
          pd = new mvaprod(entrypath, verbose); //
          if (!pd->valid()) { __XJJLOG << "!! has a bad weight for " << prefix << ", mva deactivated." << std::endl;
            pd = nullptr;
          } else if (iprod.find(pd->method()) != iprod.end()) { __XJJLOG << "!! repeated method " << pd->method() << std::endl;
            continue;
          }
          iprod[pd->method()] = pd; //
        }
      }
    }
  }
  // bad if any kinematic bin is bad
  for (auto& w1 : prods) {
    for (auto& w2 : w1) {
      if (w2.empty()) return {};
      for (auto& [_, pp] : w2) {
        if (!pp) return {};
      }
    }
  }
  // print good
  xjjc::print_tab([&prods]() {
    xjjc::array2D<std::string> result;
    for (const auto& pts : prods) {
      std::vector<std::string> rys;
      for (const auto& ys : pts) {
        rys.push_back(Form("%d methods", ys.size()));
      }
      result.push_back(rys);
    }
    return result;
  }(), 3);

  if (verbose) {
    for (const auto& pts : prods) {
      xjjc::print_tab([&pts]() {
        xjjc::array2D<std::string> result;
        for (const auto& ys : pts) {
          std::vector<std::string> rys;
          for (const auto& [_, m] : ys) {
            rys.push_back(m->briefing());
          }
          result.push_back(rys);
        }
        return result;
      }(), 0);
    }
  }
  return prods;
}

std::string mytmva::mvaprod::briefing() const {
  std::string r = "";
  if (valid_) {
    r += ("(" + method_ + ") ");
    std::string tvar = "";
    for (const auto& v : varnames_) tvar += ((tvar.empty()?"":"/") + v);
    r += tvar;
  }
  return r;
}

#endif
