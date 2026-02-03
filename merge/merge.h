#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#define LOG_XJJ \
  std::cout << std::left << std::setw(11) << __FUNCTION__

namespace xjjroot
{
  class merge
  {
  public:
    merge(const std::vector<std::string>& treelist, std::string filelist, std::string outputname, int ntotal = -1);
    void CloneTree();
    Long64_t GetEntries() const;
    void GetEntry(Long64_t i);
    void Fill();
    void Write() { outf_->Write(); outf_->Close(); }
    TChain* GetTree(const std::string& key);

    std::vector<std::string> GetTreeList() const;

  private: 
    TFile* outf_;
    std::vector<std::string> files_;
    std::map<std::string, TChain*> trs_;
    std::map<std::string, TDirectory*> dirs_;
    std::map<std::string, TTree*> newtrs_;
    std::map<std::string, bool> status_;

    void set_files(std::string filelist, int ntotal);
  };
}

xjjroot::merge::merge(const std::vector<std::string>& treelist, std::string filelist,
                      std::string outputname, int ntotal) {
  //
  TTree::SetMaxTreeSize(1LL * 1024 * 1024 * 1024 * 1024);
  set_files(filelist, ntotal);

  LOG_XJJ << ": -- Chain files" << std::endl;
  for (const auto& t : treelist)
    trs_[t] = new TChain(t.c_str()); //
  for (auto& i : files_) {
    std::cout << "\e[2m" << i << "\e[0m" << std::endl;
    for (auto& t : trs_) t.second->Add(i.c_str());
  }
  LOG_XJJ << ": >> Merged \e[31;1m" << files_.size() << "\e[0m files." << std::endl;

  outf_ = xjjroot::newfile(outputname);
  for (const auto& [t, _] : trs_) {
    dirs_[t] = outf_; //
    auto parse = xjjc::str_divide(t, "/");
    if (!parse.empty()) parse.pop_back();
    std::reverse(parse.begin(), parse.end());
    while (!parse.empty()) {
      dirs_[t] = dirs_.at(t)->mkdir(parse.back().c_str(), "");
      parse.pop_back();
    }
  }

}

void xjjroot::merge::CloneTree() {
  LOG_XJJ << ": -- Clone trees" << std::endl;
  for (const auto& [t, _] : trs_) {
    dirs_.at(t)->cd();
    newtrs_[t] = trs_.at(t)->CloneTree(0); //
    status_[t] = false;
  }
  outf_->cd();  
}

Long64_t xjjroot::merge::GetEntries() const {
  if (trs_.empty()) return -1;
  auto n = trs_.begin()->second->GetEntries();
  for (const auto& t : trs_) {
    auto diff_n = t.second->GetEntries() != n;
    LOG_XJJ << ": " << (diff_n ? "!!" : ">>") << " " << t.first << " (" << t.second->GetEntries() << ")" << std::endl;
    n = std::max(n, t.second->GetEntries());
  }
  return n;
}

void xjjroot::merge::GetEntry(Long64_t i) {
  for(auto& t : trs_) {
    status_[t.first] = false;
    if (i < t.second->GetEntries()) {
      t.second->GetEntry(i);
      status_[t.first] = true;
    }
  }
}

void xjjroot::merge::Fill() {
  for (const auto& [t, _] : newtrs_) {
    if (!status_.at(t)) continue;
    dirs_.at(t)->cd();
    newtrs_.at(t)->Fill();
  }
  outf_->cd();
}

std::vector<std::string> xjjroot::merge::GetTreeList() const {
  std::vector<std::string> keys;
  keys.reserve(trs_.size());
  for (const auto& [t, _] : trs_) {
    keys.push_back(t);
  }
  return keys;
}

TChain* xjjroot::merge::GetTree(const std::string& key) {
  auto it = trs_.find(key);
  if (it == trs_.end()) {
    LOG_XJJ << ": !! bad key " << key << std::endl;
    return nullptr;
  }
  return it->second;
}

void xjjroot::merge::set_files(std::string filelist, int ntotal = -1) {
  std::ifstream getfile(filelist.c_str());
  for (int i=0; (i<ntotal || ntotal<0); i++) {
    std::string filename;
    getfile >> filename;
    if (getfile.eof()) { break; }
    files_.push_back(filename);
  }
}
