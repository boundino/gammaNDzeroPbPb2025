#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "merge.h"

enum class skimPreset { BRANCH = 0, ZDCOR, ZDCXORJET, DSIZE };

namespace globals {
  const std::vector<std::string> trees_regexp =
    { // std::regex pattern
      ".*",
      // "hiEvtAnalyzer/HiTree",
    };
  std::vector<int> e_skim = {};
  const std::map<skimPreset, std::string> skim_map = {
    { skimPreset::BRANCH, "BRANCH" },
    { skimPreset::ZDCOR, "ZDCOR" },
    { skimPreset::ZDCXORJET, "ZDCXORJET" },
    { skimPreset::DSIZE, "DSIZE" }
  };
}
void recursive_enter(TDirectory*, const std::string&, std::vector<std::string>&, const std::regex&);
std::vector<std::string> get_trees_regexp(TFile* inf);
TFile* dump_file(const std::string& filelist);
bool has_skim(const skimPreset& l) {
  return std::ranges::find(globals::e_skim, static_cast<int>(l)) != globals::e_skim.end();
}
std::vector<std::string> str_skim() {
  std::vector<std::string> r;
  for (const auto& [s, str] : globals::skim_map) {
    if (has_skim(s)) r.push_back(str);
  }
  return r;
}

int macro(std::string outputname, std::string filelist, int ntotal = -1) {
  //
  LOG_XJJ << ": -- Skim presets selected" << std::endl;
  xjjc::print_vec_v(str_skim(), 0);

  LOG_XJJ << ": -- Finalize tree list" << std::endl;
  auto trees = get_trees_regexp(dump_file(filelist));
  if (trees.empty()) {
    LOG_XJJ << ": !! No valid trees. End." << std::endl;
    return 0;
  }

  xjjroot::merge m(trees, filelist, outputname, ntotal);

  // SetBranchAddress and SetBranchStatus
  auto* nt = m.GetTree("Tree");
  if (has_skim(skimPreset::BRANCH)) {
    nt->SetBranchStatus("G*", 0);
    nt->SetBranchStatus("Dtrk*Score", 0);
    nt->SetBranchStatus("Dtrk*dedx", 0);
    nt->SetBranchStatus("gammaN", 0);
    nt->SetBranchStatus("Ngamma", 0);
  }
  int Dsize; nt->SetBranchAddress("Dsize", &Dsize);
  bool isL1ZDCOr; nt->SetBranchAddress("isL1ZDCOr", &isL1ZDCOr);
  bool isL1ZDCXORJet8; nt->SetBranchAddress("isL1ZDCXORJet8", &isL1ZDCXORJet8);

  auto has_skim_ZDCOR = has_skim(skimPreset::ZDCOR),
    has_skim_ZDCXORJET = has_skim(skimPreset::ZDCXORJET),
    has_skim_DSIZE = has_skim(skimPreset::DSIZE);
  
  m.CloneTree();
  const auto nentries = m.GetEntries();
  LOG_XJJ << ": -- Process events" << std::endl;
  for (Long64_t i=0; i<nentries; i++) {
    xjjc::progressbar(i, nentries, 10000);

    m.GetEntry(i);

    if (has_skim_ZDCOR && !(isL1ZDCOr) ) continue;
    if (has_skim_ZDCXORJET && !(isL1ZDCXORJet8) ) continue;
    if (has_skim_DSIZE && !(Dsize > 0) ) continue;

    m.Fill();
  }
  xjjc::progressbar_summary(nentries);
      
  m.Write();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc > 3) {
    globals::e_skim = xjjc::str_convert_vector<int>(argv[3]);
    std::erase_if(globals::e_skim, [](int x) {
      return globals::skim_map.find(static_cast<skimPreset>(x)) == globals::skim_map.end();
    });
  }

  if (argc==5) { return macro(argv[1], argv[2], atoi(argv[4])); }
  if (argc==3 || argc==4) { return macro(argv[1], argv[2]); }
  LOG_XJJ << ": ./macro.exe [outputname] [filelist] ([skimPreset]) ([number of files])" << std::endl;
  return 1;
}

void recursive_enter(TDirectory* dir, const std::string& path_current,
                     std::vector<std::string>& out, const std::regex& re) {
  TIter next(dir->GetListOfKeys());
  TKey* key = nullptr;

  while ((key = static_cast<TKey*>(next()))) {
    TObject* obj = key->ReadObj();
    std::string path_full = path_current + (path_current.empty() ? "" : "/") + obj->GetName();

    if (obj->InheritsFrom(TDirectory::Class())) {
      recursive_enter(static_cast<TDirectory*>(obj), path_full, out, re);
    } else if (obj->InheritsFrom(TTree::Class())) {
      if (std::regex_match(path_full, re)) {
        out.push_back(path_full);
      }
    }
  }
}
std::vector<std::string> get_trees_regexp(TFile* inf) {
  std::vector<std::string> result;
  if (!inf) return result;

  for (const auto &pattern : globals::trees_regexp) {
    std::regex re(pattern);
    std::vector<std::string> strs;
    recursive_enter(inf, "", strs, re);
    for (const auto& istr : strs) {
      if (std::find(result.begin(), result.end(), istr) == result.end())
        result.push_back(istr);
    }
  }
  xjjc::print_vec_v(result, 0);
  
  return result;
}

TFile* dump_file(const std::string& filelist) {
  TFile* inf = nullptr;
  std::ifstream getfile(filelist.c_str());
  while (!inf) {
    std::string filename;
    getfile >> filename;
    inf = TFile::Open(filename.c_str());
    if (getfile.eof()) { break; }
  }
  return inf;
}
