#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "xjjcuti.h"
#include "merge.h"

namespace globals {
  std::vector<std::string> trees_regexp =
    {
      // ".*",
      "Tree",
      // ".*/HltTree",
      // "hiEvtAnalyzer/HiTree",
    };
}
void recursive_enter(TDirectory*, const std::string&, std::vector<std::string>&, const std::regex&);
std::vector<std::string> get_trees_regexp(TFile* infile);
TFile* dump_file(const std::string& filelist);

int macro(std::string outputname, std::string filelist, bool isskim = false, int ntotal = -1) {
  //
  std::cout<<__FUNCTION__<<": -- Finalized tree list"<<std::endl;
  auto trees = get_trees_regexp(dump_file(filelist));
  if (trees.empty()) {
    std::cout<<__FUNCTION__<<": !! No valid trees. End."<<std::endl;
    return 0;
  }

  xjjroot::merge m(trees, filelist, outputname, ntotal);
  auto* nt = m.GetTree("Tree");
  int Dsize; nt->SetBranchAddress("Dsize", &Dsize);
  bool isL1ZDCOr; nt->SetBranchAddress("isL1ZDCOr", &isL1ZDCOr);

  const auto nentries = m.GetEntries();
  std::cout<<__FUNCTION__<<": -- Event reading"<<std::endl;
  for (Long64_t i=0; i<nentries; i++) {
    xjjc::progressbar(i, nentries, 10000);

    m.GetEntry(i);

    // if (isskim && !( true )) continue; // skim
    if (isskim && !( Dsize > 0 && isL1ZDCOr)) continue;

    m.Fill();
  }
  xjjc::progressbar_summary(nentries);
      
  std::cout<<__FUNCTION__<<": -- Writing new trees done"<<std::endl;
  m.Write();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc==5) { return macro(argv[1], argv[2], atoi(argv[3]), atoi(argv[4])); }
  if (argc==4) { return macro(argv[1], argv[2], atoi(argv[3])); }
  if (argc==3) { return macro(argv[1], argv[2]); }
  std::cout<<__FUNCTION__<<": ./macro.exe [outputname] [filelist] ([number of files])";
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
std::vector<std::string> get_trees_regexp(TFile* infile) {
  std::vector<std::string> result;
  if (!infile) return result;

  for (const auto &pattern : globals::trees_regexp) {
    std::regex re(pattern);
    std::vector<std::string> strs;
    recursive_enter(infile, "", strs, re);
    for (const auto& istr : strs) {
      if (std::find(result.begin(), result.end(), istr) == result.end())
        result.push_back(istr);
    }
  }
  xjjc::print_tab(result);
  
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
