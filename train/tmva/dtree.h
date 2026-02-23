#ifndef __PHOD_DTREE__
#define __PHOD_DTREE__

#include <TTree.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2F.h>

#include <vector>
#include <map>
#include <string>
#include <type_traits>

#include "xjjcuti.h"

namespace hfupc
{
  class dtree
  {
  public:
    dtree(TTree* nt);
    TTree* nt() { return nt_; }

    // read
    template<typename T = float> T val(const std::string& br, int j);
    bool has(const std::string& br);
    int GetEntries() { return nt_->GetEntries(); }
    void GetEntry(int i) { nt_->GetEntry(i); }

    // tools
    int Dsize() { return Dsize_; }

  private:
    //Members set in constructor
    TTree* nt_;

    void setbranchaddress();
    //List of non-array variables, all datatypes
    std::vector<std::string> tbfio_ =
      {
	"Dsize"
      };
    //TTree variables, non-array, all datatypes
    int Dsize_;

    //List of variables, arrays tied to Dsize -> floats, ints, bools
    std::vector<std::string> tbvf_ =
      {
	"Dmass",
	"Dpt",
	"Deta",
	"Dphi",
	"Dy",
	"Dchi2cl",
	"Ddtheta",
	"Dalpha",
	"DsvpvDistance",
	"DsvpvDisErr",
	"DsvpvDistance_2D",
	"DsvpvDisErr_2D",
	"Dtrk1Pt",
	"Dtrk2Pt",
	"Dtrk1PtErr",
	"Dtrk2PtErr",
	"Dtrk1Eta",
	"Dtrk2Eta",
	"Dtrk1Dz1",
	"Dtrk2Dz1",
	"Dtrk1DzError1",
	"Dtrk2DzError1",
	"Dtrk1Dxy1",
	"Dtrk2Dxy1",
	"Dtrk1DxyError1",
	"Dtrk2DxyError1",
	"Dtrk1PixelHit",
	"Dtrk2PixelHit",
	"Dtrk1StripHit",
	"Dtrk2StripHit",
	"Dtrk1nStripLayer",
	"Dtrk2nStripLayer",
	"Dtrk1nPixelLayer",
	"Dtrk2nPixelLayer",
	"Dtrk1Chi2ndf",
	"Dtrk2Chi2ndf",
	"BDT",
	"Dgenpt",
	"Dgeneta",
	"Dgenphi",
	"Dgeny"
      };
    std::vector<std::string> tbvi_ =
      {
	"Dgen",
	"DgencollisionId",
      };
    std::vector<std::string> tbvo_ =
      {
      };

    std::map<std::string, std::vector<float>*> bvf_;
    std::map<std::string, std::vector<int>*> bvi_;
    std::map<std::string, std::vector<bool>*> bvo_;

    // status
    std::map<std::string, bool> bvs_; //

  };
}

hfupc::dtree::dtree(TTree* nt) :
  nt_(nt) {
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  for (auto& i : tbvf_) { bvf_[i] = nullptr; }
  for (auto& i : tbvi_) { bvi_[i] = nullptr; }
  for (auto& i : tbvo_) { bvo_[i] = nullptr; }

  nt_->SetBranchStatus("*", 0);

  // std::vector<std::vector<std::string> > tbSuperVect = { tbfio_, tbvf_, tbvi_, tbvo_ };
  for (auto& tv : { tbfio_, tbvf_, tbvi_, tbvo_ }) {
    for (auto& b : tv) {
      bvs_[b] = false;
      if (nt_->FindBranch(b.c_str())) {
        nt_->SetBranchStatus(b.c_str(), 1);
        bvs_[b] = true;
      }
    }
  }

  xjjc::print_tab(bvs_);
  setbranchaddress();
}

bool hfupc::dtree::has(const std::string& br) {
  bool r = false;
  if (bvs_.find(br) != bvs_.end()) {
    r = bvs_.at(br);
  }
  return r;
}

void hfupc::dtree::setbranchaddress() {
  if (has("Dsize")) nt_->SetBranchAddress("Dsize", &Dsize_);
  for (auto& b : tbvf_) { if (has(b)) { nt_->SetBranchAddress(b.c_str(), &(bvf_.at(b))); } }
  for (auto& b : tbvi_) { if (has(b)) { nt_->SetBranchAddress(b.c_str(), &(bvi_.at(b))); } }
  for (auto& b : tbvo_) { if (has(b)) { nt_->SetBranchAddress(b.c_str(), &(bvo_.at(b))); } }
}

template<typename T> T hfupc::dtree::val(const std::string& br, int j) {
  if (!has(br)) {
    std::cout<<__FUNCTION__<<" error: bad branch ["<<br<<"]."<<std::endl;
    return static_cast<T>(0);
  } 
  if (std::is_same<T, float>::value) { return bvf_.at(br)->at(j); }
  if (std::is_same<T, int>::value) { return bvi_.at(br)->at(j); }
  if (std::is_same<T, bool>::value) { return bvo_.at(br)->at(j); }

  std::cout<<__FUNCTION__<<" error: bad type for branch ["<<br<<"]."<<std::endl; 
  return static_cast<T>(0);
}

#endif
