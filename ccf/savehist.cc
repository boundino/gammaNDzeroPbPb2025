#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "xjjconfig.h"

#include "clusComp.h"

bool filter(const double& quality, const double& nPxlHits);
int macro(std::string configfile) {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto* nt = xjjana::chain_files(conf.get_vec("Inputs"), "Tree");
  auto outputname = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);
  auto cut = conf.has("Cut") ? conf.get("Cut") : "(1)";
  auto cut_pass = cut + " && " + globals::makestrcut("clusComp_nPixHits", "clusComp_quality");
  auto cut_filter = cut + " && ClusterCompatibilityFilter";

  nt->SetBranchStatus("D*", 0);
  nt->SetBranchStatus("G*", 0);
  
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");

#define SET_AND_PROJECT_HIST(q)                                         \
  auto* h_qual_nhit##q = new TH2F("h_qual_nhit" #q, ";Number of Pixel Hits;Cluster-Vertex Compatibility", 50, 0, 1500, 50, 0, 10); \
  std::cout<<cut##q<<std::endl;                                         \
  nt->Project(h_qual_nhit##q->GetName(), "clusComp_quality:clusComp_nPixHits", cut##q.c_str()); \
  xjjroot::writehist(h_qual_nhit##q);

  SET_AND_PROJECT_HIST();
  SET_AND_PROJECT_HIST(_pass);
  SET_AND_PROJECT_HIST(_filter);
  
  outf->Close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

// bool filter(const double& quality, const double& nPxlHits) {
//   bool accept = true;
//   // construct polynomial cut on cluster vertex quality vs. npixelhits
//   double polyCut = 0;
//   for (unsigned int i = 0; i < globals::clusterPars_.size(); i++) {
//     polyCut += globals::clusterPars_[i] * std::pow((double)nPxlHits, (int)i);
//   }
//   if (nPxlHits < globals::nhitsTrunc_)
//     polyCut = 0;  // don't use cut below nhitsTrunc_ pixel hits
//   if (polyCut > globals::clusterTrunc_ && globals::clusterTrunc_ > 0)
//     polyCut = globals::clusterTrunc_;  // no cut above clusterTrunc_

//   if (quality < polyCut)
//     accept = false;

//   // return with final filter decision
//   return accept;
// }

