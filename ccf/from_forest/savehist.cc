#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "config.h"

#include "clusComp.h"

double determineQuality(const clusComp& cc, double minZ, double maxZ);
bool filter(const double& quality, const double& nPxlHits);

int macro(std::string configfile) {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto* nt = xjjana::chain_files(conf.get_vec("Inputs"), "hiEvtAnalyzer/HiTree");
  auto outputname = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);

  clusComp cc;
  cc.nPixHits = 0; nt->SetBranchAddress("clusComp_nPixHits", &cc.nPixHits);
  cc.z0 = nullptr; nt->SetBranchAddress("clusComp_z0", &cc.z0);
  cc.nHit = nullptr; nt->SetBranchAddress("clusComp_nHit", &cc.nHit);
  cc.chi = nullptr; nt->SetBranchAddress("clusComp_chi", &cc.chi);

  auto* h_qual_nhit = new TH2F("h_qual_nhit", ";Number of pixel hits;Cluster-vertex compatibility", 50, 0, 1500, 50, 0, 10);
  auto* h_qual_nhit_pass = new TH2F("h_qual_nhit_pass", ";Number of pixel hits;Cluster-vertex compatibility", 50, 0, 1500, 50, 0, 10);
  
  const auto& nentries = nt->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    nt->GetEntry(i);

    if (!i) {
      std::map<std::string, int> size = {
        { "nPixHits", cc.nPixHits },
        { "z0.size", cc.z0->size() },
        { "nHit.size", cc.nHit->size() },
        { "chi.size", cc.chi->size() },
      }; xjjc::print_tab(size);
    } 
    
    xjjc::progressslide(i, nentries);
    auto quality = determineQuality(cc, -20.0, 20.05);
    h_qual_nhit->Fill(cc.nPixHits, quality);
    if (filter(quality, cc.nPixHits)) {
      h_qual_nhit_pass->Fill(cc.nPixHits, quality);
    }
  }
  xjjc::progressbar_summary(nentries);

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(h_qual_nhit);
  xjjroot::writehist(h_qual_nhit_pass);

  delete cc.z0;
  delete cc.nHit;
  delete cc.chi;

  outf->Write();
  outf->Close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

bool filter(const double& quality, const double& nPxlHits) {
  bool accept = true;
  // construct polynomial cut on cluster vertex quality vs. npixelhits
  double polyCut = 0;
  for (unsigned int i = 0; i < globals::clusterPars_.size(); i++) {
    polyCut += globals::clusterPars_[i] * std::pow((double)nPxlHits, (int)i);
  }
  if (nPxlHits < globals::nhitsTrunc_)
    polyCut = 0;  // don't use cut below nhitsTrunc_ pixel hits
  if (polyCut > globals::clusterTrunc_ && globals::clusterTrunc_ > 0)
    polyCut = globals::clusterTrunc_;  // no cut above clusterTrunc_

  if (quality < polyCut)
    accept = false;

  // return with final filter decision
  return accept;
}

// https://github.com/CmsHI/cmssw/blob/forest_CMSSW_15_1_X/HeavyIonsAnalysis/EventAnalysis/plugins/HIClusterCompatibilityFilter.cc
double determineQuality(const clusComp& cc, double minZ, double maxZ) {
  // will compare cluster compatibility at a determined best
  // z position to + and - 10 cm from the best position

  float best_z = 0.;
  int best_n = 0., low_n = 0., high_n = 0.;

  // look for best vertex z position within zMin to zMax range
  // best position is determined by maximum nHit with
  // chi used for breaking a tie
  int nhits_max = 0;
  double chi_max = 1e+9;
  for (int i = 0; i < cc.z0->size(); i++) {
    if (cc.z0->at(i) > maxZ || cc.z0->at(i) < minZ)
      continue;
    if (cc.nHit->at(i) == 0)
      continue;
    if (cc.nHit->at(i) > nhits_max) {
      chi_max = 1e+9;
      nhits_max = cc.nHit->at(i);
    }
    if (cc.nHit->at(i) >= nhits_max && cc.chi->at(i) < chi_max) {
      chi_max = cc.chi->at(i);
      best_z = cc.z0->at(i);
      best_n = cc.nHit->at(i);
    }
  }

  // find compatible clusters at + or - 10 cm of the best,
  // or get as close as possible in terms of z position.
  double low_target = best_z - 10.0;
  double high_target = best_z + 10.0;
  double low_match = 1000., high_match = 1000.;
  for (int i = 0; i < cc.z0->size(); i++) {
    if (fabs(cc.z0->at(i) - low_target) < low_match) {
      low_n = cc.nHit->at(i);
      low_match = fabs(cc.z0->at(i) - low_target);
    }
    if (fabs(cc.z0->at(i) - high_target) < high_match) {
      high_n = cc.nHit->at(i);
      high_match = fabs(cc.z0->at(i) - high_target);
    }
  }
  
  // determine vertex compatibility quality score
  double clusVtxQual = 0.0;
  if ((low_n + high_n) > 0)
    clusVtxQual = (2.0 * best_n) / (low_n + high_n);  // A/B
  else if (best_n > 0)
    clusVtxQual = 1000.0;  // A/0 (set to arbitrarily large number)
  else
    clusVtxQual = 0;

  return clusVtxQual;
}
