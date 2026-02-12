#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "xjjconfig.h"

int macro(std::string configfile) {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto* nt = xjjana::chain_files(conf.get_vec("Inputs"), "Tree");
  auto name = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);
  std::cout<<"\e[1m[ "<<name<<" ]\e[0m"<<std::endl;
  const auto& isgammaN = conf.get<int>("IsgammaN");
  const std::string& v_gamma = isgammaN ? "Plus" : "Minus",
    v_N = isgammaN ? "Minus" : "Plus";
  const auto& cut = conf.has("Cut") ? conf.get("Cut") : "(1)";
  const auto& cut_gap = cut + " && HFEMax" + v_gamma + "%s < " + (isgammaN ? "9.2" : "8.6");
  const auto& cut_gaploose = cut + " && HFEMax" + v_gamma + "%s < 15";

  nt->SetBranchStatus("D*", 0);
  nt->SetBranchStatus("G*", 0);
  nt->SetBranchStatus("V*", 0);
  
  auto* outf = xjjroot::newfile("rootfiles/" + name + ".root");

#define SET_AND_PROJECT_HIST(g, var, c)                                 \
  std::cout<<"\e[1m"<<"hHFEMax" #g #var #c<<"\e[0m \e[2m<< "<<Form("HFEMax%s" #var, v##g.c_str())<<"\e[0m"<<std::endl; \
  auto* hHFEMax##g##var##c = new TH1F("hHFEMax" #g #var #c, Form(";HFEMax%s;Events", v##g.c_str()), 100, 0, 25); \
  std::cout<<xjjc::str_replaceall(cut##c, "%s", #var)<<std::endl;       \
  nt->Project(hHFEMax##g##var##c->GetName(), Form("HFEMax%s" #var, v##g.c_str()), xjjc::str_replaceall(cut##c, "%s", #var).c_str()); \
  xjjroot::writehist(hHFEMax##g##var##c);
  
  SET_AND_PROJECT_HIST(_gamma, , );
  SET_AND_PROJECT_HIST(_gamma, _forest, );
  SET_AND_PROJECT_HIST(_gamma, _eta5, );
  SET_AND_PROJECT_HIST(_N, , );
  SET_AND_PROJECT_HIST(_N, _forest, );
  SET_AND_PROJECT_HIST(_N, _eta5, );
  SET_AND_PROJECT_HIST(_gamma, , _gap);
  SET_AND_PROJECT_HIST(_gamma, _forest, _gap);
  SET_AND_PROJECT_HIST(_gamma, _eta5, _gap);
  SET_AND_PROJECT_HIST(_gamma, , _gaploose);
  SET_AND_PROJECT_HIST(_gamma, _forest, _gaploose);
  SET_AND_PROJECT_HIST(_gamma, _eta5, _gaploose);
  
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

