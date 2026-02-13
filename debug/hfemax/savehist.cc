#include "xjjrootuti.h"
#include "xjjanauti.h"
#include "xjjconfig.h"

int macro(std::string configfile,
          std::string var = "") {
  //
  std::cout<<std::endl;
  xjjc::config conf(configfile);
  auto* nt = xjjana::chain_files(conf.get_vec("Inputs"), "Tree");
  auto name = conf.has("Name") ? conf.get("Name") : xjjc::str_gettag_from_file(configfile);
  std::cout<<"\e[1m[ "<<name<<" ]\e[0m"<<std::endl;
  const auto& isgammaN = conf.get<int>("IsgammaN");
  const std::string& v_gamma = isgammaN ? "Plus" : "Minus",
    v_N = isgammaN ? "Minus" : "Plus",
    var_gamma = Form("HFEMax%s%s", v_gamma.c_str(), var.c_str()),
    var_N = Form("HFEMax%s%s", v_N.c_str(), var.c_str());
  const auto& cut = conf.has("Cut") ? conf.get("Cut") : "(1)";
  const auto& cut_gap = cut + " && " + var_gamma + " < " + (isgammaN ? "9.2" : "8.6");
  const auto& cut_gaploose = cut + " && " + var_gamma + " < 15";

  nt->SetBranchStatus("D*", 0);
  nt->SetBranchStatus("G*", 0);
  nt->SetBranchStatus("V*", 0);
  
  if (!nt->FindBranch(var_gamma.c_str()) ||
      !nt->FindBranch(var_N.c_str())) {
    return 2;
  }

  auto* outf = xjjroot::newfile("rootfiles/" + name + var + ".root");

#define SET_AND_PROJECT_HIST(g, c)                                      \
  auto* hHFEMax##g##c = new TH1F("hHFEMax" #g #c, Form(";%s;Events", var##g.c_str()), 100, 0, 25); \
  std::cout<<"\e[1m"<<hHFEMax##g##c->GetName()<<"\e[0m \e[2m<< "<<var##g<<"\e[0m"<<std::endl; \
  std::cout<<cut##c<<std::endl;                                         \
  nt->Project(hHFEMax##g##c->GetName(), var##g.c_str(), cut##c.c_str()); \
  xjjroot::writehist(hHFEMax##g##c);
  
  SET_AND_PROJECT_HIST(_gamma, );
  SET_AND_PROJECT_HIST(_N, );
  SET_AND_PROJECT_HIST(_gamma, _gap);
  SET_AND_PROJECT_HIST(_gamma, _gaploose);
  
  outf->Close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3)
    return macro(argv[1], argv[2]);
  if (argc == 2)
    return macro(argv[1]);
  return 1;
}
