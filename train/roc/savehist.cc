#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"
#include "xjjmypdf.h"

namespace globals {
  const std::vector<std::pair<std::string, int>> mask_br = {
    { "DpassCut23*", 0 },
    { "DpassCut23PAS", 1 },
    { "Dtrk*P", 0 },
    { "Dtrk*MassHypo", 0 },
    { "V*", 0 },
    { "clusComp*", 0 },
  };
}

void rocs(TDirectory*, std::vector<TH1D*>&, std::string); // "rejBvsS"
int macro(std::string configfile_train, std::string configfile_cut) {
  std::cout<<std::endl;
  xjjc::config conf_t(configfile_train), conf_c(configfile_cut);

  auto dir  = conf_t.has("Name") ? conf_t.get("Name") : xjjc::str_tag_from_file(configfile_train);
  auto rootfname = conf_t.get("Rootoutput");
  auto* rootf = TFile::Open(rootfname.c_str());
  if (!rootf) { std::cout<<"bad mva root file."<<std::endl; return 2; }
  auto* dataset = (TDirectory*)rootf->Get("dataset");
  if (!dataset) { std::cout<<"no dir dataset."<<std::endl; return 2; }

  std::vector<TH1D*> hrocs;
  rocs(dataset, hrocs, "rejBvsS");
  xjjc::print_vec_v([&] {
    std::vector<std::string> hnames;
    std::transform(hrocs.begin(), hrocs.end(), std::back_inserter(hnames),
                   [](TH1D* h) { return std::string(h->GetName()); });
    return hnames;
  }(), 0);

  auto* rinfo = (TTree*)rootf->Get("dataset/tmvainfo");
  if (!rinfo) { std::cout<<"no tree dataset/tmvainfo."<<std::endl; return 2; }
  auto inputSname = conf_t.get("Inputs"), inputBname = conf_t.get("Inputb"); // can save in root while training
  TString *pcuts = nullptr, *pcutb = nullptr;
  rinfo->SetBranchAddress("cuts", &pcuts);
  rinfo->SetBranchAddress("cutb", &pcutb);
  rinfo->GetEntry(0);
  std::string cuts(*pcuts), cutb(*pcutb);
  
  auto cut = conf_c.get("Cut");
  auto name = conf_c.has("Name") ? conf_c.get("Name") : xjjc::str_tag_from_file(configfile_cut);
  std::cout<<"\e[1m"<<dir<<" / "<<name<<"\e[0m"<<std::endl;

  auto* trS = xjjana::chain_files({ inputSname }, "Tree");
  if (!trS) { std::cout<<"error: bad input file "<<conf_t.get("InputsS")<<"."<<std::endl; return 2; }
  for (const auto& m : globals::mask_br)
    trS->SetBranchStatus(m.first.c_str(), m.second);
  auto* trB = xjjana::chain_files({ inputBname }, "Tree");
  if (!trB) { std::cout<<"error: bad input file "<<conf_t.get("InputsB")<<"."<<std::endl; return 2; }
  for (const auto& m : globals::mask_br)
    trB->SetBranchStatus(m.first.c_str(), m.second);

  std::map<std::string, long long int> nentries;
  auto count = [&nentries](TTree* tr, std::string type, std::string icut) {
    std::cout<<type<<std::endl;
    std::cout<<"\e[2m"<<icut<<"\e[0m"<<std::endl;
    nentries[type] = tr->Draw("", icut.c_str(), "goff");
    std::cout<<"\e[36m"<<nentries[type]<<"\e[0m"<<std::endl;
  };

  count(trS, "ns_den", cuts);
  count(trS, "ns_num", Form("%s && %s", cuts.c_str(), cut.c_str()));
  std::cout<<nentries.at("ns_num")*1. / nentries.at("ns_den")<<std::endl;
  count(trB, "nb_den", cutb);
  count(trB, "nb_num", Form("%s && !(%s)", cutb.c_str(), cut.c_str()));
  std::cout<<nentries.at("nb_num")*1. / nentries.at("nb_den")<<std::endl;
  auto* grcut = xjjroot::drawpoint(nentries.at("ns_num")*1. / nentries.at("ns_den"),
                                   nentries.at("nb_num")*1. / nentries.at("nb_den"),
                                   kBlack, 47, 4, "goff");
  grcut->SetName("grcut");

  auto* outf = xjjroot::newfile("rootfiles/" + dir + "/" + name + ".root");
  for (auto& h : hrocs)
    xjjroot::writehist(h);
  xjjroot::writehist(grcut);
  auto* t = new TTree("info", "");
  t->Branch("cuts", &cuts);
  t->Branch("cutb", &cutb);
  t->Branch("cut", &cut);
  t->Branch("inputSname", &inputSname);
  t->Branch("inputBname", &inputBname);
  t->Branch("rootfname", &rootfname);
  t->Branch("dir", &dir);
  t->Branch("name", &name);
  t->Fill();
  t->Write();

  outf->Close();
  
  return 0;
  
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}

// https://github.com/root-project/root/blob/master/tmva/tmvagui/src/efficiencies.cxx#L116-L158
// https://github.com/root-project/root/blob/master/tmva/tmvagui/src/tmvaglob.cxx#L590
void rocs(TDirectory* dir, std::vector<TH1D*>& hists, std::string refname/* = "rejBvsS"*/) {
  TIter next(dir->GetListOfKeys());
  TKey* key;

  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class())) {
      rocs((TDirectory*)obj, hists, refname); //
    } else if (obj->InheritsFrom(TH1::Class())) {
      std::string name = obj->GetName();
      if (xjjc::str_contains(name, "MVA_") &&
          xjjc::str_contains(name, refname))
        hists.push_back((TH1D*)obj);
    }
  }
}

