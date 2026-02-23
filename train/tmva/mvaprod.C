#include <string>

#include "TMVAClassification.h"
#include "mvaprod.h"

void mvaprob_main(std::string inputname, std::string treename, std::string outputtag, std::string outputfiledir,
                  std::string mymethod = "", std::string stage = "0,1,2,3,4,5,6,7,8,9,10") {

  outputtag = mytmva::mkname(outputtag, mymethod, stage);
  
  std::vector<std::map<std::string, mytmva::mvaprod*>> prods(mytmva::nptbins);
  for (int i=0; i<mytmva::nptbins; i++) {
    auto rootfname = outputtag + "_" + mytmva::mkname(mytmva::ptbins[i], mytmva::ptbins[i+1]);
    auto weightdir = "dataset/weights/" + xjjc::str_tag_from_file(rootfname);
    std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": found weight files:"<<mytmva::nocolor<<std::endl;
    for (const auto & entry : fs::directory_iterator(weightdir)) {
      std::string entrypath(entry.path());
      if (!xjjc::str_contains(entrypath, ".weights.xml")) continue;
      std::cout<<entrypath<<std::endl;
      auto* pd = new mytmva::mvaprod(entrypath);
      if (prods.find(pd->method()) != prods.end()) {
        std::cout<<"error: repeated method "<<pd->method()<<", skip."<<std::endl;
        continue;
      }
      pd->getnote(rootfname);
      prods[i][pd->method()] = pd;
    }
  }

  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": input file:"<<mytmva::nocolor<<std::endl<<inputname<<mytmva::nocolor<<std::endl;
  std::string outfname = outputfiledir + "/" + xjjc::str_tag_from_file(inputname)
    + "_" + xjjc::str_tag_from_file(outputtag);
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": output file:"<<mytmva::nocolor<<std::endl<<outfname<<mytmva::nocolor<<std::endl;
  if (std::experimental::filesystem::exists(outfname)) {
    std::cout<<mytmva::errorcolor<<"==> "<<__FUNCTION__<<": warning: output file already exists."<<mytmva::nocolor<<std::endl; }
  gSystem->Exec(Form("rsync --progress %s %s", inputname.c_str(), outfname.c_str()));
  
  auto* outf = TFile::Open(outfname.c_str(), "update");
  auto* dir = outf->mkdir("dataset");
  dir->cd();
  auto* info = new TTree("tmvainfo", "TMVA info");
  for (int i=0; i<mytmva::nptbins; i++) {
    prods[i].begin()->second->writenote(info, "_" + mytmva::mkname(mytmva::ptbins[i], mytmva::ptbins[i+1]));
  }
  info->Fill();
  info->Write("", TObject::kOverwrite);

  dir->cd();
  auto* mvatree = new TTree("mva", "TMVA value");
  std::map<std::string, std::vector<float>*> __mvaval;
  for (const auto& [mm, _] : prods.front()) {
    __mvaval[mm] = new std::vector<float>(); //
    mvatree->Branch(mm.c_str(), &(__mvaval[mm]));
  }
  
  auto* inf = TFile::Open(inputname.c_str());
  auto* value = new mytmva::varval( (TTree*)inf->Get("Tree") );
  auto* dnt = value->getnt();

  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Filling mva values:"<<mytmva::nocolor<<std::endl;
  outf->cd();
  int nentries = dnt->GetEntries();
  for(int i=0; i<nentries; i++) {
    xjjc::progressslide(i, nentries, 1000);
    dnt->GetEntry(i);

    for (auto& [mm, vval] : __mvaval) {
      __mvaval[mm]->clear();
    }
    for (int j=0; j<dnt->Dsize(); j++) {
      int idxpt = mytmva::whichbin(dnt->val("Dpt", j));
      if (idxpt < 0 || idxpt >= prods.size()) {
        std::cout<<"error: bad idxpt: "<<idxpt<<std::endl;
      }
      for (auto& [mm, vval] : __mvaval) {
        auto val = prods[idxpt][mm]->evalmva(value, j);
        vval->push_back(val);
      }
    }
    mvatree->Fill();
  }
  xjjc::progressbar_summary(nentries);

  dir->cd();
  mvatree->Write("", TObject::kOverwrite);
  outf->cd();
  outf->Close();
}

int main(int argc, char* argv[]) {
  if (argc==7) {
    mvaprob_main(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
    return 0;
  }
  return 1;
}
