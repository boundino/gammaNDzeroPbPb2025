#include <string>

#include "TMVAClassification.h"
#include "mvaprod.h"

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

void mvaprob_main(std::string inputname, std::string outputtag, std::string mymethod, std::string stage,
                  std::string outputfiledir) {

  outputtag = mytmva::mkname(outputtag, mymethod, stage);
  
  std::vector<std::map<std::string, mytmva::mvaprod*>> prods(mytmva::nptbins);
  for (int i=0; i<mytmva::nptbins; i++) {
    auto rootfname = outputtag + "_" + mytmva::mkname(mytmva::ptbins.at(i), mytmva::ptbins.at(i+1)) + ".root";
    auto weightdir = "dataset/weights/" + xjjc::str_tag_from_file(rootfname);
    std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": found weight files:"<<mytmva::nocolor<<std::endl;
    for (const auto & entry : fs::directory_iterator(weightdir)) {
      std::string entrypath(entry.path());
      if (!xjjc::str_contains(entrypath, ".weights.xml")) continue;
      std::cout<<entrypath<<std::endl;
      auto* pd = new mytmva::mvaprod(entrypath); //
      if (prods.at(i).find(pd->method()) != prods.at(i).end()) {
        std::cout<<"error: repeated method "<<pd->method()<<", skip."<<std::endl;
        continue;
      }
      pd->getnote(rootfname);
      prods.at(i)[pd->method()] = pd; //
    }
  }

  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": input file:"<<mytmva::nocolor<<std::endl<<inputname<<mytmva::nocolor<<std::endl;
  std::string outfname = outputfiledir + "/" + xjjc::str_tag_from_file(inputname)
    + "_" + xjjc::str_tag_from_file(outputtag) + ".root";
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": output file:"<<mytmva::nocolor<<std::endl<<outfname<<mytmva::nocolor<<std::endl;
  if (fs::exists(outfname)) {
    std::cout<<mytmva::errorcolor<<"==> "<<__FUNCTION__<<": warning: output file already exists."<<mytmva::nocolor<<std::endl; }
  gSystem->Exec(Form("rsync --progress %s %s", inputname.c_str(), outfname.c_str()));
  
  auto* outf = new TFile(outfname.c_str(), "update");
  auto* dir = outf->mkdir("dataset");
  dir->cd();
  auto* info = new TTree("tmvainfo", "TMVA info");
  for (int i=0; i<mytmva::nptbins; i++) {
    prods.at(i).begin()->second->writenote(info, "_" + mytmva::mkname(mytmva::ptbins.at(i), mytmva::ptbins.at(i+1)));
  }
  info->Fill();
  info->Write("", TObject::kOverwrite);

  std::cout<<"start mva value"<<std::endl;
  dir->cd();
  auto* mvatree = new TTree("mva", "TMVA value");
  std::map<std::string, std::vector<float>*> __mvaval;
  for (const auto& [mm, _] : prods.front()) {
    __mvaval[mm] = new std::vector<float>(); //
    mvatree->Branch(mm.c_str(), &(__mvaval.at(mm)));
  }
  
  auto* inf = TFile::Open(inputname.c_str());
  auto* tr = (TTree*)inf->Get("Tree");
  auto* value = new mytmva::varval(tr);
  auto* dnt = value->getnt();

  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Filling mva values:"<<mytmva::nocolor<<std::endl;
  outf->cd();
  auto nentries = dnt->GetEntries();
  for(long long int i=0; i<nentries; i++) {
    xjjc::progressslide(i, nentries, 1000);
    dnt->GetEntry(i);

    for (auto& [_, vval] : __mvaval) {
      vval->clear();
    }
    for (int j=0; j<dnt->Dsize(); j++) {
      int idxpt = mytmva::whichbin(dnt->val("Dpt", j));
      if (idxpt < 0 || idxpt >= prods.size()) {
        std::cout<<"error: bad idxpt: "<<idxpt<<std::endl;
      }
      for (auto& [mm, vval] : __mvaval) {
        auto val = prods.at(idxpt).at(mm)->evalmva(value, j);
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
  if (argc==6) {
    mvaprob_main(argv[1], argv[2], argv[3], argv[4], argv[5]);
    return 0;
  }
  return 1;
}
