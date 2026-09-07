#include "xjjanauti.h"

int macro(std::string input)
{
  auto* inf = TFile::Open(input.c_str());
  if (!inf) {
    __XJJLOG << "warning: bad input " << input << ", abort."<< std::endl;
    return 2;
  }

  // 
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  info["input_tag"] = xjjc::str_eraseall(xjjc::str_tag_from_file(input), "save_");
  xjjc::print_tab(info, -1);

  //
  auto h3s = xjjana::getobj_regexp<TH3D>(inf, "h3_.+_.+_.+");
  if (h3s.empty() || h3s.size() != 1) {
    __XJJLOG << "!! no good h3 in the file, abort." << std::endl;
    return 2;
  }
  auto* h3 = h3s.front();
  auto axisvars = xjjc::str_divide_trim(xjjc::str_eraseall(h3->GetName(), "h3_"), "_");

  //
  std::map<std::string, TH1D*> h1s;
  std::map<std::string, std::vector<TH1D*>> h1vs;
  h1s["var_overflow"] = (TH1D*)h3->ProjectionY("h1_var_overflow", 0, -1, 0, -1, "e");
  h1s["var"] = (TH1D*)h3->ProjectionY("h1_var", 1, h3->GetNbinsX(), 1, h3->GetNbinsZ(), "e");

  for (int i=0; i<h3->GetXaxis()->GetNbins(); i++) {
    auto* h = (TH1D*)h3->ProjectionY(Form("h1_var__%s-%d", axisvars[0].c_str(), i),
                                     i+1, i+1, 1, h3->GetNbinsZ(), "e");
    h1vs["var"].push_back(h);
  }

  auto* outf = xjjroot::newfile(xjjc::str_replaceall(input, "savehist", "calchist"));
  xjjroot::writehist(h3);
  for (auto& [_, h] : h1s)
    xjjroot::writehist(h);
  for (auto& [_, hh] : h1vs)
    for (auto& h : hh)
      xjjroot::writehist(h);

  auto* tinfo = xjjana::write_info(info);
  tinfo->Fill();
  tinfo->Write();
  outf->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc==2) {
    return macro(argv[1]);
  }
  return 1;
}
