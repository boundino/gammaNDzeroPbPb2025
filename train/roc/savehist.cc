#include <TH3F.h>
#include <TH2F.h>

#include "xjjanauti.h"
#include "xjjconfig.h"

#define __BINS_PTY_EQ__
#define __BINS_MASS__
#include "../../include/bins.h"

#include "mvaroot.h"
#include "save_calc.h"

namespace globals {
  int nbin_high = 0; double minbin_high = 1., maxbin_high = -1.;
}

std::vector<double> interval_to_bins(std::vector<std::pair<double, double>> intervals) {
  if (intervals.empty()) return {};
  std::vector<double> result;
  std::sort(intervals.begin(), intervals.end(),
            [](const auto& a, const auto& b) {
              return a.first < b.first;
            });
  result.push_back(intervals[0].first);
  for (size_t i = 0; i < intervals.size(); ++i) {
    if (i > 0 && intervals[i].first != intervals[i-1].second)
      return {};  // cannot connect
    result.push_back(intervals[i].second);
  }
  return result;
}

template<typename T>
std::vector<double> uniform_bin_vec(int nbin, T binmin, T binmax) {
  std::vector<double> bins;
  bins.reserve(nbin + 1);
  const auto step = static_cast<double>(binmax - binmin) / nbin;
  for (int i = 0; i <= nbin; ++i) {
    bins.push_back(static_cast<double>(binmin) + i * step);
  }
  return bins;
}

int macro(std::string config_train, std::string config_sig) {
  std::cout<<std::endl;
  xjjc::config conf_t(config_train), conf_s(config_sig);

  auto dir = conf_t.has("Name") ? conf_t.get("Name") : xjjc::str_tag_from_file(config_train);
  auto method = conf_s.get("Method"),
    name = conf_s.has("Name") ? conf_s.get("Name") : xjjc::str_tag_from_file(config_sig);
  std::cout<<"\e[1m"<<dir<<" / "<<name<<" [ "<<method<<" ]\e[0m"<<std::endl;

  // train
  std::vector<mytmva::mvaroot*> rmva;
  for (const auto& rfname : conf_t.get_vec("Rootoutputs")) {
    __XJJLOG << ">> " << rfname << std::endl;
    auto* rr = new mytmva::mvaroot(rfname);
    if (!rr->good()) {
      __XJJLOG << "!! bad mva input " << rfname << ", abort." << std::endl;
      return 2;
    }
    if (!( rr->has_info<float>("ptmin") && rr->has_info<float>("ptmax") && rr->has_info<float>("ymin") && rr->has_info<float>("ymax") )) {
      __XJJLOG << "?? " << rfname << " doesn't have [ptmin, ptmax, ymin, ymax], set to inclusive" << std::endl;
      rr->add_info<float>("ptmin", 2); rr->add_info<float>("ptmax", 5); rr->add_info<float>("ymin", -2); rr->add_info<float>("ymax", 2);
    } 
    // rr->print();
    rmva.push_back(rr);
  }
  if (rmva.empty()) {
    __XJJLOG << "!! no good mva root files, abort." << std::endl;
    return 2;
  }

  // reorganize binning and cut
  std::map<std::string, std::string> cut;
  std::vector<std::pair<double, double>> yintervals;
  for (const auto& rr : rmva) {
    // cut
    cut["S"] = (cut.find("S")==cut.end() ? rr->info("cutS") : save::common_cut(cut.at("S"), rr->info("cutS")));
    cut["B"] = (cut.find("B")==cut.end() ? rr->info("cutB") : save::common_cut(cut.at("B"), rr->info("cutB")));
    auto rrcutMass = (rr->has_info("cutData") ? rr->info("cutData") : xjjc::str_eraseall(rr->info("cutB"), " && TMath::Abs(Dmass-1.8648) > 0.05 && TMath::Abs(Dmass-1.8648) < 0.12"));
    cut["Mass"] = (cut.find("Mass")==cut.end() ? rrcutMass : save::common_cut(cut.at("Mass"), rrcutMass));
    // y binning 
    yintervals.push_back({ static_cast<double>(rr->info<float>("ymin")), static_cast<double>(rr->info<float>("ymax")) });
    // mva binning
    auto* high_bin = rr->hroc("MVA_" + method + "_S_high")->GetXaxis();
    // __XJJLOG << ">> " << high_bin->GetNbins() << " , " << high_bin->GetXmin() << " , " << high_bin->GetXmax() << std::endl;
    globals::nbin_high = std::max(globals::nbin_high, high_bin->GetNbins());
    globals::minbin_high = std::min(globals::minbin_high, high_bin->GetXmin());
    globals::maxbin_high = std::max(globals::maxbin_high, high_bin->GetXmax());    
  }
  // cut
  cut["Swap"] = xjjc::str_replaceall(cut.at("S"), "23333", "23344"); // !!
  xjjc::print_tab(cut, -1);
  // y binning
  auto ybins_train = interval_to_bins(yintervals);
  xjjc::print_vec_h(ybins_train, 0);
  if (ybins_train.empty()) {
    __XJJLOG << "!! bad ybins_train, abort." << std::endl; return 3;
  }
  // mva binning
  xjjc::print_vec_h(std::vector<std::string>{ xjjc::to_string(globals::nbin_high), xjjc::to_string(globals::minbin_high), xjjc::to_string(globals::maxbin_high) }, 0);
  if (globals::nbin_high == 0) {
    __XJJLOG << "!! bad universal binning, abort." << std::endl; return 3; 
  }

  // significance
  std::map<std::string, TChain*> trs; std::map<std::string, std::string> inputstr;
  auto prepare_tree = [&trs, &conf_s, &inputstr](const std::string &key, const std::string &Input) {
    trs[key] = xjjana::chain_files(conf_s.get_vec(Input), "Tree");
    if (!trs[key]) { __XJJLOG<<"!! bad input file "<<conf_s.get(Input)<<std::endl; return 2; }
    save::mask_branch(trs[key]);
    inputstr[key] = conf_s.get(Input); // 
    return 0;
  };
  if (prepare_tree("S", "InputsS")) return 2;
  if (prepare_tree("B", "InputsB")) return 2;
  if (prepare_tree("Mass", conf_s.has("InputsMass") ? "InputsMass" : "InputsB")) return 2;

  auto* outf = xjjroot::newfile("rootfiles/" + dir + "/" + name + ".root");

  std::map<std::string, TH3F*> h3_meths;
  auto make_h3 = [&h3_meths, &method, &trs](const std::string &key, TChain* tr, const std::string &icut,
                                            std::vector<double> ybins = uniform_bin_vec(bins::ny, bins::miny, bins::maxy),
                                            std::vector<double> mbins = uniform_bin_vec(globals::nbin_high, globals::minbin_high, globals::maxbin_high)) { // can add href
    if (h3_meths.find(key+"_high") == h3_meths.end()) {
      h3_meths[key+"_high"] = new TH3F(Form("h3_%s_high", key.c_str()), Form(";y;m_{K#pi} (GeV);%s", method.c_str()),
                                       ybins.size()-1, ybins.data(),
                                       bins::nmass, uniform_bin_vec(bins::nmass, bins::minmass, bins::maxmass).data(),
                                       mbins.size()-1, mbins.data());
    }
    __XJJLOG << ">> " << key << " \e[2m[ " << icut << " ]\e[0m" << std::endl;
    tr->Project(h3_meths[key+"_high"]->GetName(), Form("Dmva_%s:Dmass:Dy", method.c_str()), icut.c_str());
    xjjroot::writehist(h3_meths[key+"_high"]);
  };
  // analysis binning
  make_h3("S", trs["S"], cut.at("S"));
  make_h3("B", trs["B"], cut.at("B"));
  make_h3("Swap", trs["S"], cut.at("Swap"));
  make_h3("Mass", trs["Mass"], cut.at("Mass"));
  // training binning
  for (int i=0; i<rmva.size(); i++) {
    auto* rr = rmva.at(i);
    auto xaxis = rr->hroc("MVA_" + method + "_S_high")->GetXaxis();
    make_h3(Form("S_trainbin__y-%d", i), trs["S"], cut.at("S"), ybins_train, uniform_bin_vec(xaxis->GetNbins(), xaxis->GetXmin(), xaxis->GetXmax()));
    make_h3(Form("B_trainbin__y-%d", i), trs["B"], cut.at("B"), ybins_train, uniform_bin_vec(xaxis->GetNbins(), xaxis->GetXmin(), xaxis->GetXmax()));
  }
  for (auto& rr : rmva) {
    auto* dr = outf->mkdir(Form("pt-%s-%s_y-%s-%s",
                                xjjc::number_to_string(rr->info<float>("ptmin")).c_str(),
                                xjjc::number_to_string(rr->info<float>("ptmax")).c_str(),
                                xjjc::number_to_string(rr->info<float>("ymin")).c_str(),
                                xjjc::number_to_string(rr->info<float>("ymax")).c_str()
                                ));
    dr->mkdir("dataset")->cd();
    rr->write_hrocs(".*_" + method + "_.*");
    auto* t = new TTree("tmvainfo", "");
    rr->write_info(t);
    t->Fill();
    t->Write();
    outf->cd();
  }

  auto* t = new TTree("info", "");
  auto branchrval = [&t](const std::string& bname, const std::string& content) {
    auto* tr = new std::string(content);
    t->Branch(bname.c_str(), tr);
  };
  for (auto& [k, c] : cut)
    t->Branch(Form("cut%s", k.c_str()), &c);
  for (auto& [k, c] : inputstr)
    t->Branch(Form("input%sname", k.c_str()), &c);
  t->Branch("method", &method);
  branchrval("label", conf_t.get("Label"));
  branchrval("label2", conf_t.get("Label2"));
  branchrval("scaleto25", conf_s.get("Scale_to_25"));
  t->Branch("config_train", &config_train);
  t->Branch("config_sig", &config_sig);
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
