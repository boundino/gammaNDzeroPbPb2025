
namespace mytmva {
  class mvaroot {
  public:
    mvaroot(const std::string &inputname);

    // hrocs
    std::vector<TH1D*> hrocs(const std::string &pattern) const;
    TH1D* hroc(const std::string &key) const;
    // info
    bool has_info(const std::string &key) const { return info_.find(key) != info_.end(); }
    std::string info(const std::string &key) const { return has_info(key) ? info_.at(key) : ""; }

    // tool
    bool good() const { return good_; }
    void print() const;
    void write_hrocs(const std::string &pattern);
    
  private:
    bool good_;
    std::vector<TH1D*> hrocs_;
    std::map<std::string, std::string> info_;
    
    std::map<std::string, int> methods_;

    void hrocs_from_file(TDirectory* dir, const std::string &pattern);
  };
}

mytmva::mvaroot::mvaroot(const std::string &inputname) :
  good_(false) {
  auto* inf = TFile::Open(inputname.c_str());
  if (inf) {
    auto* dataset_ = (TDirectory*)inf->Get("dataset");
    if (dataset_) {
      hrocs_from_file(dataset_, "MVA_[^_]+_.+");
      auto* rinfo = (TTree*)inf->Get("dataset/tmvainfo");
      if (rinfo) {
        info_ = xjjana::getstr_regexp(rinfo, ".*");
        info_["rootfile"] = inputname;
        good_ = true;
      }
    }
  }
}

std::vector<TH1D*> mytmva::mvaroot::hrocs(const std::string &pattern) const {
  std::regex re(pattern);
  std::vector<TH1D*> result;
  for (const auto& h : hrocs_) {
    if (std::regex_match(h->GetName(), re)) {
      result.push_back(h);
    }
  }
  return result;
}

TH1D* mytmva::mvaroot::hroc(const std::string &key) const {
  TH1D* result = nullptr;
  auto hs = hrocs(key);
  return hs.size() > 0 ? hs.front() : nullptr;
}

// https://github.com/root-project/root/blob/master/tmva/tmvagui/src/efficiencies.cxx#L116-L158
// https://github.com/root-project/root/blob/master/tmva/tmvagui/src/tmvaglob.cxx#L590
void mytmva::mvaroot::hrocs_from_file(TDirectory* dir, const std::string &pattern) {
  TIter next(dir->GetListOfKeys());
  TKey* key;
  std::regex re(pattern);
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class())) {
      hrocs_from_file((TDirectory*)obj, pattern); //
    } else if (obj->InheritsFrom(TH1::Class())) {
      if (std::regex_match(obj->GetName(), re))
	hrocs_.push_back((TH1D*)obj);
    }
  }  
  for (const auto& h : hrocs_) {
    const auto m = xjjc::str_divide(h->GetName(), "_").at(1); // MVA_[method]_*
    if (methods_.find(m) == methods_.end()) methods_[m] = 0; //
    else if (xjjc::str_contains(h->GetName(), "_S_high")) methods_[m]++;
  }
}

void mytmva::mvaroot::print() const {
  xjjc::print_vec_v([&] {
    std::vector<std::string> hnames;
    std::transform(hrocs_.begin(), hrocs_.end(), std::back_inserter(hnames),
                   [](TH1D* h) { return std::string(h->GetName()); });
    return hnames;
  }(), 0);

  xjjc::print_tab(info_, -1);
}

void mytmva::mvaroot::write_hrocs(const std::string &pattern) {
  std::regex re(pattern);
  for (auto& h : hrocs_) {
    if (std::regex_match(h->GetName(), re))
      xjjroot::writehist(h);
  }
}
