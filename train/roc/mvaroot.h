
namespace mytmva {
  class mvaroot {
  public:
    mvaroot(TFile*);
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
    void print_info() const { xjjc::print_tab(info_, -1); }
    void write_hrocs(const std::string &pattern);
    void write_info(TTree*);
    auto methods() const { return methods_; }
    
  private:
    bool good_;
    std::vector<TH1D*> hrocs_;
    std::map<std::string, std::string> info_;
    std::vector<std::string> methods_;
  };
}

mytmva::mvaroot::mvaroot(TFile* inf)
  : good_(false) {
  if (inf) {
    auto* dataset = (TDirectory*)inf->Get("dataset");
    if (dataset) {
      hrocs_ = xjjana::getobj_regexp_recur<TH1D>(dataset, "MVA_[^_]+_.+");
      for (const auto h : hrocs_) {
        auto meth = xjjc::str_divide(h->GetName(), "_").at(1);
        if (std::ranges::find(methods_, meth) == methods_.end())
          methods_.push_back(meth);
      }
    }
    auto* tinfo = (TTree*)inf->Get("dataset/tmvainfo");
    if (!tinfo) tinfo = (TTree*)inf->Get("info");
    if (tinfo) {
      info_ = xjjana::getstr_regexp(tinfo, ".*");
    }
    good_ = dataset && tinfo;
  }
}

mytmva::mvaroot::mvaroot(const std::string &inputname)
  : mvaroot(TFile::Open(inputname.c_str())) {
  if (!has_info("rootfile"))
    info_["rootfile"] = inputname;
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

void mytmva::mvaroot::print() const {
  xjjc::print_vec_v([&] {
    std::vector<std::string> hnames;
    std::transform(hrocs_.begin(), hrocs_.end(), std::back_inserter(hnames),
                   [](TH1D* h) { return std::string(h->GetName()); });
    return hnames;
  }(), 0);

}

void mytmva::mvaroot::write_hrocs(const std::string &pattern) {
  std::regex re(pattern);
  for (auto& h : hrocs_) {
    if (std::regex_match(h->GetName(), re))
      xjjroot::writehist(h);
  }
}

void mytmva::mvaroot::write_info(TTree* t) {
  for (auto& [key, content] : info_) {
    t->Branch(key.c_str(), &content);
  }
}
