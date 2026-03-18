
template<typename T>
concept _supported_typename =
  std::same_as<T,float>
  || std::same_as<T,int>
  || std::same_as<T,std::string>;

namespace mytmva {
  class mvaroot {
  public:
    mvaroot(TDirectory*);
    mvaroot(const std::string &);

    // hrocs
    std::vector<TH1D*> hrocs(const std::string &pattern = ".+") const;
    TH1D* hroc(const std::string &key) const;
    // info
    template<_supported_typename T = std::string> bool has_info(const std::string &key) const;
    template<_supported_typename T = std::string> T info(const std::string &key) const;
    template<_supported_typename T1 = std::string, typename T2> void add_info(const std::string &key, const T2& value);
    template<_supported_typename T = std::string> void rm_info(const std::string &key);

    // tool
    bool good() const { return good_; }
    void print() const;
    void print_info() const { xjjc::print_tab(info_, -1); xjjc::print_tab(info_f_, -1); xjjc::print_tab(info_i_, -1); }
    void write_hrocs(const std::string &pattern = ".+");
    void write_info(TTree*);
    auto methods() const { return methods_; }
    
  private:
    bool good_;
    std::vector<TH1D*> hrocs_;
    std::map<std::string, std::string> info_;
    std::map<std::string, float> info_f_;
    std::map<std::string, int> info_i_;
    std::vector<std::string> methods_;
  };

  std::vector<mvaroot*> read_rmvas(TDirectory* dir, const std::string &pattern);
}

mytmva::mvaroot::mvaroot(TDirectory* inf)
  : good_(false) {
  if (inf) {
    auto* dataset = (TDirectory*)inf->Get("dataset");
    if (dataset) {
      hrocs_ = xjjana::getobj_regexp_recur<TH1D>(dataset, "MVA_[^_]+_.+", "", false);
      for (const auto h : hrocs_) {
        auto meth = xjjc::str_divide(h->GetName(), "_").at(1);
        if (std::ranges::find(methods_, meth) == methods_.end())
          methods_.push_back(meth);
      }
    }
    auto* tinfo = (TTree*)inf->Get("dataset/tmvainfo");
    if (!tinfo) tinfo = (TTree*)inf->Get("info");
    if (tinfo) {
      info_ = xjjana::getval_regexp(tinfo, ".*");
      info_f_ = xjjana::getval_regexp<float>(tinfo, ".*");
      info_i_ = xjjana::getval_regexp<int>(tinfo, ".*");
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
  if (hs.size() > 1)
    __XJJLOG << "?? more than one hist matching key " << key << std::endl;
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

template<_supported_typename T>
bool mytmva::mvaroot::has_info(const std::string &key) const {
  if constexpr (std::is_same_v<T, std::string>) {
    return info_.find(key) != info_.end();
  } else if constexpr (std::is_same_v<T, float>) {
    return info_f_.find(key) != info_f_.end();
  } else if constexpr (std::is_same_v<T, int>) {
    return info_i_.find(key) != info_i_.end();
  } else {
    __XJJLOG << "!! bad type other than string, float and int." << std::endl;
  }
  return false;
}

template<_supported_typename T>
T mytmva::mvaroot::info(const std::string &key) const {
  if constexpr (std::is_same_v<T, std::string>) {
    return has_info<T>(key) ? info_.at(key) : "";
  } else if constexpr (std::is_same_v<T, float>) {
    return has_info<T>(key) ? info_f_.at(key) : (float)0.;
  } else if constexpr (std::is_same_v<T, int>) {
    return has_info<T>(key) ? info_i_.at(key) : (int)0.;
  } else {
    __XJJLOG << "!! bad type other than string, float and int." << std::endl;
  }
}

template<_supported_typename T1, typename T2>
void mytmva::mvaroot::add_info(const std::string &key, const T2 &value) {
  if constexpr (std::is_same_v<T1, std::string>) {
    if (!has_info<T1>(key)) info_[key] = static_cast<T1>(value); return;
  } else if constexpr (std::is_same_v<T1, float>) {
    if (!has_info<T1>(key)) info_[key] = static_cast<T1>(value); return;
  } else if constexpr (std::is_same_v<T1, int>) {
    if (!has_info<T1>(key)) info_[key] = static_cast<T1>(value); return;
  } else {
    __XJJLOG << "!! bad type other than string, float and int." << std::endl; return;
  }
  __XJJLOG << "!! key " << key << " already exists, abort." << std::endl;
  return;
}

template<_supported_typename T>
void mytmva::mvaroot::rm_info(const std::string &key) {
  if constexpr (std::is_same_v<T, std::string>) {
    if (has_info<T>(key)) info_.erase(key); return;
  } else if constexpr (std::is_same_v<T, float>) {
    if (has_info<T>(key)) info_f_.erase(key); return;
  } else if constexpr (std::is_same_v<T, int>) {
    if (has_info<T>(key)) info_i_.erase(key); return;
  } else {
    __XJJLOG << "!! bad type other than string, float and int." << std::endl;
  }
  __XJJLOG << "!! key " << key << " does not exist, abort." << std::endl;
  return;
}

void mytmva::mvaroot::write_info(TTree* t) {
  for (auto& [key, content] : info_) {
    t->Branch(key.c_str(), &content);
  }
  for (auto& [key, content] : info_f_) {
    t->Branch(key.c_str(), &content);
  }
  for (auto& [key, content] : info_i_) {
    t->Branch(key.c_str(), &content);
  }
}


std::vector<mytmva::mvaroot*> mytmva::read_rmvas(TDirectory* dir, const std::string &pattern) {
  std::vector<mytmva::mvaroot*> rmva;
  TIter next(dir->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    auto* obj = key->ReadObj();
    if (!obj->InheritsFrom(TDirectory::Class())) continue;
    if (!std::regex_match(obj->GetName(), std::regex(pattern))) continue;
    __XJJLOG << ">> " << obj->GetName() << " >> " << static_cast<TDirectory*>(obj)->GetPath() << std::endl;
    auto* rr = new mytmva::mvaroot(static_cast<TDirectory*>(obj));
    if (!rr->good()) continue;
    rr->add_info("dir", obj->GetName());
    rmva.push_back(rr);
    // rr->print();
    // rr->print_info();
  }
  return rmva;
}
