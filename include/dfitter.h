#ifndef _XJJROOT_DFITTER_H_
#define _XJJROOT_DFITTER_H_

/******************************************************************************************
 * Class : xjjroot::dfitter                                                               * 
 * This class provides tools fitting Dzero invariant mass spectra.                        * 
 * The object to be used can be declared by                                               * 
 *                                                                                        * 
 *    xjjroot::dfitter obj(options);                                                      * 
 *                                                                                        * 
 * Options supported are listed below                                                     * 
 *                                                                                        * 
 *   "3" : Using 3-Gaussian function to model signal (default is 2-Gaussian function)     * 
 *   "S" : Draw significance info and lines at signal region                              * 
 *   "V" : Switch off Quiet mode of fitting                                               * 
 *                                                                                        * 
 * The core function of this class is                                                     * 
 *                                                                                        * 
 *    TF1* fit(TH1*, TH1*, TH1*, std::string, std::string, std::vector<std::string>)                  * 
 *                                                                                        * 
 ******************************************************************************************/

#include <TString.h>
#include <TMath.h>
#include <TFitResult.h>

namespace xjjroot {
  const std::map<std::string, thgrstyle> fstyle = {
    { "h", thgrstyle(kBlack, 20, 1.3, kBlack, 1, 1) },
    { "f", thgrstyle(-1, -1, -1, 2, 1, 3) },
    { "mass", thgrstyle(-1, -1, -1, kOrange-3, 2, 3, kOrange-3, 0.4, 1001) },
    { "swap", thgrstyle(-1, -1, -1, kGreen+4, 1, 3, kGreen+4, 1, 3005) },
    { "background", thgrstyle(-1, -1, -1, 4, 2, 3) },
    { "notmass", thgrstyle(-1, -1, -1, kGray, 2, 3) },
  };

  class dfitter
  {
  public:
    dfitter(Option_t* option = "") : option_(option) { parse_opt(); reset(); }
    ~dfitter() {};

    void fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped);
    bool fitted() const { return fitted_; }
    static void set_hist(TH1* h);

    std::vector<std::string> draw_result(float x = 0.25, float y = 0.86, float tsize = 0.035, float lspacescale = 1.15) const;
    void draw_params(float x = 0.25, float y = 0.86, float tsize = 0.035, float lspacescale = 1.15) const;
    void draw_fmc() const { fun_mc_swap_->Draw("same"); for (const auto& f : vfun_mc_mass_) { f->Draw("same"); } fun_mc_mass_->Draw("same"); }
    void draw_leg(float x1 = 0.65, float y2 = 0.88) { xjjroot::moveleg_n_draw(leg_, x1, y2); }
    
    double S() const { return S_; }
    double B() const { return B_; }
    double yield() const { return yield_; }
    double yieldErr() const { return yieldErr_; }
    double chi2() const { return 2.*r_->MinFcnValue(); }
    double ndf() const { return fun_f_->GetNDF(); }
    double chi2prob() const { return TMath::Prob(chi2(), ndf()); }

    TF1* f_f(const std::string& name) const;
    TF1* f_mass(const std::string& name = "") const;
    TF1* f_swap(const std::string& name = "") const;
    TF1* f_background(const std::string& name = "") const;
    TF1* f_notmass(const std::string& name = "") const;

  private:
    double S_;
    double B_;
    double yield_;
    double yieldErr_;

    TF1* fun_f_;
    TF1* fun_mc_mass_;
    TF1* fun_mc_swap_;
    std::vector<TF1*> vfun_mc_mass_;
    TLegend* leg_;
    
    TFitResultPtr r_;

    std::string option_;
    bool opt_3gaus_;
    bool opt_sig_;
    bool opt_verbose_;

    // status
    bool fitted_;

    double xmin_ = 0, xmax_ = 0, binwidth_ = 0;
    double signal_region_l_ = 1.8649 - 0.045;
    double signal_region_h_ = 1.8649 + 0.045;
    
    void reset();
    void calculate_SnB();
    void parse_opt();
    void parse_fmc();

    TF1* clone_fun(const TF1* fun, const std::string& fun_name) const;
  };
}

void xjjroot::dfitter::parse_opt() {
  opt_3gaus_ = xjjc::str_contains(option_, "3");
  opt_sig_ = xjjc::str_contains(option_, "S");
  opt_verbose_ = xjjc::str_contains(option_, "V");
}

void xjjroot::dfitter::reset() {
  fitted_ = false;
  // for (auto& s : { S_, B_, yield_, yieldErr_ }) s = -1;
  // for (auto& f : { fun_f_, fun_mc_mass_, fun_mc_swap_ }) f = nullptr;
  S_ = -1;
  B_ = -1;
  yield_ = -1;
  yieldErr_ = -1;
  fun_f_ = nullptr;
  fun_mc_mass_ = nullptr;
  fun_mc_swap_ = nullptr;
  vfun_mc_mass_.clear();
  float tsize = 0.035;
  leg_ = new TLegend(0.6, 0.86-tsize*1.25*5, 0.85, 0.86);
  xjjroot::setleg(leg_, 0.04);
  xjjroot::addentrybystyle(leg_, "Data", "pl", fstyle.at("h"));
  xjjroot::addentrybystyle(leg_, "Fit", "l", fstyle.at("f"));
  xjjroot::addentrybystyle(leg_, "D^{0}+#bar{D^{#lower[0.2]{0}}} Signal", "f", fstyle.at("mass"));
  xjjroot::addentrybystyle(leg_, "K-#pi swapped", "f", fstyle.at("swap"));
  xjjroot::addentrybystyle(leg_, "Combinatorial", "l", fstyle.at("background"));
  leg_->Draw();
}

void xjjroot::dfitter::fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped) {

  reset();

  if (!hmass || !hmassMCSignal || !hmassMCSwapped) {
    __XJJLOG << "!! bad histograms" << std::endl;
    return;
  }
  xmin_ = hmass->GetXaxis()->GetXmin(); //
  xmax_ = hmass->GetXaxis()->GetXmax(); //
  binwidth_ = (xmax_ - xmin_) / hmass->GetXaxis()->GetNbins();
  if (fabs(binwidth_-hmass->GetBinWidth(1))/hmass->GetBinWidth(1) > 1.e-7) {
    __XJJLOG << "!! bad bin width: xmin = " << xmin_
             << ", xmax = " << xmax_
             << ", \e[1mbinwidth = (xmax_ - xmin_) / hmass->GetXaxis()->GetNbins() = " << binwidth_ << "\e[0m"
             << ", hmass->GetBinWidth(1) = " << hmass->GetBinWidth(1)
             << std::endl;
    if (binwidth_ == 0) binwidth_ = hmass->GetBinWidth(1);
  }

  fitted_ = true;
  
  std::string str_fun_f = opt_3gaus_ ?
    "[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*([12]*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11]))+(1-[12])*TMath::Gaus(x,[1],[13]*(1+[11]))/(sqrt(2*3.14159)*[13]*(1+[11]))))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x" :
    "[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11])))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x";
  fun_f_ = new TF1(Form("f_%s", xjjc::unique_str().c_str()), str_fun_f.c_str(), xmin_, xmax_);
  fun_f_->SetNpx(2000);
  xjjroot::setthgrstyle(fun_f_, fstyle.at("f"));
  
  auto* h = (TH1F*)hmass->Clone(Form("h_%s", xjjc::unique_str().c_str()));
  set_hist(h);
  auto* hMCSignal = (TH1F*)hmassMCSignal->Clone(Form("hMCSignal_%s", xjjc::unique_str().c_str()));
  set_hist(hMCSignal);
  auto* hMCSwapped = (TH1F*)hmassMCSwapped->Clone(Form("hMCSwapped_%s", xjjc::unique_str().c_str()));
  set_hist(hMCSwapped);

  const char* fitopt = opt_verbose_?"L m":"L m q";
  
  const double param_init_0 = 100.,
    param_init_1 = 1.8649,
    param_init_2 = 0.03,
    param_init_10 = 0.005,
    param_init_13 = 0.002,
    param_init_8 = 0.1,
    param_init_9 = 0.1,
    param_init_12 = 0.5;

  fun_f_->SetParLimits(4,  -1000, 1000);
  fun_f_->SetParLimits(10, 0.001, 0.05);
  fun_f_->SetParLimits(2,  0.01,  0.5);
  fun_f_->SetParLimits(8,  0.02,  0.2);
  if (opt_3gaus_) fun_f_->SetParLimits(13, 0.002, 0.1);
  fun_f_->SetParLimits(7,  0,     1);
  fun_f_->SetParLimits(9,  0,     1);
  if (opt_3gaus_) fun_f_->SetParLimits(12, 0, 1);
    
  // -- fit MC
  fun_f_->FixParameter(3, 0);
  fun_f_->FixParameter(4, 0);
  fun_f_->FixParameter(5, 0);
  fun_f_->FixParameter(6, 0);

  //  - fit signal
  fun_f_->FixParameter(11, 0);
  fun_f_->FixParameter(1,  param_init_1);
  fun_f_->FixParameter(7,  1);
  fun_f_->FixParameter(8,  1);
  fun_f_->SetParameter(0,  param_init_0);
  fun_f_->SetParameter(1,  param_init_1);
  fun_f_->SetParameter(2,  param_init_2);
  fun_f_->SetParameter(10, param_init_10);
  if(opt_3gaus_) fun_f_->SetParameter(13, param_init_13);
  fun_f_->SetParameter(9,  param_init_9);
  if(opt_3gaus_) fun_f_->SetParameter(12, param_init_12);

  hMCSignal->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  hMCSignal->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  fun_f_->ReleaseParameter(1);
  hMCSignal->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSignal->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSignal->Fit(fun_f_->GetName(), fitopt, "", xmin_, xmax_);

  fun_f_->FixParameter(1, fun_f_->GetParameter(1));
  fun_f_->FixParameter(2, fun_f_->GetParameter(2));
  fun_f_->FixParameter(10, fun_f_->GetParameter(10));
  if (opt_3gaus_) fun_f_->FixParameter(13, fun_f_->GetParameter(13));
  fun_f_->FixParameter(9, fun_f_->GetParameter(9));
  if (opt_3gaus_) fun_f_->FixParameter(12, fun_f_->GetParameter(12));

  fun_mc_mass_ = f_mass(Form("%s_mass_mc", fun_f_->GetName()));
  xjjroot::setthgrstyle(fun_mc_mass_, -1, -1, -1, fstyle.at("mass").lcolor, 1, fstyle.at("mass").lwidth, 0, 0, 0);
  const auto fixparam7 = fun_f_->GetParameter(0);

  //   - fit swapped
  fun_f_->FixParameter(7, 0);
  fun_f_->ReleaseParameter(8);
  fun_f_->SetParameter(8, param_init_8);
  
  hMCSwapped->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSwapped->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSwapped->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSwapped->Fit(fun_f_->GetName(), fitopt,"", xmin_, xmax_);
  
  fun_f_->FixParameter(8, fun_f_->GetParameter(8));

  fun_mc_swap_ = f_swap(Form("%s_mc_swap", fun_f_->GetName()));
  xjjroot::setthgrstyle(fun_mc_swap_, -1, -1, -1, fstyle.at("swap").lcolor, 1, fstyle.at("swap").lwidth, 0, 0, 0);

  parse_fmc();
  
  //  -- fit data
  fun_f_->FixParameter(7, fixparam7/(fun_f_->GetParameter(0)+fixparam7));
  fun_f_->ReleaseParameter(3);
  fun_f_->ReleaseParameter(4);
  fun_f_->ReleaseParameter(5);
  fun_f_->ReleaseParameter(6);
  
  h->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  h->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  fun_f_->ReleaseParameter(1);
  fun_f_->SetParLimits(1, 1.86, 1.87);
  // fun_f_->ReleaseParameter(11);
  // fun_f_->SetParLimits(11, -0.5, 0.5);
  h->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  h->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  h->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  r_ = h->Fit(fun_f_->GetName(), Form("%s S", fitopt),"", xmin_, xmax_);

  auto* fun_background = f_background();
  auto* fun_mass = f_mass();
  auto* fun_swap = f_swap();
  auto* fun_notmass = f_notmass();

  yield_ = fun_mass->Integral(xmin_, xmax_)/binwidth_;
  yieldErr_ = fun_mass->Integral(xmin_, xmax_)/binwidth_*fun_mass->GetParError(0)/fun_mass->GetParameter(0);

  calculate_SnB();
  
  h->Draw("e");
  fun_mass->Draw("same");
  fun_background->Draw("same");
  fun_swap->Draw("same");
  if (opt_sig_) {
    fun_notmass->SetRange(signal_region_l_, signal_region_h_);
    fun_notmass->Draw("same");
    xjjroot::drawline(signal_region_l_, 0, signal_region_l_, fun_f_->Eval(signal_region_l_), fun_notmass->GetLineColor(), fun_notmass->GetLineStyle(), fun_notmass->GetLineWidth());
    xjjroot::drawline(signal_region_h_, 0, signal_region_h_, fun_f_->Eval(signal_region_h_), fun_notmass->GetLineColor(), fun_notmass->GetLineStyle(), fun_notmass->GetLineWidth());
  }
  fun_f_->Draw("same");
}

std::vector<std::string> xjjroot::dfitter::draw_result(float x, float y, float tsize, float lspacescale) const {
  std::vector<std::string> text = {
    Form("N = %.0f#scale[0.5]{ }#pm %.0f", yield_, yieldErr_),
    Form("#chi^{2} / ndf = %.1f / %.0f", chi2(), ndf()),
    Form("Prob = %.2f%s", chi2prob()*100, "%"),
  };
  std::vector<Color_t> color(text.size(), kBlack);
  if (opt_sig_) {
    text.push_back(Form("S = %.0f, B = %.0f", S_, B_));
    text.push_back(Form("S/#sqrt{S+B} = %.1f", S_/TMath::Sqrt(S_ + B_)));
    xjjc::vec_append(color, std::vector<Color_t>(2, fstyle.at("notmass").lcolor));
  }
  xjjroot::drawtexgroup(x, y, text, tsize, 13, 42, lspacescale,
                        1, 0.2, color);

  // if (fdrawdetail) {
  //   drawtex(texxpos, texypos=(texypos-texdypos), Form("N#scale[0.6]{#lower[0.7]{sig}}/(N#scale[0.6]{#lower[0.7]{sig}}+N#scale[0.6]{#lower[0.7]{swap}}) = %.2f",fun_f_->GetParameter(7)));
  // }
  return text;
}

void xjjroot::dfitter::calculate_SnB() {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return;
  }
  auto* fun_mass = f_mass(Form("fun_mass_%s", xjjc::unique_str().c_str()));
  S_ = fun_mass->Integral(signal_region_l_, signal_region_h_)/binwidth_;
  delete fun_mass;
  auto* fun_notmass = f_notmass(Form("fun_notmass_%s", xjjc::unique_str().c_str()));
  B_ = fun_notmass->Integral(signal_region_l_, signal_region_h_)/binwidth_;
  delete fun_notmass;
}

TF1* xjjroot::dfitter::f_mass(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_mass", fun_f_->GetName()) : name;
  std::string str_fun_mass = "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*([7]*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[6]))/(sqrt(2*3.14159)*[8]*(1+[6])))))";
  auto* fun = new TF1(fname.c_str(), str_fun_mass.c_str(), fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 7 }, { 4, 9 }, { 5, 10 }, { 6, 11 }
  };
  if (opt_3gaus_) { params[7] = 12; params[8] = 13; }
  else {
    fun->SetParameter(7, 1); fun->SetParError(7, 0);
    fun->SetParameter(8, 0.1); fun->SetParError(8, 0);
  }
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  fun->SetNpx(2000);
  xjjroot::setthgrstyle(fun, fstyle.at("mass"));

  return fun;
}

TF1* xjjroot::dfitter::f_swap(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_swap", fun_f_->GetName()) : name;
  auto* fun = new TF1(fname.c_str(), "[0]*(1-[2])*TMath::Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 1 }, { 2, 7 }, { 3, 8 }, { 4, 11 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("swap"));

  return fun;
}

TF1* xjjroot::dfitter::f_background(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_background", fun_f_->GetName()) : name;
  auto* fun = new TF1(fname.c_str(), "[0]+[1]*x+[2]*x*x+[3]*x*x*x", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 3 }, { 1, 4 }, { 2, 5 }, { 3, 6 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("background"));

  return fun;
}

TF1* xjjroot::dfitter::f_notmass(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_notmass", fun_f_->GetName()) : name;
  auto* fun = new TF1(name.c_str(), "[0]*(1-[2])*TMath::Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))+[5]+[6]*x+[7]*x*x+[8]*x*x*x", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 1 }, { 2, 7 }, { 3, 8 }, { 4, 11 },
    { 5, 3 }, { 6, 4 }, { 7, 5 }, { 8, 6 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("notmass"));

  return fun;
}

TF1* xjjroot::dfitter::clone_fun(const TF1* fun, const std::string& name) const {
  auto* newfun = new TF1(*fun);
  newfun->SetName(name.c_str());
  return newfun;
}

void xjjroot::dfitter::set_hist(TH1* h) {
  h->SetXTitle("m_{#piK} (GeV/c^{2})");
  h->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", h->GetBinWidth(1)*1.e+3));
  xjjroot::sethempty(h, 0.045);
  xjjroot::setthgrstyle(h, fstyle.at("h"));
  h->SetMaximum(-1111);
  h->SetAxisRange(0, h->GetMaximum()*1.4*1.2, "Y");
}

void xjjroot::dfitter::parse_fmc() {
  std::string str_fun_mass = "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*([7]*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[6]))/(sqrt(2*3.14159)*[8]*(1+[6])))))";
  std::vector<std::vector<double>> params = {
    { fun_mc_mass_->GetParameter(0)*fun_mc_mass_->GetParameter(3)*fun_mc_mass_->GetParameter(4), fun_mc_mass_->GetParameter(1), fun_mc_mass_->GetParameter(2) },
    { fun_mc_mass_->GetParameter(0)*fun_mc_mass_->GetParameter(3)*(1-fun_mc_mass_->GetParameter(4))*fun_mc_mass_->GetParameter(7), fun_mc_mass_->GetParameter(1), fun_mc_mass_->GetParameter(5) },    
  };
  if (opt_3gaus_) {
    params.push_back( { fun_mc_mass_->GetParameter(0)*fun_mc_mass_->GetParameter(3)*(1-fun_mc_mass_->GetParameter(4))*(1-fun_mc_mass_->GetParameter(7)), fun_mc_mass_->GetParameter(1), fun_mc_mass_->GetParameter(8) } );
  }
  for (int i=0; i<params.size(); i++) {
    auto* f1 = new TF1(Form("%s-%d", fun_mc_mass_->GetName(), i), "[0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])", xmin_, xmax_);
    for (int j=0; j<params.at(i).size(); j++) {
      f1->SetParameter(j, params.at(i).at(j));
      f1->SetParError(j, 0);
    }
    xjjroot::setlinestyle(f1, fun_mc_mass_->GetLineColor(), 2, fun_mc_mass_->GetLineWidth());
    vfun_mc_mass_.push_back(f1);
  }
}

void xjjroot::dfitter::draw_params(float x, float y, float tsize, float lspacescale) const {
  std::string sigma_mass, frac_mass;
  float norm = 0;
  for (const auto& f : vfun_mc_mass_) {
    sigma_mass += std::string(Form("%s%.3f", (sigma_mass.empty() ? "" : ", "), f->GetParameter(2)));
    norm += f->GetParameter(0);
  }
  for (const auto& f : vfun_mc_mass_) {
    frac_mass += std::string(Form("%s%.2f", (frac_mass.empty() ? "" : ", "), f->GetParameter(0) / norm));
  }
  std::vector<std::string> rtex = {
    "Mean (data) = " + std::string(Form("%.3f#pm%.3f", fun_f_->GetParameter(1), fun_f_->GetParError(1))),
    "Signal#scale[0.5]{ }#sigma in MC = " + sigma_mass,
    "#sigma_{data}/#sigma_{MC} - 1 = " + std::string(Form(" %.2f#pm%.2f", fun_f_->GetParameter(11), fun_f_->GetParError(11))),
    "Fraction of each gaus = " + frac_mass,
    Form("N_{sig}/(N_{sig}+N_{swap}) = %.2f", fun_f_->GetParameter(7)),
  };
  xjjroot::drawtexgroup(x, y, rtex, tsize, 13, 42, lspacescale);
}

#endif
