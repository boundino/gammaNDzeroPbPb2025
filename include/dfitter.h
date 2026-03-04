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
    { "f", thgrstyle(-1, -1, -1, 2, 1, 3) },
    { "mass", thgrstyle(-1, -1, -1, kOrange-3, 2, 3, kOrange-3, 0.4, 1001) },
    { "swap", thgrstyle(-1, -1, -1, kGreen+4, 1, 3, kGreen+4, 1, 3005) },
    { "background", thgrstyle(-1, -1, -1, 4, 2, 3) },
    { "notmass", thgrstyle(-1, -1, -1, kGray, 2, 3) },
  };

  class dfitter
  {
  public:
    dfitter(Option_t* option = "") : option_(option) { parseopt(); reset(); }
    ~dfitter() {};

    void fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped);
    bool fitted() const { return fitted_; }
    static void sethist(TH1* h);

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

    TF1* f_mc_mass() const { return fun_mass_mc_; }
    TF1* f_mc_swap() const { return fun_swap_mc_; }
    
    std::vector<std::string> draw_result(float x = 0.22, float y = 0.80, float tsize = 0.04, float lspacescale = 1.15) const;
    
  private:
    double S_;
    double B_;
    double yield_;
    double yieldErr_;

    TF1* fun_f_;
    TF1* fun_mass_mc_;
    TF1* fun_swap_mc_;
    
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
    void parseopt();

    TF1* clonefun(const TF1* fun, const std::string& fun_name) const;
    void drawleg(TH1* h) const;
  };
}

void xjjroot::dfitter::parseopt() {
  opt_3gaus_ = xjjc::str_contains(option_, "3");
  opt_sig_ = xjjc::str_contains(option_, "S");
  opt_verbose_ = xjjc::str_contains(option_, "V");
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
  if (binwidth_ != hmass->GetBinWidth(1)) {
    __XJJLOG << "!! bad bin width" << std::endl;
    return;
  }

  std::string str_fun_f = opt_3gaus_ ?
    "[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*([12]*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11]))+(1-[12])*TMath::Gaus(x,[1],[13]*(1+[11]))/(sqrt(2*3.14159)*[13]*(1+[11]))))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x" :
    "[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11])))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x";
  fun_f_ = new TF1(Form("f_%s", xjjc::unique_str().c_str()), str_fun_f.c_str(), xmin_, xmax_);
  fun_f_->SetNpx(2000);
  xjjroot::setthgrstyle(fun_f_, fstyle.at("f"));
  
  auto* h = (TH1F*)hmass->Clone(Form("h_%s", xjjc::unique_str().c_str()));
  sethist(h);
  auto* hMCSignal = (TH1F*)hmassMCSignal->Clone(Form("hMCSignal_%s", xjjc::unique_str().c_str()));
  sethist(hMCSignal);
  auto* hMCSwapped = (TH1F*)hmassMCSwapped->Clone(Form("hMCSwapped_%s", xjjc::unique_str().c_str()));
  sethist(hMCSwapped);

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

  fitted_ = true;
  
  fun_mass_mc_ = f_mass(Form("%s_mass_mc", fun_f_->GetName()));
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

  fun_swap_mc_ = f_swap(Form("%s_swap_mc", fun_f_->GetName()));

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

  drawleg(h);
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

void xjjroot::dfitter::reset() {
  S_ = -1;
  B_ = -1;
  yield_ = -1;
  yieldErr_ = -1;
  fitted_ = false;
  fun_f_ = nullptr;
  fun_mass_mc_ = nullptr;
  fun_swap_mc_ = nullptr;
  // if (fun_f_) delete fun_f_;
  // if (fun_mass_mc_) delete fun_mass_mc_;
  // if (fun_swap_mc_) delete fun_swap_mc_;
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
  std::string str_fun_mass = opt_3gaus_ ?
    "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*([7]*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[6]))/(sqrt(2*3.14159)*[8]*(1+[6])))))" :
    "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))))";
  auto* fun = new TF1(fname.c_str(), str_fun_mass.c_str(), fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 7 }, { 4, 9 }, { 5, 10 }, { 6, 11 }
  };
  if (opt_3gaus_) {
    params[7] = 12;
    params[8] = 13;
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
  fun->SetNpx(2000);
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
  // fun->SetNpx(2000);
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
  // fun->SetNpx(2000);
  xjjroot::setthgrstyle(fun, fstyle.at("notmass"));

  return fun;
}

TF1* xjjroot::dfitter::clonefun(const TF1* fun, const std::string& name) const {
  auto* newfun = new TF1(*fun);
  newfun->SetName(name.c_str());
  return newfun;
}

void xjjroot::dfitter::sethist(TH1* h) {
  h->SetXTitle("m_{#piK} (GeV/c^{2})");
  h->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", h->GetBinWidth(1)*1.e+3));
  xjjroot::sethempty(h, 0.045);
  xjjroot::setthgrstyle(h, kBlack, 20, 1.3, kBlack, 1, 1);
  h->SetMaximum(-1111);
  h->SetAxisRange(0, h->GetMaximum()*1.4*1.2, "Y");
}

void xjjroot::dfitter::drawleg(TH1* h) const {
  auto* leg = new TLegend(0.65, 0.58, 0.85, 0.88);
  xjjroot::setleg(leg, 0.04);
  leg->AddEntry(h, "Data", "pl");
  leg->AddEntry(fun_f_, "Fit", "l");
  xjjroot::addentrybystyle(leg, "D^{0}+#bar{D^{#lower[0.2]{0}}} Signal", "f", fstyle.at("mass"));
  xjjroot::addentrybystyle(leg, "K-#pi swapped", "f", fstyle.at("swap"));
  xjjroot::addentrybystyle(leg, "Combinatorial", "l", fstyle.at("background"));
  leg->Draw("same");
}

#endif
