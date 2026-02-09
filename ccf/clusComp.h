#include <TF1.h>

struct clusComp {
  int nPixHits;
  std::vector<float>* z0;
  std::vector<int>* nHit;
  std::vector<float>* chi;
};

namespace globals {
  std::vector<double> clusterPars_ = { 0.0, 0.0045 };
  int nhitsTrunc_ = 150;
  double clusterTrunc_ = 2.;
  TF1* makefcut(std::string name = "f");
  std::string makestrcut(const std::string& v_nhits, const std::string& v_quality);
  double val_y(double x);
  double val_x(double y, double xmin, double xmax);
  TF1* drawline(TH2F* hempty);
}

TF1* globals::makefcut(std::string name) {
  auto* f = new TF1(name.c_str(), Form("pol%d", clusterPars_.size()-1), 0, 1500);
  for (int i=0; i<clusterPars_.size(); i++) {
    f->SetParameter(i, clusterPars_.at(i));
  }
  return f;
}

std::string globals::makestrcut(const std::string& v_nhits, const std::string& v_quality) {
  std::string result;
  for (int i=0; i<clusterPars_.size(); i++) {
    if (!result.empty()) result += "+";
    result += Form("%f", clusterPars_.at(i));
    for (int j=0; j<i; j++) {
      result += Form("*%s", v_nhits.c_str());
    }
  }
  result = result.empty() ? "" : Form(" || %s >= %s", v_quality.c_str(), result.c_str());
  result = Form("((%s < %d && %s >= 0) || %s >= %.0f%s)", v_nhits.c_str(), nhitsTrunc_, v_quality.c_str(), v_quality.c_str(), clusterTrunc_, result.c_str());
  return result;
}

double globals::val_y(double x) {
  auto* f = makefcut();
  auto result = f->Eval(x);
  delete f;
  return result;
}

double globals::val_x(double y, double xmin, double xmax) {
  auto* f = makefcut();
  double x = f->GetX(y, xmin, xmax);
  delete f;
  return x;
}

TF1* globals::drawline(TH2F* hempty) {
  xjjroot::drawline(nhitsTrunc_, 0 , nhitsTrunc_, val_y(nhitsTrunc_), kRed, 1, 2);
  xjjroot::drawline(val_x(clusterTrunc_, hempty->GetXaxis()->GetXmin(), hempty->GetXaxis()->GetXmax()), clusterTrunc_, hempty->GetXaxis()->GetXmax(), clusterTrunc_, kRed, 1, 2);
  auto* f = makefcut(xjjc::unique_str());
  f->SetRange(nhitsTrunc_, val_x(clusterTrunc_, hempty->GetXaxis()->GetXmin(), hempty->GetXaxis()->GetXmax()));
  xjjroot::setlinestyle(f, kRed, 1, 2);
  f->Draw("same");
  return f;
}
