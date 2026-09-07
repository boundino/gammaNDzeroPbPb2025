#include "xjjcuti.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace draw {
  class bintex {
  public:
    explicit bintex(TH1* h = nullptr, int xyz_y = 0,  int xyz_pt = 0)
      : h_bins(h), axis_pt(decide_axis(xyz_pt)), axis_y(decide_axis(xyz_y)) { }

    void seth(TH1* h = nullptr, int xyz_y = 0,  int xyz_pt = 0) {
      h_bins = h;
      axis_pt = decide_axis(xyz_pt);
      axis_y = decide_axis(xyz_y);
    }
    bool valid() { return h_bins != nullptr; }
    
    std::string label_y(int i = -1) const {
      const auto ymin = (i >= 0 ? axis_y->GetBinLowEdge(i+1) : axis_y->GetBinLowEdge(1));
      const auto ymax = (i >= 0 ? axis_y->GetBinUpEdge(i+1) : axis_y->GetBinUpEdge(axis_y->GetNbins()));
      return xjjc::str_replaceall(xjjc::number_range_string(ymin, ymax, "y", -1.e1), " #", "#scale[0.4]{ }#");
    }
    int ny() const { return axis_y->GetNbins(); }
    double binwidth_y(int i = -1) const {
      const auto ymin = (i >= 0 ? axis_y->GetBinLowEdge(i+1) : axis_y->GetBinLowEdge(1));
      const auto ymax = (i >= 0 ? axis_y->GetBinUpEdge(i+1) : axis_y->GetBinUpEdge(axis_y->GetNbins()));
      return ymax-ymin;
    }
    std::string label_pt(int i = -1) const {
      const auto ptmin = (i >= 0 ? axis_pt->GetBinLowEdge(i+1) : axis_pt->GetBinLowEdge(1));
      const auto ptmax = (i >= 0 ? axis_pt->GetBinUpEdge(i+1) : axis_pt->GetBinUpEdge(axis_pt->GetNbins()));
      return xjjc::str_replaceall(xjjc::number_range_string(ptmin, ptmax, "#it{p}_{T}"), " #", "#scale[0.4]{ }#") + " GeV";
    }
    double binwidth_pt(int i = -1) const {
      const auto ptmin = (i >= 0 ? axis_pt->GetBinLowEdge(i+1) : axis_pt->GetBinLowEdge(1));
      const auto ptmax = (i >= 0 ? axis_pt->GetBinUpEdge(i+1) : axis_pt->GetBinUpEdge(axis_pt->GetNbins()));
      return ptmax-ptmin;
    }
    int npt() const { return axis_pt->GetNbins(); }
    
    template<class T> T* make_h1_y(std::string name) const {
      auto* h = new T(name.c_str(), Form(";%s;", axis_y->GetTitle()), axis_y->GetNbins(), axis_y->GetXbins()->GetArray());
      return h;
    }
    
  private:
    TH1* h_bins;
    TAxis *axis_pt, *axis_y;
    TAxis* decide_axis(int xyz) {
      TAxis* axis = nullptr;
      if (!h_bins) return axis;
      if (xyz == 0) {
        axis = h_bins->GetXaxis();
      }
      else if (xyz == 1) {
        if (!dynamic_cast<TH2*>(h_bins))
          throw std::runtime_error("Y axis requested but histogram is not TH2/TH3");
        axis = h_bins->GetYaxis();
      }
      else if (xyz == 2) {
        if (!dynamic_cast<TH3*>(h_bins))
          throw std::runtime_error("Z axis requested but histogram is not TH3");
        axis = h_bins->GetZaxis();
      }
      else {
        throw std::runtime_error("xyz must be 0, 1, or 2");
      }
      return axis;
    }
  };

}
