#ifndef _TMVACLASSIFICATION_H_
#define _TMVACLASSIFICATION_H_

#include <string>
#include <vector>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>

#include "xjjcuti.h"
#include "dtree.h"

namespace mytmva
{
  std::vector<float> ptbins({2., 5.});
  int nptbins = ptbins.size()-1;
  struct tmvavar
  {
    std::string varname;
    std::string vartex;
    std::string var;
    std::string cutsign;
    float varmin;
    float varmax;
    tmvavar(const std::string& varname_, const std::string& var_, const std::string& cutsign_, const std::string& vartex_, const float& varmin_, const float& varmax_) 
      : varname(varname_), var(var_), cutsign(cutsign_), vartex(vartex_), varmin(varmin_), varmax(varmax_) { ; }
  };

  const std::vector<mytmva::tmvavar> varlist = {
    /*0: */ mytmva::tmvavar("Dchi2cl"    , "Dchi2cl"                                                                                        , "FMax", "vertex #chi^{2} prob"                   , 0   , 1)  ,
    /*1: */ mytmva::tmvavar("Dalpha"     , "Dalpha"                                                                                         , "FMin", "#alpha"                                 , 0   , 3.2),
    /*2: */ mytmva::tmvavar("costheta"   , "costheta := TMath::Cos(Ddtheta)"                                                                , "FMax", "cos(#theta)"                            , -1  , 1)  ,
    /*3: */ mytmva::tmvavar("dls3D"      , "dls3D := TMath::Abs(DsvpvDistance/DsvpvDisErr)"                                                 , "FMax", "l_{xyz}/#sigma(l_{xyz})"                , 0   , 10) ,
    /*4: */ mytmva::tmvavar("dls2D"      , "dls2D := TMath::Abs(DsvpvDistance_2D/DsvpvDisErr_2D)"                                           , "FMax", "l_{xy}/#sigma(l_{xy})"                  , 0   , 10) ,
    /*5: */ mytmva::tmvavar("Dtrk1Pt"    , "Dtrk1Pt"                                                                                        , "FMax", "#pi_{1} p_{T} (GeV/c)"                  , 0   , 10) ,
    /*6: */ mytmva::tmvavar("Dtrk2Pt"    , "Dtrk2Pt"                                                                                        , "FMax", "#pi_{2} p_{T} (GeV/c)"                  , 0   , 10) ,
    /*7: */ mytmva::tmvavar("trkptimba"  , "trkptimba := TMath::Abs((Dtrk1Pt-Dtrk2Pt) / (Dtrk1Pt+Dtrk2Pt))"                                 , "FMax", "|p_{T}(1)-p_{T}(2)| / p_{T}(1)+p_{T}(2)", 0   , 1)  ,
    /*8: */ mytmva::tmvavar("Dy"         , "Dy"                                                                                             , ""    , "y"                                      , -2.4, 2.4),
    /*9: */ mytmva::tmvavar("Dmass"      , "Dmass"                                                                                          , ""    , "m_{K#pi}"                               , 1.7 , 2.0),
    /*10:*/ mytmva::tmvavar("Dtrk1Eta"   , "Dtrk1Eta"                                                                                       , ""    , "#pi_{1} #eta"                           , -2  , 2)  ,
    /*11:*/ mytmva::tmvavar("Dtrk2Eta"   , "Dtrk2Eta"                                                                                       , ""    , "#pi_{2} #eta"                           , -2  , 2)  ,
    /*12:*/ mytmva::tmvavar("Dtrk1DxySig", "Dtrk1DxySig := TMath::Abs(Dtrk1Dxy1/Dtrk1DxyError1)"                                            , ""    , "#pi_{1} |D_{xy}/#sigma(D_{xy})|"        , 0   , 4)  ,
    /*13:*/ mytmva::tmvavar("Dtrk2DxySig", "Dtrk2DxySig := TMath::Abs(Dtrk2Dxy1/Dtrk2DxyError1)"                                            , ""    , "#pi_{2} |D_{xy}/#sigma(D_{xy})|"        , 0   , 4)  ,
  };
  const mytmva::tmvavar* findvar(std::string varlabel);

  class varval
  {
  public:
    varval(TTree* nttree) : fnt(nullptr), fvalid(true), rr(new TRandom()) { fnt = new hfupc::dtree(nttree); fval.clear(); fvalid = checkvarlist(); }
    float getval(std::string varname, int j) { refreshval(j); if (fval.find(varname) == fval.end()) { std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<varname<<std::endl; return 0; } ; return fval[varname]; }
    hfupc::dtree* getnt() { return fnt; }
    bool isvalid() { return fvalid; }

  private:
    bool fvalid;
    std::map<std::string, float> fval;
    hfupc::dtree* fnt; //~
    TRandom* rr;

    void refreshval(int j) {
      fval["Dmass"]       = j<0?0:fnt->val("Dmass", j);
      fval["costheta"]    = j<0?0:TMath::Cos(fnt->val("Ddtheta", j));
      fval["dls3D"]       = j<0?0:TMath::Abs(fnt->val("DsvpvDistance", j)/fnt->val("DsvpvDisErr", j));
      fval["Dchi2cl"]     = j<0?0:fnt->val("Dchi2cl", j);
      fval["Dalpha"]      = j<0?0:fnt->val("Dalpha", j);
      fval["dls2D"]       = j<0?0:TMath::Abs(fnt->val("DsvpvDistance_2D", j)/fnt->val("DsvpvDisErr_2D", j));
      fval["trkptimba"]   = j<0?0:TMath::Abs((fnt->val("Dtrk1Pt", j)-fnt->val("Dtrk2Pt", j)) / (fnt->val("Dtrk1Pt", j)+fnt->val("Dtrk2Pt", j)));
      fval["Dy"]          = j<0?0:fnt->val("Dy", j);
      fval["Dtrk1Pt"]     = j<0?0:fnt->val("Dtrk1Pt", j);
      fval["Dtrk2Pt"]     = j<0?0:fnt->val("Dtrk2Pt", j);
      fval["Dtrk1Eta"]    = j<0?0:fnt->val("Dtrk1Eta", j);
      fval["Dtrk2Eta"]    = j<0?0:fnt->val("Dtrk2Eta", j);
      fval["Dtrk1DxySig"] = j<0?0:TMath::Abs(fnt->val("Dtrk1Dxy1", j)/fnt->val("Dtrk1DxyError1", j));
      fval["Dtrk2DxySig"] = j<0?0:TMath::Abs(fnt->val("Dtrk2Dxy1", j)/fnt->val("Dtrk2DxyError1", j));
    }
    bool checkvarlist() {
      refreshval(-1);
      for (auto& vn : varlist) {
        if (fval.find(vn.varname) == fval.end()) {
          std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<vn.varname<<std::endl;
          return false;
        }
      }
      return true;
    }
  };

  std::vector<std::string> argmethods; std::vector<int> argstages;
  std::string mkname(std::string outputname, std::string mymethod, std::string stage, float ptmin, float ptmax);
  std::string mkname(std::string outputname, std::string mymethod, std::string stage);
  std::string mkname(float ptmin, float ptmax);
}

const mytmva::tmvavar* mytmva::findvar(std::string varlabel) {
  for (auto& vv : varlist) {
    if (vv.varname == varlabel)
      return &vv;
  }
  return 0;
}

std::string mytmva::mkname(std::string outputname, std::string mymethod, std::string stage) {
  mymethod = xjjc::str_replaceall(mymethod, " ", "");
  stage = xjjc::str_replaceall(stage, " ", "");
  auto outfname = outputname
    + "_" + xjjc::str_replaceall(mymethod, ",", "-")
    + "_" + xjjc::str_replaceall(stage, ",", "-");
  return outfname;
}

std::string mytmva::mkname(float ptmin, float ptmax) {
  auto outfname = "pt-" + xjjc::number_to_string(ptmin) + "-" + (ptmax<0?"inf":xjjc::number_to_string(ptmax));
  return outfname;
}

std::string mytmva::mkname(std::string outputname, std::string mymethod, std::string stage, float ptmin, float ptmax) {
  mymethod = xjjc::str_replaceall(mymethod, " ", "");
  stage = xjjc::str_replaceall(stage, " ", "");
  auto outfname = mkname(outputname, mymethod, stage)
    + "_" + mkname(ptmin, ptmax)
    + ".root";
  // std::string outfname(Form("%s_%s_%s-%s_%s.root", outputname.c_str(), xjjc::str_replaceall(mymethod, ",", "-").c_str(),
  //                           xjjc::number_to_string(ptmin).c_str(), (ptmax<0?"inf":xjjc::number_to_string(ptmax).c_str()),
  //                           xjjc::str_replaceall(stage, ",", "-").c_str()));
  return outfname;
}


#endif
