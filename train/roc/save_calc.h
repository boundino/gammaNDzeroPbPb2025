namespace globals {
  const std::vector<std::pair<std::string, int>> mask_br = {
    { "DpassCut23*", 0 },
    { "DpassCut23PAS", 1 },
    { "Dtrk*P", 0 },
    { "Dtrk*MassHypo", 0 },
    { "V*", 0 },
    { "clusComp*", 0 },
  };
}

namespace save {
  void mask_branch(TChain* tr) {
    for (const auto& m : globals::mask_br)
      tr->SetBranchStatus(m.first.c_str(), m.second);
  }
}

namespace calc {
  float eff_cut(TH3F* h3) {
    // if (Z nbins is not 2) ...
    auto* h1 = (TH1F*)h3->ProjectionZ(Form("%s_pass", xjjc::str_replaceall(h3->GetName(), "h3_", "h1_").c_str()),
                                      0, h3->GetXaxis()->GetNbins()+1, // y
                                      0, h3->GetYaxis()->GetNbins()+1); // mass
    float eff = h1->GetBinContent(2) / h1->Integral(1, h1->GetNbinsX()); // passed
    delete h1;
    return eff;
  }
}
