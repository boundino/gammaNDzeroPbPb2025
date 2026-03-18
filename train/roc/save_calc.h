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

  std::string common_cut(const std::string& cut_a, const std::string& cut_b) {
    auto cuts_a = xjjc::str_divide_trim(cut_a, "&&"), cuts_b = xjjc::str_divide_trim(cut_b, "&&");
    std::string r;
    for (const auto& cc : cuts_a) {
      if (std::ranges::find(cuts_b, cc) != cuts_b.end())
        r += ((r.empty()?"":" && ") + cc);
    }
    return r;
  }
}

namespace calc {
  float eff_cut(TH3F* h3, int ymin = 0, int ymax = -1) {
    // if (Z nbins is not 2) ...
    ymax = (ymax < 0 ? h3->GetXaxis()->GetNbins()+1 : ymax);
    auto* h1 = (TH1F*)h3->ProjectionZ(Form("%s_pass", xjjc::str_replaceall(h3->GetName(), "h3_", "h1_").c_str()),
                                      ymin, ymax, // y
                                      0, h3->GetYaxis()->GetNbins()+1); // mass
    float eff = h1->GetBinContent(2) / h1->Integral(1, h1->GetNbinsX()); // passed
    delete h1;
    return eff;
  }
}
