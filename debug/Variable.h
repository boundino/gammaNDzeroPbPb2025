
struct Variable {
  std::string name;
  std::string var;
  std::string latex;
  int nbin;
  float minbin;
  float maxbin;
};

namespace globals {
  const std::vector<Variable> variables = {
    { .name = "Dmass", .var = "Dmass", .latex = "m_{K#pi} (GeV/c^{2})", .nbin = 48, .minbin = 1.66 , .maxbin = 2.16 },
    { .name = "Dchi2cl", .var = "Dchi2cl", .latex = "Secondary vertex prob", .nbin = 50, .minbin = 0. , .maxbin = 1. },
    { .name = "Dalpha", .var = "Dalpha", .latex = "3D pointing angle", .nbin = 50, .minbin = 0. , .maxbin = 1. },
    { .name = "Ddtheta", .var = "Ddtheta", .latex = "2D pointing angle", .nbin = 50, .minbin = 0. , .maxbin = 1. },
    { .name = "dls", .var = "DsvpvDistance/DsvpvDisErr", .latex = "3D decay length significance", .nbin = 40, .minbin = 0. , .maxbin = 20. },
    { .name = "Dtrk1Pt", .var = "Dtrk1Pt", .latex = "Track 1 p_{T} (GeV/c)", .nbin = 50, .minbin = 0. , .maxbin = 5. },
    { .name = "Dtrk2Pt", .var = "Dtrk2Pt", .latex = "Track 1 p_{T} (GeV/c)", .nbin = 50, .minbin = 0. , .maxbin = 5. },
  };
}


int i_variables(std::string name,
                const std::vector<Variable> &vars = globals::variables) {
  int count = 0;
  for (const auto &v : vars) {
    if (name == v.name) {
      return count;
    }
    ++count;
  }
  return -1;
}
