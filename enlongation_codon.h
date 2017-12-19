#ifndef SIMULATIONS_ENLONGATION_CODON_H
#define SIMULATIONS_ENLONGATION_CODON_H

#include "mrnaelement.h"
#include "ribosomesimulator.h"

namespace Simulations {

class EnlongationCodon : public mRNAElement {
 public:
  EnlongationCodon();
  std::string concentrations_file_name;
  RibosomeSimulator ribosome;
  void setCodon(const std::string&);
  void loadConcentrations(const std::string&);
  void getAlphas(std::vector<double>&, std::vector<int>&) override;
  void executeReaction(int) override;
  int getState() override;
  void setState(int) override;
  void updateAlphas() override;
};
}  // namespace Simulations

#endif  // SIMULATIONS_ENLONGATION_CODON_H
