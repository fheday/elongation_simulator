#ifndef SIMULATIONS_ELONGATION_CODON_H
#define SIMULATIONS_ELONGATION_CODON_H

/*
 * @file  elongation_codon.h
 * 
 * @brief definition of general representation of codon
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "mrnaelement.h"
#include "ribosomesimulator.h"

namespace Simulations {

class ElongationCodon : public mRNAElement {
 public:
  ElongationCodon();
  void setCodon(const std::string&);
  void loadConcentrations(const std::string&);
  void loadConcentrationsFromString(const std::string&);
  void executeReaction(int) override;
  int getState() override;
  void setState(int) override;
  void updateAlphas() override;

  void setPropensities(std::map<std::string, double> prop) override;
  void setNoNonCognate(bool noNonCog) override;
  std::map<std::string, double> getPropensities() override;

 private:
  // std::string concentrations_file_name;
  RibosomeSimulator ribosome;
};
}  // namespace Simulations

#endif  // SIMULATIONS_ELONGATION_CODON_H
