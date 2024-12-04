#ifndef INITIATIONTERMINATIONCODON_H
#define INITIATIONTERMINATIONCODON_H

/*
 * @file  initiationterminationcodon.h
 * 
 * @brief general definition of a non-elongation codon.
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include <cmath>
#include "mrnaelement.h"

namespace Simulations {
class InitiationTerminationCodon final : public mRNAElement {
 public:
  InitiationTerminationCodon(float, bool);
  void executeReaction(int r) override;
  int getState() override;
  void setState(int s) override;
  void updateAlphas() override;

 private:
  float propensity;
  float a0;
  int state = 0;
  bool is_initiation;
};
}  // namespace Simulations
#endif  // INITIATIONTERMINATIONCODON_H
