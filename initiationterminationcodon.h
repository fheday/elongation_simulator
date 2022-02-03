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

#include <math.h>
#include <random>
#include "mrnaelement.h"

namespace Simulations {
class InitiationTerminationCodon : public mRNAElement {
 public:
  InitiationTerminationCodon(double, bool);
  void executeReaction(int r) override;
  int getState() override;
  void setState(int s) override;
  void updateAlphas() override;

 private:
  double propensity;
  double a0;
  int state = 0;
  bool is_initiation;
};
}  // namespace Simulations
#endif  // INITIATIONTERMINATIONCODON_H
