#ifndef INITIATIONTERMINATIONCODON_H
#define INITIATIONTERMINATIONCODON_H

#include <math.h>
#include <random>
#include "mrnaelement.h"

namespace Simulations {
class InitiationTerminationCodon : public mRNAElement {
 public:
  InitiationTerminationCodon(double, bool);
  void getAlphas(std::vector<double>& as,
                 std::vector<int>& reactions_index) override;
  void executeReaction(int r) override;
  int getState() override;
  void setState(int s) override;
  void updateAlphas() override;
  double propensity;
  double a0;

 private:
  int state = 0;
  bool is_initiation;
};
}  // namespace Simulations
#endif  // INITIATIONTERMINATIONCODON_H
