#include "initiationterminationcodon.h"

Simulations::InitiationTerminationCodon::InitiationTerminationCodon(double prop,
                                                                    bool init) {
  propensity = prop;
  is_initiation =
      init;  // boolean to mark if the codon  is initiation or termination.
  alphas = std::vector<double>(1);
  alphas[0] = propensity;
  reactions_index = std::vector<int>(1);
  reactions_index[0] = 0;
  a0 = 0;
}

void Simulations::InitiationTerminationCodon::getAlphas(std::vector<double>& as,
                                                        std::vector<int>& r_i) {
  as = alphas;
  r_i = reactions_index;
}

int Simulations::InitiationTerminationCodon::getState() { return state; }

void Simulations::InitiationTerminationCodon::setState(int s) {
  state = s;
  if (state == 0) {
    alphas = std::vector<double>(1);
    alphas[0] = propensity;
    reactions_index = std::vector<int>(1);
    reactions_index[0] = 23;
  } else if (state == 23 &&
             (!is_initiation ||
              (is_initiation && next_mRNA_element.lock()->isAvailable()))) {
    alphas = std::vector<double>(1);
    alphas[0] = 10000;  // verify
    reactions_index = std::vector<int>(1);
    reactions_index[0] = 31;
  } else {
    alphas = std::vector<double>(0);
    reactions_index = std::vector<int>(0);
  }
}

void Simulations::InitiationTerminationCodon::updateAlphas() {
  setState(state);
}

void Simulations::InitiationTerminationCodon::executeReaction(int r) {
  (void)r;  // just to silence unused parameter warning.
  if (state == 0) {
    setState(23);
  } else if (state == 23) {
    setState(31);
  }
}
