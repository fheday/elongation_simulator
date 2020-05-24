#include "elongation_codon.h"
#include <fstream>

Simulations::ElongationCodon::ElongationCodon() {}

void Simulations::ElongationCodon::loadConcentrations(
    const std::string& file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    // concentrations_file_name = file_name;
    // when setting the concentrations file name, we can also
    // initialize the RibosomeSimulator object.
    ribosome.loadConcentrations(file_name);
  }
}

void Simulations::ElongationCodon::setPropensities(
    std::array<double, 40> prop) {
  ribosome.setPropensities(prop);
  updateAlphas();
}

void Simulations::ElongationCodon::setNoNonCognate(bool noNonCog) {
  if (noNonCog) {
    ribosome.setNonCognate(0.0);
  }
  updateAlphas();
}

std::map<std::string, double> Simulations::ElongationCodon::getPropensities() {
  return ribosome.getPropensities();
}

void Simulations::ElongationCodon::setCodon(const std::string& cdn) {
  ribosome.setCodonForSimulation(cdn);
  // update reactions.
  ribosome.getAlphas(alphas, reactions_index);
}

void Simulations::ElongationCodon::getAlphas(std::vector<double>& as,
                                              std::vector<int>& r_i) {
  as = alphas;
  r_i = reactions_index;
}

void Simulations::ElongationCodon::executeReaction(int r) {
  // execute reaction.
  ribosome.setState(r);
  updateAlphas();
}

int Simulations::ElongationCodon::getState() { return ribosome.getState(); }

void Simulations::ElongationCodon::setState(int s) {
  ribosome.setState(s);
  updateAlphas();
}

void Simulations::ElongationCodon::updateAlphas() {
  if (next_mRNA_element->isAvailable()) {
    ribosome.getAlphas(alphas, reactions_index);
  } else {
    ribosome.getDecodingAlphas(alphas, reactions_index);
  }
}
