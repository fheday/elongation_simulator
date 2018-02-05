#include "enlongation_codon.h"
#include <fstream>

Simulations::EnlongationCodon::EnlongationCodon() {}

void Simulations::EnlongationCodon::loadConcentrations(
    const std::string& file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    concentrations_file_name = file_name;
    // when setting the concentrations file name, we can also
    // initialize the RibosomeSimulator object.
    ribosome.loadConcentrations(file_name);
  }
}

void Simulations::EnlongationCodon::setWCPropensities(
    std::array<double, 10> prop) {
  ribosome.setWCPropensities(prop);
  updateAlphas();
}

void Simulations::EnlongationCodon::setWooblePropensities(
    std::array<double, 10> prop) {
  ribosome.setWooblePropensities(prop);
  updateAlphas();
}

void Simulations::EnlongationCodon::setNearCognatePropensities(
    std::array<double, 10> prop) {
  ribosome.setNearCognatePropensities(prop);
  updateAlphas();
}

void Simulations::EnlongationCodon::setNonCogPropensities(
    std::array<double, 2> prop) {
  ribosome.setNonCogPropensities(prop);
  updateAlphas();
}

void Simulations::EnlongationCodon::setTranslocationPropensities(
    std::array<double, 9> prop) {
  ribosome.setTranslocationPropensities(prop);
  updateAlphas();
}

std::map<std::string, double> Simulations::EnlongationCodon::getPropensities() {
  return ribosome.getPropensities();
}

void Simulations::EnlongationCodon::setCodon(const std::string& cdn) {
  ribosome.setCodonForSimulation(cdn);
  // update reactions.
  ribosome.getAlphas(alphas, reactions_index);
}

void Simulations::EnlongationCodon::getAlphas(std::vector<double>& as,
                                              std::vector<int>& r_i) {
  as = alphas;
  r_i = reactions_index;
}

void Simulations::EnlongationCodon::executeReaction(int r) {
  // execute reaction.
  ribosome.setState(r);
  updateAlphas();
}

int Simulations::EnlongationCodon::getState() { return ribosome.getState(); }

void Simulations::EnlongationCodon::setState(int s) {
  ribosome.setState(s);
  updateAlphas();
}

void Simulations::EnlongationCodon::updateAlphas() {
  if (next_mRNA_element->isAvailable()) {
    ribosome.getAlphas(alphas, reactions_index);
  } else {
    ribosome.getDecodingAlphas(alphas, reactions_index);
  }
}
