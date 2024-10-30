/*
 * @file  mrnaelement.cpp
 * 
 * @brief class to represent codons in a mRNA
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "mrnaelement.h"
#include <cstddef>

Simulations::mRNAElement::mRNAElement() {
  index = -1;
  alphas.resize(2);
  reactions_index.resize(2);
  previous_mRNA_element = nullptr;
  next_mRNA_element = nullptr;
}

void Simulations::mRNAElement::setAvailable(bool avail) {
  is_available = avail;
  if (avail && previous_mRNA_element != nullptr) {
    // update the next codon.
    previous_mRNA_element->updateAlphas();
  }
}
void Simulations::mRNAElement::setOccupied(bool occup) {
  is_occupied = occup;
  if (occup) {
    updateAlphas();
  }
}

void Simulations::mRNAElement::setNextCodon(mRNAElement* n_c) {
  next_mRNA_element = n_c;
}

void Simulations::mRNAElement::setPreviousCodon(mRNAElement* p_c) {
  previous_mRNA_element = p_c;
}

bool Simulations::mRNAElement::isAvailable() const { return is_available; }

bool Simulations::mRNAElement::isOccupied() const { return is_occupied; }

void Simulations::mRNAElement::addReactionToHistory(std::size_t state, float dt) {
  state_history.push_back(state);
  dt_history.push_back(dt);
}

std::pair<std::vector<std::size_t>, std::vector<float>>
Simulations::mRNAElement::getHistory() {
  return {state_history, dt_history};
}

Simulations::mRNAElement::~mRNAElement() {}
