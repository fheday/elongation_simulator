#include "mrnaelement.h"

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

bool Simulations::mRNAElement::isAvailable() { return is_available; }

bool Simulations::mRNAElement::isOccupied() { return is_occupied; }

Simulations::mRNAElement::~mRNAElement() {}
