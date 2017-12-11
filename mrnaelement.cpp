#include "mrnaelement.h"

void Simulations::mRNAElement::setAvailable(bool avail)
{
    is_available = avail;
    if (avail) {
        //update the next codon.
        if (auto tmp = previous_mRNA_element.lock()){
            tmp->updateAlphas();
        }
    }
}
void Simulations::mRNAElement::setOccupied(bool occup)
{
    is_occupied = occup;
    if (occup) updateAlphas();
}

void Simulations::mRNAElement::setNextCodon(std::shared_ptr<mRNAElement> n_c)
{
    next_mRNA_element = n_c;
}

void Simulations::mRNAElement::setPreviousCodon(std::shared_ptr<mRNAElement> p_c)
{
    previous_mRNA_element = p_c;
}

bool Simulations::mRNAElement::isAvailable()
{
    return is_available;
}

bool Simulations::mRNAElement::isOccupied()
{
    return is_occupied;
}
