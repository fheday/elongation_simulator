#include "mrnaelement.h"

void Simulations::mRNAElement::setAvailable(bool avail)
{
    is_available = avail;
    if (avail) {
        //update the next codon.
        if (previousMRNAElement.get() != nullptr){
            previousMRNAElement->updateAlphas();
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
    nextMRNAElement = n_c;
}

void Simulations::mRNAElement::setPreviousCodon(std::shared_ptr<mRNAElement> p_c)
{
    previousMRNAElement = p_c;
}

bool Simulations::mRNAElement::isAvailable()
{
    return is_available;
}

bool Simulations::mRNAElement::isOccupied()
{
    return is_occupied;
}
