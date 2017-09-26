#include "initiationterminationcodon.h"

Simulations::InitiationTerminationCodon::InitiationTerminationCodon(float prop, bool init)
{
    propensity = prop;
    is_initiation = init; // boolean to mark if the codon  is initiation or termination.
}

void Simulations::InitiationTerminationCodon::getAlphas(std::vector<double>& as, std::vector<int>& reactions_index)
{
    if ((is_initiation && state == 0 && isAvailable)  || (!is_initiation && state == 0 && isOccupied)) {
        as = std::vector<double>(1);
        as[0] = propensity;
        reactions_index = std::vector<int>(1);
        reactions_index[0] = 0;
    } else if (state==23) {
        as = std::vector<double>(1);
        as[0] =1000; // verify
        reactions_index = std::vector<int>(1);
        reactions_index[0] = 0;
    } else {
        as = std::vector<double>(0);
        reactions_index = std::vector<int>(0);
    }
}

int Simulations::InitiationTerminationCodon::getState()
{
    return state;
}

void Simulations::InitiationTerminationCodon::setState(int s)
{
    state = s;
}


void Simulations::InitiationTerminationCodon::executeReaction(int r)
{
    if (state == 0) {
        state = 23;
        isOccupied = true;
        isAvailable = false;
    } else if (state == 23) {
        state = 31;
    }
}
