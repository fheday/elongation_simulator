#include "initiationterminationcodon.h"

Simulations::InitiationTerminationCodon::InitiationTerminationCodon(float prop, bool init)
{
    propensity = prop;
    is_initiation = init; // boolean to mark if the codon  is initiation or termination.
}

void Simulations::InitiationTerminationCodon::getAlphas(Eigen::VectorXd& as, Eigen::VectorXi& reactions_index)
{
    if ((is_initiation && state == 0 && isAvailable)  || (!is_initiation && state == 0 && isOccupied)) {
        as = Eigen::VectorXd(1);
        as[0] = propensity;
        reactions_index = Eigen::VectorXi(1);
        reactions_index[0] = 0;
    } else if (state==23) {
        as = Eigen::VectorXd(1);
        as[0] = 1000; // verify
        reactions_index = Eigen::VectorXi(1);
        reactions_index[0] = 0;
    } else {
        as = Eigen::VectorXd(0);
        reactions_index = Eigen::VectorXi(0);
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
