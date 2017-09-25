#include "enlongation_codon.h"
#include <fstream>

using namespace Simulations;

EnlongationCodon::EnlongationCodon()
{
    
}

void EnlongationCodon::loadConcentrations(std::string file_name)
{
    std::ifstream ist{file_name};

    if (!ist) {
        throw std::runtime_error("can't open input file: "+ file_name);
    } else {
        concentrationsFileName = file_name;
        // when setting the concentrations file name, we can also
        // initialize the RibosomeSimulator object.
        ribosome.loadConcentrations(file_name);
        ribosome.setIterationLimit(2000);
        ribosome.setNumberOfRibosomes(1);
    }
}

void EnlongationCodon::setCodon(std::string cdn)
{
    ribosome.setCodonForSimulation(cdn);
}

void EnlongationCodon::getAlphas(std::vector<double>& as, std::vector<int>& reactions_index)
{
    ribosome.getAlphas(as, reactions_index);
}

void EnlongationCodon::executeReaction(int r)
{
    //execute reaction.
    ribosome.setState(r);
}

int Simulations::EnlongationCodon::getState()
{
    return ribosome.getState();
}

void Simulations::EnlongationCodon::setState(int s)
{
    ribosome.setState(s);
}

Simulations::EnlongationCodon::~EnlongationCodon()
{
    //we are not creating any naked pointers, therefore it is fine to have an empty destructor.
}


