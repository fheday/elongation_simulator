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
    }
}

void EnlongationCodon::setCodon(std::string cdn)
{
    ribosome.setCodonForSimulation(cdn);
    //update reactions.
    ribosome.getAlphas(alphas, reactions_index);
}

void EnlongationCodon::getAlphas(std::vector<double>& as, std::vector<int>& r_i)
{
    as = alphas;
    r_i = reactions_index;
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


void Simulations::EnlongationCodon::updateAlphas(bool is_next_codon_available)
{
    if (is_next_codon_available){
        ribosome.getAlphas(alphas, reactions_index);
    } else {
        ribosome.getDecodingAlphas(alphas, reactions_index);
    }
}



Simulations::EnlongationCodon::~EnlongationCodon()
{
    //we are not creating any naked pointers, therefore it is fine to have an empty destructor.
}


