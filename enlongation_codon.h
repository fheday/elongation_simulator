#ifndef SIMULATIONS_ENLONGATION_CODON_H
#define SIMULATIONS_ENLONGATION_CODON_H

#include "ribosomesimulator.h"
#include "mrnaelement.h"

namespace Simulations {
    
    class EnlongationCodon: public mRNAElement
    {
    public:
        EnlongationCodon();
        ~EnlongationCodon();
        std::string concentrationsFileName;
        RibosomeSimulator ribosome;
        void setCodon(std::string);
        void loadConcentrations(std::string);
        void getAlphas(std::vector<double>&, std::vector<int>&) override;
        void executeReaction(int) override;
        int getState() override;
        void setState(int) override;
    };
}

#endif // SIMULATIONS_ENLONGATION_CODON_H
