#ifndef INITIATIONTERMINATIONCODON_H
#define INITIATIONTERMINATIONCODON_H

#include "mrnaelement.h"
#include <random>
#include <math.h>

namespace Simulations {
    class InitiationTerminationCodon : public mRNAElement
    {
    public:
        InitiationTerminationCodon(float, bool);
        void getAlphas(Eigen::VectorXd& as, Eigen::VectorXi& reactions_index) override;
        void executeReaction(int) override;
        int getState() override;
        void setState(int ) override;
        float propensity;
        float a0;
    private:
        std::mt19937 gen;
        std::uniform_real_distribution<> dis;
        int state = 0;
        bool is_initiation;
    };
}
#endif // INITIATIONTERMINATIONCODON_H
