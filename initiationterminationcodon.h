#ifndef INITIATIONTERMINATIONCODON_H
#define INITIATIONTERMINATIONCODON_H

#include "mrnaelement.h"
#include <random>
#include <math.h>
#include <eigen3/Eigen/Dense>

namespace Simulations {
    class InitiationTerminationCodon : public mRNAElement
    {
    public:
        InitiationTerminationCodon(float, bool);
        void getAlphas(std::vector<double>& as, std::vector<int>& reactions_index) override;
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
