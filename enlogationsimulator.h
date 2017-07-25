#ifndef SIMULATIONS_ENLOGATIONSIMULATOR_H
#define SIMULATIONS_ENLOGATIONSIMULATOR_H

#include "gillespie.h"
#include "enlogationsimulator.h"
#include "ribosomesimulator.h"
#include "mrna_reader.h"

namespace Simulations {
    
    class EnlogationSimulator : public Simulations::Gillespie
    {
    public:
        EnlogationSimulator();
        void set_termination_rate(double);
        void set_initiation_rate(double);
        void set_mRna_file_name(std::string);
        void set_concentrations_file_name(std::string);
        void setAverageTimesFileName(std::string);
        double getReactionTime(double, double, std::string) override;
        void updateRibosomeHistory();
        double termination_rate = -1;
        double initiation_rate = -1;
        std::string mRNAFileName;
        std::string concentrationsFileName;
        std::string average_times_file_name;
        mRNA_utils::mRNAReader mrna_reader;
        RibosomeSimulator ribosome_simulator;
        std::vector<std::vector<int>> ribosome_positions_history;
    private:
        void intializeMRNAReader();
        std::vector<double> translocation_times;
    };
}

#endif // SIMULATIONS_ENLOGATIONSIMULATOR_H
