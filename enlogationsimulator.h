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
        void setTerminationRate(double);
        void setInitiationRate(double);
        void setMRnaFileName(std::string);
        void setConcentrationsFileName(std::string);
        void setAverageTimesFileName(std::string);
        double getReactionTime(double&, double&, std::string&) override;
        void updateRibosomeHistory(bool=false);
        void calculateAverageTimes();
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
        std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"}; // list of stop codons.
    };
}

#endif // SIMULATIONS_ENLOGATIONSIMULATOR_H
