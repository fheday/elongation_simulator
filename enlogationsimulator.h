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
        std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> getEnlogationDuration();
        void calculateAverageTimes();
        double termination_rate = -1;
        double initiation_rate = -1;
        std::string mRNA_file_name;
        std::string concentrations_file_name;
        std::string average_times_file_name;
        mRNA_utils::mRNAReader mrna_reader;
        RibosomeSimulator ribosome_simulator;
        std::vector<std::vector<int>> ribosome_positions_history;
        // array with the total times the ribosomes spent in the codons
        std::vector<double> total_time;
        // number of times a codon was occupied
        std::vector<int> n_times_occupied;
        // average occupation time
        std::vector<double> codons_average_occupation_time;
    private:
        void intializeMRNAReader();
        std::vector<double> translocation_times;
        std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"}; // list of stop codons.
    };
}

#endif // SIMULATIONS_ENLOGATIONSIMULATOR_H
