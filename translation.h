#ifndef TRANSLATION_H
#define TRANSLATION_H

#include <eigen3/Eigen/Dense>
#include <vector>
#include "enlongation_codon.h"

namespace Simulations{
    
    class Translation
    {
    public:
        void loadMRNA(std::string);
        void loadConcentrations(std::string);
        
        void setInitiationRate(double);
        void setTerminationRate(double);

        void setIterationLimit(int);
        void setTimeLimit(double);

        void getAlphas();
        void run();
        
        void calculateAverageTimes();
        std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> getEnlogationDuration();
        
        
        double termination_rate = -1;
        double initiation_rate = -1;
        int iteration_limit = -1;
        double time_limit = -1;
        

        std::vector<double>  alphas; // reactions alphas - all available ones.
        std::vector<int>  codon_index; // indexes of the codon where the alpha belongs to.
        std::vector<int> reaction_index; // in the codon, the index of the reaction.
        std::vector<std::unique_ptr<Simulations::mRNAElement>> codons_vector;
        std::string mrna_file_name;
        std::string concentrations_file_name;
        
        std::vector<double> dt_history;
        std::vector<std::vector<int>> ribosome_positions_history;
        
        // array with the total times the ribosomes spent in the codons
        std::vector<double> total_time;
        // number of times a codon was occupied
        std::vector<int> n_times_occupied;
        // average occupation time
        std::vector<double> codons_average_occupation_time;
    private:
        void initializeMRNAReader();
    };
}
#endif // TRANSLATION_H
