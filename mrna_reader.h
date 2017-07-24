#ifndef MRNA_UTILS_MRNA_READER_H
#define MRNA_UTILS_MRNA_READER_H

#include <string>
#include <eigen3/Eigen/Dense>
#include "ratecalculator.h"
#include "reactionsset.h"

namespace mRNA_utils {
    
    class mRNAReader
    {
    public:
        mRNAReader();
        void loadRateCalculatorFile(std::string);
        void loadmRNAFile(std::string);
        void generateInitialPopulation();
        void generateReactions();
        void setInitiationRate(double);
        void setTerminationRate(double);
        std::string mRNA_sequence;
        double termination_rate;
        double initiation_rate;
        std::string mRNA_file_name;
        csv_utils::RateCalculator rate_calculator;
        Eigen::MatrixXi initial_population;
        Simulations::ReactionsSet reactions_set;
    };
}

#endif // MRNA_UTILS_MRNA_READER_H
