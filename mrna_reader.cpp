#include "mrna_reader.h"

using namespace mRNA_utils;

mRNAReader::mRNAReader()
{

}
void mRNAReader::loadRateCalculatorFile(std::string file_name)
{
    rate_calculator.loadRates(file_name);
}

void mRNAReader::loadmRNAFile(std::string mRNA_file_name)
{
    std::ifstream ist{mRNA_file_name};

    if (!ist) {
        throw std::runtime_error("can't open input file: "+ mRNA_file_name);
    }
    mRNA_sequence.clear();
    std::string line;
    while(ist.good()) {
        std::getline ( ist, line );
        // some file formats start with a '>' symbol on the first line.
        // we need to skip that line.
        if (line[0]=='>') continue;
        mRNA_sequence.append(line);
    }
    // replace all T's for U's.
    std::size_t found = mRNA_sequence.find('T');
    while (found != std::string::npos) {
        mRNA_sequence[found] = 'U';
        found = mRNA_sequence.find('T', found+1);
    }
}

void mRNAReader::setInitiationRate(double val)
{
    if (val > 0) {
        initiation_rate = val;
    } else {
        throw std::runtime_error("invalid initiation rate.");
    }
}

void mRNAReader::setTerminationRate(double val)
{
    if (val > 0) {
        termination_rate = val;
    } else {
        throw std::runtime_error("invalid termination rate.");
    }
}

void mRNAReader::generateInitialPopulation()
{
    int n_codons = mRNA_sequence.size()/3;
    initial_population = Eigen::MatrixXi(4, n_codons);
    initial_population.fill(0);
    initial_population.row(0).setOnes();
    initial_population.row(1).setOnes();
}

void mRNAReader::generateReactions()
{
    int n_codons = mRNA_sequence.size()/3;
    std::string codon;
    Eigen::MatrixXi matrix;
    
    // initiation
    matrix = Eigen::MatrixXi(4, n_codons);
    matrix.fill(0);
    matrix.col(0) << 0, -1, 1, 0;
    std::string start_codon = mRNA_sequence.substr(0, 3);
    reactions_set.addReaction(matrix, initiation_rate, start_codon);

    for(int i =0; i< n_codons; i++) {
        codon = mRNA_sequence.substr(i * 3, 3);
        // Decoding
        matrix = Eigen::MatrixXi(4, n_codons);
        matrix.fill(0);
        matrix.col(i) << -1, 0, -1, 1;
        reactions_set.addReaction(matrix, rate_calculator.codon_rates[codon], codon);
        if (i < n_codons - 1){
            // Translocating reactions  - to be optimised
            matrix = Eigen::MatrixXi(4, n_codons);
            matrix.fill(0);
            matrix(3, i) = -1;
            matrix.col(i+1) << 0, -1, +1, 0;
            if (i > 9 && i <= n_codons - 2){
                // we need to replenish the mRNA slots
                matrix(0, i - 10) = 1;
                matrix(1, i - 10) = 1;
            }
            reactions_set.addReaction(matrix, rate_calculator.codon_rates["tra"], "tra");
        }
    }
    // termination
    std::string termination_codon = mRNA_sequence.substr((n_codons-1)*3, (n_codons-1)*3 +2);
    matrix = Eigen::MatrixXi(4, n_codons);
    matrix.fill(0);
    matrix.col(n_codons - 1) << 0, 0, 0, -1;
    // we need to replenish the mRNA slots
    for (int j = n_codons - 10; j < n_codons; j++) {
        matrix(0, j) = 1;
        matrix(1, j) = 1;
    }
    reactions_set.addReaction(matrix, termination_rate, termination_codon);
}
