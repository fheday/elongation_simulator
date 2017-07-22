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
}
