#include "ratecalculator.h"
#include <vector>
#include <fstream>
#include <error.h>

using namespace csv_utils;

RateCalculator::RateCalculator()
{

}


/**
 * @brief Load a file with the average times and calculates the rates of the codons.
 * 
 * @param file_name string with the location of the file to be open.
 */
void RateCalculator::loadRates(std::string file_name)
{
    std::ifstream ist{file_name};

    if (!ist) {
        throw std::runtime_error("can't open file: "+ file_name);
    }
    codon_rates.clear();
    std::string codon;
    double decoding_time;
    std::string tmp_str;
    std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"}; // list of stop codons.
    bool header = true;
    while(ist.good()) {
        if (!std::getline ( ist, codon, ',' )){
            break;
        }

        std::getline ( ist, tmp_str, '\n' );
        decoding_time = std::atof(tmp_str.c_str());
        codon.erase(std::remove(codon.begin(), codon.end(), '\"'), codon.end()); // remove \" from string.
        if (!header){
            auto result = std::find(stop_codons.begin(), stop_codons.end(), codon);
            //only add if not a stop codon.
            if (result == end(stop_codons)) {
                codon_rates[codon] = 1/decoding_time;
            } else {
                //it is a stop codon. It has a fixed decoding rate.
                codon_rates[codon] = 1.0;
            }
        } else {
            header = false;
        }
    }
}
