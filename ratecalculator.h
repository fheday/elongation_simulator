#ifndef CSV_UTILS_RATECALCULATOR_H
#define CSV_UTILS_RATECALCULATOR_H

#include <string>
#include <map>
namespace csv_utils {

class RateCalculator
{
public:
RateCalculator();
void loadRates(std::string);
    std::string rates_file_name;
    std::map<std::string, double> codon_rates;
};
}

#endif // CSV_UTILS_RATECALCULATOR_H
