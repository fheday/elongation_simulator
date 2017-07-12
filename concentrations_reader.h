#ifndef CONCENTRATIONS_READER_H
#define CONCENTRATIONS_READER_H

#include <vector>
namespace csv_utils {
    struct concentration_entry{
        std::string codon;
        std::string three_letter;
        float wc_cognate_conc;
        float wobblecognate_conc;
        float nearcognate_conc;
    };
    
    class concentrations_reader
    {
    public:
        concentrations_reader(std::string);
        concentrations_reader() : concentrations_reader("../../../Projects/RSim/data_with_times/concentrations.csv") {};
        void get_contents(std::vector<concentration_entry>&);
        std::vector<concentration_entry> contents;
    };
}
#endif // CONCENTRATIONS_READER_H
