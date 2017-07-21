#ifndef CONCENTRATIONS_READER_H
#define CONCENTRATIONS_READER_H

#include <vector>
namespace csv_utils {
    struct concentration_entry{
        std::string codon;
        std::string three_letter;
        double wc_cognate_conc;
        double wobblecognate_conc;
        double nearcognate_conc;
    };
    
    class concentrations_reader
    {
    public:
        concentrations_reader();
        void load_concentrations(std::string);
        void get_contents(std::vector<concentration_entry>&);
        void get_codons_vector(std::vector<std::string>&);
        std::vector<concentration_entry> contents;
    };
}
#endif // CONCENTRATIONS_READER_H
