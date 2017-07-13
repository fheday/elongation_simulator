#include "concentrations_reader.h"
#include <fstream>
#include <error.h>
using namespace csv_utils;



concentrations_reader::concentrations_reader(std::string file_name)
{
    std::ifstream ist{file_name};
    
    if (!ist) {
        std::cout<<"can't open input file "<< file_name;
        throw -1;
    }
    contents.clear();
    std::string codon;
    std::string three_letter;
    float wc_cognate_conc;
    float wobblecognate_conc;
    float nearcognate_conc;
    std::string tmp_str;
    while(ist.good()) {
        if (!std::getline ( ist, codon, ',' )){
            break;
        }
        std::getline ( ist, three_letter, ',' );
        std::getline ( ist, tmp_str, ',' );
        wc_cognate_conc = std::atof(tmp_str.c_str());
        std::getline ( ist, tmp_str, ',' );
        wobblecognate_conc = std::atof(tmp_str.c_str());
        std::getline ( ist, tmp_str, '\n' );
        nearcognate_conc = std::atof(tmp_str.c_str());
        codon.erase(std::remove(codon.begin(), codon.end(), '\"'), codon.end()); // remove \" from string.
        three_letter.erase(std::remove(three_letter.begin(), three_letter.end(), '\"'), three_letter.end()); //remove \" from string.
        contents.push_back(concentration_entry{codon, three_letter, wc_cognate_conc, wobblecognate_conc, nearcognate_conc});
    }
}


void concentrations_reader::get_contents(std::vector<concentration_entry>& result)
{
    result = contents;
}
