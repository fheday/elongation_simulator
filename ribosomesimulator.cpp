#include "ribosomesimulator.h"
#include <eigen3/Eigen/Dense>
#include "concentrations_reader.h"
using namespace Simulations;

RibosomeSimulator::RibosomeSimulator(csv_utils::concentrations_reader& cr)
{
    concentrations_reader = cr;
    std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"};
    std::vector<csv_utils::concentration_entry> codons_concentrations;
    cr.get_contents(codons_concentrations);
    for (csv_utils::concentration_entry entry:codons_concentrations) {
        auto result = std::find(stop_codons.begin(), stop_codons.end(), entry.codon);
        if (result == end(stop_codons)) {
            //Not a stop codon. Proceed.
            ReactionsSet rs = createReactionSet(entry);
            reactions_map[entry.codon]=(rs);
        }
    }
    //print all reactions.
    for (std::pair<std::string, ReactionsSet> element:reactions_map){
        std::cout<<"Starting new reactions set: "<< element.first <<"\n ++++++++++++\n";
        for (Eigen::MatrixXi m: element.second.reactions_vector) {
            std::cout<<"Reaction Matrix = \n"<< m<<"\n---------\n";
        }
        std::cout<<"Finished reactions set\n ++++++++++++\n";
    }
}

ReactionsSet RibosomeSimulator::createReactionSet(const csv_utils::concentration_entry& codon)
{
    float totalconc = 1.9e-4;
    float nonconc = totalconc - codon.wc_cognate_conc - codon.wobblecognate_conc - codon.nearcognate_conc;
    // based on yeast value of 226000 molecules per cell as determined
    // in von der Haar 2008 (PMID 18925958)
    float eEF2conc = 1.36e-5;
    // constants for WCcognate interaction in 1/sec
    float WC1f = 1.4e8*codon.wc_cognate_conc;
    float WC1r = 85;
    float WC2f = 190;
    float WC2r = 0.23;
    float WC3f = 260;
    float WC4f = 1000;
    float WC5f = 1000;
    float WCdiss = 60;
    float WC6f = 1000;
    float WC7f = 200;

    // constants for wobblecognate interaction in 1/sec
    float wobble1f = 1.4e8*codon.wobblecognate_conc;
    float wobble1r = 85;
    float wobble2f = 190;
    float wobble2r = 1;
    float wobble3f = 25;
    float wobble4f = 1000;
    float wobble5f = 1000;
    float wobblediss = 1.1;
    float wobble6f = 1.6;
    float wobble7f = 200;

    // constants for nearcognate interaction in 1/sec
    float near1f = 1.4e8*codon.nearcognate_conc;
    float near1r = 85;
    float near2f = 190;
    float near2r = 80;
    float near3f = 0.4;
    float near4f = 1000;
    float near5f = 1000;
    float neardiss = 1000;
    float near6f = 60;
    float near7f = 200;

    // constants for noncognate interaction in 1/sec.
    // Non-cognates are assumed to not undergo any significant
    // interaction but to simply dissociate quickly.
    float non1f = 1.4e8*nonconc;
    float non1r = 1e5;

    // constants for translocation in 1/sec
    // 150 uM-1 s-1 = is from Fluitt et al 2007 (PMID 17897886)
    float trans1f = eEF2conc*1.5e8;
    float trans1r = 140;
    float trans2 = 250;
    float trans3 = 350;
    float trans4 = 1000;
    float trans5 = 1000;
    float trans6 = 1000;
    float trans7 = 1000;
    float trans8 = 1000;
    float trans9 = 1000;
    
    std::vector<float> ks = {non1f, near1f, wobble1f, WC1f, non1r, near1r, near2f, near2r, near3f, near4f, near5f, neardiss, near6f, near4f, trans1f, wobble1r, wobble2f, wobble2r, wobble3f, wobble4f, wobble5f, wobblediss, wobble6f, wobble7f, trans1f, WC1r, WC2f, WC2r, WC3f, WC4f, WC5f, WCdiss, WC6f, WC7f, trans1f, trans1r, trans2, trans3, trans4, trans5, trans6, trans7, trans8, trans9};

    std::vector<Eigen::MatrixXi> vector_of_reactios;
    // build the vector of reactions.
    // [] x=0 -> non1f:(x'=1);
    Eigen::MatrixXi* reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(0,0) = -1;
    (*reactionMatrix)(1,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=0 -> near1f:(x'=2);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(0,0) = -1;
    (*reactionMatrix)(2,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=0 -> wobble1f:(x'=9);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(0,0) = -1;
    (*reactionMatrix)(9,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);

    // [] x=0 -> WC1f:(x'=16);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(0,0) = -1;
    (*reactionMatrix)(16,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);

    // [] x=1 -> non1r:(x'=0);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(1,0) = -1;
    (*reactionMatrix)(0,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=2 -> near1r:(x'=0);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(2,0) = -1;
    (*reactionMatrix)(0,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=2 -> near2f:(x'=3);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(2,0) = -1;
    (*reactionMatrix)(3,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=3 -> near2r:(x'=2);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(3,0) = -1;
    (*reactionMatrix)(2,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=3 -> near3f:(x'=4);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(3,0) = -1;
    (*reactionMatrix)(4,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=4 -> near4f:(x'=5);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(4,0) = -1;
    (*reactionMatrix)(5,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=5 -> near5f:(x'=6);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(5,0) = -1;
    (*reactionMatrix)(6,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=6 -> neardiss:(x'=0);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(6,0) = -1;
    (*reactionMatrix)(0,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=6 -> near6f:(x'=7);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(6,0) = -1;
    (*reactionMatrix)(7,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=7 -> near7f:(x'=8);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(7,0) = -1;
    (*reactionMatrix)(8,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=8 -> trans1f:(x'=23);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(8,0) = -1;
    (*reactionMatrix)(23,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=9 -> wobble1r:(x'=0);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(9,0) = -1;
    (*reactionMatrix)(0,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=9 -> wobble2f:(x'=10);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(9,0) = -1;
    (*reactionMatrix)(10,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=10 -> wobble2r:(x'=9);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(10,0) = -1;
    (*reactionMatrix)(9,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);

    // [] x=10 -> wobble3f:(x'=11);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(10,0) = -1;
    (*reactionMatrix)(11,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);

    // [] x=11 -> wobble4f:(x'=12);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(11,0) = -1;
    (*reactionMatrix)(12,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);

    // [] x=12 -> wobble5f:(x'=13);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(12,0) = -1;
    (*reactionMatrix)(13,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=13 -> wobblediss:(x'=0);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(13,0) = -1;
    (*reactionMatrix)(0,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=13 -> wobble6f:(x'=14);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(13,0) = -1;
    (*reactionMatrix)(14,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=14 -> wobble7f:(x'=15);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(14,0) = -1;
    (*reactionMatrix)(15,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=15 -> trans1f:(x'=23);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(15,0) = -1;
    (*reactionMatrix)(23,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=16 -> WC1r:(x'=0);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(16,0) = -1;
    (*reactionMatrix)(0,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=16 -> WC2f:(x'=17);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(16,0) = -1;
    (*reactionMatrix)(17,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=17 -> WC2r:(x'=16);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(17,0) = -1;
    (*reactionMatrix)(16,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=17 -> WC3f:(x'=18);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(17,0) = -1;
    (*reactionMatrix)(18,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=18 -> WC4f:(x'=19);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(18,0) = -1;
    (*reactionMatrix)(19,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=19 -> WC5f:(x'=20);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(19,0) = -1;
    (*reactionMatrix)(20,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=20 -> WCdiss:(x'=0);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(20,0) = -1;
    (*reactionMatrix)(0,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=20 -> WC6f:(x'=21);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(20,0) = -1;
    (*reactionMatrix)(21,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=21 -> WC7f:(x'=22);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(21,0) = -1;
    (*reactionMatrix)(22,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=22 -> trans1f:(x'=23);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(22,0) = -1;
    (*reactionMatrix)(23,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=23 -> trans1r:(x'=22);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(23,0) = -1;
    (*reactionMatrix)(22,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=23 -> trans2:(x'=24);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(23,0) = -1;
    (*reactionMatrix)(24,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=24 -> trans3:(x'=25);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(24,0) = -1;
    (*reactionMatrix)(25,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=25 -> trans4:(x'=26);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(25,0) = -1;
    (*reactionMatrix)(26,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=26 -> trans5:(x'=27);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(26,0) = -1;
    (*reactionMatrix)(27,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=27 -> trans6:(x'=28);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(27,0) = -1;
    (*reactionMatrix)(28,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=28 -> trans7:(x'=29);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(28,0) = -1;
    (*reactionMatrix)(29,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=29 -> trans8:(x'=30);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(29,0) = -1;
    (*reactionMatrix)(30,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);
    
    // [] x=30 -> trans9:(x'=31);
    reactionMatrix = new Eigen::MatrixXi(32, 1);
    reactionMatrix->fill(0);
    (*reactionMatrix)(30,0) = -1;
    (*reactionMatrix)(31,0) = 1;
    vector_of_reactios.push_back(*reactionMatrix);

    ReactionsSet rs;
    int ii = 0;
    for (Eigen::MatrixXi m:vector_of_reactios) {
        rs.addReaction(m, ks.at(ii++));
    }
    return rs;
}

