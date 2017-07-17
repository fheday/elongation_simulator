/*
<%
cfg['dependencies'] = ['reactionsset.h', 'gillespie.h']

cfg['sources'] = ['reactionsset.cpp', 'gillespie.cpp', 'concentrations_reader.cpp']

setup_pybind11(cfg)
%>
*/
//cfg['include_dirs'] = ['/opt/anaconda/include/', '/opt/anaconda/include/eigen3']
//cfg['compiler_args'] = ['-std=c++11', '-stdlib=libc++', '-std=c++14', '-shared-libgcc', '-static-libstdc++']
//cfg['compiler_args'] = ['-std=c++14']
//cfg['parallel'] = False
//

#ifndef CMAKE_BUILD
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#endif

#include "ribosomesimulator.h"
#include <eigen3/Eigen/Dense>
#include "concentrations_reader.h"
#include <numeric>
using namespace Simulations;

#ifndef CMAKE_BUILD
PYBIND11_PLUGIN(ribosomesimulator){
    pybind11::module mod("ribosomesimulator", "auto-compiled c++ extension");

    py::class_<Gillespie> (mod, "gillespie")
    .def(py::init<>()) //constructor
    .def("setIterationLimit", &Gillespie::setIterationLimit)
    .def("run", &Gillespie::run);

    py::class_<RibosomeSimulator, Gillespie> (mod, "ribosomesimulator")
    .def(py::init<std::string&>()) //constructor
    .def("setNumberOfRibosomes", &RibosomeSimulator::setNumberOfRibosomes)
    .def("setCodonForSimulation", &RibosomeSimulator::setCodonForSimulation)
    .def("run_and_get_times", [](RibosomeSimulator &rs) {float d=0.0; float t=0.0; rs.run_and_get_times(d, t); return std::make_tuple(d, t); });
    
    return mod.ptr();

}
#endif


RibosomeSimulator::RibosomeSimulator(const std::string file_name)
{
    csv_utils::concentrations_reader cr;
    cr.load_concentrations(file_name);
    std::vector<csv_utils::concentration_entry> concentrations_vector;
    concentrations_reader = cr;
    // copied code from RibosomeSimulator(csv_utils::concentrations_reader& cr)
    //TODO: NEEDS TO IMPROVE software engineering here.
    std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"};
    std::vector<csv_utils::concentration_entry> codons_concentrations;
    cr.get_contents(codons_concentrations);
    for (csv_utils::concentration_entry entry:codons_concentrations) {
        auto result = std::find(stop_codons.begin(), stop_codons.end(), entry.codon);
        if (result == end(stop_codons)) {
            //Not a stop codon. Proceed.
            ReactionsSet rs = createReactionSet(entry);
            reactions_map[entry.codon] = rs;
        }
    }
}

void RibosomeSimulator::setNumberOfRibosomes(int nrib)
{
    Eigen::MatrixXi population(32, 1);
    population.fill(0);
    population(0,0) = nrib;
    Gillespie::setInitialPopulation(population);
}

void RibosomeSimulator::setCodonForSimulation(const std::string& codon)
{
    setReactionsSet(reactions_map.at(codon));
}

void RibosomeSimulator::run_and_get_times(float& decoding_time, float& translocation_time)
{
    Gillespie::run();
    getDecodingAndTranslocationTimes(decoding_time, translocation_time);
}

void RibosomeSimulator::getDecodingAndTranslocationTimes(float& decoding_time, float& translocation_time)
{
    int translocation_index = 0, i = 0;
    decoding_time = 0; // avoid propagating errors.
    translocation_time = 0; // avoid propagating errors.
    for (Eigen::MatrixXi population:population_history) {
        if (population(24,0) == 1) translocation_index = i;
        i++;
    }
    for (int ii = 0; ii < translocation_index; ii++) {
        decoding_time += dt_history[ii];
    }
    for (unsigned int ii = translocation_index; ii < dt_history.size(); ii++) translocation_time += dt_history[ii];
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
//    float near7f = 200;

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

    Eigen::MatrixXi reactionMatrix[44];
    // build the vector of reactions.
    // [] x=0 -> non1f:(x'=1);
    reactionMatrix[0].resize(32, 1);
    reactionMatrix[0].fill(0);
    reactionMatrix[0](0,0) = -1;
    reactionMatrix[0](1,0) = 1;
    
    // [] x=0 -> near1f:(x'=2);
    reactionMatrix[1].resize(32, 1);
    reactionMatrix[1].fill(0);
    reactionMatrix[1](0,0) = -1;
    reactionMatrix[1](2,0) = 1;
    
    // [] x=0 -> wobble1f:(x'=9);
    reactionMatrix[2].resize(32, 1);
    reactionMatrix[2].fill(0);
    reactionMatrix[2](0,0) = -1;
    reactionMatrix[2](9,0) = 1;

    // [] x=0 -> WC1f:(x'=16);
    reactionMatrix[3].resize(32, 1);
    reactionMatrix[3].fill(0);
    reactionMatrix[3](0,0) = -1;
    reactionMatrix[3](16,0) = 1;

    // [] x=1 -> non1r:(x'=0);
    reactionMatrix[4].resize(32, 1);
    reactionMatrix[4].fill(0);
    reactionMatrix[4](1,0) = -1;
    reactionMatrix[4](0,0) = 1;
    
    // [] x=2 -> near1r:(x'=0);
    reactionMatrix[5].resize(32, 1);
    reactionMatrix[5].fill(0);
    reactionMatrix[5](2,0) = -1;
    reactionMatrix[5](0,0) = 1;
    
    // [] x=2 -> near2f:(x'=3);
    reactionMatrix[6].resize(32, 1);
    reactionMatrix[6].fill(0);
    reactionMatrix[6](2,0) = -1;
    reactionMatrix[6](3,0) = 1;
    
    // [] x=3 -> near2r:(x'=2);
    reactionMatrix[7].resize(32, 1);
    reactionMatrix[7].fill(0);
    reactionMatrix[7](3,0) = -1;
    reactionMatrix[7](2,0) = 1;
    
    // [] x=3 -> near3f:(x'=4);
    reactionMatrix[8].resize(32, 1);
    reactionMatrix[8].fill(0);
    reactionMatrix[8](3,0) = -1;
    reactionMatrix[8](4,0) = 1;
    
    // [] x=4 -> near4f:(x'=5);
    reactionMatrix[9].resize(32, 1);
    reactionMatrix[9].fill(0);
    reactionMatrix[9](4,0) = -1;
    reactionMatrix[9](5,0) = 1;
    
    // [] x=5 -> near5f:(x'=6);
    reactionMatrix[10].resize(32, 1);
    reactionMatrix[10].fill(0);
    reactionMatrix[10](5,0) = -1;
    reactionMatrix[10](6,0) = 1;
    
    // [] x=6 -> neardiss:(x'=0);
    reactionMatrix[11].resize(32, 1);
    reactionMatrix[11].fill(0);
    reactionMatrix[11](6,0) = -1;
    reactionMatrix[11](0,0) = 1;
    
    // [] x=6 -> near6f:(x'=7);
    reactionMatrix[12].resize(32, 1);
    reactionMatrix[12].fill(0);
    reactionMatrix[12](6,0) = -1;
    reactionMatrix[12](7,0) = 1;
    
    // [] x=7 -> near7f:(x'=8);
    reactionMatrix[13].resize(32, 1);
    reactionMatrix[13].fill(0);
    reactionMatrix[13](7,0) = -1;
    reactionMatrix[13](8,0) = 1;
    
    // [] x=8 -> trans1f:(x'=23);
    reactionMatrix[14].resize(32, 1);
    reactionMatrix[14].fill(0);
    reactionMatrix[14](8,0) = -1;
    reactionMatrix[14](23,0) = 1;
    
    // [] x=9 -> wobble1r:(x'=0);
    reactionMatrix[15].resize(32, 1);
    reactionMatrix[15].fill(0);
    reactionMatrix[15](9,0) = -1;
    reactionMatrix[15](0,0) = 1;
    
    // [] x=9 -> wobble2f:(x'=10);
    reactionMatrix[16].resize(32, 1);
    reactionMatrix[16].fill(0);
    reactionMatrix[16](9,0) = -1;
    reactionMatrix[16](10,0) = 1;
    
    // [] x=10 -> wobble2r:(x'=9);
    reactionMatrix[17].resize(32, 1);
    reactionMatrix[17].fill(0);
    reactionMatrix[17](10,0) = -1;
    reactionMatrix[17](9,0) = 1;

    // [] x=10 -> wobble3f:(x'=11);
    reactionMatrix[18].resize(32, 1);
    reactionMatrix[18].fill(0);
    reactionMatrix[18](10,0) = -1;
    reactionMatrix[18](11,0) = 1;

    // [] x=11 -> wobble4f:(x'=12);
    reactionMatrix[19].resize(32, 1);
    reactionMatrix[19].fill(0);
    reactionMatrix[19](11,0) = -1;
    reactionMatrix[19](12,0) = 1;

    // [] x=12 -> wobble5f:(x'=13);
    reactionMatrix[20].resize(32, 1);
    reactionMatrix[20].fill(0);
    reactionMatrix[20](12,0) = -1;
    reactionMatrix[20](13,0) = 1;
    
    // [] x=13 -> wobblediss:(x'=0);
    reactionMatrix[21].resize(32, 1);
    reactionMatrix[21].fill(0);
    reactionMatrix[21](13,0) = -1;
    reactionMatrix[21](0,0) = 1;
    
    // [] x=13 -> wobble6f:(x'=14);
    reactionMatrix[22].resize(32, 1);
    reactionMatrix[22].fill(0);
    reactionMatrix[22](13,0) = -1;
    reactionMatrix[22](14,0) = 1;
    
    // [] x=14 -> wobble7f:(x'=15);
    reactionMatrix[23].resize(32, 1);
    reactionMatrix[23].fill(0);
    reactionMatrix[23](14,0) = -1;
    reactionMatrix[23](15,0) = 1;
    
    // [] x=15 -> trans1f:(x'=23);
    reactionMatrix[24].resize(32, 1);
    reactionMatrix[24].fill(0);
    reactionMatrix[24](15,0) = -1;
    reactionMatrix[24](23,0) = 1;
    
    // [] x=16 -> WC1r:(x'=0);
    reactionMatrix[25].resize(32, 1);
    reactionMatrix[25].fill(0);
    reactionMatrix[25](16,0) = -1;
    reactionMatrix[25](0,0) = 1;
    
    // [] x=16 -> WC2f:(x'=17);
    reactionMatrix[26].resize(32, 1);
    reactionMatrix[26].fill(0);
    reactionMatrix[26](16,0) = -1;
    reactionMatrix[26](17,0) = 1;
    
    // [] x=17 -> WC2r:(x'=16);
    reactionMatrix[27].resize(32, 1);
    reactionMatrix[27].fill(0);
    reactionMatrix[27](17,0) = -1;
    reactionMatrix[27](16,0) = 1;
    
    // [] x=17 -> WC3f:(x'=18);
    reactionMatrix[28].resize(32, 1);
    reactionMatrix[28].fill(0);
    reactionMatrix[28](17,0) = -1;
    reactionMatrix[28](18,0) = 1;
    
    // [] x=18 -> WC4f:(x'=19);
    reactionMatrix[29].resize(32, 1);
    reactionMatrix[29].fill(0);
    reactionMatrix[29](18,0) = -1;
    reactionMatrix[29](19,0) = 1;
    
    // [] x=19 -> WC5f:(x'=20);
    reactionMatrix[30].resize(32, 1);
    reactionMatrix[30].fill(0);
    reactionMatrix[30](19,0) = -1;
    reactionMatrix[30](20,0) = 1;
    
    // [] x=20 -> WCdiss:(x'=0);
    reactionMatrix[31].resize(32, 1);
    reactionMatrix[31].fill(0);
    reactionMatrix[31](20,0) = -1;
    reactionMatrix[31](0,0) = 1;
    
    // [] x=20 -> WC6f:(x'=21);
    reactionMatrix[32].resize(32, 1);
    reactionMatrix[32].fill(0);
    reactionMatrix[32](20,0) = -1;
    reactionMatrix[32](21,0) = 1;
    
    // [] x=21 -> WC7f:(x'=22);
    reactionMatrix[33].resize(32, 1);
    reactionMatrix[33].fill(0);
    reactionMatrix[33](21,0) = -1;
    reactionMatrix[33](22,0) = 1;
    
    // [] x=22 -> trans1f:(x'=23);
    reactionMatrix[34].resize(32, 1);
    reactionMatrix[34].fill(0);
    reactionMatrix[34](22,0) = -1;
    reactionMatrix[34](23,0) = 1;
    
    // [] x=23 -> trans1r:(x'=22);
    reactionMatrix[35].resize(32, 1);
    reactionMatrix[35].fill(0);
    reactionMatrix[35](23,0) = -1;
    reactionMatrix[35](22,0) = 1;
    
    // [] x=23 -> trans2:(x'=24);
    reactionMatrix[36].resize(32, 1);
    reactionMatrix[36].fill(0);
    reactionMatrix[36](23,0) = -1;
    reactionMatrix[36](24,0) = 1;
    
    // [] x=24 -> trans3:(x'=25);
    reactionMatrix[37].resize(32, 1);
    reactionMatrix[37].fill(0);
    reactionMatrix[37](24,0) = -1;
    reactionMatrix[37](25,0) = 1;
    
    // [] x=25 -> trans4:(x'=26);
    reactionMatrix[38].resize(32, 1);
    reactionMatrix[38].fill(0);
    reactionMatrix[38](25,0) = -1;
    reactionMatrix[38](26,0) = 1;
    
    // [] x=26 -> trans5:(x'=27);
    reactionMatrix[39].resize(32, 1);
    reactionMatrix[39].fill(0);
    reactionMatrix[39](26,0) = -1;
    reactionMatrix[39](27,0) = 1;
    
    // [] x=27 -> trans6:(x'=28);
    reactionMatrix[40].resize(32, 1);
    reactionMatrix[40].fill(0);
    reactionMatrix[40](27,0) = -1;
    reactionMatrix[40](28,0) = 1;
    
    // [] x=28 -> trans7:(x'=29);
    reactionMatrix[41].resize(32, 1);
    reactionMatrix[41].fill(0);
    reactionMatrix[41](28,0) = -1;
    reactionMatrix[41](29,0) = 1;
    
    // [] x=29 -> trans8:(x'=30);
    reactionMatrix[42].resize(32, 1);
    reactionMatrix[42].fill(0);
    reactionMatrix[42](29,0) = -1;
    reactionMatrix[42](30,0) = 1;
    
    // [] x=30 -> trans9:(x'=31);
    reactionMatrix[43].resize(32, 1);
    reactionMatrix[43].fill(0);
    reactionMatrix[43](30,0) = -1;
    reactionMatrix[43](31,0) = 1;

    ReactionsSet rs;
    int ii = 0;
    for (Eigen::MatrixXi m:reactionMatrix) {
        rs.addReaction(m, ks.at(ii++));
    }
    return rs;
}

