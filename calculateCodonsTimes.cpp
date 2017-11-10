#include <iostream>
#include <fstream>
#include <iomanip>
#include "gillespie.h"
#include "reactionsset.h"
#include "concentrationsreader.h"
#include "enlongation_codon.h"


/**
  * @brief use the concentrations informed in concentrations_file_name, execute the number of informed iterations and then calculates the average decoding and translocating times, writing the output as a csv file into average_times_file_name.
  * This function is usually called from the function calculate_codons_propensities and used by it.
  * 
  * This procedure should usually be used only for initializing values for the EnlogationSimulator class.
  * 
  * @param concentrations_file_name string containing the path to the csv file containing the concentrations in the cell.
  * @param iterations number of iterations to run per codon base.
  * @param average_times_file_name string containing the path to write the average times calculated by the algorithm.
  * @param translocating_times boolean if true, all codons have only decoding times and the translocating time is represented by the codon 'tra'. If false, all the codons times are decoding + translocating.
  * @return std::map< std::__cxx11::string, double > a map with codons and average decoding times. Average Translocating time is given by entry 'tra'
  */
 std::map<std::string, double> calculate_codons_times(std::string concentrations_file_name, int iterations, std::string average_times_file_name, std::string times_vector_file_name, bool translocating_times)
 {
     Simulations::EnlongationCodon enlongating_ribosome;
     enlongating_ribosome.loadConcentrations(concentrations_file_name);
     csv_utils::ConcentrationsReader cr;
     cr.loadConcentrations(concentrations_file_name);
     double decoding, translocating;
     std::vector<std::string> codons;
     std::map<std::string, double> codons_times;
     cr.getCodonsVector(codons);
     
     double total_translocating = 0, total_decoding=0, n = 0;
     std::ofstream averageTimesFile;
     std::ofstream timesVectorFile;
     //set numbers precision in the files.
     averageTimesFile<<std::setprecision(15);
     timesVectorFile<<std::setprecision(15);
     //open the files for writing.
     averageTimesFile.open (average_times_file_name);
     timesVectorFile.open (times_vector_file_name);
     //create header line.
     averageTimesFile<<"codon, time\n";
     timesVectorFile<<"codon";
     for (int i = 0; i < 2 * iterations; i++) timesVectorFile<<", V"<<i;
     timesVectorFile<<"\n";
     //calculate times and generate the vectors.
     std::vector<double> vector (2 * iterations, 0);
     int codon_total_translocating = 0;
     //double codon_time = 0;
     for (std::string codon:codons){
         total_decoding=0;
         codon_total_translocating = 0;
         enlongating_ribosome.setCodon(codon);
         averageTimesFile<<"\"" <<codon<<"\"";
         std::cout<< "Starting codon: "<< codon;
         for (int i = 0 ; i < iterations; i++){
             enlongating_ribosome.ribosome.setState(0);
             enlongating_ribosome.ribosome.run_and_get_times(decoding, translocating);
             if (decoding == 0 || translocating == 0) throw std::runtime_error("decoding nor translocation cannot be zero.");
             total_decoding += decoding;
             total_translocating += translocating;
             codon_total_translocating += translocating;
             n++;
             //save vector.
             vector[i] = decoding;
             vector[iterations + i] = translocating;
         }
         //write times and vector to files.
         //codon_time = (translocating_times) ? (total_decoding)/iterations : (total_decoding + codon_total_translocating)/iterations;
         averageTimesFile<<", "<<(total_decoding)/iterations<<"\n";
         for (int j=0; j < (2 * iterations) - 1; j++) timesVectorFile<<vector[j]<<",";
         timesVectorFile<<vector[(2 * iterations) - 1]<<vector[(2 * iterations) - 1]<<"\n";
         codons_times[codon] = decoding;
         std::cout<<". Finished. Average time: "<< ((total_decoding + codon_total_translocating)/iterations) <<"\n";
     }
     //save translocation times.
     if (translocating_times){
        averageTimesFile << "tra, "<<(total_translocating/n)<<"\n";
        codons_times["tra"] = (total_translocating/n);
        std::cout<<"Average translocating time: "<< (total_translocating/n) << "\n";
     }
     // close files.
     averageTimesFile.close();
     timesVectorFile.close();
     return codons_times;
 }
 
 int main(int argc, char **argv) {
     if (argc < 6){
         std::cout<<"Wrong number of parameters informed.\n";
         std::cout<<"Usage: calculateCodonsTimes concentrations_file_name iterations average_times_file_name translocating_times\n";
         std::cout<<"concentrations_file_name  - The path to the csv file containing the concentrations in the cell.\n";
         std::cout<<"iterations - Number of iterations to run per codon base.\n";
         std::cout<<"average_times_file_name - The path to write the average times calculated by the algorithm.\n";
         std::cout<<"translocating_times - if true or 1, all codons have only decoding times and the translocating time is represented by the codon 'tra'. Otherwise, all the codons times are decoding + translocating.\n";
         return 0;
     }
     calculate_codons_times(argv[1], std::stoi(argv[2]), argv[3], argv[4], std::string(argv[5])=="true" || std::string(argv[5])=="1");

     return 0;
 }