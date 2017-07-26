 #include <iostream>
 #include <fstream>
 #include <iomanip>
 #include <fstream>
 #include "gillespie.cpp"
 #include "reactionsset.cpp"
 #include "concentrationsreader.cpp"
 #include <vector>
 #include "ribosomesimulator.cpp"
 #include "enlogationsimulator.cpp"
 #include "ratecalculator.cpp"
 #include "mrna_reader.cpp"
 
 
 
 void testGillespie()
 {
     // create one reaction matrix
     Eigen::MatrixXi reactions(4,10);
     reactions.fill(0);
     reactions(1,1) = -1;
     reactions(2,1) = 1;
     // create a ReactionSet object and add the reaction there.
     Simulations::ReactionsSet reactions_set;
     // the reaction propensity is the second argument.
     reactions_set.addReaction(reactions, .1);
     // create a matrix with the initial species population.
     Eigen::MatrixXi population(4,10);
     population.fill(0);
     population(1,1) = 800;
     std::cout<<"initial population = "<<population<<"\n------\n";
     // create the Gillespie simulator object with the population and reactions
     Simulations::Gillespie simulation(300, population, reactions_set);
     //run the simulation.
     simulation.run();
 }
 
 /**
  * @brief use the concentrations informed in concentrations_file_name, execute the number of informed iterations and then calculates the average decoding and translocating times, writing the output as a csv file into average_times_file_name.
  * This function is usually called from the function calculate_codons_propensities and used by it.
  * 
  * This procedure should usually be used only for initializing values for the EnlogationSimulator class.
  * 
  * @param concentrations_file_name string containing the path to the csv file containing the concentrations in the cell.
  * @param iterations number of iterations to run per codon base.
  * @param average_times_file_name string containing the path to write the average times calculated by the algorithm.
  * @return std::map< std::__cxx11::string, double > a map with codons and average decoding times. Average Translocating time is given by entry 'tra'
  */
 std::map<std::string, double> calculate_codons_times(std::string concentrations_file_name, int iterations, std::string average_times_file_name, std::string times_vector_file_name)
 {
     csv_utils::ConcentrationsReader cr;
     cr.loadConcentrations(concentrations_file_name);
     Simulations::RibosomeSimulator rs;
     rs.loadConcentrations(concentrations_file_name);
     rs.setIterationLimit(2000);
     // create a matrix with the initial species population.
     rs.setNumberOfRibosomes(1);
     
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
     for (std::string codon:codons){
         total_decoding=0;
         rs.setCodonForSimulation(codon);
         averageTimesFile<<"\"" <<codon<<"\"";
         for (int i = 0 ; i < iterations; i++){
             rs.run_and_get_times(decoding, translocating);
             total_decoding += decoding;
             total_translocating += translocating;
             n++;
             //save vector.
             vector[i] = decoding;
             vector[iterations + i] = translocating;
         }
         //write times and vector to files.
         averageTimesFile<<", "<<(total_decoding)/iterations<<"\n";
         for (int j=0; j < (2 * iterations) - 1; j++) timesVectorFile<<vector[j]<<",";
         timesVectorFile<<vector[(2 * iterations) - 1]<<vector[(2 * iterations) - 1]<<"\n";
         codons_times[codon] = decoding;
     }
     //save translocation times.
     averageTimesFile << "tra, "<<(total_translocating/n)<<"\n";
     codons_times["tra"] = (total_translocating/n);
     // close files.
     averageTimesFile.close();
     timesVectorFile.close();
     return codons_times;
 }
 
 std::map<std::string, double> load_average_times_file(std::string average_times_file_name)
 {
     std::ifstream ist{average_times_file_name};
     std::map<std::string, double> codons_times;
     
     if (!ist) {
         throw std::runtime_error("can't open input file: "+ average_times_file_name);
     }
     bool header = true;
     std::string codon, tmp_str;
     while(ist.good()) {
         if (!std::getline ( ist, codon, ',' )){
             break;
         }
         std::getline ( ist, codon, ',' );
         std::getline ( ist, tmp_str, ',' );
         if (!header){
             codons_times[codon] = std::atof(tmp_str.c_str());
         } else {
             header = false;
         }
     }
     return codons_times;
 }
 
 void test_ribosome_simulator(std::string concentrations_file_name="../../RSim/data_with_times/concentrations.csv")
 {
     csv_utils::ConcentrationsReader cr;
     cr.loadConcentrations(concentrations_file_name);
     Simulations::RibosomeSimulator rs;
     rs.loadConcentrations(concentrations_file_name);
     rs.setIterationLimit(2000);
     // create a matrix with the initial species population.
     rs.setNumberOfRibosomes(1);
     rs.setCodonForSimulation("AAA");
     rs.run();
     for (float dt:rs.dt_history) {
         std::cout<<"dt = "<<dt<<"\n";
     }
     double decoding, translocating;
     rs.run_and_get_times(decoding, translocating);
     std::cout<< "decoding : "<< decoding<<", translocating = "<<translocating<<"\n";
 }
 
 void test_mRNA_reader()
 {
     mRNAReader mrr;
     mrr.loadRateCalculatorFile("../data/codons/average_time.csv");
     mrr.loadmRNAFile("../../PolisomeSimulator/mRNA/S288C_YOR271C_FSF1_coding.fsa");
     std::cout<<"mRNA = "<<"\n"<<mrr.mRNA_sequence;
     mrr.generateInitialPopulation();
     std::cout<<"initial population = \n";
     std::cout<<mrr.initial_population<<"\n";
     mrr.generateReactions();
 }
 
 void test_enlogation_simulator(std::string average_times_file_name="../data/codons/average_time.csv", std::string concentrations_file_name="../../RSim/data_with_times/concentrations.csv", int init_rate = 1, int term_rate = 10, int iterations = 100000)
 {
     //run enlogation simulator.
     Simulations::EnlogationSimulator es;
     es.setAverageTimesFileName(average_times_file_name);
     es.setConcentrationsFileName(concentrations_file_name);
     es.setMRnaFileName("../../PolisomeSimulator/mRNA/mRNA_sample_MinCFLuc.txt");
     es.setInitiationRate(init_rate);
     es.setTerminationRate(term_rate);
     es.setIterationLimit(iterations);
     es.run();
     es.updateRibosomeHistory();
     int i = 0;
     for (std::vector<int> rib_pos_vec:es.ribosome_positions_history) {
         std::cout<<"iteration = "<< i <<"positions = ";
         for (int pos:rib_pos_vec) std::cout<<" "<<pos;
         std::cout<<"\n";
         i++;
     }
 }
 
 int main(int argc, char **argv) {
     std::cout << "Hello, world!" << std::endl;
     test_enlogation_simulator();
     
     return 0;
 }
 
 
