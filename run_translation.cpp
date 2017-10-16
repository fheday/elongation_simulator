#include <iostream>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <vector>
#include <execinfo.h>
#include <unistd.h>
#include <signal.h>
#include <getopt.h>
#include "translation.h"
#include "concentrationsreader.h"

#define EIGEN_NO_DEBUG //disables Eigen's assertions.


void handler(int sig) {
    void *array[10];
    size_t size;
    
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

void execute_translation(std::string concentrations_file, std::string mrna_file, float initiation_rate, float termination_rate, int time_limit, int number_iterations, int number_ribosomes,  std::string output_file_name ) {
    //separate the path from the file name.
    std::size_t found = output_file_name.find_last_of("/\\");
    std::string path = "./"; // current path.
    std::string file_name = output_file_name;
    if (found != std::string::npos) {
        //there is a path.
        path = output_file_name.substr(0,found + 1);
        file_name = output_file_name.substr(found + 1);
    }

    //prepare and run the simulation.
    Simulations::Translation ts;
    ts.loadConcentrations(concentrations_file);
    ts.loadMRNA(mrna_file);
    ts.setInitiationRate(initiation_rate);
    ts.setTerminationRate(termination_rate);
    if (time_limit > 0){
        ts.setTimeLimit(time_limit);
    } else if (number_iterations > 0) {
        ts.setIterationLimit(number_iterations);
    } else if (number_ribosomes > 0) {
        ts.setFinishedRibosomes(number_ribosomes);
    }
    
    ts.setPrepopulate(true); // simulations pre-populate the mRNA by default. This can be changed in the future.
    ts.run();
    ts.calculateAverageTimes();

    std::vector<double> enlongation_duration;
    std::vector<int> iteration_initiation, iteration_termination;
    
    std::tie(enlongation_duration, iteration_initiation, iteration_termination) = ts.getEnlogationDuration();
    //save enlongation data into csv file.
    std::vector<double> clock, clock_at_initiation;
    double c = 0;
    for(auto dt:ts.dt_history) {
        c += dt;
        clock.push_back(c);
    }
    //get the clock at the initiation of each terminating ribosome.
    for (int iteration:iteration_initiation) clock_at_initiation.push_back(clock[iteration]);
    //now we save the clock_at_initiation and enlongation_duration.
    std::ofstream clock_and_enlongation_csv_file;
    clock_and_enlongation_csv_file.open(output_file_name);
    //header
    clock_and_enlongation_csv_file<<"Clock at initiation, Ribosome enlongation duration\n";
    //data
    for (int i = 0; i <clock_at_initiation.size(); i++) clock_and_enlongation_csv_file<<std::fixed<<std::setprecision(10)<<clock_at_initiation[i]<<", "<<enlongation_duration[i]<<"\n";
    clock_and_enlongation_csv_file.close();

    std::ofstream codon_average_time_file;
    codon_average_time_file.open(path + "codon_average_time_" + file_name);
    //header
    codon_average_time_file<<"codon average time\n";
    //data
    for (auto average_occupation_time:ts.codons_average_occupation_time) codon_average_time_file<<std::fixed<<std::setprecision(10)<<average_occupation_time<<"\n";
    codon_average_time_file.close();
}

void printHelp() {
    std::cout<<"Wrong number of parameters informed.\n";
    std::cout<<"run_translation concentrations mRNA Initiation Termination Time output\n";
    std::cout<<"Concentrations = path to the file containing the concentrations to be used in the simulation.\n";
    std::cout<<"mRNA = path to the file with the mRNA to be used.\n";
    std::cout<<"Initiation = value to be used as the initiation factor.\n";
    std::cout<<"Termination = value to be used as the termination factor.\n";
    std::cout<<"Time limit = time limit for when the simulation should stop. This is in yeast time, not in real life time\n";
    std::cout<<"output = file to be created with the simulation results.\n";
}

int main(int argc, char **argv) {
    signal(SIGSEGV, handler);   // install our handler
    const char* const short_opts = "c:m:i:t:y:r:l:o:h";
    const option long_opts[] = {
        {"concentration", 1, nullptr, 'w'},
        {"mrna", 1, nullptr, 'w'},
        {"initiation", 1, nullptr, 's'},
        {"termination", 1, nullptr, 's'},
        {"yeasttime", 1, nullptr, 's'},
        {"ribosomes", 1, nullptr, 's'},
        {"iterations", 1, nullptr, 's'},
        {"output", 1, nullptr, 'w'},
        {"help", 0, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };
    
    std::string concentration_file, mrna_file, output_file;
    double initiation, termination, yeast_time, ribosomes, iterations;
    bool stop_condition_passed = false;
    yeast_time = ribosomes = iterations = -1;
    
    std::string halting_condition_error = "only one of the following halting options can be used: yeast time, terminating ribosomes, or iteration limit\n";
    while (optind < argc) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
        if (opt != -1) {
            // Option argument
            switch (opt) {
                case 'c': {
                    concentration_file = std::string(optarg);
                    break;
                }
                case 'm': {
                    mrna_file = std::string(optarg);
                    break;
                }
                case 'i':
                    initiation = std::stof(optarg);
                    break;
                case 't':
                    termination = std::stof(optarg);
                    break;
                case 'y':
                    if (!stop_condition_passed){
                        yeast_time = std::stoi(optarg);
                        stop_condition_passed = true;
                    } else {
                        std::cout<<halting_condition_error;
                        return -1;
                    }
                    
                    break;
                case 'r':
                    if (!stop_condition_passed){
                        ribosomes = std::stoi(optarg);
                        stop_condition_passed = true;
                    } else {
                        std::cout<<halting_condition_error;
                        return -1;
                    }
                    break;
                case 'l':
                    if (!stop_condition_passed) {
                        iterations = std::stoi(optarg);
                        stop_condition_passed = true;
                    } else {
                        std::cout<<halting_condition_error;
                        return -1;
                    }
                    break;
                case 'o':
                    output_file = std::string(optarg);
                    break;
                case 'h': // -h or --help
                case '?': // Unrecognized option
                default:
                    printHelp();
                    break;
            }
        } else {
            break;
        }
        
    }
    
    if (optind == 1) {
        // Regular argument
        int index = optind;
        while (index < argc){
            switch(index){
                case 1:
                    concentration_file = argv[index];
                    break;
                case 2:
                    mrna_file = argv[index];
                    break;
                case 3:
                    initiation = std::stof(argv[index]);
                    break;
                case 4:
                    termination = std::stof(argv[index]);
                    break;
                case 5:
                    yeast_time = std::stoi(argv[index]);
                    break;
                case 6:
                    output_file = argv[index];
                    break;
                    
            }
//             std::cout<<": "<<argv[index]<<"\n";
            index++;  // Skip to the next argument
        }
    }
    
    
    
    //     
        execute_translation(concentration_file, mrna_file, initiation, termination, yeast_time, iterations, ribosomes, output_file);
    
    return 0;
}
