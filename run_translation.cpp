#include <iostream>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <vector>
#include <execinfo.h>
#include <unistd.h>
#include <signal.h>
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

void execute_translation(std::string concentrations_file, std::string mrna_file, float initiation_rate, float termination_rate, int time_limit, std::string output_file_name ) {
    Simulations::Translation ts;
    ts.loadConcentrations(concentrations_file);
    ts.loadMRNA(mrna_file);
    ts.setInitiationRate(initiation_rate);
    ts.setTerminationRate(termination_rate);
    ts.setTimeLimit(time_limit);
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
    //TODO: save codon average time into csv file.
    std::ofstream codon_average_time_file;
    codon_average_time_file.open("codon_average_time_"+output_file_name);
    //header
    codon_average_time_file<<"codon average time\n";
    //data
    for (auto average_occupation_time:ts.codons_average_occupation_time) codon_average_time_file<<std::fixed<<std::setprecision(10)<<average_occupation_time<<"\n";
    codon_average_time_file.close();
}

int main(int argc, char **argv) {
    signal(SIGSEGV, handler);   // install our handler
    if (argc != 7){
        std::cout<<"Wrong number of parameters informed.\n";
        std::cout<<"run_translation concentrations mRNA Initiation Termination Time output\n";
        std::cout<<"Concentrations = path to the file containing the concentrations to be used in the simulation.\n";
        std::cout<<"mRNA = path to the file with the mRNA to be used.\n";
        std::cout<<"Initiation = value to be used as the initiation factor.\n";
        std::cout<<"Termination = value to be used as the termination factor.\n";
        std::cout<<"Time limit = time limit for when the simulation should stop. This is in simulation time, not in real life time\n";
        std::cout<<"output = file to be created with the simulation results.\n";
        return 0;
    }
    
    execute_translation(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stoi(argv[5]), argv[6]);
    
    return 0;
}
