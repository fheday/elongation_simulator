#include <gtest/gtest.h>
#include <stdlib.h>
#include <string>
#include "../elongation_simulation_manager.h"
#include "../sequence_simulator.h"


TEST(SequenceSimulatorManagerTester, openjsonconfig_file01)
{
    std::string conf_file_path = "../../data/configurations/sim_test_01.json";
    Elongation_manager::SimulationManager sim_man(conf_file_path);
    ASSERT_EQ(sim_man.get_configuration_file_path(), conf_file_path);
    ASSERT_TRUE(sim_man.get_concentration_file_path().find("concentrations.csv") != std::string::npos);
    ASSERT_EQ(sim_man.get_pre_populate(), false);
    ASSERT_EQ(sim_man.get_stop_condition_type(), Elongation_manager::stop_condition_enum::TIME);
    ASSERT_EQ(sim_man.get_stop_condition_value(), 60.0);
    auto simulations_configurations = sim_man.get_simulations_configurations();
    ASSERT_EQ(simulations_configurations.size(), 0);
}

TEST(SequenceSimulatorManagerTester, openjsonconfig_file02)
{
    std::string conf_file_path = "../../data/configurations/sim_test_02.json";
    Elongation_manager::SimulationManager sim_man(conf_file_path);
    ASSERT_EQ(sim_man.get_configuration_file_path(), conf_file_path);
    std::cout<<sim_man.get_concentration_file_path();
    ASSERT_TRUE(sim_man.get_concentration_file_path().find("../../concentrations/Saccharomyces_cerevisiae.csv") != std::string::npos);
    ASSERT_EQ(sim_man.get_pre_populate(), false);
    ASSERT_EQ(sim_man.get_stop_condition_type(), Elongation_manager::stop_condition_enum::TIME);
    // ASSERT_EQ(sim_man.get_stop_condition_value(), 60.0);
    auto simulations_configurations = sim_man.get_simulations_configurations();
    ASSERT_EQ(simulations_configurations.size(), 19);
    std::vector<std::string> genes{
        "YAL001C",
        "YAL002W",
        "YAL003W",
        "YAL005C",
        "YAL007C",
        "YAL008W",
        "YAL009W",
        "YAL010C",
        "YAL011W",
        "YAL012W",
        "YAL013W",
        "YAL014C",
        "YAL015C",
        "YAL016W",
        "YAL017W",
        "YAL018C",
        "YAL019W",
        "YAL020C",
        "YAL021C",
    };
    std::vector<float> init_rates{
        0.2154329023064577,
        0.09514953185201884,
        2.808042992925658,
        0.0,
        0.6743205952990898,
        1.1747337607652857,
        0.3721681442743443,
        0.5751902380783647,
        0.3225249102645953,
        1.3902447184711295,
        0.7298960396622051,
        0.35694734139399675,
        0.4447596657036942,
        0.7988189573203942,
        0.26250030813645553,
        1.1242319173622863,
        0.0,
        0.6726033765125892,
        0.5831518888157773,
    };

    std::vector<float> copy_numbers{
        4.5748826822239215,
        3.4008593956858126,
        117.12649114831306,
        46.59278180819215,
        8.743500440522217,
        2.0358152806345093,
        2.283805139523944,
        1.9460999389449092,
        2.1264098112282053,
        26.00494089214123,
        2.978069539423594,
        2.943713707607156,
        2.3996130091556167,
        12.327135682446167,
        5.4415462409526905,
        0.023552614479173062,
        7.665040692561605,
        2.6149421432521502,
        6.1667358845246065
    };
    int i = 0;
    for (auto [fasta_file_path, gene_name, initiation_rate, termination_rate, gene_copy_number] : simulations_configurations)
    {
        ASSERT_TRUE(fasta_file_path.find("Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa") != std::string::npos);
        ASSERT_EQ(gene_name, genes[i]);
        ASSERT_EQ(initiation_rate, init_rates[i]);
        ASSERT_EQ(termination_rate, 10);
        ASSERT_EQ(gene_copy_number, copy_numbers[i]);
        i++;
    }
}
TEST(SequenceSimulatorManagerTester, parallel_simulation_file02) {
    std::string conf_file_path = "../../data/configurations/sim_test_02.json";
    Elongation_manager::SimulationManager sim_man(conf_file_path);
    sim_man.start(true, 4);
}

TEST(SequenceSimulatorManagerTester, propensity_change) {
    std::string conf_file_path = "../../data/configurations/sim_test_03.json";
    Elongation_manager::SimulationManager sim_man(conf_file_path);
    for (const auto &item : sim_man.get_reactions_modifiers()) {
        std::cout<<item.first<<":"<<item.second<<"\n";
    }
    sim_man.start();
}