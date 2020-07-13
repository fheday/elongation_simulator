#include <gtest/gtest.h>
#include "../elongation_simulation_manager.h"
TEST(ElongationSimulatorTester, openjsonconfig_file)
{
    std::string conf_file_path = "Unit-Tests/sim_test_02.json";
  Elongation_manager::SimulationManager sim_man(conf_file_path);
}
