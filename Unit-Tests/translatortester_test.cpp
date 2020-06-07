#include <gtest/gtest.h>
#include <stdlib.h>
#include <numeric>
#include "../translation.h"

TEST(TranslatorTester, simulateAAAx100xlowInitxHighTerm)
{
  Simulations::Translation ts;
  char const *home_dir;
  home_dir = getenv("HOME");
  std::string home_path(home_dir);
  ts.loadConcentrations(home_path +
                        "/Projects/RSim/data_2011/concentrations.csv");
  ts.loadMRNA(home_path +
              "/Projects/translation_results/initiation_frequency_tests/"
              "20-09-2017/genes/AAA.txt");
  ts.setInitiationRate(0.001);
  ts.setTerminationRate(1000);
  ts.setFinishedRibosomes(1000);

  ts.setPrepopulate(false);
  ts.run();
  ts.calculateAverageTimes();

  std::vector<double> enlongation_duration;
  std::vector<int> iteration_initiation;
  std::tie(enlongation_duration, iteration_initiation) =
      ts.getEnlogationDuration();
  ts.getInitiationElongationTermination();
  double total = 0;
  for (auto dur : ts.elongations_durations)
    total += dur;
  double average = total / ts.elongations_durations.size();
  ASSERT_LE(average / 100, 0.052 * 1.1);
  ASSERT_GE(average / 100, 0.052 * 0.9);
  std::cerr << "\naverage enlongation: " << average;
  total = 0;
  for (unsigned int i = 1; i < ts.codons_average_occupation_time.size() - 1;
       i++)
    total += ts.codons_average_occupation_time[i];
  average = total / ts.codons_average_occupation_time.size();
  // 0.052+-10% tolerance is acceptable
  ASSERT_LE(average, 0.052 * 1.1);
  ASSERT_GE(average, 0.052 * 0.9);
  std::cerr << "\naverage codon_usage: " << average << "\n";
}

TEST(TranslatorTester, checkSpaceBetweenRibosomes)
{
  Simulations::Translation ts;
  char const *home_dir;
  home_dir = getenv("HOME");
  std::string home_path(home_dir);
  ts.loadConcentrations(home_path +
                        "/Projects/RSim/data_2011/concentrations.csv");
  ts.loadMRNA(home_path +
              "/Projects/translation_results/initiation_frequency_tests/"
              "20-09-2017/genes/AAA.txt");
  ts.setInitiationRate(100);
  ts.setTerminationRate(1000);
  ts.setFinishedRibosomes(1000);

  ts.setPrepopulate(true);
  ts.run();
  for (auto ribosomes_positions : ts.ribosome_positions_history)
  {
    for (std::size_t i = 1; i < ribosomes_positions.size(); i++)
    {
      ASSERT_GE(ribosomes_positions[i] - ribosomes_positions[i - 1], 10) << "\nOverlapping ribosomomes detected at positions: " << i - 1 << " and " << i <<", codons: "<<ribosomes_positions[i]<<" "<<ribosomes_positions[i-1]<< "\n";
    }
  }
}

TEST(TranslatorTester, Initialize) { EXPECT_EQ(1, 1); }
