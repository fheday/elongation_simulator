#include <gtest/gtest.h>
#include <stdlib.h>
#include <numeric>
#include <string>
#include "../translation.h"

std::string concentrationsString = R"(codon,three.letter,WCcognate.conc,wobblecognate.conc,nearcognate.conc
AAA,Lys,4.90774907749077E-06,0,1.40221402214022E-06
AAC,Asn,7.01107011070111E-06,0,1.05166051660517E-05
AAG,Lys,9.81549815498155E-06,0,6.309963099631E-06
AAU,Asn,7.01107011070111E-06,0,1.19188191881919E-05
ACA,Thr ,2.80442804428044E-06,0,1.75276752767528E-05
ACC,Thr ,7.71217712177122E-06,0,1.12177121771218E-05
ACG,Thr ,7.01107011070111E-07,0,2.80442804428044E-06
ACU,Thr ,7.71217712177122E-06,0,1.05166051660517E-05
AGA,Arg,7.71217712177122E-06,0,0
AGC,Ser,2.10332103321033E-06,0,1.05166051660517E-05
AGG,Arg,1.40221402214022E-06,0,1.19188191881919E-05
AGU,Ser,2.10332103321033E-06,0,1.05166051660517E-05
AUA,Ile,1.40221402214022E-06,0,1.40221402214022E-05
AUC,Ile,9.11439114391144E-06,0,8.41328413284133E-06
AUG,Met,3.50553505535055E-06,0,8.41328413284133E-06
AUU,Ile,9.11439114391144E-06,0,8.41328413284133E-06
CAA,Gln,5.60885608856089E-06,0,7.01107011070111E-06
CAC,His,4.90774907749077E-06,0,1.82287822878229E-05
CAG,Gln,7.01107011070111E-07,0,1.54243542435424E-05
CAU,His,4.90774907749077E-06,0,1.2619926199262E-05
CCA,Pro,7.01107011070111E-06,0,2.38376383763838E-05
CCC,Pro,1.40221402214022E-06,0,2.24354243542435E-05
CCG,Pro,7.01107011070111E-06,0,2.80442804428044E-06
CCU,Pro,1.40221402214022E-06,0,2.24354243542435E-05
CGA,Arg,4.20664206642066E-06,0,7.71217712177122E-06
CGC,Arg,4.20664206642066E-06,0,4.90774907749077E-06
CGG,Arg,7.01107011070111E-07,0,4.90774907749077E-06
CGU,Arg,4.20664206642066E-06,0,4.90774907749077E-06
CUA,Leu,2.10332103321033E-06,0,1.54243542435424E-05
CUC,Leu,7.01107011070111E-07,0,1.82287822878229E-05
CUG,Leu,2.10332103321033E-06,0,1.05166051660517E-05
CUU,Leu,2.10332103321033E-06,0,1.68265682656827E-05
GAA,Glu,9.81549815498155E-06,0,6.309963099631E-06
GAC,Asp,1.12177121771218E-05,0,2.24354243542435E-05
GAG,Glu,1.40221402214022E-06,0,1.96309963099631E-05
GAU,Asp,1.12177121771218E-05,0,2.38376383763838E-05
GCA,Ala,3.50553505535055E-06,0,2.80442804428044E-05
GCC,Ala,7.71217712177122E-06,0,1.40221402214022E-05
GCG,Ala,3.50553505535055E-06,0,1.40221402214022E-06
GCU,Ala,7.71217712177122E-06,0,1.40221402214022E-05
GGA,Gly,2.10332103321033E-06,0,7.71217712177122E-06
GGC,Gly,1.12177121771218E-05,0,7.01107011070111E-06
GGG,Gly,1.40221402214022E-06,0,7.01107011070111E-06
GGU,Gly,1.12177121771218E-05,0,7.01107011070111E-06
GUA,Val,1.40221402214022E-06,0,2.5239852398524E-05
GUC,Val,9.81549815498155E-06,0,1.61254612546125E-05
GUG,Val,1.40221402214022E-06,0,1.05166051660517E-05
GUU,Val,9.81549815498155E-06,0,1.61254612546125E-05
UAA,Stop,1.14E-06,0,1.19188191881919E-05
UAC,Tyr,5.60885608856089E-06,0,1.75276752767528E-05
UAG,Stop,1.14E-06,0,2.24354243542435E-05
UAU,Tyr,5.60885608856089E-06,0,1.8929889298893E-05
UCA,Ser,2.10332103321033E-06,0,3.01476014760148E-05
UCC,Ser,7.71217712177122E-06,0,1.12177121771218E-05
UCG,Ser,7.01107011070111E-07,0,1.75276752767528E-05
UCU,Ser,7.71217712177122E-06,0,1.12177121771218E-05
UGA,Stop,1.14E-06,0,1.19188191881919E-05
UGC,Cys,2.80442804428044E-06,0,6.309963099631E-06
UGG,Trp,4.20664206642066E-06,0,2.54221402214022E-06
UGU,Cys,2.80442804428044E-06,0,6.309963099631E-06
UUA,Leu,4.90774907749077E-06,0,1.2619926199262E-05
UUC,Phe,7.01107011070111E-06,0,1.68265682656827E-05
UUG,Leu,7.01107011070111E-06,0,1.05166051660517E-05
UUU,Phe,7.01107011070111E-06,0,1.68265682656827E-05
)";

std::string mRNA_aaa = "AUGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUGA";

TEST(TranslatorTester, simulateAAAx100xlowInitxHighTerm)
{
  Simulations::Translation ts;
  char const *home_dir;
  home_dir = getenv("HOME");
  std::string home_path(home_dir);
  ts.loadConcentrationsFromString(concentrationsString);
  ts.inputMRNA(mRNA_aaa);
  ts.setInitiationRate(0.001);
  ts.setTerminationRate(1000);
  ts.setFinishedRibosomes(1000);

  ts.setPrepopulate(false);
  ts.run();
  ts.getAverageTimes();

  std::vector<double> enlongation_duration;
  std::vector<int> iteration_initiation;
  std::tie(enlongation_duration, iteration_initiation) =
      ts.getElongationDuration();
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
  ts.loadConcentrationsFromString(concentrationsString);
  ts.inputMRNA(mRNA_aaa);
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

TEST(TranslatorTester, loadConcentrationsFromString)
{
  Simulations::Translation ts;
  
  ts.loadConcentrationsFromString(concentrationsString);
  ts.inputMRNA(mRNA_aaa);
  
  

  ts.setInitiationRate(0.001);
  ts.setTerminationRate(1000);
  ts.setFinishedRibosomes(1000);

  ts.setPrepopulate(false);
  ts.run();
  ts.getAverageTimes();

  std::vector<double> enlongation_duration;
  std::vector<int> iteration_initiation;
  std::tie(enlongation_duration, iteration_initiation) =
      ts.getElongationDuration();
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
