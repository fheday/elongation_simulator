#include <gtest/gtest.h>
#include <cstdlib>
#include <numeric>
#include <string>
#include "../sequence_simulator.h"

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

TEST(SequenceSimulatorTester, simulateAAAx100xlowInitxHighTerm)
{
  Simulations::SequenceSimulator ts;
  char const *home_dir;
  home_dir = getenv("HOME");
  if (home_dir == nullptr) {
    ADD_FAILURE();
  }
  std::string home_path(home_dir);
  ts.loadConcentrationsFromString(concentrationsString);
  ts.inputMRNA(mRNA_aaa);
  ts.setInitiationRate(0.001);
  ts.setTerminationRate(1000);
  ts.setFinishedRibosomes(1000);

  ts.setPrepopulate(false);
  ts.run();
  ts.getAverageTimes();

  std::vector<float> enlongation_duration;
  std::vector<int> iteration_initiation;
  std::tie(enlongation_duration, iteration_initiation) =
      ts.getElongationDuration();
  ts.getInitiationElongationTermination();
  float total = 0;
  for (auto dur : ts.elongations_durations)
    total += dur;
  float average = total / static_cast<float>(ts.elongations_durations.size());
  ASSERT_LE(average / 100, 0.052 * 1.1);
  ASSERT_GE(average / 100, 0.052 * 0.9);
  std::cerr << "\naverage enlongation: " << average;
  total = 0;
  for (unsigned int i = 1; i < ts.codons_average_occupation_time.size() - 1;
       i++)
    total += ts.codons_average_occupation_time[i];
  average = total / static_cast<float>(ts.codons_average_occupation_time.size());
  // 0.052+-10% tolerance is acceptable
  ASSERT_LE(average, 0.052 * 1.1);
  ASSERT_GE(average, 0.052 * 0.9);
  std::cerr << "\naverage codon_usage: " << average << "\n";
}

TEST(SequenceSimulatorTester, checkSpaceBetweenRibosomes)
{
  Simulations::SequenceSimulator ts;
  char const *home_dir;
  home_dir = getenv("HOME");
  if (home_dir == nullptr) {
    ADD_FAILURE();
  }
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

TEST(SequenceSimulatorTester, loadConcentrationsFromString)
{
  Simulations::SequenceSimulator ts;
  
  ts.loadConcentrationsFromString(concentrationsString);
  ts.inputMRNA(mRNA_aaa);
  
  

  ts.setInitiationRate(0.001);
  ts.setTerminationRate(1000);
  ts.setFinishedRibosomes(1000);

  ts.setPrepopulate(false);
  ts.run();
  ts.getAverageTimes();

  std::vector<float> enlongation_duration;
  std::vector<int> iteration_initiation;
  std::tie(enlongation_duration, iteration_initiation) =
      ts.getElongationDuration();
  ts.getInitiationElongationTermination();
  float total = 0;
  for (auto dur : ts.elongations_durations)
    total += dur;
  float average = total / static_cast<float>(ts.elongations_durations.size());
  ASSERT_LE(average / 100, 0.052 * 1.1);
  ASSERT_GE(average / 100, 0.052 * 0.9);
  std::cerr << "\naverage enlongation: " << average;
  total = 0;
  for (unsigned int i = 1; i < ts.codons_average_occupation_time.size() - 1;
       i++)
    total += ts.codons_average_occupation_time[i];
  average = total / static_cast<float>(ts.codons_average_occupation_time.size());
  // 0.052+-10% tolerance is acceptable
  ASSERT_LE(average, 0.052 * 1.1);
  ASSERT_GE(average, 0.052 * 0.9);
  std::cerr << "\naverage codon_usage: " << average << "\n";
}

TEST(SequenceSimulatorTester, testStopAfter5kRibosomesFinished)
{
  Simulations::SequenceSimulator ts;
  std::string mRNA = "ATGAAATGGGTCACCTTTATATCCTTATTATTTTTATTTTCCTCCGCCTACTCCGTCTTTACCTTAGAGGATTTTGTCGGTGATTGGCGCCAGACCGCCGGTTACAATTTAGATCAGGTCTTAGAGCAGGGTGGTGTCTCCTCCTTATTTCAGAATTTAGGTGTCTCCGTCACCCCCATACAGCGCATAGTCTTATCCGGTGAGAATGGTTTAAAAATAGATATACATGTCATAATACCCTACGAGGGTTTATCCGGTGATCAGATGGGTCAGATAGAGAAAATATTTAAAGTCGTCTACCCCGTCGATGATCATCATTTTAAAGTCATATTACATTACGGTACCTTAGTCATAGATGGTGTCACCCCCAATATGATAGATTACTTTGGTCGCCCCTACGAGGGTATAGCCGTCTTTGATGGTAAAAAAATAACCGTCACCGGTACCTTATGGAATGGTAATAAAATAATAGATGAGCGCTTAATAAATCCCGATGGTTCCTTATTATTTCGCGTCACCATAAATGGTGTCACCGGTTGGCGCTTATGTGAGCGCATATTAGCCTAA";
  std::string SKMC_concentration= R"(codon,three.letter,WCcognate.conc,wobblecognate.conc,nearcognate.conc
AAA,Lys,2.799655413453454e-07,0.0,0.0
AAG,Lys,8.742289137825855e-06,2.799655413453454e-07,0.0
AAC,Asn,0.0,4.0427162409271773e-07,1.2602261705628357e-07
AAU,Asn,0.0,4.0427162409271773e-07,1.2602261705628357e-07
AGA,Arg,1.3001590905276952e-07,0.0,8.615848398706653e-08
AGG,Arg,4.954401906136258e-07,1.3001590905276952e-07,1.4041642881543133e-07
AGC,Ser,1.95222916111089e-06,0.0,5.976032159054159e-07
AGU,Ser,0.0,1.95222916111089e-06,5.976032159054159e-07
ACA,Thr,1.3593112881062503e-07,1.7769062748534778e-07,1.4233988185365818e-06
ACG,Thr,1.4489408458534457e-07,3.136217562959728e-07,1.8427944064050975e-06
ACC,Thr,0.0,3.136217562959728e-07,1.4233988185365818e-06
ACU,Thr,1.7769062748534778e-07,1.3593112881062503e-07,1.4233988185365818e-06
AUA,Ile,2.0685143589853748e-09,8.447957805395312e-08,4.928212722521385e-07
AUG,Met,1.8330032065783554e-06,0.0,2.6010206495582227e-06
AUC,Ile,5.965425817089019e-10,8.654809241293849e-08,8.118894584919365e-08
AUU,Ile,8.447957805395312e-08,2.6650569406942766e-09,8.118894584919365e-08
GAA,Glu,1.2030189060348958e-05,0.0,2.799655413453454e-07
GAG,Glu,2.0739284418578572e-05,1.2030189060348958e-05,9.0222546791712e-06
GAC,Asp,2.0762194427982553e-05,0.0,1.256048330149796e-05
GAU,Asp,0.0,2.0762194427982553e-05,1.256048330149796e-05
GGA,Gly,5.877465446236063e-06,0.0,2.1617439303983605e-07
GGG,Gly,2.0085334765806074e-06,5.877465446236063e-06,7.658725284818267e-07
GGC,Gly,1.3821708360519286e-05,5.877465446236063e-06,2.549832377016306e-06
GGU,Gly,0.0,1.969917380675535e-05,2.549832377016306e-06
GCA,Ala,1.3840658242456562e-06,2.694175405276694e-06,1.7370205748325546e-06
GCG,Ala,5.860853384274279e-07,1.3840658242456562e-06,2.3013102472864148e-06
GCC,Ala,0.0,4.0782412295223505e-06,1.7370205748325546e-06
GCU,Ala,0.0,4.0782412295223505e-06,1.7370205748325546e-06
GUA,Val,1.1017248100696978e-07,3.60707404910023e-07,5.79369364665077e-07
GUG,Val,4.805700533079567e-07,4.708798859169928e-07,4.4340238561365785e-06
GUC,Val,0.0,3.60707404910023e-07,1.6833358084384105e-07
GUU,Val,3.60707404910023e-07,0.0,1.6833358084384105e-07
CAA,Gln,2.8191637298206445e-07,0.0,2.799655413453454e-07
CAG,Gln,1.539943158404245e-06,2.8191637298206445e-07,9.0222546791712e-06
CAC,His,7.867549703031669e-06,0.0,5.302942411490014e-07
CAU,His,0.0,7.867549703031669e-06,5.302942411490014e-07
CGA,Arg,7.85735205935732e-08,2.804276145377727e-07,2.1617439303983605e-07
CGG,Arg,1.644928678589188e-07,7.85735205935732e-08,7.658725284818267e-07
CGC,Arg,0.0,2.804276145377727e-07,2.549832377016306e-06
CGU,Arg,0.0,2.804276145377727e-07,2.549832377016306e-06
CCA,Pro,4.5429860571364995e-06,8.673342048393995e-06,1.7370205748325546e-06
CCG,Pro,2.9675901574042905e-06,1.3216328105530494e-05,2.3013102472864148e-06
CCC,Pro,0.0,1.3216328105530494e-05,1.7370205748325546e-06
CCU,Pro,8.673342048393995e-06,4.5429860571364995e-06,1.7370205748325546e-06
CUA,Leu,1.687829190445213e-05,3.247993715165786e-05,5.79369364665077e-07
CUG,Leu,1.3790582812441081e-05,4.9358229056109994e-05,4.4340238561365785e-06
CUC,Leu,0.0,4.9358229056109994e-05,1.6833358084384105e-07
CUU,Leu,3.247993715165786e-05,1.687829190445213e-05,1.6833358084384105e-07
UAA,Stop,0.0,0.0,5.618819143274099e-07
UAG,Stop,0.0,0.0,1.0844114210557508e-05
UAC,Tyr,0.0,1.2602261705628357e-07,8.271821327124387e-06
UAU,Tyr,0.0,1.2602261705628357e-07,8.271821327124387e-06
UGA,Stop,8.615848398706653e-08,0.0,4.890170441841154e-07
UGG,Trp,5.42579448283648e-08,0.0,9.546809721059538e-07
UGC,Cys,5.114447319183494e-07,0.0,2.3188152596357293e-06
UGU,Cys,0.0,5.114447319183494e-07,2.3188152596357293e-06
UCA,Ser,5.409893906661975e-07,8.824094278703843e-07,1.3529949861826465e-05
UCG,Ser,4.1939558786851586e-07,1.4233988185365818e-06,1.66424341038161e-05
UCC,Ser,0.0,1.4233988185365818e-06,1.3529949861826465e-05
UCU,Ser,8.824094278703843e-07,5.409893906661975e-07,1.3529949861826465e-05
UUA,Leu,4.928212722521385e-07,0.0,4.944477714852293e-05
UUG,Leu,2.0216512848931454e-06,4.928212722521385e-07,6.506836316754237e-05
UUC,Phe,8.118894584919365e-08,0.0,4.944537369110464e-05
UUU,Phe,0.0,8.118894584919365e-08,4.944537369110464e-05
)";

  ts.loadConcentrationsFromString(SKMC_concentration);
  ts.inputMRNA(mRNA);
  ts.setInitiationRate(0.001);
  ts.setTerminationRate(10000);
  ts.setFinishedRibosomes(5000);
  ts.setHistorySize(6000000);

  ts.setPrepopulate(false);
  ts.run();
  ts.getAverageTimes();

  std::vector<float> enlongation_duration;
  std::vector<int> iteration_initiation;
  std::tie(enlongation_duration, iteration_initiation) =
      ts.getElongationDuration();
  ts.getInitiationElongationTermination();
  float average = std::reduce(ts.elongations_durations.begin(), ts.elongations_durations.end()) / static_cast<float>(ts.elongations_durations.size());
  
  // ASSERT_LE(average / 100, 0.052 * 1.1);
  // ASSERT_GE(average / 100, 0.052 * 0.9);
  std::cerr << "\naverage enlongation: " << average;
  average = std::reduce(ts.codons_average_occupation_time.begin() + 1, ts.codons_average_occupation_time.end() - 1) / static_cast<float>(ts.codons_average_occupation_time.size() - 2);
  
  // total = 0;
  // for (unsigned int i = 1; i < ts.codons_average_occupation_time.size() - 1;
  //      i++)
  //   total += ts.codons_average_occupation_time[i];
  // average = (total / static_cast<float>(ts.codons_average_occupation_time.size()) - 2);
  // 0.052+-10% tolerance is acceptable
  std::cerr << "\naverage codon_usage: " << average << "\n";
  ASSERT_LE(average, 2.52 * 1.1);
  ASSERT_GE(average, 2.52 * 0.9);
}