#include <gtest/gtest.h>
#include "../mrna_reader.h"

TEST(mRNAReaderTester, open_simple_fasta_file)
{
    //Test to read file with only one DNA sequence without header.
    mRNA_utils::mRNAReader reader;
    reader.loadmRNAFile("data/mRNAs/MaxCFLuc.txt");
    ASSERT_EQ(reader.getCodon(23), "GAA");
    ASSERT_EQ(reader.getCodon(21), "GCU"); // the reader converts T's to U's automatically.
}

TEST(mRNAReaderTester, open_fasta_file_one_gene_with_header) {
    //Test to read file with only one DNA sequence with header.
    mRNA_utils::mRNAReader reader;
    reader.loadmRNAFile("data/mRNAs/RTN3.txt");
    ASSERT_EQ(reader.getCodon(23), "GCG");
    ASSERT_EQ(reader.getCodon(22), "UCC"); // the reader converts T's to U's automatically.
    ASSERT_EQ(reader.getCodon(21), "CCG");
    //check if this works when a file has more than one gene.
    //In this case, we should get the first one only.
    reader.loadmRNAFile("data/mRNAs/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa");
    ASSERT_EQ(reader.getCodon(0), "AUG");
    ASSERT_EQ(reader.getCodon(1), "UUC");
    ASSERT_EQ(reader.getCodon(2), "AGC");
    ASSERT_EQ(reader.sizeInCodons(), 62);
}

TEST(mRNAReaderTester, open_fasta_file_many_genes_with_header) {
    //Test to read file with only one DNA sequence with header.
    mRNA_utils::mRNAReader reader;
    const std::string fasta_file = "data/mRNAs/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa";
    reader.loadGene(fasta_file, "YHR055C");
    ASSERT_EQ(reader.getCodon(0), "AUG");
    ASSERT_EQ(reader.getCodon(1), "UUC");
    ASSERT_EQ(reader.getCodon(2), "AGC");
    reader.loadGene(fasta_file, "YPR161C");
    ASSERT_EQ(reader.getCodon(0), "AUG");
    ASSERT_EQ(reader.getCodon(1), "AGU");
    ASSERT_EQ(reader.getCodon(2), "GAU");
    reader.loadGene(fasta_file, "YKL202W");
    ASSERT_EQ(reader.getCodon(0), "AUG");
    ASSERT_EQ(reader.getCodon(1), "CAG");
    ASSERT_EQ(reader.getCodon(3), "UAU");
    auto genes = reader.get_names_in_file(fasta_file);
    ASSERT_EQ(genes.size(), 6713);
    ASSERT_EQ(genes[0], "YHR055C");
    ASSERT_EQ(genes[2], "YOL138C");
    ASSERT_EQ(genes[3], "YDR395W");
    ASSERT_EQ(genes[6712], "YAR061W");
    ASSERT_EQ(genes[6711], "YIL168W");
}
