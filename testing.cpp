/*
  This file is make to test the performance of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  

  Author: Cole Koryto
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <random>

// Declares testing constants 
const int MAX_P = 512;
const int MAX_T = 1024;
const int MAX_SEQ_LEN = 256;
const int MAX_SEQ_NUM = 1;
const int MIN_P = 2;
const int MIN_T = 1;
const int MIN_SEQ_LEN = 1;
const int P_GROW_FACTOR = 2;
const int T_GROW_FACTOR = 2;
const int SEQ_LEN_GROW_FACTOR = 2;
const char *SEQ_LETTERS[] = {"A", "C", "T", "G"};
const char *SEQUENCES_ALL[] = {"TATA", "CCAAT", "GGCGC"}; // For findAll()
const char *SEQUENCES_FIRST[] = {"AGGGATTGAAGACTTCAGGGGAGGAAGTAACTGCAGATGTAATTTAAAAAAATCAAAAGAACAGGAAATGGAGCCTGGACATGTGACTGATT", "AACTGAAACAGTGGCCACTCTCAACCAAGGAAGAGGGAACTCAGGACAAAGAAAAGCAGCATTTAGCTGACAATCTTCAACCTTTCTTCTCT", "AAAAAAAAGAAAAAAACACCACACAAAAAAAAAAAAAATTTTTTTTTTCTTCTGTTTCTTGTTTCGTCGCACACACCAAAAAACAAGCAGTG", "XGAAAATAAATCCACCAAAACAATCATGGAAGGTTATGTGTTCAGCAAAGATCCTTCTTTGAATCTATCAGATATCTACTGAGCACTGAATA"}; // For findFirst() , //10%, 50%, 90%, NONE 
const char *SEQUENCES_GREP[] = {"", "", "", ""}; 
const std::string INPUT_FILE = "input8GB.fastq";
const std::string CONNECT_LOGIC_STR = "\n";
const std::string COMP_FILE = "comCmds.sh";
const std::string DECOMP_FILE = "decomCmds.sh";
const std::string FINDFIRST_FILE = "findFirstCmds.sh";
const std::string FINDALL_FILE = "findAllCmds.sh";
const std::string GREP_FILE = "grepCmds.sh";
const std::string FASTA_FILE = "FASTACmds.sh";
const int MAX_INSTANCES = 33000; // Observed maximum for Jigwe 


// --------------------------------------------------------------------------------------------
// Generates a random DNA sequence given a length
std::string getRanSeq(int length)
{

    // Sets up random number generation to pick sequence letters
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> randIndex(0,3); // distribution in range [0, 3]

    // Builds random sequence of given length
    std::string sequence = "";
    for (int i = 0; i < length; ++i)
        sequence += SEQ_LETTERS[randIndex(rng)];

    // Returns final sequence
    return sequence;
}

// --------------------------------------------------------------------------------------------
// Tests the performance of compression
void testCompression()
{

    // Runs through all process combinations
    std::string beginning = "mpiexec -np ";
    std::string middle = " ./main -c ";
    std::string end = " " + INPUT_FILE + " testing_out.ngsc";
    std::string combined_commands = "";
    for (int num_process = MIN_P; num_process <= MAX_P; num_process *= P_GROW_FACTOR)
    {
        // Runs through all thread combinations
        for (int num_thread = MIN_T; num_thread <= MAX_T; num_thread *= T_GROW_FACTOR)
        {
            if(num_process * num_thread >= MAX_INSTANCES)
                break;

            combined_commands += beginning + std::to_string(num_process) + middle + std::to_string(num_thread) + end + CONNECT_LOGIC_STR;
        }
    }

    // Outputs final commands 
    std::ofstream command_file;
    command_file.open(COMP_FILE);
    command_file << "#! /usr/bin/bash\n";
    command_file << combined_commands.c_str();
    command_file.close();
}

// --------------------------------------------------------------------------------------------
// Tests the performance of decompression
void testDecompression()
{

    // Runs through all process combinations
    std::string beginning = "mpiexec -np ";
    std::string middle = " ./main -d ";
    std::string end = " testing_out.ngsc decomp_testing_out.fastq";
    std::string combined_commands = "";
    for (int num_process = MIN_P; num_process <= MAX_P; num_process *= P_GROW_FACTOR)
    {
        // Runs through all thread combinations
        for (int num_thread = MIN_T; num_thread <= MAX_T; num_thread *= T_GROW_FACTOR)
        {
            if(num_process * num_thread >= MAX_INSTANCES)
                break;

            combined_commands += beginning + std::to_string(num_process) + middle + std::to_string(num_thread) + end + CONNECT_LOGIC_STR;
        }
    }

    // Outputs final commands 
    std::ofstream command_file;
    command_file.open(DECOMP_FILE);
    command_file << "#! /usr/bin/bash\n";
    command_file << combined_commands.c_str();
    command_file.close();
}

// --------------------------------------------------------------------------------------------
// Tests the performance of findFirst
void testFindFirst()
{

    // Runs through all process combinations
    std::string beginning = "mpiexec -np ";
    std::string middle = " ./main -i ";
    std::string end = " testing_out.ngsc -findfirst ";
    std::string combined_commands = "";
    for (int num_process = MIN_P; num_process <= MAX_P; num_process *= P_GROW_FACTOR)
    {
        // Runs through all thread combinations
        for (int num_thread = MIN_T; num_thread <= MAX_T; num_thread *= T_GROW_FACTOR)
        {
            // Runs through a set of the sequences of the given length
            for (int num_seq = 0; num_seq < 4; ++num_seq)
            {
                if(num_process * num_thread >= MAX_INSTANCES)
                    break;

                combined_commands += beginning + std::to_string(num_process) + middle + std::to_string(num_thread) + end + SEQUENCES_FIRST[num_seq] + CONNECT_LOGIC_STR;
            }  
        }
    }

    // Outputs final commands 
    std::ofstream command_file;
    command_file.open(FINDFIRST_FILE);
    command_file << "#! /usr/bin/bash\n";
    command_file << combined_commands.c_str();
    command_file.close();
}

// --------------------------------------------------------------------------------------------
// Tests the performance of findAll
void testFindAll()
{

    // Runs through all process combinations
    std::string beginning = "mpiexec -np ";
    std::string middle = " ./main -i ";
    std::string end = " testing_out.ngsc -findall ";
    std::string combined_commands = "";
    for (int num_process = MIN_P; num_process <= MAX_P; num_process *= P_GROW_FACTOR)
    {
        // Runs through all thread combinations
        for (int num_thread = MIN_T; num_thread <= MAX_T; num_thread *= T_GROW_FACTOR)
        {
            if(num_process * num_thread >= MAX_INSTANCES)
                break;

            // Runs through a set of the sequences of the given length
            for (int num_seq = 0; num_seq < 3; ++num_seq)
            {
                combined_commands += beginning + std::to_string(num_process) + middle + std::to_string(num_thread) + end + SEQUENCES_ALL[num_seq] + CONNECT_LOGIC_STR;
            }      
        }
    }

    // Outputs final commands 
    std::ofstream command_file;
    command_file.open(FINDALL_FILE);
    command_file << "#! /usr/bin/bash\n";
    command_file << combined_commands.c_str();
    command_file.close();
}

// --------------------------------------------------------------------------------------------
// Tests the performance of grep
void testGrep()
{

    // Runs through all process combinations
    std::string beginning = "time grep ";
    std::string combined_commands = "";
    // Runs through a set of the sequences of the given length
    for (int num_seq = 0; num_seq < 3; ++num_seq)
    {
        combined_commands += beginning + SEQUENCES_GREP[num_seq] + " " + INPUT_FILE + CONNECT_LOGIC_STR;
    }      

    // Outputs final commands 
    std::ofstream command_file;
    command_file.open(GREP_FILE);
    command_file << "#! /usr/bin/bash\n";
    command_file << combined_commands.c_str();
    command_file.close();
}

// --------------------------------------------------------------------------------------------
// Tests the performance of ToFASTA
void testToFASTA()
{

    // Runs through all process combinations
    std::string beginning = "mpiexec -np ";
    std::string middle = " ./main -i ";
    std::string end = " testing_out.ngsc fasta_testing_out.fasta -fasta";
    std::string combined_commands = "";
    for (int num_process = MIN_P; num_process <= MAX_P; num_process *= P_GROW_FACTOR)
    {
        // Runs through all thread combinations
        for (int num_thread = MIN_T; num_thread <= MAX_T; num_thread *= T_GROW_FACTOR)
        {
            if(num_process * num_thread >= MAX_INSTANCES)
                break;

            combined_commands += beginning + std::to_string(num_process) + middle + std::to_string(num_thread) + end + CONNECT_LOGIC_STR;
        }
    }

    // Outputs final commands 
    std::ofstream command_file;
    command_file.open(FASTA_FILE);
    command_file << "#! /usr/bin/bash\n";
    command_file << combined_commands.c_str();
    command_file.close();
}

// --------------------------------------------------------------------------------------------
// main()
// --------------------------------------------------------------------------------------------
int main(int argc, char ** argv)
{

    // Generates all commands to test functions
    testCompression();
    testDecompression();
    testFindFirst();
    testGrep();
    testFindAll();
    testToFASTA();
    
    // Ends program
    return 0;
}


