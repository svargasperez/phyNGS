/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
//#include <vector>
//#include <map>
//#include <math.h>
//#include <algorithm>
//#include "defs.h"
//#include "utils.h"
//#include "structures.h"
//#include "bit_stream.h"
//#include "huffman.h"
//#include "tasks.h"
#include "phyNGSC.h"
#include "phyNGSD.h"
#include "incompresso.h"

//TODO
// Remove unneeded includes
// make sure that parameter error output has correct program name (not ./main)
// What other checks need to be added?
// add error handling for incompresso parameters
// test all error types
// spell check all comments
// Test decompression still works after changing parameters
// Does pattern variable need to be ended in \0?

// incompresso
// ngsc file is input
// For find seq, do yes no search, do all match search


//QUESTIONS
// for the neq freq, do we want to handle overflow of the 20 bit number?
// What is up with delete[] rec; memory leak? (see TODOs)


// --------------------------------------------------------------------------------------------
// main()
// --------------------------------------------------------------------------------------------
int main(int argc, char ** argv)
{
    int32 g_size, p_rank, provided;
    int32 p_mode = -1, i_mode = -1;
    char *pattern = NULL;
    bool to_print = false;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);

    // Checks that there is at least one argument provided by user
    if (argc < 2)
    {
        if (p_rank == 0)
            fprintf(stderr, "\n[E] ERROR: Less than one argument provided. See program documentation for command examples.\n");
        MPI_Finalize();
        exit(1);
    }
    
    // Determines program mode
    if (strcmp(argv[1], "-c") == 0)
        p_mode = 0;
    else if (strcmp(argv[1], "-d") == 0)
        p_mode = 1;
    else if (strcmp(argv[1], "-i") == 0)
        p_mode = 2;
    else
    {
        if (p_rank == 0)
        {
            fprintf(stderr, "\n[E] ERROR: Incorrect program mode value. Options:\n");
            fprintf(stderr, "                 -c (or) -d (or) -i:                           Program mode is compress, decompress, or incompresso.\n");
        }
        MPI_Finalize();
        exit(1);
    }

    // Gets the number of threads
    int32 no_threads = atoi(argv[2]);

    // Checks that num_threads parameter is a number
    if (no_threads < 1)
    {
        if (p_rank == 0)
        {
            fprintf(stderr, "\n[E] ERROR: Invalid value for number of threads to be used per MPI process (must be an integer at least 1).\n");
        }
        MPI_Finalize();
        exit(1);
    }

    // Handles errors for compression and decompression parameters
    if (p_mode == 0 || p_mode == 1)
    {
        // Checks that there is the correct number of parameters
        if (argc != 5) 
        {
            if (p_rank == 0)
            {
                fprintf(stderr, "\n[E] ERROR: Incorrect number of arguments. Usage:\n");
                fprintf(stderr, "\n           mpiexec -np p ./main -c num_threads input_filename.fastq output_filename.ngsc.\n");
                fprintf(stderr, "                                                 or\n");
                fprintf(stderr, "           mpiexec -np p ./main -d num_threads input_filename.ngsc output_filename.fastq.\n");
                fprintf(stderr, "                 mpiexec:                              Command to run MPI applications.\n");
                fprintf(stderr, "                 -np p:                                Number of MPI processes to be used, where p is a number greater than 1.\n");
                fprintf(stderr, "                 ./main:                               main application.\n");
                fprintf(stderr, "                 -c (or) -d:                           Command to compress or decompress data.\n");
                fprintf(stderr, "                 num_threads:                          Number of threads to be used per MPI process (must be at least 1).\n");
                fprintf(stderr, "                 input_filename.fastq (or) .ngsc:      Name of FASTQ file for compression or NGSC for decompression.\n");
                fprintf(stderr, "                 output_filename.ngsc (or) .fastq:     Name of resulting NGSC file for compression or FASTQ for decompression (created).\n");
            }
            MPI_Finalize();
            exit(1);
        }
    }

    // Handles errors for incompresso
    if (p_mode == 2)
    {
        // Determines incompresso mode
        if (strcmp(argv[4], "-findfirst") == 0)
            i_mode = 0;
        else if (strcmp(argv[4], "-findall") == 0)
            i_mode = 1;
        else if (strcmp(argv[5], "-fasta") == 0)
            i_mode = 2;
        else if (strcmp(argv[4], "-freqinfo") == 0)
            i_mode = 3;
        else if (strcmp(argv[4], "-trim") == 0)
            i_mode = 4;
        else
        {
            if (p_rank == 0)
            {
                fprintf(stderr, "\n[E] ERROR: Invalid incompresso mode. Options:\n");
                fprintf(stderr, "                 -findfirst (or) -findall (or) -fasta (or) -freqinfo (or) -trim:                           Program mode is find first sequence match, final all sequence matches, convert to fasta, find nucleotide frequency information, or trim sequence.\n");
            }
            MPI_Finalize();
            exit(1);
        }

        // Checks parameters for search modes
        if ((i_mode == 0 || i_mode == 1) && (argc != 6 && argc != 7)) 
        {
            if (p_rank == 0)
            {
                fprintf(stderr, "\n[E] ERROR: Incorrect number of arguments for sequence search. Usage:\n");
                fprintf(stderr, "\n           mpiexec -np p ./main -i num_threads input_filename.ngsc -findfirst [-p] sequence_to_find.\n");
                fprintf(stderr, "                                                 or\n");
                fprintf(stderr, "\n           mpiexec -np p ./main -i num_threads input_filename.ngsc -findall [-p] sequence_to_find.\n");
                fprintf(stderr, "                 mpiexec:                              Command to run MPI applications.\n");
                fprintf(stderr, "                 -np p:                                Number of MPI processes to be used, where p is a number greater than 1.\n");
                fprintf(stderr, "                 ./main:                               main application.\n");
                fprintf(stderr, "                 -i:                                   Program mode is incompresso.\n");                
                fprintf(stderr, "                 num_threads:                          Number of threads to be used per MPI process (must be at least than 1).\n");
                fprintf(stderr, "                 input_filename.ngsc:                  Name of NGSC input file.\n");
                fprintf(stderr, "                 -findfirst (or) -findall:             Incompresso mode is find first or all sequence matches.\n");
                fprintf(stderr, "                 -p: (Optional)                        Option to print records of matches found (Off by default).\n");
                fprintf(stderr, "                 sequence_to_find:                     The sequence for the program to search for.\n");
            }
            MPI_Finalize();
            exit(1);
        }

        // Checks parameters for convert to fasta mode
        if (i_mode == 2 && argc != 6) 
        {
            if (p_rank == 0)
            {
                fprintf(stderr, "\n[E] ERROR: Incorrect number of arguments. Usage:\n");
                fprintf(stderr, "\n           mpiexec -np p ./main -i num_threads input_filename.ngsc output_filename.fasta -fasta.\n");
                fprintf(stderr, "                 mpiexec:                              Command to run MPI applications.\n");
                fprintf(stderr, "                 -np p:                                Number of MPI processes to be used, where p is a number greater than 1.\n");
                fprintf(stderr, "                 ./main:                               main application.\n");
                fprintf(stderr, "                 -i:                                   Program mode is incompresso.\n");                
                fprintf(stderr, "                 num_threads:                          Number of threads to be used per MPI process (must be at least than 1).\n");
                fprintf(stderr, "                 input_filename.ngsc:                  Name of NGSC file for decompression.\n");
                fprintf(stderr, "                 output_filename.fasta:                Name of resulting FASTQ file (created).\n");
                fprintf(stderr, "                 -fasta:                               Incompresso mode is convert to fasta.\n");
            }
            MPI_Finalize();
            exit(1);
        }

        // Checks parameters for convert to nucleotide frequency mode
        if (i_mode == 3 && argc != 5) 
        {
            if (p_rank == 0)
            {
                fprintf(stderr, "\n[E] ERROR: Incorrect number of arguments. Usage:\n");
                fprintf(stderr, "\n           mpiexec -np p ./main -i num_threads input_filename.ngsc -freqinfo.\n");
                fprintf(stderr, "                 mpiexec:                              Command to run MPI applications.\n");
                fprintf(stderr, "                 -np p:                                Number of MPI processes to be used, where p is a number greater than 1.\n");
                fprintf(stderr, "                 ./main:                               main application.\n");
                fprintf(stderr, "                 -i:                                   Program mode is incompresso.\n");                
                fprintf(stderr, "                 num_threads:                          Number of threads to be used per MPI process (must be at least than 1).\n");
                fprintf(stderr, "                 input_filename.ngsc:                  Name of NGSC input file.\n");
                fprintf(stderr, "                 -freqinfo:                            Incompresso mode is find nucleotide frequency information.\n");
            }
            MPI_Finalize();
            exit(1);
        }


        // Checks parameters for sequence trimming mode
        if (i_mode == 4 && (argc != 8 && argc != 9)) 
        {
            if (p_rank == 0)
            {
                fprintf(stderr, "\n[E] ERROR: Incorrect number of arguments. Usage:\n");
                fprintf(stderr, "\n           mpiexec -np p ./main -i num_threads input_filename.ngsc output_filename.ngsc -trim -5 [-p] trim_sequence.\n");
                fprintf(stderr, "                                                 or\n");
                fprintf(stderr, "\n           mpiexec -np p ./main -i num_threads input_filename.ngsc output_filename.ngsc -trim -3 [-p] trim_sequence.\n");
                fprintf(stderr, "                 mpiexec:                              Command to run MPI applications.\n");
                fprintf(stderr, "                 -np p:                                Number of MPI processes to be used, where p is a number greater than 1.\n");
                fprintf(stderr, "                 ./main:                               main application.\n");
                fprintf(stderr, "                 -i:                                   Program mode is incompresso.\n");                
                fprintf(stderr, "                 num_threads:                          Number of threads to be used per MPI process (must be at least than 1).\n");
                fprintf(stderr, "                 input_filename.ngsc:                  Name of NGSC file for decompression.\n");
                fprintf(stderr, "                 output_filename.ngsc:                 Name of resulting NGSC file (created).\n");
                fprintf(stderr, "                 -trim:                                Incompresso mode is trim sequence.\n");
                fprintf(stderr, "                 -5 (or) -3:                           Sequence trimming occurs at 5' end or 3' end.\n");
                fprintf(stderr, "                 -p: (Optional)                        Option to print records of first matches found (Off by default).\n");
                fprintf(stderr, "                 trim_sequence:                        The sequence for the program to search for and trim sequences based on.\n");
            }
            MPI_Finalize();
            exit(1);
        }
    }

    // Gets file input and output names
    std::string in_file = argv[3];
    std::string out_file = argv[4];

    // Finds the input and output file endings
    std::string in_file_end = in_file.substr(in_file.find_last_of(".")+1);
    std::string out_file_end = out_file.substr(out_file.find_last_of(".")+1);

    // Checks parameter file types for compression
    if (p_mode == 0 && (in_file_end.compare("fastq") != 0 || out_file_end.compare("ngsc") != 0))
    {
        if (p_rank == 0)
        {
            fprintf(stderr, "\n[E] ERROR: Incorrect input or output file endings. Correct command:\n");
            fprintf(stderr, "                 mpiexec -np p ./main -c num_threads input_filename.fastq output_filename.ngsc.\n");
        }
        MPI_Finalize();
        exit(1);
    }

    // Checks parameter order for decompression
    if (p_mode == 1 && (in_file_end.compare("ngsc") != 0 || out_file_end.compare("fastq") != 0))
    {
        if (p_rank == 0)
        {
            fprintf(stderr, "\n[E] ERROR: Incorrect input or output file endings. Correct command:\n");
            fprintf(stderr, "                 mpiexec -np p ./main -d num_threads input_filename.ngsc output_filename.fastq.\n");
        }
        MPI_Finalize();
        exit(1);
    }

    // Checks input parameter file type for sequence searches, nucleotide frequency information, and trim modes
    if ((i_mode == 0 || i_mode == 1 || i_mode == 3 || i_mode == 4) && in_file_end.compare("ngsc") != 0 )
    {
        if (p_rank == 0)
            fprintf(stderr, "\n[E] ERROR: Incorrect input file ending. Correct ending is: input_filename.ngsc.\n");
        MPI_Finalize();
        exit(1);
    }

    // Checks parameter order for fasta conversion
    if (i_mode == 2 && (in_file_end.compare("ngsc") != 0 || out_file_end.compare("fasta") != 0))
    {
        if (p_rank == 0)
        {
            fprintf(stderr, "\n[E] ERROR: Incorrect input or output file endings. Correct command:\n");
            fprintf(stderr, "                 mpiexec -np p ./main -i num_threads input_filename.ngsc output_filename.fasta -fasta.\n");
        }
        MPI_Finalize();
        exit(1);
    }

    // Gets pattern to search for from user if incompresso mode is sequence search
    if (p_mode == 2 && (i_mode == 0 || i_mode == 1))
    {
        if (strcmp(argv[5], "-p") == 0)
        {
            pattern = argv[6];
            to_print = true;
        }
        else
            pattern = argv[5];
    }

    // Gets pattern to search and trim for from user if incompresso mode is sequence trim
    if (p_mode == 2 && i_mode == 4)
    {
        if (strcmp(argv[7], "-p") == 0)
        {
            pattern = argv[8];
            to_print = true;
        }
        else
            pattern = argv[7];
    }

    // If hybrid not supported, then set to only one thread
    if (provided < MPI_THREAD_FUNNELED)
    {
        if (p_rank == 0)
            printf("\n[W] WARNING: Hybrid parallelization not supported. Setting number of threads to 1.\n");
        no_threads = 1;
    }

    // Determines if program is compressing or decompressing data
    if (p_mode == 0)
        CompressData(in_file.c_str(), out_file.c_str(), g_size, p_rank, no_threads);
    else if (p_mode == 1)
        DecompressData(in_file.c_str(), out_file.c_str(), g_size, p_rank, no_threads);
    else if (p_mode == 2)
    {
        if (i_mode == 0)
            FindFirst(in_file.c_str(), g_size, p_rank, pattern, to_print, no_threads); 
        else if (i_mode == 1)
            FindAll(in_file.c_str(), g_size, p_rank, pattern, to_print, no_threads); 
        else if (i_mode == 2)
            ToFASTA(in_file.c_str(), out_file.c_str(), g_size, p_rank, no_threads); 
        else if (i_mode == 3)
            printf("Not implemented yet.\n"); 
        else if (i_mode == 4)
            // Trim(in_file.c_str(), out_file.c_str(), g_size, p_rank, pattern, to_print, no_threads); 
            printf("Not implemented yet.\n"); 
    } 
    // TODO when delete[] rec; is in in incompresso, memory error occurs here
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}