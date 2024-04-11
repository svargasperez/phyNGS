/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.

  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
  Incompresso Implementation Author: Cole Koryto
*/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "defs.h"
#include "structures.h"
#include "bit_stream.h"
#include "tasks.h"
#include "incompresso.h"
#include <sys/stat.h>

// --------------------------------------------------------------------------------------------
bool FindFirst(const char *in_file, int32 g_size, int32 p_rank, char *pat, bool to_print, int32 no_threads)
{
  bool local_pattern_found = false;
  bool global_pattern_found = false;
  bool local_search_done = false;
  bool stop_communication = false;
  Record found_rec;

  MPI_File input_NGSC;
  Footer footer;
  BlockHeader block_header;
  std::vector<SubBlock> p_subblocks;
  std::vector<int> block_pos;
  MPI_Offset fileSize;
  uint32 r_buffer_size = 8 << 15;
  uchar *read_buffer;
  MPI_Info Lustre_info;

  BitStream footer_bit_stream;
  double p_timer_start, p_timer_end;

  MPI_Info_create(&Lustre_info);
  // Number of OST ("disks") per file
  MPI_Info_set(Lustre_info, "striping_factor", "64");
  // Set the striping unit to 8MiB
  MPI_Info_set(Lustre_info, "striping_unit", "8388608");

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_open(MPI_COMM_WORLD, in_file, MPI_MODE_RDONLY, Lustre_info, &input_NGSC);

  MPI_Info_free(&Lustre_info);

  if (p_rank == 0)
    printf("[I] INFO: Processing NGSC footer.\n\n");

  p_timer_start = MPI_Wtime();
  read_buffer = (uchar *)malloc(r_buffer_size * sizeof(uchar));
  MPI_File_get_size(input_NGSC, &fileSize);

  // try to use read at all or something collective
  MPI_File_read_at(input_NGSC, fileSize - r_buffer_size, read_buffer, r_buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);
  footer_bit_stream.Create(1);

  uint32 footer_len = read_buffer[r_buffer_size - 2];
  footer_len = (footer_len << 8) + read_buffer[r_buffer_size - 1];

  uint64 pos = fileSize - footer_len - 2;

  if (footer_len + 2 > r_buffer_size)
  {
    free(read_buffer);
    read_buffer = (uchar *)malloc((footer_len + 1) * sizeof(uchar));
    MPI_File_read_at(input_NGSC, pos, read_buffer, footer_len, MPI_CHAR, MPI_STATUS_IGNORE);
    r_buffer_size = footer_len;
    read_buffer[footer_len] = '\0';
  }

  footer_bit_stream.SetIO_Buffer(read_buffer, r_buffer_size, r_buffer_size - footer_len - 2);
  free(read_buffer);

  ReadFooter(footer_bit_stream, footer, p_subblocks, input_NGSC, g_size, p_rank);
  footer_bit_stream.Close();

  // Creates lps[] for pattern search that will hold the longest prefix suffix values for pattern
  int32 M = strlen(pat);
  int32 *lps = (int32 *)malloc(M * sizeof(int32));

  // Preprocess the pattern (calculate lps[] array)
  ComputeLPSArray(pat, M, lps);

  uchar untrans_amb_codes[256];
  std::fill_n(untrans_amb_codes + 128 + 0 * 8, 8, 'Y');
  std::fill_n(untrans_amb_codes + 128 + 1 * 8, 8, 'R');
  std::fill_n(untrans_amb_codes + 128 + 2 * 8, 8, 'W');
  std::fill_n(untrans_amb_codes + 128 + 3 * 8, 8, 'S');
  std::fill_n(untrans_amb_codes + 128 + 4 * 8, 8, 'K');
  std::fill_n(untrans_amb_codes + 128 + 5 * 8, 8, 'M');
  std::fill_n(untrans_amb_codes + 128 + 6 * 8, 8, 'D');
  std::fill_n(untrans_amb_codes + 128 + 7 * 8, 8, 'V');
  std::fill_n(untrans_amb_codes + 128 + 8 * 8, 8, 'H');
  std::fill_n(untrans_amb_codes + 128 + 9 * 8, 8, 'B');
  std::fill_n(untrans_amb_codes + 128 + 10 * 8, 8, 'N');
  std::fill_n(untrans_amb_codes + 128 + 11 * 8, 8, 'X');
  std::fill_n(untrans_amb_codes + 128 + 12 * 8, 8, 'U');
  std::fill_n(untrans_amb_codes + 128 + 13 * 8, 8, '.');
  std::fill_n(untrans_amb_codes + 128 + 14 * 8, 8, '-');

  // Makes sure that all processes enter the subblock loop so that blocking communication can be used
  bool not_looping = false;
  if (p_subblocks.begin() == p_subblocks.end())
    not_looping = true;
  MPI_Allreduce(&not_looping, &stop_communication, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

  // Make calculations for initial position of the decompressed block in FASTQ file
  bool initial_FASTQ_pos = true;
  for (std::vector<SubBlock>::iterator it = p_subblocks.begin(); it != p_subblocks.end(); ++it)
  {
    SubBlock sb = *it;
    BitStream sb_header_stream, sb_title_stream, sb_qua_stream;
    std::vector<Field> fields;
    std::vector<uint32> dna_occ(256, 0);
    uint32 fastq_flags;
    std::vector<uchar> symbols;
    std::vector<uchar> qualities;
    std::vector<uint32> no_ambiguity;
    uint32 no_records = 0;
    uint32 no_symbols;
    uint32 no_qualities;
    uint32 quality_stats_mode;
    uint32 max_quality_length;
    uint32 global_max_sequence_length;
    int64 prev_title_len = 0, prev_qua_len = 0;

    sb_header_stream.Create(p_rank);
    sb_title_stream.Create(p_rank);
    sb_qua_stream.Create(p_rank);

    uint32 total_sb_len = sb.is_splitted ? (sb.sb_length + sb.sb_cont_length) : sb.sb_length;

    read_buffer = (uchar *)malloc(total_sb_len * sizeof(uchar));

    if (!read_buffer)
    {
      printf("\n[E] ERROR: p_Rank %d can't allocate memory for the read_buffer.\n", p_rank);
      printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
      MPI_Finalize();
      exit(3);
    }

    MPI_File_read_at(input_NGSC, sb.sb_start_pos, read_buffer, total_sb_len, MPI_CHAR, MPI_STATUS_IGNORE);

    if (sb.is_splitted)
    {
      MPI_File_read_at(input_NGSC, sb.sb_cont_pos, read_buffer + sb.sb_length, sb.sb_cont_length, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    sb_header_stream.SetIO_Buffer(read_buffer, (total_sb_len / 2));

    sb_header_stream.GetWord(no_records);
    sb_header_stream.GetWord(max_quality_length);
    sb_header_stream.GetWord(global_max_sequence_length);
    sb_header_stream.GetByte(no_symbols);
    my_assert(no_symbols <= 256);

    sb_header_stream.GetByte(quality_stats_mode);
    my_assert(quality_stats_mode < 4);

    sb_header_stream.GetByte(no_qualities);
    my_assert(no_qualities <= 256);

    sb_header_stream.GetWord(fastq_flags);
    my_assert(fastq_flags < 1 << 8);

    sb_header_stream.FlushInputWordBuffer();

    int32 quality_len_bits = BitStream::BitLength(max_quality_length);
    Record *rec = new Record[no_records];

    if ((fastq_flags & FLAG_VARIABLE_LENGTH) != 0)
    {
      uint32 tmp;
      for (uint32 i = 0; i < no_records; ++i)
      {
        sb_header_stream.GetBits(tmp, quality_len_bits);
        rec[i].qua_len = tmp;
        rec[i].seq_len = tmp;
        rec[i].dna_seq = new uchar[tmp + 1];
        rec[i].quality = new uchar[tmp + 1];
      }
      sb_header_stream.FlushInputWordBuffer();
    }
    else
    {
      for (uint32 i = 0; i < no_records; ++i)
      {
        rec[i].qua_len = max_quality_length;
        rec[i].seq_len = global_max_sequence_length;
        rec[i].dna_seq = new uchar[global_max_sequence_length + 1];
        rec[i].quality = new uchar[max_quality_length + 1];
      }
    }

    uint64 sb_offset_pos;
    sb_header_stream.GetDWord(sb_offset_pos);

    if (initial_FASTQ_pos)
    {
      initial_FASTQ_pos = false;
    }

    uint32 title_ln = 0, qua_ln = 0;
    sb_header_stream.GetWord(title_ln);
    sb_header_stream.GetWord(qua_ln);

    int32 curr_pos = sb_header_stream.GetIO_Buffer_Pos();
    sb_header_stream.Close();

    title_ln += curr_pos;
    qua_ln += title_ln;

    sb_title_stream.SetIO_Buffer(read_buffer, title_ln, curr_pos);
    sb_qua_stream.SetIO_Buffer(read_buffer, total_sb_len, title_ln);
    // This is for when you separte DNA and Quality
    // sb_qua_stream.SetIO_Buffer(read_buffer, qua_ln, title_ln);

    free(read_buffer);

    // OpenMP Parallel Region: Threads read and decompress records in subblock.
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel sections num_threads(no_threads)
    {
      // Fetch tiltes from the subblock
      #pragma omp section
      {
        FetchTitleHeader(sb_title_stream, fields);
        FetchTitleBody(sb_title_stream, fields, rec, no_records, fastq_flags, prev_title_len);
        sb_title_stream.Close();
      }

      // Fetch DNA sequences and Qualities from the subblock
      #pragma omp section
      {
        FetchQuality(sb_qua_stream, qualities, no_qualities, rec, no_records, fastq_flags, no_ambiguity, max_quality_length, prev_qua_len);
        FetchDNA(sb_qua_stream, symbols, no_symbols, fastq_flags, rec, no_records, no_ambiguity);
        sb_qua_stream.Close();
      }
    }

    #pragma omp parallel num_threads(no_threads)
    {
      // Remove ambiguity and calculate positions
      #pragma omp for
      for (uint32 i = 0; i < no_records; ++i)
      {
        if (rec[i].seq_len == rec[i].qua_len)
          continue;

        int32 k = rec[i].seq_len;
        int32 j = rec[i].qua_len;

        for (; j >= 0; --j)
        {
          if (rec[i].quality[j] >= 128)
          {
            rec[i].dna_seq[j] = untrans_amb_codes[rec[i].quality[j]];
            rec[i].quality[j] = 33 + (rec[i].quality[j] & 7);
          }
          else
          {
            rec[i].dna_seq[j] = rec[i].dna_seq[k--];
          }
        }
        rec[i].seq_len = rec[i].qua_len;
      }
    } 

    // Searches for pattern
    #pragma omp parallel for num_threads(no_threads) reduction(||:local_pattern_found)
    for (uint32 i = 0; i < no_records; ++i)
    {
      // This is needed because random data was being included after seq end
      char * trimmed_seq = (char *)malloc(rec[i].seq_len+1);
      strncpy(trimmed_seq, (char *)rec[i].dna_seq, rec[i].seq_len);
      trimmed_seq[rec[i].seq_len] = '\0';

      // If the pattern is found then set variable to true
      if (KMPSearch(pat, trimmed_seq, lps) != -1)
      {
        if (to_print && !local_pattern_found)
          #pragma omp critical(matched)
            printf("p%d:%.*sp%d:%.*s\np%d:%.*s\n", p_rank, rec[i].title_len, rec[i].title, p_rank, rec[i].seq_len, rec[i].dna_seq, p_rank, rec[i].qua_len, rec[i].quality);
        local_pattern_found = true;
      }
    }

    // Makes sure that all processes still have more subblocks to loop through before using blocking Allreduce function
    if (!stop_communication)
    {
      if (it == p_subblocks.end()-1)
        local_search_done = true;
      MPI_Allreduce(&local_search_done, &stop_communication, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

      // Uses all reduce if possible to determine if any process has found match, if so, end program
      MPI_Allreduce(&local_pattern_found, &global_pattern_found, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
      if (global_pattern_found)
      {
        delete[] rec;
        no_records = 0;
        break;
      }
    }

    delete[] rec;
    no_records = 0;
  }

  // If communication was stopped before the loop, then check if any processes in loop found a match
  if(stop_communication)
    MPI_Allreduce(&local_pattern_found, &global_pattern_found, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

  // Stop timer
  p_timer_end = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);

  if (p_rank == 0)
  {
    if (global_pattern_found)
      printf("Sequence \"%s\" found in file.\n", pat);
    else
      printf("Sequence not found in file.\n");
  }

  if (p_rank == 0)
    printf("\nRANK\tDECO_TIME\n----------------------------------------------\n");

  MPI_Barrier(MPI_COMM_WORLD);
  printf("%03d\t%f\t%ld\n", p_rank, p_timer_end - p_timer_start, p_subblocks.size());

  // TODO Temp for performance testing
  if (p_rank == 0)
  {
    std::string folder_name = "./performanceOutput-Findfirst";

    // makes output directory if needed
    if (mkdir(folder_name.c_str(), 0777) != -1)
      printf("Performance output directory made.\n");

    std::string file_path = "./" + folder_name + "/" + std::to_string(g_size) + "," + std::to_string(no_threads) + "," + pat + ".txt";
    std::ofstream test_file;
    test_file.open(file_path);
    test_file << p_timer_end-p_timer_start;
  }

  MPI_File_close(&input_NGSC);
  return global_pattern_found;
}


//--------------------------------------------------------------------------------------------
void FindAll(const char *in_file, int32 g_size, int32 p_rank, char *pat, bool to_print, int32 no_threads)
{
  int32 global_num_matches = 0;
  std::string found_out_buffer;
  int32 local_num_matches = 0;

  MPI_File input_NGSC;
  Footer footer;
  BlockHeader block_header;
  std::vector<SubBlock> p_subblocks;
  std::vector<int> block_pos;
  MPI_Offset fileSize;
  uint32 r_buffer_size = 8 << 15;
  uchar *read_buffer;
  MPI_Info Lustre_info;

  BitStream footer_bit_stream;
  double p_timer_start, p_timer_end;

  MPI_Info_create(&Lustre_info);
  // Number of OST ("disks") per file
  MPI_Info_set(Lustre_info, "striping_factor", "64");
  // Set the striping unit to 8MiB
  MPI_Info_set(Lustre_info, "striping_unit", "8388608");

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_open(MPI_COMM_WORLD, in_file, MPI_MODE_RDONLY, Lustre_info, &input_NGSC);

  MPI_Info_free(&Lustre_info);

  if (p_rank == 0)
    printf("[I] INFO: Processing NGSC footer.\n\n");

  p_timer_start = MPI_Wtime();
  read_buffer = (uchar *)malloc(r_buffer_size * sizeof(uchar));
  MPI_File_get_size(input_NGSC, &fileSize);

  // try to use read at all or something collective
  MPI_File_read_at(input_NGSC, fileSize - r_buffer_size, read_buffer, r_buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);
  footer_bit_stream.Create(1);

  uint32 footer_len = read_buffer[r_buffer_size - 2];
  footer_len = (footer_len << 8) + read_buffer[r_buffer_size - 1];

  uint64 pos = fileSize - footer_len - 2;

  if (footer_len + 2 > r_buffer_size)
  {
    free(read_buffer);
    read_buffer = (uchar *)malloc((footer_len + 1) * sizeof(uchar));
    MPI_File_read_at(input_NGSC, pos, read_buffer, footer_len, MPI_CHAR, MPI_STATUS_IGNORE);
    r_buffer_size = footer_len;
    read_buffer[footer_len] = '\0';
  }

  footer_bit_stream.SetIO_Buffer(read_buffer, r_buffer_size, r_buffer_size - footer_len - 2);
  free(read_buffer);

  ReadFooter(footer_bit_stream, footer, p_subblocks, input_NGSC, g_size, p_rank);
  footer_bit_stream.Close();

  // Creates lps[] for pattern search that will hold the longest prefix suffix values for pattern
  int32 M = strlen(pat);
  int32 *lps = (int32 *)malloc(M * sizeof(int32));

  // Preprocess the pattern (calculate lps[] array)
  ComputeLPSArray(pat, M, lps);

  uchar untrans_amb_codes[256];
  std::fill_n(untrans_amb_codes + 128 + 0 * 8, 8, 'Y');
  std::fill_n(untrans_amb_codes + 128 + 1 * 8, 8, 'R');
  std::fill_n(untrans_amb_codes + 128 + 2 * 8, 8, 'W');
  std::fill_n(untrans_amb_codes + 128 + 3 * 8, 8, 'S');
  std::fill_n(untrans_amb_codes + 128 + 4 * 8, 8, 'K');
  std::fill_n(untrans_amb_codes + 128 + 5 * 8, 8, 'M');
  std::fill_n(untrans_amb_codes + 128 + 6 * 8, 8, 'D');
  std::fill_n(untrans_amb_codes + 128 + 7 * 8, 8, 'V');
  std::fill_n(untrans_amb_codes + 128 + 8 * 8, 8, 'H');
  std::fill_n(untrans_amb_codes + 128 + 9 * 8, 8, 'B');
  std::fill_n(untrans_amb_codes + 128 + 10 * 8, 8, 'N');
  std::fill_n(untrans_amb_codes + 128 + 11 * 8, 8, 'X');
  std::fill_n(untrans_amb_codes + 128 + 12 * 8, 8, 'U');
  std::fill_n(untrans_amb_codes + 128 + 13 * 8, 8, '.');
  std::fill_n(untrans_amb_codes + 128 + 14 * 8, 8, '-');

  // Make calculations for initial position of the decompressed block in FASTQ file
  bool initial_FASTQ_pos = true;
  for (std::vector<SubBlock>::iterator it = p_subblocks.begin(); it != p_subblocks.end(); ++it)
  {
    SubBlock sb = *it;
    BitStream sb_header_stream, sb_title_stream, sb_qua_stream;
    std::vector<Field> fields;
    std::vector<uint32> dna_occ(256, 0);
    uint32 fastq_flags;
    std::vector<uchar> symbols;
    std::vector<uchar> qualities;
    std::vector<uint32> no_ambiguity;
    uint32 no_records = 0;
    uint32 no_symbols;
    uint32 no_qualities;
    uint32 quality_stats_mode;
    uint32 max_quality_length;
    uint32 global_max_sequence_length;
    int64 prev_title_len = 0, prev_qua_len = 0;

    sb_header_stream.Create(p_rank);
    sb_title_stream.Create(p_rank);
    sb_qua_stream.Create(p_rank);

    uint32 total_sb_len = sb.is_splitted ? (sb.sb_length + sb.sb_cont_length) : sb.sb_length;

    read_buffer = (uchar *)malloc(total_sb_len * sizeof(uchar));

    if (!read_buffer)
    {
      printf("\n[E] ERROR: p_Rank %d can't allocate memory for the read_buffer.\n", p_rank);
      printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
      MPI_Finalize();
      exit(3);
    }

    MPI_File_read_at(input_NGSC, sb.sb_start_pos, read_buffer, total_sb_len, MPI_CHAR, MPI_STATUS_IGNORE);

    if (sb.is_splitted)
    {
      MPI_File_read_at(input_NGSC, sb.sb_cont_pos, read_buffer + sb.sb_length, sb.sb_cont_length, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    sb_header_stream.SetIO_Buffer(read_buffer, (total_sb_len / 2));

    sb_header_stream.GetWord(no_records);
    sb_header_stream.GetWord(max_quality_length);
    sb_header_stream.GetWord(global_max_sequence_length);
    sb_header_stream.GetByte(no_symbols);
    my_assert(no_symbols <= 256);

    sb_header_stream.GetByte(quality_stats_mode);
    my_assert(quality_stats_mode < 4);

    sb_header_stream.GetByte(no_qualities);
    my_assert(no_qualities <= 256);

    sb_header_stream.GetWord(fastq_flags);
    my_assert(fastq_flags < 1 << 8);

    sb_header_stream.FlushInputWordBuffer();

    int32 quality_len_bits = BitStream::BitLength(max_quality_length);
    Record *rec = new Record[no_records];

    if ((fastq_flags & FLAG_VARIABLE_LENGTH) != 0)
    {
      uint32 tmp;
      for (uint32 i = 0; i < no_records; ++i)
      {
        sb_header_stream.GetBits(tmp, quality_len_bits);
        rec[i].qua_len = tmp;
        rec[i].seq_len = tmp;
        rec[i].dna_seq = new uchar[tmp + 1];
        rec[i].quality = new uchar[tmp + 1];
      }
      sb_header_stream.FlushInputWordBuffer();
    }
    else
    {
      for (uint32 i = 0; i < no_records; ++i)
      {
        rec[i].qua_len = max_quality_length;
        rec[i].seq_len = global_max_sequence_length;
        rec[i].dna_seq = new uchar[global_max_sequence_length + 1];
        rec[i].quality = new uchar[max_quality_length + 1];
      }
    }

    uint64 sb_offset_pos;
    sb_header_stream.GetDWord(sb_offset_pos);

    if (initial_FASTQ_pos)
    {
      initial_FASTQ_pos = false;
    }

    uint32 title_ln = 0, qua_ln = 0;
    sb_header_stream.GetWord(title_ln);
    sb_header_stream.GetWord(qua_ln);

    int32 curr_pos = sb_header_stream.GetIO_Buffer_Pos();
    sb_header_stream.Close();

    title_ln += curr_pos;
    qua_ln += title_ln;

    sb_title_stream.SetIO_Buffer(read_buffer, title_ln, curr_pos);
    sb_qua_stream.SetIO_Buffer(read_buffer, total_sb_len, title_ln);
    // This is for when you separte DNA and Quality
    // sb_qua_stream.SetIO_Buffer(read_buffer, qua_ln, title_ln);

    free(read_buffer);

    // OpenMP Parallel Region: Threads read and decompress records in subblock.
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel sections num_threads(no_threads)
    {
      // Fetch tiltes from the subblock
      #pragma omp section
      {
        FetchTitleHeader(sb_title_stream, fields);
        FetchTitleBody(sb_title_stream, fields, rec, no_records, fastq_flags, prev_title_len);
        sb_title_stream.Close();
      }

      // Fetch DNA sequences and Qualities from the subblock
      #pragma omp section
      {
        FetchQuality(sb_qua_stream, qualities, no_qualities, rec, no_records, fastq_flags, no_ambiguity, max_quality_length, prev_qua_len);
        FetchDNA(sb_qua_stream, symbols, no_symbols, fastq_flags, rec, no_records, no_ambiguity);
        sb_qua_stream.Close();
      }
    }

    // Removes ambiguity symbols from quality and adds "N" in place of missing dna letters
    #pragma omp parallel num_threads(no_threads)
    {
      // Remove ambiguity and calculate positions
      #pragma omp for
      for (uint32 i = 0; i < no_records; ++i)
      {
        if (rec[i].seq_len == rec[i].qua_len)
          continue;

        int32 k = rec[i].seq_len;
        int32 j = rec[i].qua_len;

        for (; j >= 0; --j)
        {
          if (rec[i].quality[j] >= 128)
          {
            rec[i].dna_seq[j] = untrans_amb_codes[rec[i].quality[j]];
            rec[i].quality[j] = 33 + (rec[i].quality[j] & 7);
          }
          else
          {
            rec[i].dna_seq[j] = rec[i].dna_seq[k--];
          }
        }
        rec[i].seq_len = rec[i].qua_len;
      }
    }

    // Searches for pattern
    #pragma omp parallel for num_threads(no_threads) reduction(+:local_num_matches)
    for (uint32 i = 0; i < no_records; ++i)
    {
      // TODO this is needed because a bunch of random data was being included after seq end
      char * trimmed_seq = (char *)malloc(rec[i].seq_len+1);
      strncpy(trimmed_seq, (char *)rec[i].dna_seq, rec[i].seq_len);
      trimmed_seq[rec[i].seq_len] = '\0';
      
      if (KMPSearch(pat, trimmed_seq, lps) != -1)
      {

        // Prints out record match was found
        if (to_print)
          #pragma omp critical(all_print)
            printf("%.*s%.*s\n%.*s\n", rec[i].title_len, rec[i].title, rec[i].seq_len, rec[i].dna_seq, rec[i].qua_len, rec[i].quality);

        // Adds to match counter 
        local_num_matches++;
      }
    }   

    delete[] rec;
    no_records = 0;
  }

  // Uses reduce to find total number of matches
  MPI_Reduce(&local_num_matches, &global_num_matches, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  // Stop timer
  p_timer_end = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (p_rank == 0)
    printf("Number of sequence matches in file: %d.\n", global_num_matches);
    
  if (p_rank == 0)
    printf("\nRANK\tDECO_TIME\n----------------------------------------------\n");

  MPI_Barrier(MPI_COMM_WORLD);
  printf("%03d\t%f\t%ld\n", p_rank, p_timer_end - p_timer_start, p_subblocks.size());


  // TODO Temp for performance testing
  if (p_rank == 0)
  {
    std::string folder_name = "./performanceOutput-Findall";

    // makes output directory if needed
    if (mkdir(folder_name.c_str(), 0777) != -1)
      printf("Performance output directory made.\n");

    std::string file_path = "./" + folder_name + "/" + std::to_string(g_size) + "," + std::to_string(no_threads) + "," + pat + ".txt";
    std::ofstream test_file;
    test_file.open(file_path);
    test_file << p_timer_end-p_timer_start;
  }


  MPI_File_close(&input_NGSC);
}

// --------------------------------------------------------------------------------------------
// Prints occurrences of seq[] in pat[]
int32 KMPSearch(char *pat, char *seq, int32 *lps)
{
  int32 M = strlen(pat);
  int32 N = strlen(seq);

  int32 i = 0; // index for seq[]
  int32 j = 0; // index for pat[]
  while ((N - i) >= (M - j))
  {

    // Checks if pattern matches sequence 
    if (pat[j] == seq[i])
    {
      j++;
      i++;
    }

    // Returns when end of pattern match complete
    if (j == M)
    {
      return (i - j);
    }

    // Handles mismatch after j matches
    else if (i < N && pat[j] != seq[i])
    {
      // Do not match lps[0..lps[j-1]] characters since they will match anyway
      if (j != 0)
        j = lps[j - 1];
      else
        i = i + 1;
    }
  }

  // Pattern not found
  return (-1);
}

// --------------------------------------------------------------------------------------------
// Fills lps[] for given pattern pat[0..M-1]
void ComputeLPSArray(char *pat, int32 M, int32 *lps)
{
  int32 len = 0;
  lps[0] = 0; 

  // Calculates lps[i] for all remaining i
  int32 i = 1;
  while (i < M)
  {
    if (pat[i] == pat[len])
    {
      len++;
      lps[i] = len;
      i++;
    }
    else
    {
      if (len != 0)
      {
        len = lps[len - 1];
      }
      else
      {
        lps[i] = 0;
        i++;
      }
    }
  }
}

// --------------------------------------------------------------------------------------------
void ToFASTA(const char *in_file, const char *out_file, int32 g_size, int32 p_rank, int32 no_threads)
{
  MPI_File input_NGSC, output_FASTA;
  Footer footer;
  BlockHeader block_header;
  std::vector<SubBlock> p_subblocks;
  std::vector<int> block_pos;
  MPI_Offset fileSize;
  uint32 w_buffer_size, r_buffer_size = 8 << 15;
  uchar *read_buffer, *write_buffer = NULL;
  MPI_Offset p_working_region = 0, p_curr_start_pos = 0;
  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Info Lustre_info;

  BitStream footer_bit_stream;
  double p_timer_start, p_timer_end;

  MPI_Info_create(&Lustre_info);
  // Number of OST ("disks") per file
  MPI_Info_set(Lustre_info, "striping_factor", "64");
  // Set the striping unit to 8MiB
  MPI_Info_set(Lustre_info, "striping_unit", "8388608");

  // Deletes output file if it exits
  if (p_rank == 0)
      MPI_File_delete(out_file, MPI_INFO_NULL);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_open(MPI_COMM_WORLD, in_file, MPI_MODE_RDONLY, Lustre_info, &input_NGSC);
  MPI_File_open(MPI_COMM_WORLD, out_file, MPI_MODE_CREATE | MPI_MODE_WRONLY, Lustre_info, &output_FASTA);

  MPI_Info_free(&Lustre_info);

  if (p_rank == 0)
    printf("[I] INFO: Processing NGSC footer.\n\n");

  p_timer_start = MPI_Wtime();
  read_buffer = (uchar *)malloc(r_buffer_size * sizeof(uchar));
  MPI_File_get_size(input_NGSC, &fileSize);

  // try to use read at all or something collective
  MPI_File_read_at(input_NGSC, fileSize - r_buffer_size, read_buffer, r_buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);
  footer_bit_stream.Create(1);

  uint32 footer_len = read_buffer[r_buffer_size - 2];
  footer_len = (footer_len << 8) + read_buffer[r_buffer_size - 1];

  uint64 pos = fileSize - footer_len - 2;

  if (footer_len + 2 > r_buffer_size)
  {
    free(read_buffer);
    read_buffer = (uchar *)malloc((footer_len + 1) * sizeof(uchar));
    MPI_File_read_at(input_NGSC, pos, read_buffer, footer_len, MPI_CHAR, MPI_STATUS_IGNORE);
    r_buffer_size = footer_len;
    read_buffer[footer_len] = '\0';
  }

  footer_bit_stream.SetIO_Buffer(read_buffer, r_buffer_size, r_buffer_size - footer_len - 2);
  free(read_buffer);

  ReadFooter(footer_bit_stream, footer, p_subblocks, input_NGSC, g_size, p_rank);
  footer_bit_stream.Close();

  // Start Decompression
  // --------------------------------------------------------------------------------------------
  if (p_subblocks.size())
  {
    int32 wr_id = p_subblocks.at(0).wr_id;
    p_working_region = footer.FS / footer.PS;
    p_curr_start_pos = (wr_id * p_working_region);
  }

  uchar untrans_amb_codes[256];
  std::fill_n(untrans_amb_codes + 128 + 0 * 8, 8, 'Y');
  std::fill_n(untrans_amb_codes + 128 + 1 * 8, 8, 'R');
  std::fill_n(untrans_amb_codes + 128 + 2 * 8, 8, 'W');
  std::fill_n(untrans_amb_codes + 128 + 3 * 8, 8, 'S');
  std::fill_n(untrans_amb_codes + 128 + 4 * 8, 8, 'K');
  std::fill_n(untrans_amb_codes + 128 + 5 * 8, 8, 'M');
  std::fill_n(untrans_amb_codes + 128 + 6 * 8, 8, 'D');
  std::fill_n(untrans_amb_codes + 128 + 7 * 8, 8, 'V');
  std::fill_n(untrans_amb_codes + 128 + 8 * 8, 8, 'H');
  std::fill_n(untrans_amb_codes + 128 + 9 * 8, 8, 'B');
  std::fill_n(untrans_amb_codes + 128 + 10 * 8, 8, 'N');
  std::fill_n(untrans_amb_codes + 128 + 11 * 8, 8, 'X');
  std::fill_n(untrans_amb_codes + 128 + 12 * 8, 8, 'U');
  std::fill_n(untrans_amb_codes + 128 + 13 * 8, 8, '.');
  std::fill_n(untrans_amb_codes + 128 + 14 * 8, 8, '-');

  // Make calculations for initial position of the decompressed block in FASTA file
  bool initial_FASTA_pos = true;

  for (std::vector<SubBlock>::iterator it = p_subblocks.begin(); it != p_subblocks.end(); ++it)
  {
    SubBlock sb = *it;
    BitStream sb_header_stream, sb_title_stream, sb_qua_stream;
    std::vector<Field> fields;
    std::vector<uint32> dna_occ(256, 0);
    uint32 fasta_flags;
    std::vector<uchar> symbols;
    std::vector<uchar> qualities;
    std::vector<uint32> no_ambiguity;
    uint32 no_records = 0;
    uint32 no_symbols;
    uint32 no_qualities;
    uint32 quality_stats_mode;
    uint32 max_quality_length;
    uint32 global_max_sequence_length;
    int64 prev_title_len = 0, prev_qua_len = 0, prev_seq_len = 0;

    sb_header_stream.Create(p_rank);
    sb_title_stream.Create(p_rank);
    sb_qua_stream.Create(p_rank);

    uint32 total_sb_len = sb.is_splitted ? (sb.sb_length + sb.sb_cont_length) : sb.sb_length;

    read_buffer = (uchar *)malloc(total_sb_len * sizeof(uchar));

    if (!read_buffer)
    {
      printf("\n[E] ERROR: p_Rank %d can't allocate memory for the read_buffer.\n", p_rank);
      printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
      MPI_Finalize();
      exit(3);
    }

    MPI_File_read_at(input_NGSC, sb.sb_start_pos, read_buffer, total_sb_len, MPI_CHAR, MPI_STATUS_IGNORE);

    if (sb.is_splitted)
    {
      MPI_File_read_at(input_NGSC, sb.sb_cont_pos, read_buffer + sb.sb_length, sb.sb_cont_length, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    sb_header_stream.SetIO_Buffer(read_buffer, (total_sb_len / 2));

    sb_header_stream.GetWord(no_records);
    sb_header_stream.GetWord(max_quality_length);
    sb_header_stream.GetWord(global_max_sequence_length);
    sb_header_stream.GetByte(no_symbols);
    my_assert(no_symbols <= 256);

    sb_header_stream.GetByte(quality_stats_mode);
    my_assert(quality_stats_mode < 4);

    sb_header_stream.GetByte(no_qualities);
    my_assert(no_qualities <= 256);

    sb_header_stream.GetWord(fasta_flags);
    my_assert(fasta_flags < 1 << 8);

    sb_header_stream.FlushInputWordBuffer();

    int32 quality_len_bits = BitStream::BitLength(max_quality_length);
    Record *rec = new Record[no_records];

    if ((fasta_flags & FLAG_VARIABLE_LENGTH) != 0)
    {
      uint32 tmp;
      for (uint32 i = 0; i < no_records; ++i)
      {
        sb_header_stream.GetBits(tmp, quality_len_bits);
        rec[i].qua_len = tmp;
        rec[i].seq_len = tmp;
        rec[i].dna_seq = new uchar[tmp + 1];
        rec[i].quality = new uchar[tmp + 1];
      }
      sb_header_stream.FlushInputWordBuffer();
    }
    else
    {
      for (uint32 i = 0; i < no_records; ++i)
      {
        rec[i].qua_len = max_quality_length;
        rec[i].seq_len = global_max_sequence_length;
        rec[i].dna_seq = new uchar[global_max_sequence_length + 1];
        rec[i].quality = new uchar[max_quality_length + 1];
      }
    }

    uint64 sb_offset_pos;
    sb_header_stream.GetDWord(sb_offset_pos);

    if (initial_FASTA_pos)
    {
      p_curr_start_pos = sb_offset_pos;
      initial_FASTA_pos = false;
    }

    uint32 title_ln = 0, qua_ln = 0;
    sb_header_stream.GetWord(title_ln);
    sb_header_stream.GetWord(qua_ln);

    int32 curr_pos = sb_header_stream.GetIO_Buffer_Pos();
    sb_header_stream.Close();

    title_ln += curr_pos;
    qua_ln += title_ln;

    sb_title_stream.SetIO_Buffer(read_buffer, title_ln, curr_pos);
    sb_qua_stream.SetIO_Buffer(read_buffer, total_sb_len, title_ln);
    // This is for when you separte DNA and Quality
    // sb_qua_stream.SetIO_Buffer(read_buffer, qua_ln, title_ln);

    free(read_buffer);

    // OpenMP Parallel Region: Threads read and decompress records in subblock.
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel sections num_threads(no_threads)
    {

      // Fetch titles from the subblock
      #pragma omp section
      {
        FetchTitleHeader(sb_title_stream, fields);
        FetchTitleBody(sb_title_stream, fields, rec, no_records, fasta_flags, prev_title_len);
        sb_title_stream.Close();
      }
      
      // Fetch DNA sequences and Qualities from the subblock
      #pragma omp section
      {
        FetchQuality(sb_qua_stream, qualities, no_qualities, rec, no_records, fasta_flags, no_ambiguity, max_quality_length, prev_qua_len);
        FetchDNA(sb_qua_stream, symbols, no_symbols, fasta_flags, rec, no_records, no_ambiguity);
        sb_qua_stream.Close();
      }
    }
    

    #pragma omp parallel num_threads(no_threads)
    {
    // Remove ambiguity and calculate positions
    #pragma omp for
      for (uint32 i = 0; i < no_records; ++i)
      {
        if (rec[i].seq_len == rec[i].qua_len)
          continue;

        int32 k = rec[i].seq_len;
        int32 j = rec[i].qua_len;

        for (; j >= 0; --j)
        {
          if (rec[i].quality[j] >= 128)
          {
            rec[i].dna_seq[j] = untrans_amb_codes[rec[i].quality[j]];
            rec[i].quality[j] = 33 + (rec[i].quality[j] & 7);
          }
          else
          {
            rec[i].dna_seq[j] = rec[i].dna_seq[k--];
          }
        }
        rec[i].seq_len = rec[i].qua_len;
      }
    }

    // Including mpi_wait to be able to implement non-blocking asynchronous writing
    MPI_Wait(&request, MPI_STATUS_IGNORE);

    if (write_buffer)
      free(write_buffer);

    // Calculates sum of previous sequence lengths plus newline
    int64 new_len = 0;
    for (uint32 i = 0; i < no_records; ++i)
    {
      rec[i].prev_seq_len = prev_seq_len;
      prev_seq_len += rec[i].seq_len + 1;

      // new_len += rec[i].prev_seq_qua_len == 0 ? ;
      // if(prev_seq_len != new_len)
      // {
      // printf("not equal: %lld != %lld, accumulation = %lld\n", prev_seq_len, new_len, rec[i].prev_seq_qua_len);
      // }
    }

    w_buffer_size = prev_title_len + prev_seq_len;
    write_buffer = (uchar *)malloc(w_buffer_size * sizeof(uchar));

    if (!write_buffer)
    {
      printf("[E] ERROR: p_Rank %d can't allocate memory for the write_buffer.\n", p_rank);
      printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
      MPI_Finalize();
      exit(3);
    }

    // Copy rec to the write_buffer
    uchar FASTATitleChar[] = ">";
    int32 writen = 0;
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel default(shared) num_threads(no_threads)
    {
      #pragma omp for
      for (uint32 i = 0; i < no_records; ++i)
      {
        // Adjusts rec_start_pos to fit FASTA format
        int64 rec_start_pos = rec[i].prev_title_len + rec[i].prev_seq_len;

        // Creates length variables
        int32 t_len = rec[i].title_len;
        int32 d_len = rec[i].seq_len + 1; // +1 added since \n is not included in seq_len

        std::copy(FASTATitleChar, FASTATitleChar + 1, write_buffer + rec_start_pos);
        rec_start_pos += 1;
        writen += 1;

        std::copy(rec[i].title + 1, rec[i].title + t_len, write_buffer + rec_start_pos);
        rec_start_pos += t_len - 1; // - 1 here because first char of title was written in previous statement
        writen += t_len - 1 ;

        std::copy(rec[i].dna_seq, rec[i].dna_seq + d_len, write_buffer + rec_start_pos);
        writen += d_len;
      }
    }

    delete[] rec;
    no_records = 0;

    MPI_File_iwrite_at(output_FASTA, p_curr_start_pos, write_buffer, w_buffer_size, MPI_CHAR, &request);
    p_curr_start_pos += w_buffer_size;
    // free(write_buffer); DO NOT FREE THE BUFFER HERE IF USING NON-Blocking
  }

  // Wait for the Writing operation to finish before closing the file
  if (p_subblocks.size())
  {
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    free(write_buffer);
  }

  // Stop timer
  p_timer_end = MPI_Wtime();

  if (p_rank == 0)
    printf("\nRANK\tDECO_TIME\n----------------------------------------------\n");

  MPI_Barrier(MPI_COMM_WORLD);
  printf("%03d\t%f\t%ld\n", p_rank, p_timer_end - p_timer_start, p_subblocks.size());


  // TODO Temp for performance testing
  if (p_rank == 0)
  {
    std::string folder_name = "./performanceOutput-FASTA";

    // makes output directory if needed
    if (mkdir(folder_name.c_str(), 0777) != -1)
      printf("Performance output directory made.\n");

    std::string file_path = "./" + folder_name + "/" + std::to_string(g_size) + "," + std::to_string(no_threads) + ".txt";
    std::ofstream test_file;
    test_file.open(file_path);
    test_file << p_timer_end-p_timer_start;
  }

  MPI_File_close(&input_NGSC);
  MPI_File_close(&output_FASTA);
}