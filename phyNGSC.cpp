/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>
#include "defs.h"
#include "utils.h"
#include "structures.h"
#include "bit_stream.h"
#include "huffman.h"
#include "tasks.h"
#include <sys/stat.h>

// --------------------------------------------------------------------------------------------
void CompressData(const char * in_file, const char * out_file, int32 g_size, int32 p_rank, int32 no_threads)
{
  MPI_Status status;
  int32 p_subblock_count = 0;

  MPI_File input_FASTQ, output_NGSC;
  int32 err_read, err_write;

  double p_timer_start, p_timer_end;
  
  MPI_Offset FASTQ_size;

  MPI_Offset p_working_region, r_buffer_size = READ_BUFFER_SIZE, w_buffer_size = WRITE_BUFFER_SIZE;
  MPI_Offset p_bytes_read = 0, p_bytes_written = 2, p_bytes_to_copy = 0;
  MPI_Offset p_wr_start, p_wr_end, r_buffer_curr_pos = 0;
  uchar *read_buffer, *write_buffer, *copy_buffer;
  std::string separators = ". _:/=#,\n";

  int32 overlap = 500;
  int64 no_records = 0, rec_start_pos = 0, rec_end_pos = 0;
  int64 records_per_th = 100000;

  ProcessCompressionInfo p_compress_info, *all_info = NULL;
  std::vector<BlockHeader> p_blocks_data;
  BlockHeader p_block;
  std::vector<double> timestamps;
  MPI_Info Lustre_info;

  MPI_Info_create(&Lustre_info);
  // Number of OST ("disks") per file
  MPI_Info_set(Lustre_info, "striping_factor", "64");
  // Set the striping unit to 8MiB
  MPI_Info_set(Lustre_info, "striping_unit", "8388608");

  // Open input/output files
  err_read  = MPI_File_open(MPI_COMM_WORLD, in_file, MPI_MODE_RDONLY, Lustre_info, &input_FASTQ);
  err_write = MPI_File_open(MPI_COMM_WORLD, out_file, MPI_MODE_CREATE | MPI_MODE_RDWR, Lustre_info, &output_NGSC);

  MPI_Info_free(&Lustre_info);
  
  if (no_threads)
    records_per_th /= no_threads;
  else
  {
    if (p_rank == 0)
      fprintf(stderr, "\n[E] ERROR: Number of threads is less than 1.\n");
    MPI_Finalize();
    exit(1);
  }

  if (g_size < 2)
  {
    if (p_rank == 0)
      fprintf(stderr, "\n[E] ERROR: Number of MPI Processes (-np) is less than 2.\n");
    MPI_Finalize();
    exit(1);
  }

  if (err_read || err_write) {

    if (p_rank == 0)
      fprintf(stderr, "\n[E] ERROR: MPI failed to open file <<%s>>.\n", err_read ? out_file : in_file);
    MPI_Finalize();
    exit(2);
  }

  if(p_rank == 0)
    printf("\n[I] INFO: Processing <<%s>> with %d MPI processes and %d threads per process.\n", in_file, g_size, no_threads);

  // Start timer
  p_timer_start = MPI_Wtime();

  // MPI_Bcast(&p_timer_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_File_get_size(input_FASTQ, &FASTQ_size);

  p_working_region = FASTQ_size/g_size;
  p_wr_start = p_rank * p_working_region;
  p_wr_end = p_wr_start + p_working_region + overlap - 1;

  if (p_rank == g_size - 1)
    p_wr_end = FASTQ_size - 1;

  // If the p_working_region is smaller than 8Mb, then read that amount
  if (p_working_region < r_buffer_size)
    r_buffer_size = p_wr_end - p_wr_start + 1;

  write_buffer = (uchar*) malloc(w_buffer_size * sizeof(uchar));
  read_buffer  = (uchar*) malloc(r_buffer_size * sizeof(uchar));
  MPI_File_read_at(input_FASTQ, p_wr_start, read_buffer, r_buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);

  // Each p_rank (except for 0) locate the start of the first complete record in their working region
  if (p_rank != 0) 
    r_buffer_curr_pos = utils::get_first_record(read_buffer, 0);

  rec_start_pos = r_buffer_curr_pos;
  p_block.BCSS = 0;

  // printf("[%d]{OV::%lld}\n", p_rank, r_buffer_curr_pos);
  // Begin processing FASTQ file with each p_rank's corresponding working region
  // --------------------------------------------------------------------------------------------
  while(p_bytes_read < p_working_region)
  {
    std::vector<Field> fields;
    uint32 n_fields = 0;
    std::vector<uint32> dna_occ(256,0);
    uint32 fastq_flags;
    std::vector<uchar> symbols;
    std::vector<uchar> qualities;
    uchar sym_code[256]={0};
    uchar qua_code[256]={0};
    uint32 *sym_stats = NULL;
    uint32 **quality_stats = NULL;
    uint32 *raw_quality_stats = NULL;
    HuffmanEncoder::Code **qua_huf_codes = NULL;
    HuffmanEncoder::Code *raw_qua_huf_codes = NULL;

    bool last_read = (p_bytes_read + r_buffer_size >= p_working_region) ? 
                     true : false;

    uchar trans_amb_codes[256];
    bool alphabet_qua[256] = {false};

    std::fill_n(trans_amb_codes, 256, 0);
    trans_amb_codes['Y'] = 2;
    trans_amb_codes['R'] = 3;
    trans_amb_codes['W'] = 4;
    trans_amb_codes['S'] = 5;
    trans_amb_codes['K'] = 6;
    trans_amb_codes['M'] = 7;
    trans_amb_codes['D'] = 8;
    trans_amb_codes['V'] = 9;
    trans_amb_codes['H'] = 10;
    trans_amb_codes['B'] = 11;
    trans_amb_codes['N'] = 12;
    trans_amb_codes['X'] = 13;
    trans_amb_codes['U'] = 14;
    trans_amb_codes['.'] = 15;
    trans_amb_codes['-'] = 16;
    trans_amb_codes['A'] = 1;
    trans_amb_codes['C'] = 1;
    trans_amb_codes['T'] = 1;
    trans_amb_codes['G'] = 1;

    // Moved outside while loop, std::string now
    // const char *c_separators = " ._,=:/-#\n";
    // const std::vector<uchar> separators(c_separators, c_separators+strlen(c_separators));

    int64 *tmp_title_pos = new int64[no_threads*records_per_th]();
    int64 *tmp_title_end = new int64[no_threads*records_per_th]();
    int64 *tmp_seq_end   = new int64[no_threads*records_per_th]();
    int64 *th_rec_count  = new int64[no_threads]();

    ++p_subblock_count;

    BitStream p_write_buff_bit_stream;
    BitStream title_bit_stream;
    BitStream dna_bit_stream;
    BitStream quality_bit_stream;

    // Variables to analize quality and dna seq
    uint32 max_quality_length = 0;
    uint32 min_quality_length = (uint32) -1;
    uint32 max_sequence_length = 0;
    uint32 min_sequence_length = (uint32) -1;
    uint32 global_max_sequence_length = 0;

    if(p_bytes_read != 0)
    {
      read_buffer = (uchar*) malloc(r_buffer_size * sizeof(uchar));
      MPI_File_read_at(input_FASTQ, r_buffer_curr_pos, read_buffer, r_buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    if (last_read)
    {
      if (p_rank == g_size - 1)
        rec_end_pos = r_buffer_size - 1;
      else
        rec_end_pos = utils::get_first_record(read_buffer, r_buffer_size - overlap) - 1;
    }
    else
      rec_end_pos = utils::get_last_record(read_buffer, r_buffer_size - 1);
   
    // OpenMP Parallel Region: Threads locate the start of each record in read_buffer
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel num_threads(no_threads) reduction(+ : no_records)
    {
      // Obtain and print thread id
      int32 th_id = omp_get_thread_num();
      int32 i, th_i_start , th_i_end;
      int32 th_rec_no = th_id * records_per_th;

      th_i_start = th_id * (rec_end_pos) / no_threads;
      th_i_end = (th_id + 1) * (rec_end_pos) / no_threads;

      if(th_id == 0)
        th_i_start = rec_start_pos;
      else 
        th_i_start = utils::get_first_record(read_buffer, th_i_start);

      if (th_id == no_threads - 1)
        th_i_end = rec_end_pos;

      for (i = th_i_start; i < th_i_end; ++i)
      {
        if (read_buffer[i] != '@')
          printf("\n[1CK] WARNING: Thread %d in p_Rank %d will not start good in buffer[%d]=%c.SubBlock %d.\n", th_id, p_rank,i,read_buffer[i], p_subblock_count);
     
        tmp_title_pos[th_rec_no] = i;
        
        while(read_buffer[++i] != '\n')
          continue;

        tmp_title_end[th_rec_no] = i;

        while(read_buffer[++i] != '\n')
          continue;

        tmp_seq_end[th_rec_no] = i;

        i +=  (tmp_seq_end[th_rec_no] - tmp_title_end[th_rec_no]) + 2;
        ++th_rec_no;

        // Check that threads keep using their memory space
        if (th_rec_no > ((th_id * records_per_th) + records_per_th))
        {
          printf("\n[!] WARNING: Thread %d in p_Rank %d can't allocate more records.\n", th_id, p_rank);
          printf("             Thread will exit loop with only %lld records\n", (th_rec_no - (th_id * records_per_th)));
          i = th_i_end ;
        }
      }

      th_rec_count[th_id] = th_rec_no - (th_id * records_per_th);
      // Total number of records in the team of threads
      no_records  += th_rec_count[th_id];
    }

    Record *records = new Record[no_records];
    int64 f_start = rec_start_pos;

    // Initializing vector Fields to analize record titles
    for (uint32 i = rec_start_pos; i <= tmp_title_end[0]; ++i)
    {
      if (separators.find(read_buffer[i]) == std::string::npos)
        continue;

      fields.push_back(Field());

      fields[n_fields].data            = new uchar[i-f_start+1];
      std::copy(read_buffer+f_start, read_buffer+i, fields[n_fields].data);
      fields[n_fields].data[i-f_start] = '\0';
      fields[n_fields].len             = i - f_start;
      fields[n_fields].max_len         = fields[n_fields].len;
      fields[n_fields].min_len         = fields[n_fields].len;
      fields[n_fields].sep             = read_buffer[i];
      fields[n_fields].is_constant     = true;
      fields[n_fields].is_len_constant = true;
      fields[n_fields].is_numeric      = fields[n_fields].IsNum();
      fields[n_fields].Ham_mask        = new bool[fields[n_fields].len];

      if (fields[n_fields].is_numeric)
      {
        fields[n_fields].min_value = fields[n_fields].ToNum();
        fields[n_fields].max_value = fields[n_fields].min_value;
        fields[n_fields].num_values[fields[n_fields].min_value]++;
      }

      for (uint32 k = 0; k < fields[n_fields].len; ++k)
      {
        fields[n_fields].Ham_mask[k] = true;
      }
      fields[n_fields].block_desc.clear();

      f_start = i+1;
      n_fields++;
    }

    // OpenMP Parallel Region: Threads populate the Records structure
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel num_threads(no_threads)
    {
      int32 th_id       = omp_get_thread_num();
      int64 th_cur_rec  = 0;
      int64 th_f_start  = 0;
      int64 cur_pos     = th_id * records_per_th;

      for (int32 r = th_id; r > 0; --r)
        th_cur_rec += th_rec_count[r - 1];

      for (int64 r = th_cur_rec; r < (th_cur_rec + th_rec_count[th_id]); ++r)
      {
        records[r].title_end_pos  = tmp_title_end[cur_pos];
        records[r].seq_end_pos    = tmp_seq_end[cur_pos];
        th_f_start                = tmp_title_pos[cur_pos];

        ++cur_pos;

        if (read_buffer[th_f_start] != '@')
        {
          printf("\n[!] WARNING: Thread %d in p_Rank %d REC[%lld] val is %c. SubBlock %d\n", th_id, p_rank, r, read_buffer[th_f_start], p_subblock_count);
        }
        //my_assert(read_buffer[th_f_start] != '@');

        for (int64 f = th_f_start; f <= records[r].title_end_pos; ++f)
        {
          if (separators.find(read_buffer[f]) == std::string::npos)
            continue;

          records[r].fields_end_pos.push_back(f);
        }

        if(records[r].fields_end_pos.size() != n_fields)
        {
          printf("\n[!] WARNING: Thread %d in p_Rank %d has diff num_fields: %lu/%u. SubBlock %d\n", th_id, p_rank, records[r].fields_end_pos.size(), n_fields, p_subblock_count);
        }
      }
    }

    // Release temporary variables memory
    delete[] tmp_title_pos;
    delete[] tmp_title_end;
    delete[] tmp_seq_end;
    delete[] th_rec_count;
    tmp_title_pos = NULL;
    tmp_title_end = NULL;
    tmp_seq_end   = NULL;
    th_rec_count  = NULL;

    fastq_flags = FLAG_PLUS_ONLY | FLAG_DNA_PLAIN | FLAG_CONST_NUM_FIELDS;
    fastq_flags &= ~(FLAG_LINE_BREAKS | FLAG_VARIABLE_LENGTH | FLAG_TRY_LZ | FLAG_USE_DELTA | 
    FLAG_DELTA_CONSTANT | FLAG_DELTA_NO_BEGIN_NUC);
    
    // OpenMP Parallel Region: Threads analize, process and compress records in buffer.
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel default(shared) num_threads(no_threads)
    {
      std::vector<uint32> th_dna_occ(256, 0);
      bool th_alphabet_qua[256] = {false};
      uint32 th_max_qua_len  = 0, th_max_seq_len = 0;
      uint32 th_min_qua_len  = (uint32) -1, th_min_seq_len = (uint32) -1;

      bool is_delta             = false;
      bool is_delta_constant    = true;
      bool has_no_begin_nuc     = false;
      bool is_first_th_rec      = true;

      uchar sequence_start      = 0;
      uchar quality_start       = 0;

      // Analize fields data in title portion of records
      #pragma omp for nowait
      for (uint32 i = 0; i < n_fields; ++i)
      {
        AnalyzeTitleFields(read_buffer, rec_start_pos, i, fields, records, no_records);
      }
      
      // Analize ambiguity in DNA and transfer to Quality
      #pragma omp for nowait
      for (int64 i = 0; i < no_records; ++i)
      { 
        uint32 j, k, q, tmp;
        uint32 seq_start = records[i].title_end_pos + 1;
        uint32 seq_end   = records[i].seq_end_pos;
        uint32 seq_len   = seq_end - seq_start;
        uint32 qua_start = records[i].seq_end_pos + 3; 
        uint32 qua_end   = qua_start + seq_len;
        uint32 qua_len   = seq_len;
        
        if (is_first_th_rec)
        {
          sequence_start = read_buffer[seq_start];
          quality_start  = read_buffer[qua_start];

          is_delta = read_buffer[seq_start+1] >= '0' && read_buffer[seq_start+1] <= '3';

          if (sequence_start <= '3' && sequence_start >= '0')
          {
            has_no_begin_nuc = true;
            is_delta_constant = false;
          }

          is_first_th_rec = false;
        }
        
        if (is_delta && !has_no_begin_nuc && is_delta_constant)
        {
          is_delta_constant &= read_buffer[seq_start] == sequence_start;
          is_delta_constant &= read_buffer[qua_start] == quality_start;  
        }

        if(is_delta)
        {
          static const char delta_A[] = {'N', 'N', 'A', 'C', 'G', 'T'};
          static const char delta_C[] = {'N', 'N', 'C', 'A', 'T', 'G'};
          static const char delta_G[] = {'N', 'N', 'G', 'T', 'A', 'C'};
          static const char delta_T[] = {'N', 'N', 'T', 'G', 'C', 'A'};

          const uint32 translation = (uint32)is_delta_constant;
          const char* last_matrix = delta_A;

          uchar symbol;
          if (!has_no_begin_nuc)
          {
            symbol = read_buffer[seq_start]; 
            th_dna_occ[symbol]++;
          }
          else
          {
            symbol = 'A';
            my_assert(is_delta_constant == false);
          }
          
          uint32 n = 1 - (uint32)has_no_begin_nuc;
          for ( ; n < seq_len; ++n)
          {
            switch (symbol)
            {
              case 'A': last_matrix = delta_A; break;
              case 'C': last_matrix = delta_C; break;
              case 'G': last_matrix = delta_G; break;
              case 'T': last_matrix = delta_T; break;
              case 'N': 
              default: break;
            }
            symbol = last_matrix[read_buffer[seq_start+n]-'.'];
            my_assert(symbol =='A' || symbol == 'C' || symbol == 'G' || symbol == 'T' || symbol == 'N');
            th_dna_occ[symbol]++;

            read_buffer[seq_start+n-translation] = symbol;
            read_buffer[qua_start+n-translation] = read_buffer[seq_start+n];
          }
          
          for ( ; n < qua_len; ++n)
          {
            read_buffer[qua_start+n-translation] = read_buffer[seq_start+n];
          }

          read_buffer[seq_end-translation] = 0;
          read_buffer[qua_end-translation] = 0;
          seq_end -= translation;
          qua_len -= translation;
          qua_end -= translation;
        }

        // Check whether make transfer or not
        bool make_transfer = false;
        bool possible_transfer = true;
        for (j = seq_start, q = qua_start; j < seq_end; ++j, ++q)
        {
          if (trans_amb_codes[read_buffer[j]] == 1)
            continue;

          if (trans_amb_codes[read_buffer[j]] == 0)
          {
            possible_transfer = false;
            break;
          }
          else
          {
            if (read_buffer[q] < 33 || read_buffer[q] > 40)
            {
              possible_transfer = false;
              break;
            }
            make_transfer = true;
          }
        }

        if (make_transfer && possible_transfer)
        {
          for (j = k = seq_start; j < seq_end; ++j)
          {
            if ((tmp = trans_amb_codes[read_buffer[j]]) > 1)
            {
              read_buffer[j + seq_len + 3] = (uchar)(128 + (tmp << 3) - 16 + (read_buffer[j + seq_len + 3] - 33));
            }
            else
            {
              read_buffer[k++] = read_buffer[j];
            }
          }
          read_buffer[k] = '\0';
          seq_end = k;
        }

        // Recalculate seq_len after amb
        seq_len = seq_end - seq_start;

        if (!is_delta)
        {
          for (k = seq_start; k < seq_end; ++k)
          {
            ++th_dna_occ[read_buffer[k]];
          }
        }

        if (qua_len > th_max_qua_len)
          th_max_qua_len = qua_len;

        if (qua_len < th_min_qua_len)
          th_min_qua_len = qua_len;

        if (seq_len > th_max_seq_len)
          th_max_seq_len = seq_len;

        if (seq_len < th_min_seq_len)
          th_min_seq_len = seq_len;

        // Identify quality alphabet
        for (j = qua_start; j < qua_end; ++j)
          th_alphabet_qua[read_buffer[j]] = true;

        records[i].seq_len = seq_len;
        records[i].qua_len = qua_len;
      }

      // Reduce each thread's record analysis result
      #pragma omp critical
      { 
        for (uint32 k = 0; k < 256; ++k)
        {
          dna_occ[k] += th_dna_occ[k];
          alphabet_qua[k] |= th_alphabet_qua[k];
        }

        if (th_max_qua_len > max_quality_length)
          max_quality_length = th_max_qua_len;

        if (th_min_qua_len > min_quality_length)
          min_quality_length = th_min_qua_len;

        if (th_max_seq_len > max_sequence_length)
          max_sequence_length = th_max_seq_len;

        if (th_min_seq_len > min_sequence_length)
          min_sequence_length = th_min_seq_len;

        if (max_sequence_length > global_max_sequence_length)
          global_max_sequence_length = max_sequence_length;

        if (is_delta)
          fastq_flags |= FLAG_USE_DELTA;

        if (is_delta_constant)
          fastq_flags |= FLAG_DELTA_CONSTANT;

        if (has_no_begin_nuc)
          fastq_flags |= FLAG_DELTA_NO_BEGIN_NUC;
      }
    }

    bool is_length_variable = min_quality_length != max_quality_length;
    is_length_variable &= min_sequence_length != max_sequence_length;

    if (is_length_variable)
    {
      fastq_flags |= FLAG_VARIABLE_LENGTH;
    }

    symbols.clear();
    qualities.clear();

    int32 sym_idx = 0;
    int32 qua_idx = 0;
    for (uint32 i = 0; i < 256; ++i)
    {
      if (dna_occ[i])
      {
        symbols.push_back((uchar) i);
        sym_code[i] = (uchar) sym_idx++;
      }
      if (alphabet_qua[i])
      {
        qualities.push_back((uchar) i);
        qua_code[i] = (uchar) qua_idx++;
      }
    }  

    // Analize, compress and store record information
    #pragma omp parallel sections num_threads(no_threads)
    {
      // Store title header information
      #pragma omp section
      {
        title_bit_stream.Create(1);
        StoreTitle(title_bit_stream, fields, read_buffer, records, rec_start_pos, no_records);
      }

      // Analize DNA data and store information
      #pragma omp section
      {
        AnalyzeDNA(fastq_flags, dna_occ, symbols, sym_stats);
        dna_bit_stream.Create(1);
        StoreDNA(dna_bit_stream, read_buffer, records, no_records, fastq_flags, symbols, sym_stats, sym_code);
        dna_occ.clear();
      }

      // Analyze quality data and store information
      #pragma omp section
      {
        AnalyzeQuality(read_buffer, records, no_records, max_quality_length, qualities, qua_code, quality_stats, raw_quality_stats);
        quality_bit_stream.Create(1);
        StoreQuality(quality_bit_stream, read_buffer, records, no_records, max_quality_length, qua_code, qualities, quality_stats, qua_huf_codes, raw_qua_huf_codes);
      }

      // Store meta information about the current data to be written to output file
      #pragma omp section
      {
        p_write_buff_bit_stream.Create(1);
        p_write_buff_bit_stream.PutWord(no_records);
        p_write_buff_bit_stream.PutWord(max_quality_length);
        p_write_buff_bit_stream.PutWord(global_max_sequence_length);
        p_write_buff_bit_stream.PutByte((uchar) symbols.size());
        p_write_buff_bit_stream.PutByte((uchar) QUALITY_PLAIN);
        p_write_buff_bit_stream.PutByte((uchar) qualities.size());
      }
    }

    p_write_buff_bit_stream.PutWord(fastq_flags);
    p_write_buff_bit_stream.FlushPartialWordBuffer();

    uint32 quality_len_bits = BitStream::BitLength(max_quality_length);
    if ((fastq_flags & FLAG_VARIABLE_LENGTH) != 0)
    {
      for (uint32 i = 0; i < no_records; ++i)
      {
        uint32 qua_len = records[i].qua_len;
        p_write_buff_bit_stream.PutBits(qua_len, quality_len_bits);
      }
      p_write_buff_bit_stream.FlushPartialWordBuffer();
    }

    uint64 sb_offset = p_subblock_count != 1 ? r_buffer_curr_pos : (p_wr_start + r_buffer_curr_pos) ;
    p_write_buff_bit_stream.PutDWord(sb_offset);


    // printf("%d\t%llu\n", p_rank, sb_offset);


    // Prepare for next read from the FASTQ working region
    p_bytes_read += rec_end_pos + 1;
    r_buffer_curr_pos = p_bytes_read + p_wr_start;

    if ((r_buffer_curr_pos + r_buffer_size) > p_wr_end)
      r_buffer_size = p_wr_end - r_buffer_curr_pos + 1;

    // Release read_buffer allocated memory
    free(read_buffer);

    // Release records allocated memory
    delete[] records;
    records = NULL;
    no_records = 0;
    rec_start_pos = 0;

    // Release fields allocated memory
    fields.clear();
    n_fields = 0;

    // Release DNA and quality stats allocated memory
    if (sym_stats)
    {
      delete[] sym_stats;
      sym_stats = NULL;
    }
    if (quality_stats)
    {
      delete[] quality_stats;
      delete[] raw_quality_stats;
      quality_stats = NULL;
      raw_quality_stats = NULL;
    }

    // Release quality Huffman codes allocated memory
    if (qua_huf_codes)
    {
      delete[] qua_huf_codes;
      delete[] raw_qua_huf_codes;
      qua_huf_codes = NULL;
      raw_qua_huf_codes = NULL;
    }

    uint32 w_title_len = title_bit_stream.GetIO_Buffer_Pos();
    uint32 w_dna_len = dna_bit_stream.GetIO_Buffer_Pos();
    uint32 w_quality_len = quality_bit_stream.GetIO_Buffer_Pos();

    p_write_buff_bit_stream.PutWord(w_title_len);
    p_write_buff_bit_stream.PutWord(w_quality_len);

    uint32 w_info_buff_len = p_write_buff_bit_stream.GetIO_Buffer_Pos();
    
    // Allocate an specific ammount of memory for the copy buffer
    p_bytes_to_copy = w_title_len + w_dna_len + w_quality_len + w_info_buff_len;
    copy_buffer = (uchar*) malloc((p_bytes_to_copy) * sizeof(uchar));

    // Copy SubBlock's data into copy buffer
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel sections shared(copy_buffer) num_threads(no_threads)
    {
      // Copy general information about the compressed chunk to copy buffer
      #pragma omp section
      {
        std::vector<uchar> tem_buffer(p_write_buff_bit_stream.GetIO_Buffer());
        p_write_buff_bit_stream.Close();
        std::copy(tem_buffer.begin(), tem_buffer.end(), copy_buffer);
      }

      // Copy title information to copy buffer
      #pragma omp section
      {
        std::vector<uchar> tem_buffer(title_bit_stream.GetIO_Buffer());
        title_bit_stream.Close();
        uint32 bytes_in_w_buff = w_info_buff_len;
        std::copy(tem_buffer.begin(), tem_buffer.end(), copy_buffer+bytes_in_w_buff);
      }

      // Copy quality information to copy buffer
      #pragma omp section
      {
        std::vector<uchar> tem_buffer(quality_bit_stream.GetIO_Buffer());
        quality_bit_stream.Close();
        uint32 bytes_in_w_buff = w_info_buff_len + w_title_len;
        std::copy(tem_buffer.begin(), tem_buffer.end(), copy_buffer+bytes_in_w_buff);
      }

      // Copy DNA information to copy buffer
      #pragma omp section
      {
        std::vector<uchar> tem_buffer(dna_bit_stream.GetIO_Buffer());
        dna_bit_stream.Close();
        uint32 bytes_in_w_buff = w_info_buff_len + w_title_len + w_quality_len;
        std::copy(tem_buffer.begin(), tem_buffer.end(), copy_buffer+bytes_in_w_buff);
      }
    }

    // Add information to the header for the current block
    p_block.SBOL.push_back(p_bytes_to_copy);

    // If the write buffer will be full with current subblock, 
    // then build block with splitted subblock and wirte it to NGSC file
    // else copy subblock to write buffer
    // --------------------------------------------------------------------------------------------
    if ((p_bytes_written+p_bytes_to_copy) >= w_buffer_size)
    {
      p_block.BCSS |= LSBS;
      // Get the number of bytes needed to fill the write_buffer.
      uint32 bytes_fill_buffer = w_buffer_size - p_bytes_written;
      uint32 bytes_remaining   = p_bytes_to_copy - bytes_fill_buffer;

      // Copy bytes_fill_buffer bytes to write_buffer, to fill it up.
      std::copy(copy_buffer, copy_buffer+bytes_fill_buffer, write_buffer+p_bytes_written);

      // Add identifier to Block using 2 bytes
      write_buffer[0] = p_rank >> 8;
      write_buffer[1] = p_rank & 0xFF;

      p_block.SBOL.pop_back();
      p_block.SBOL.push_back(bytes_fill_buffer);

      MPI_File_write_shared(output_NGSC, write_buffer, w_buffer_size, MPI_CHAR, &status);
      // Take the time of completion of the writing operation
      timestamps.push_back(MPI_Wtime());

      // Save this Block's information
      p_blocks_data.push_back(p_block);
      
      free(write_buffer);
      write_buffer = (uchar*) malloc((w_buffer_size) * sizeof(uchar));

      if(!write_buffer)
      {
        printf("[E] ERROR: p_Rank %d can't allocate memory for the write_buffer.\n", p_rank);
        printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
        MPI_Finalize();
        exit(3);
      }

      // Copy rest of the bytes in copy_buffer to write_buffer, leaving the first 2 bytes
      // reserved for Block identifier
      std::copy(copy_buffer+bytes_fill_buffer, copy_buffer+p_bytes_to_copy, write_buffer+2);
      p_bytes_written = bytes_remaining + 2;
      p_block.BCSS = FSBS;
      p_block.SBOL.clear();
      p_block.SBOL.push_back(bytes_remaining);
    }
    else
    { 
      // Copy subblock in copy_buffer to write_buffer
      std::copy(copy_buffer, copy_buffer+p_bytes_to_copy, write_buffer+p_bytes_written);
      p_bytes_written += p_bytes_to_copy;
    }

    free(copy_buffer);
  } 

  // If the write_buffer is not empty
  // --------------------------------------------------------------------------------------------
  if(p_bytes_written > 0)
  {
    // Add identifier to Block using 2 bytes
    write_buffer[0] = p_rank >> 8;
    write_buffer[1] = p_rank & 0xFF;

    // Save this Block's information
    p_blocks_data.push_back(p_block);

    MPI_File_write_shared(output_NGSC, write_buffer, p_bytes_written, MPI_CHAR, &status);
    timestamps.push_back(MPI_Wtime());
  }

  free(write_buffer);

  BitStream header_bit_stream;
  MakeHeader(header_bit_stream, p_blocks_data);

  std::vector<uchar> block_buffer(header_bit_stream.GetIO_Buffer());
  int32 b_data_size = header_bit_stream.GetIO_Buffer_Pos();
  block_buffer.resize(b_data_size);

  // Gather information to build footer
  // --------------------------------------------------------------------------------------------
  if (p_rank == 0) 
  {
    all_info = (ProcessCompressionInfo*) malloc(g_size * sizeof(ProcessCompressionInfo));
  } 
  p_compress_info.n_blocks         = timestamps.size();
  p_compress_info.n_subblocks      = p_subblock_count;
  p_compress_info.last_block_size  = p_bytes_written;
  p_compress_info.blocks_data_size = b_data_size;

  // create an MPI type for struct ProcessCompressionInfo
  const int32 nitems = 4;
  int32 blocklengths[4] = {1,1,1,1};
  MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_UNSIGNED, MPI_INT};
  MPI_Datatype mpi_p_info;
  MPI_Aint offsets[4];
  MPI_Aint addr[5];

  MPI_Get_address(&p_compress_info, &addr[0]);
  MPI_Get_address(&p_compress_info.n_blocks, &addr[1]);
  MPI_Get_address(&p_compress_info.n_subblocks, &addr[2]);
  MPI_Get_address(&p_compress_info.last_block_size, &addr[3]);
  MPI_Get_address(&p_compress_info.blocks_data_size, &addr[4]);

  offsets[0] = addr[1] - addr[0];
  offsets[1] = addr[2] - addr[0];
  offsets[2] = addr[3] - addr[0];
  offsets[3] = addr[4] - addr[0];

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_p_info);
  MPI_Type_commit(&mpi_p_info);

  // Gather information from all processes
  MPI_Gather(&p_compress_info, 1, mpi_p_info, &(all_info[0]), 1, mpi_p_info, 0, MPI_COMM_WORLD);

  int32  n_total_blocks     = 0;    // Total number of Blocks in NGSC.
  int32  n_total_subblocks  = 0;    // Total number of SubBlocks in NGSC.
  int32  total_b_data_size  = 0;    // Total number of bytes in blocks data.
  int32  *displs_ts         = NULL; // Displacements for Timestamp gathering.
  int32  *displs_bd         = NULL; // Displacements for Block data
  int32  *timestamps_counts = NULL; // Number of timestamps in the toc list of each p_rank.
  double *TOCW_buffer       = NULL; // Buffer to save all timestamps gathered from each p_rank.
  uchar  *b_data_buffer     = NULL; // Buffer to save all block data gathered from each p_rank.
  uint32 *lb_sizes          = NULL; // Size of each p_rank's last Block written.
  int32  *all_b_data_sizes  = NULL; // Block data sizes from all P's
  
  if (p_rank == 0)
  {
    displs_ts         = (int32*) malloc(g_size * sizeof(int32));
    displs_bd         = (int32*) malloc(g_size * sizeof(int32));
    timestamps_counts = (int32*) malloc(g_size * sizeof(int32));
    all_b_data_sizes  = (int32*) malloc(g_size * sizeof(int32));
    lb_sizes          = (uint32*) malloc(g_size * sizeof(uint32));

    for (int32 i = 0; i < g_size; ++i)
    {
      displs_ts[i]         = n_total_blocks;
      displs_bd[i]         = total_b_data_size;
      n_total_blocks      += all_info[i].n_blocks;
      n_total_subblocks   += all_info[i].n_subblocks;
      total_b_data_size   += all_info[i].blocks_data_size;
      timestamps_counts[i] = all_info[i].n_blocks;
      all_b_data_sizes[i]  = all_info[i].blocks_data_size;
      lb_sizes[i]          = all_info[i].last_block_size;
    }

    TOCW_buffer   = (double*) malloc(n_total_blocks * sizeof(double));
    b_data_buffer = (uchar*) malloc(total_b_data_size * sizeof(uchar));
  }
  // Gather time of completion lists from all processes
  MPI_Gatherv(timestamps.data(), timestamps.size(), MPI_DOUBLE, TOCW_buffer, timestamps_counts, displs_ts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // Gather block data from all processes
  MPI_Gatherv(block_buffer.data(), block_buffer.size(), MPI_UNSIGNED_CHAR, b_data_buffer, all_b_data_sizes, displs_bd, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  // Sort list of time of completion to stablish order of writing
  if (p_rank == 0)
  {
    std::multimap<double,int32> blocks_order;

    for (int32 i = 0; i < g_size; ++i)
    {
      for (int32 j = 0; j < timestamps_counts[i]; ++j)
      {
        int32 ts_pos = j + displs_ts[i];

        blocks_order.emplace(TOCW_buffer[ts_pos], i);
      }
    }
    
    BitStream footer_bit_stream;

    // Build footer
    MakeFooter(footer_bit_stream, g_size, FASTQ_size, n_total_blocks, n_total_subblocks, blocks_order, 
               lb_sizes, b_data_buffer, total_b_data_size);

    std::vector<uchar> footer_buffer(footer_bit_stream.GetIO_Buffer());
    int32 footer_size = footer_bit_stream.GetIO_Buffer_Pos();
    
    // Write footer to NGSC file
    MPI_File_write_shared(output_NGSC, &footer_buffer[0], footer_size, MPI_CHAR, &status);
    // Close stream
    footer_bit_stream.Close();

    // Release memory
    free(all_info);
    free(displs_ts);
    free(displs_bd);
    free(timestamps_counts);
    free(TOCW_buffer);
    free(b_data_buffer);
    free(lb_sizes);
    free(all_b_data_sizes); 
  }

  // Stop timer
  p_timer_end = MPI_Wtime();

  if (p_rank == 0)
    printf("\nRANK\tCOMP_TIME\tN_BLOCK\tN_SUBBLOCKS\n----------------------------------------------\n");

  MPI_Barrier(MPI_COMM_WORLD);
  printf("%03d\t%f\t%ld\t%d\n", p_rank, p_timer_end-p_timer_start, timestamps.size(), p_subblock_count);


  // TODO Temp for performance testing
  if (p_rank == 0)
  {
    std::string folder_name = "./performanceOutput-Com";

    // makes output directory if needed
    if (mkdir(folder_name.c_str(), 0777) != -1)
      printf("Performance output directory made.\n");

    std::string file_path = "./" + folder_name + "/" + std::to_string(g_size) + "," + std::to_string(no_threads) + ".txt";
    std::ofstream test_file;
    test_file.open(file_path);
    test_file << p_timer_end-p_timer_start;
  }


  // Wait for all ranks, before closing the I/O files.
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_close(&input_FASTQ);
  MPI_File_close(&output_NGSC);
}