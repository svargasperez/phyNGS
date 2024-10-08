/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#include <mpi.h>
#include <vector>
#include <stdio.h>
#include <map>
#include <math.h>
#include <algorithm>
#include "defs.h"
#include "utils.h"
#include "structures.h"
#include "bit_stream.h"
#include "huffman.h"

//todo temp is this needed
#include "tasks.h"

// --------------------------------------------------------------------------------------------
void AnalyzeTitleFields(uchar *read_buffer, uint32 field_0_rec_0_pos, int32 c_field, 
  std::vector<Field> &fields, Record *&records, int64 no_records)
{
  const uint32 MAX_FIELD_STAT_LEN  = 128;
  const uint32 DEFAULT_B_SIZE      = 32;
  const uint32 MAX_NUM_VAL_HUF     = 512;

  int32 prev_value = 0;
  for (uint32 i = 0; i < no_records; ++i)
  {
    uint32 start_pos = 0; 
    uint32 end_pos = records[i].fields_end_pos[c_field];

    if (c_field) // Other than the field 0 assing
    {
      start_pos = records[i].fields_end_pos[c_field-1] + 1;
    }
    else
    {
      if (i)
      {
        int32 rec_const_lines = (records[i-1].seq_end_pos - records[i-1].title_end_pos) * 2;
        start_pos = records[i-1].fields_end_pos[fields.size()-1] + rec_const_lines + 3;
      }
      else
      {
        start_pos = field_0_rec_0_pos;
      }
    }

    uint32 field_len = end_pos - start_pos;

    if (field_len > fields[c_field].max_len)
    {
      fields[c_field].max_len = field_len;
    }
    else if (field_len < fields[c_field].min_len)
    {
      fields[c_field].min_len = field_len;
    }

    // Check whether field in a block is constant
    if (i % DEFAULT_B_SIZE > 0)
    {
      if (fields[c_field].block_str_len != field_len)
      {
        fields[c_field].block_desc.back().is_block_constant = false;
      }
      else if (fields[c_field].block_desc.back().is_block_constant)
      {
        fields[c_field].block_desc.back().is_block_constant = std::equal(read_buffer+start_pos, read_buffer+end_pos,
          read_buffer+fields[c_field].block_str_start);
      }
    }
    else
    {
      fields[c_field].block_str_start = start_pos;
      fields[c_field].block_str_len   = field_len;
      fields[c_field].block_desc.push_back(Field::BlockDesc(true, true, true, 0));
    }

    fields[c_field].chars.resize(fields[c_field].max_len);
    uint32 chars_len = MIN(MAX_FIELD_STAT_LEN, field_len);

    for (uint32 x = 0; x < chars_len; ++x)
    {
      fields[c_field].chars[x][read_buffer[start_pos+x]]++;
    }
    for (uint32 x = MAX_FIELD_STAT_LEN; x < field_len; ++x)
    {
      fields[c_field].chars[MAX_FIELD_STAT_LEN][read_buffer[start_pos+x]]++;
    }
    
    if (fields[c_field].is_constant)
    {
      if (field_len != fields[c_field].len)
      {
        fields[c_field].is_constant = false;
      }
      else
      {
        fields[c_field].is_constant = std::equal(fields[c_field].data, fields[c_field].data + fields[c_field].len, read_buffer+start_pos);
      }
    }

    if (fields[c_field].is_len_constant)
    {
      fields[c_field].is_len_constant = fields[c_field].len == field_len;
    }
    
    if (fields[c_field].is_numeric)
    {
      fields[c_field].is_numeric = utils::is_num(read_buffer+start_pos, field_len);
      if (fields[c_field].is_numeric)
      {
        int32 value = utils::to_num(read_buffer+start_pos, field_len);
        if (value < fields[c_field].min_value)
        {
          fields[c_field].min_value = value;
        }
        else if (value > fields[c_field].max_value)
        {
          fields[c_field].max_value = value;
        }

        if (i % DEFAULT_B_SIZE > 0)
        {
          if (fields[c_field].block_value != value)
            fields[c_field].block_desc.back().is_block_value_constant = false;

          // when there are more than 2 records on the new block
          if (i % DEFAULT_B_SIZE > 1)
          {
            if (value - prev_value != fields[c_field].block_delta)
              fields[c_field].block_desc.back().is_block_delta_constant = false;
          }
          else // when there are 2 records on the new block
          {
            fields[c_field].block_delta = value - prev_value;
            fields[c_field].block_desc.back().block_delta_constant = value - prev_value;
          }
        }
        else
        {
          fields[c_field].block_value = value;
        }

        if (i >= 1)
        {
          if (i > 1)
          {
            if (value - prev_value > fields[c_field].max_delta)
            {
              fields[c_field].max_delta = value - prev_value;
            }
            if (value - prev_value < fields[c_field].min_delta)
            {
              fields[c_field].min_delta = value - prev_value;
            }
          }
          else // i == 1
          {
            fields[c_field].max_delta = value - prev_value;
            fields[c_field].min_delta = value - prev_value;
          }

          fields[c_field].delta_values[value - prev_value]++;
          if (fields[c_field].delta_values.size() > MAX_NUM_VAL_HUF)
          {
            fields[c_field].delta_values.clear();
          }
        }

        if (fields[c_field].num_values.size())
        {
          fields[c_field].num_values[value]++;
          if (fields[c_field].num_values.size() > MAX_NUM_VAL_HUF)
          {
            fields[c_field].num_values.clear();
          }
        }

        prev_value = value;
      }
    }
    if (!fields[c_field].is_constant)
    {
      for (uint32 p = 0; p < field_len && p < fields[c_field].len; ++p)
      {
        fields[c_field].Ham_mask[p] &= fields[c_field].data[p] == read_buffer[p+start_pos];
      }
    }
  }

  // Find better encoding of numeric values
  if (!fields[c_field].is_numeric)
  {
    if (!fields[c_field].is_constant)
    {
      fields[c_field].chars.resize(MIN(fields[c_field].max_len, MAX_FIELD_STAT_LEN+1));
      fields[c_field].no_of_bits_per_len = BitStream::BitLength(fields[c_field].max_len - fields[c_field].min_len);
    }
  } 
  else
  {
    int32 diff;
    if (fields[c_field].max_value - fields[c_field].min_value < fields[c_field].max_delta - fields[c_field].min_delta)
    {
      fields[c_field].is_delta_coding = false;
      diff = fields[c_field].max_value - fields[c_field].min_value;
    }
    else
    {
      fields[c_field].is_delta_coding = true;
      diff = fields[c_field].max_delta - fields[c_field].min_delta;
    }

    fields[c_field].no_of_bits_per_num = BitStream::BitLength(diff);
    diff = fields[c_field].max_value - fields[c_field].min_value;
    fields[c_field].no_of_bits_per_value = BitStream::BitLength(diff);
  }
}

// --------------------------------------------------------------------------------------------
void AnalyzeDNA(uint32 &fastq_flags, std::vector<uint32> &dna_occ, std::vector<uchar> &symbols,
  uint32 *&sym_stats)
{
  if (sym_stats)
  {
    delete[] sym_stats;
  }
  sym_stats = new uint32[symbols.size()];
  for (uint32 i = 0; i < symbols.size(); ++i)
  {
    sym_stats[i] = dna_occ[symbols[i]];
  }

  if (symbols.size() > 4)
  {
    fastq_flags &= ~FLAG_DNA_PLAIN;
  }
  else
  {
    std::vector<uint32> sym_tmp(dna_occ.begin(), dna_occ.end());
    std::sort(sym_tmp.begin(), sym_tmp.end());

    if (sym_tmp[0] > sym_tmp[2] + sym_tmp[3])
    {
      fastq_flags &= ~FLAG_DNA_PLAIN;
    }
    else
    {
      fastq_flags |= FLAG_DNA_PLAIN;
    }
  }
}

// --------------------------------------------------------------------------------------------
void AnalyzeQuality(uchar *read_buffer, Record *&records, int64 no_records, uint32 max_quality_length, 
  std::vector<uchar> &qualities, uchar *qua_code, uint32 **&quality_stats, uint32 *&raw_quality_stats)
{
  quality_stats = new uint32*[max_quality_length+1];
  raw_quality_stats = new uint32[(max_quality_length+1)*qualities.size()]();
  
  for (uint32 i = 0; i <= max_quality_length; ++i)
  {
    quality_stats[i] = raw_quality_stats+i*qualities.size();
  }
  for (uint32 i = 0; i < (max_quality_length+1)*qualities.size(); ++i)
  {
    raw_quality_stats[i] = 0;
  }
  
  for (uint32 i = 0; i < no_records; ++i)
  {
    int64 qua_start = records[i].seq_end_pos + 3; 
    int64 qua_end = qua_start + (records[i].seq_end_pos - records[i].title_end_pos) - 1;

    for (uint32 j = 0; qua_start < qua_end && read_buffer[qua_start] != '\0'; ++j, ++qua_start)
    {
      quality_stats[j+1][qua_code[read_buffer[qua_start]]]++;
      quality_stats[0][qua_code[read_buffer[qua_start]]]++;
    }
  }    
}

// --------------------------------------------------------------------------------------------
void StoreTitle(BitStream &title_bit_stream, std::vector<Field> &fields, uchar *read_buffer,
  Record *&records, uint32 field_0_rec_0_pos, int64 no_records)
{
  const int32  MAX_FIELD_STAT_LEN  = 128;
  const int32  DEFAULT_B_SIZE      = 32;
  const uint32 MAX_NUM_VAL_HUF     = 512;

  int32 n_blocks = (no_records + DEFAULT_B_SIZE -1) / DEFAULT_B_SIZE;
  int32 n_fields = fields.size();
  std::vector<int32> prev_value(n_fields);
  int32 rec_count = 0;
  int32 first_rec_in_block; 

  title_bit_stream.PutWord(n_fields);

  for (int32 i = 0; i < n_fields; ++i)
  {
    title_bit_stream.PutByte(fields[i].sep);
    title_bit_stream.PutByte(fields[i].is_constant);
    if (fields[i].is_constant)
    {
      title_bit_stream.PutWord(fields[i].len);
      title_bit_stream.PutBytes(fields[i].data, fields[i].len);
      continue;
    }

    title_bit_stream.PutByte(fields[i].is_numeric);
    if (fields[i].is_numeric)
    {
      title_bit_stream.PutWord(fields[i].min_value);
      title_bit_stream.PutWord(fields[i].max_value);
      title_bit_stream.PutWord(fields[i].min_delta);
      title_bit_stream.PutWord(fields[i].max_delta);

      int32 diff, base;
      std::map<int32, int32> &map_stats = fields[i].num_values;
      if (fields[i].max_value-fields[i].min_value < fields[i].max_delta-fields[i].min_delta)
      {
        diff = fields[i].max_value - fields[i].min_value;
        base = fields[i].min_value;
      }
      else
      {
        diff = fields[i].max_delta - fields[i].min_delta;
        base = fields[i].min_delta;
        map_stats = fields[i].delta_values;
      }

      diff++;
      if (diff <= (int32)MAX_NUM_VAL_HUF && map_stats.size())     // few values, so use Huffman for them
      {
        HuffmanEncoder* huf = fields[i].Huffman_global = new HuffmanEncoder(Field::HUF_GLOBAL_SIZE);
        for (uint32 j = 0; j < (uint32)diff; ++j)
        {
          huf->Insert(map_stats[base+j]);
        }
        huf->Complete();
        HuffmanEncoder::StoreTree(title_bit_stream, *huf);
      }

      continue;
    }

    title_bit_stream.PutByte(fields[i].is_len_constant);
    title_bit_stream.PutWord(fields[i].len);
    title_bit_stream.PutWord(fields[i].max_len);
    title_bit_stream.PutWord(fields[i].min_len);
    title_bit_stream.PutBytes(fields[i].data, fields[i].len);

    for (uint32 j = 0; j < fields[i].len; ++j)
    {
      title_bit_stream.PutBit(fields[i].Ham_mask[j]);
    }
    fields[i].Huffman_local.resize(MIN(fields[i].max_len+1, MAX_FIELD_STAT_LEN+1));

    for (uint32 j = 0; j < MIN(fields[i].max_len, MAX_FIELD_STAT_LEN); ++j)
    {
      fields[i].Huffman_local[j] = NULL;
      if (j >= fields[i].len || !fields[i].Ham_mask[j])
      {
        HuffmanEncoder* huf = fields[i].Huffman_local[j] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
        for (uint32 k = 0; k < Field::HUF_LOCAL_SIZE; ++k)
        {
          huf->Insert(fields[i].chars[j][k]);
        }
        huf->Complete(true);
        HuffmanEncoder::StoreTree(title_bit_stream, *huf);
      }
    } 
    if (fields[i].max_len >= MAX_FIELD_STAT_LEN)
    {
      HuffmanEncoder* huf = fields[i].Huffman_local[MAX_FIELD_STAT_LEN] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
      for (uint32 k = 0; k < Field::HUF_LOCAL_SIZE; ++k)
      {
        huf->Insert(fields[i].chars[(uchar) MAX_FIELD_STAT_LEN][k]);
      }
      huf->Complete(true);
      HuffmanEncoder::StoreTree(title_bit_stream, *huf);
    }

    title_bit_stream.FlushPartialWordBuffer();
  }

  // Insert specific data about non-constant records
  for (int32 block_no = 0; block_no < n_blocks; ++block_no)
  {
    first_rec_in_block = rec_count;
    rec_count += DEFAULT_B_SIZE;

    if (rec_count > no_records)
      rec_count = no_records;

    for (int32 i = 0; i < n_fields; ++i)
    {
      if (fields[i].is_constant)
        continue;

      Field::BlockDesc& block_desc = fields[i].block_desc[block_no];
      prev_value[i] = 0;
      if (!fields[i].is_numeric)
      {
        title_bit_stream.PutBit(block_desc.is_block_constant);
      }

      if (fields[i].is_numeric)
      {
        block_desc.is_block_delta_constant &= (int32)block_desc.block_delta_constant == fields[i].min_delta;
        if (fields[i].is_delta_coding)
        {
          title_bit_stream.PutBit(block_desc.is_block_delta_constant);
        }
        else
        {
          title_bit_stream.PutBit(block_desc.is_block_value_constant);
        }
      }   
    }

    for (int32 i = first_rec_in_block; i < rec_count; ++i)
    {
      for (int32 c_field = 0; c_field < n_fields; ++c_field)
      {
        Field &cur_field = fields[c_field];
        uint32 start_pos = 0;
        uint32 end_pos = records[i].fields_end_pos[c_field];
        
        if (c_field) // Other than the field 0 assing
          start_pos = records[i].fields_end_pos[c_field-1] + 1;
        else
        {
          if (i)
          {
            int32 rec_const_lines = (records[i-1].seq_end_pos - records[i-1].title_end_pos) * 2;
            start_pos = records[i-1].fields_end_pos[fields.size()-1] + rec_const_lines + 3;
          }
          else
            start_pos = field_0_rec_0_pos;
        }

        if (cur_field.is_constant)
          continue;

        if (cur_field.is_numeric)
        {
          int32 value = utils::to_num(read_buffer+start_pos, end_pos-start_pos);

          if ((i % DEFAULT_B_SIZE) == 0) // Is the first record of Block?
          {
            title_bit_stream.PutBits(value-cur_field.min_value, cur_field.no_of_bits_per_value);
          }
          else if ((cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_delta_constant) ||
            (!cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_value_constant))
          {
            int32 to_store;
            if (cur_field.is_delta_coding)
            {
              to_store = value - prev_value[c_field] - cur_field.min_delta;
            }
            else
            {
              to_store = value - cur_field.min_value;
            }

            if (cur_field.Huffman_global)
            {
              const HuffmanEncoder::Code* codes = cur_field.Huffman_global->GetCodes();
              title_bit_stream.PutBits(codes[to_store].code, codes[to_store].len);
            }
            else
            {
              title_bit_stream.PutBits(to_store, cur_field.no_of_bits_per_num);
            }
          }

          prev_value[c_field] = value;
          continue;
        }

        if ((i % DEFAULT_B_SIZE) > 0 && cur_field.block_desc[block_no].is_block_constant)
        {
          continue;
        }

        if (!cur_field.is_len_constant)
        {
          title_bit_stream.PutBits(end_pos-start_pos - cur_field.min_len, cur_field.no_of_bits_per_len);
        }
        
        for (uint32 j = 0; j < end_pos-start_pos; ++j)
        {
          if (j >= cur_field.len || !cur_field.Ham_mask[j])
          {
            uchar c = read_buffer[start_pos+j];
            const HuffmanEncoder::Code* codes = cur_field.Huffman_local[MIN(j, MAX_FIELD_STAT_LEN)]->GetCodes();
            title_bit_stream.PutBits(codes[c].code, codes[c].len);
          }
        }
      }
    }
    title_bit_stream.FlushPartialWordBuffer();
  }
}

// --------------------------------------------------------------------------------------------
void StoreDNA(BitStream &dna_bit_stream, uchar *read_buffer, Record *&records, int64 no_records, 
  uint32 fastq_flags, std::vector<uchar> &symbols, uint32 *&sym_stats, uchar *sym_code)
{
  HuffmanEncoder::Code *sym_huf_codes = NULL;
  HuffmanEncoder *Huffman_sym = NULL;

  for (uint32 i = 0; i < symbols.size(); ++i)
  {
    dna_bit_stream.PutByte(symbols[i]);
  }
  dna_bit_stream.FlushPartialWordBuffer();

  if ((fastq_flags & FLAG_DNA_PLAIN) == 0)
  {
    sym_huf_codes = new HuffmanEncoder::Code[symbols.size()];
    Huffman_sym   = new HuffmanEncoder((uint32) symbols.size());

    for (uint32 i = 0; i < symbols.size(); ++i)
    {
      Huffman_sym->Insert(sym_stats[i]);
    }

    HuffmanEncoder::Code* codes = Huffman_sym->Complete();
    for (uint32 i = 0; i < symbols.size(); ++i)
    {
      sym_huf_codes[i] = codes[i];
    }

    HuffmanEncoder::StoreTree(dna_bit_stream, *Huffman_sym);
  }

  for (int64 i = 0; i < no_records; ++i)
  {
    uint32 seq_start = records[i].title_end_pos + 1; 
    uint32 seq_end   = records[i].seq_end_pos;

    for (uint32 j = seq_start; read_buffer[j] != '\0' && j < seq_end; ++j)
    {
      uchar c = read_buffer[j];
      if ((fastq_flags & FLAG_DNA_PLAIN) != 0) 
        dna_bit_stream.Put2Bits(sym_code[c]);
      else 
        dna_bit_stream.PutBits(sym_huf_codes[sym_code[c]].code, sym_huf_codes[sym_code[c]].len);
    }
  }

  dna_bit_stream.FlushPartialWordBuffer();

  // Release Huffman trees - DNA symbols
  if (Huffman_sym)
  {
    delete Huffman_sym;
    Huffman_sym = NULL;
    delete[] sym_huf_codes;
    sym_huf_codes = NULL;
  }
}

// --------------------------------------------------------------------------------------------
void StoreQuality(BitStream &quality_bit_stream, uchar *read_buffer, Record *&records, int64 no_records, 
  uint32 max_quality_length, uchar *qua_code, std::vector<uchar> &qualities, uint32 **&quality_stats, 
  HuffmanEncoder::Code **&qua_huf_codes, HuffmanEncoder::Code *&raw_qua_huf_codes)
{
  // Store qualities
  for (uint32 i = 0; i < qualities.size(); ++i)
  {
    quality_bit_stream.PutByte(qualities[i]);
  }
  quality_bit_stream.FlushPartialWordBuffer();

  qua_huf_codes = new HuffmanEncoder::Code*[max_quality_length+1];
  raw_qua_huf_codes = new HuffmanEncoder::Code[(max_quality_length+1)*qualities.size()];
  for (uint32 i = 0; i <= max_quality_length; ++i)
  {
    qua_huf_codes[i] = &raw_qua_huf_codes[i*qualities.size()];
  }

  for (uint32 i = 0; i <= max_quality_length; ++i)
  {
    HuffmanEncoder *huf = new HuffmanEncoder((uint32) qualities.size());

    for (uint32 j = 0; j < qualities.size(); ++j)
    {
      huf->Insert(quality_stats[i][j]);
    }

    HuffmanEncoder::Code* codes = huf->Complete();
    for (uint32 j = 0; j < qualities.size(); ++j)
    {
      qua_huf_codes[i][j] = codes[j];
    }
    HuffmanEncoder::StoreTree(quality_bit_stream, *huf);
  }

  quality_bit_stream.FlushPartialWordBuffer();
  
  for (int64 i = 0; i < no_records; ++i)
  {
    uint32 qua_start = records[i].seq_end_pos + 3; 
    uint32 qua_end = qua_start + (records[i].seq_end_pos - records[i].title_end_pos) - 1;

    for (uint32 j = qua_start, k = 0; j < qua_end && read_buffer[j] != '\0'; ++j, ++k)
    {
      int32 qua = qua_code[read_buffer[j]];
      quality_bit_stream.PutBits(qua_huf_codes[k+1][qua].code, qua_huf_codes[k+1][qua].len);
    }
  }

  quality_bit_stream.FlushPartialWordBuffer();
}

// --------------------------------------------------------------------------------------------
void FetchTitleHeader(BitStream &sb_a_bit_stream, std::vector<Field> &fields)
{
  uint32 tmp;
  uint32 n_fields = 0;
  const uint32 MAX_NUM_VAL_HUF     = 512;
  const uint32 MAX_FIELD_STAT_LEN  = 128;

  sb_a_bit_stream.GetWord(n_fields);
  fields.resize(n_fields);

  for (uint32 i = 0; i < n_fields; ++i)
  {
    sb_a_bit_stream.GetByte(tmp);
    my_assert(tmp < 256);

    fields[i].sep = (uchar) tmp;
    sb_a_bit_stream.GetByte(tmp);
    fields[i].is_constant = tmp != 0;

    if (fields[i].is_constant)
    {
      sb_a_bit_stream.GetWord(tmp);
      fields[i].len = tmp;
      fields[i].data = new uchar[fields[i].len+1];
      sb_a_bit_stream.GetBytes(fields[i].data, fields[i].len);
      continue;
    }

    sb_a_bit_stream.GetByte(tmp);
    fields[i].is_numeric = tmp != 0;
    if (fields[i].is_numeric)
    {
      sb_a_bit_stream.GetWord(tmp);
      fields[i].min_value = tmp;
      sb_a_bit_stream.GetWord(tmp);
      fields[i].max_value = tmp;
      my_assert(fields[i].min_value <= fields[i].max_value);

      sb_a_bit_stream.GetWord(tmp);
      fields[i].min_delta = (int32) tmp;
      sb_a_bit_stream.GetWord(tmp);
      fields[i].max_delta = tmp;
      my_assert(fields[i].min_delta <= fields[i].max_delta);
 
      int32 diff;
      // std::map<int32, int32> &map_stats = fields[i].num_values;
      if (fields[i].max_value - fields[i].min_value < fields[i].max_delta - fields[i].min_delta)
      {
        fields[i].is_delta_coding = false;
        diff = fields[i].max_value - fields[i].min_value;
      }
      else
      {
        fields[i].is_delta_coding = true;
        diff = fields[i].max_delta - fields[i].min_delta;
        // map_stats = fields[i].delta_values;
      }

      fields[i].no_of_bits_per_num = BitStream::BitLength(diff);
      int32 v_diff = fields[i].max_value - fields[i].min_value;
      fields[i].no_of_bits_per_value = BitStream::BitLength(v_diff);

      diff++;
      if (diff <= (int32)MAX_NUM_VAL_HUF)     // a few values, so use Huffman for them
      { 
        fields[i].Huffman_global = new HuffmanEncoder();
        HuffmanEncoder::LoadTree(sb_a_bit_stream, *fields[i].Huffman_global);
        sb_a_bit_stream.FlushInputWordBuffer();
      }

      continue;
    }

    sb_a_bit_stream.GetByte(tmp);
    fields[i].is_len_constant = tmp != 0;
    sb_a_bit_stream.GetWord(tmp);
    fields[i].len = tmp;
    sb_a_bit_stream.GetWord(tmp);
    fields[i].max_len = tmp;
    sb_a_bit_stream.GetWord(tmp);
    fields[i].min_len = tmp;
    fields[i].no_of_bits_per_len = BitStream::BitLength(fields[i].max_len - fields[i].min_len);
    fields[i].data = new uchar[fields[i].len+1];
    sb_a_bit_stream.GetBytes(fields[i].data, fields[i].len);
    fields[i].Ham_mask = new bool[fields[i].len+1];

    for (uint32 j = 0; j < fields[i].len; ++j)
    {
      sb_a_bit_stream.GetBits(tmp, 1);
      fields[i].Ham_mask[j] = tmp != 0;
    }
    fields[i].Huffman_local.resize(MIN(fields[i].max_len, MAX_FIELD_STAT_LEN+1));

    for (uint32 j = 0; j < MIN(fields[i].max_len, MAX_FIELD_STAT_LEN); ++j)
    {
      fields[i].Huffman_local[j] = NULL;
      if (j >= fields[i].len || !fields[i].Ham_mask[j])
      {
        fields[i].Huffman_local[j] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
        HuffmanEncoder::LoadTree(sb_a_bit_stream, *fields[i].Huffman_local[j]);
      }
    }
    if (fields[i].max_len >= MAX_FIELD_STAT_LEN)
    {
      fields[i].Huffman_local[MAX_FIELD_STAT_LEN] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
      HuffmanEncoder::LoadTree(sb_a_bit_stream, *fields[i].Huffman_local[MAX_FIELD_STAT_LEN]);   
    }
    sb_a_bit_stream.FlushInputWordBuffer();
  }
}

// --------------------------------------------------------------------------------------------
void FetchTitleBody(BitStream &sb_a_bit_stream, std::vector<Field> &fields, Record *&records, 
  int64 no_records, uint32 fastq_flags, int64 &prev_title_len)
{
  const int32 MAX_FIELD_STAT_LEN  = 128;
  const int32 DEFAULT_B_SIZE      = 32;

  int32 n_blocks = (no_records + DEFAULT_B_SIZE -1) / DEFAULT_B_SIZE;
  int32 n_fields = fields.size();
  std::vector<int32> prev_value(n_fields);
  int32 rec_count = 0;
  int32 first_rec_in_block;
  bool is_num_fields_constant = (fastq_flags & FLAG_CONST_NUM_FIELDS) != 0;

  for (int32 i = 0; i < n_fields; ++i)
  {
    fields[i].block_desc.resize(n_blocks);
  }

  uint32 n_fields_bits = BitStream::BitLength(n_fields);

  for (int32 block_no = 0; block_no < n_blocks; ++block_no)
  {
    first_rec_in_block = rec_count;
    rec_count += DEFAULT_B_SIZE;

    if (rec_count > no_records)
      rec_count = no_records;

    uint32 tmp = 0;
    for (int32 i = 0; i < n_fields; ++i)
    {
      if (fields[i].is_constant)
        continue;

      prev_value[i] = 0;
      if (!fields[i].is_numeric)
      {
        sb_a_bit_stream.GetBit(tmp);
        fields[i].block_desc[block_no].is_block_constant = tmp != 0;
      }
      else
      {
        sb_a_bit_stream.GetBit(tmp);
        if (fields[i].is_delta_coding)
        {
          fields[i].block_desc[block_no].is_block_delta_constant = tmp != 0;
        }
        else
        {
          fields[i].block_desc[block_no].is_block_value_constant = tmp != 0;
        }
      }   
    }

    for (int32 i = first_rec_in_block; i < rec_count; ++i)
    {
      uint32 cn_fields = n_fields;
      if (!is_num_fields_constant)
      {
        sb_a_bit_stream.GetBits(tmp, n_fields_bits);
        cn_fields = tmp;
      }

      for (uint32 j = 0; j < cn_fields; ++j)
      {
        if (fields[j].is_constant)
        {
          records[i].AppendData(fields[j].data, fields[j].len, 'T');
          records[i].AppendChar(fields[j].sep, 'T');
          continue;
        }
        if (fields[j].is_numeric)
        {
          uint32 num_val = 0;

          if ((i % DEFAULT_B_SIZE) == 0) // Is the first record of Block?
          {
            sb_a_bit_stream.GetBits(num_val, fields[j].no_of_bits_per_value);
            num_val += fields[j].min_value;
            uchar tmp_str[10];
            int32 tmp_len = 0;
            tmp_len = utils::to_string(tmp_str, num_val);

            records[i].AppendData(tmp_str, tmp_len, 'T');
            prev_value[j] = num_val;
          }
          else
          {
            if ((fields[j].is_delta_coding && !fields[j].block_desc[block_no].is_block_delta_constant) ||
              (!fields[j].is_delta_coding && !fields[j].block_desc[block_no].is_block_value_constant))
            {
              if (fields[j].no_of_bits_per_num > 0)
              {
                if (fields[j].Huffman_global)
                {
                  uint32 bit;
                  int32 h_tmp;

                  sb_a_bit_stream.GetBits(bit, fields[j].Huffman_global->GetMinLen());
                  h_tmp = fields[j].Huffman_global->DecodeFast(bit);

                  while (h_tmp < 0)
                  {
                    sb_a_bit_stream.GetBit(bit);
                    h_tmp = fields[j].Huffman_global->Decode(bit);
                  };

                  num_val = h_tmp;
                }
                else
                {
                  sb_a_bit_stream.GetBits(num_val, fields[j].no_of_bits_per_num);
                }
              }
              else
              {
                num_val = 0;
              }
            }
            else
            {
              if (fields[j].is_delta_coding)
              {
                num_val = 0;
              }
              else
              {
                num_val = prev_value[j] - fields[j].min_value;
              }
            }

            if (fields[j].is_delta_coding)
            {
              num_val += prev_value[j] + fields[j].min_delta;
            }
            else
            {
              num_val += fields[j].min_value;
            }

            uchar tmp_str[10];
            int32 tmp_len = 0;
            tmp_len = utils::to_string(tmp_str, num_val);

            records[i].AppendData(tmp_str, tmp_len, 'T');
            prev_value[j] = num_val;

          }
          records[i].AppendChar(fields[j].sep, 'T');

          continue;
        }

        if ((i % DEFAULT_B_SIZE) > 0 && fields[j].block_desc[block_no].is_block_constant)
        {
          uchar *tmp_str = new uchar[fields[j].block_str_len]; 
          std::copy(records[first_rec_in_block].title+fields[j].block_str_start, 
            records[first_rec_in_block].title+fields[j].block_str_start+fields[j].block_str_len, tmp_str);

          records[i].AppendData(tmp_str, fields[j].block_str_len, 'T');
          records[i].AppendChar(fields[j].sep, 'T');
          delete[] tmp_str;

          continue;
        }

        uint32 field_len;
        if (!fields[j].is_len_constant)
        {
          sb_a_bit_stream.GetBits(field_len, fields[j].no_of_bits_per_len);      
          field_len += fields[j].min_len;
        }
        else
        {
          field_len = fields[j].len;
        }
        
        for (uint32 k = 0; k < field_len; ++k)
        {
          if (k < fields[j].len && fields[j].Ham_mask[k])
          {
            records[i].AppendChar(fields[j].data[k], 'T');
          }
          else
          {
            uint32 bit;
            int32 h_tmp;
            HuffmanEncoder *cur_huf = fields[j].Huffman_local[MIN(k, MAX_FIELD_STAT_LEN)]; 

            sb_a_bit_stream.GetBits(bit, cur_huf->GetMinLen());
            h_tmp = cur_huf->DecodeFast(bit);

            while (h_tmp < 0)
            {
              sb_a_bit_stream.GetBit(bit);
              h_tmp = cur_huf->Decode(bit);
            };

            tmp = h_tmp;
            records[i].AppendChar((uchar) tmp, 'T');
          }
        }
        if ((i % DEFAULT_B_SIZE) == 0 && fields[j].block_desc[block_no].is_block_constant)
        {
          fields[j].block_str_start = records[i].title_len - field_len;
          fields[j].block_str_len = field_len;
        }

        records[i].AppendChar(fields[j].sep, 'T');
      }

      records[i].prev_title_len = prev_title_len;
      prev_title_len += records[i].title_len;
    }
    sb_a_bit_stream.FlushInputWordBuffer();
  }
}

// --------------------------------------------------------------------------------------------
void FetchDNA(BitStream &sb_b_bit_stream, std::vector<uchar> &symbols, uint32 no_symbols, uint32 fastq_flags,
  Record *&records, int64 no_records, std::vector<uint32> &no_ambiguity)
{
  uint32 tmp;
  HuffmanEncoder *Huffman_sym = NULL;

  // DNA symbols
  symbols.resize(no_symbols);
  for (uint32 i = 0; i < no_symbols; ++i)
  {
    sb_b_bit_stream.GetByte(tmp);
    my_assert(tmp < 256);
    symbols[i] = (uchar) tmp;
  }

  if ((fastq_flags & FLAG_DNA_PLAIN) == 0)
  {
    Huffman_sym = new HuffmanEncoder((uint32) symbols.size());
    HuffmanEncoder::LoadTree(sb_b_bit_stream, *Huffman_sym); 
  }

  sb_b_bit_stream.FlushInputWordBuffer();

  for (uint32 i = 0; i < no_records; ++i)
  {
    records[i].seq_len = records[i].qua_len - no_ambiguity[i];
  }

  if ((fastq_flags & FLAG_DNA_PLAIN) != 0)
  {
    for (uint32 i = 0; i < no_records; ++i)
    {
      uint32 cur_sequence_len = records[i].seq_len;
      uchar *cur_sequence = new uchar[cur_sequence_len];

      for (uint32 j = 0; j < cur_sequence_len; ++j)
      {
        uint32 tmp(0);      // reduce warning...
        sb_b_bit_stream.Get2Bits(tmp);
        cur_sequence[j] = symbols[tmp];
      }
      //printf("if: %s", cur_sequence);
      //getchar();
      records[i].AppendData(cur_sequence, cur_sequence_len, 'D');
      records[i].AppendChar('\n','D');
      delete[] cur_sequence;
    }
  }
  // todo
  // seems to never go to else
  else
  {
    for (uint32 i = 0; i < no_records; ++i)
    {
      uint32 cur_sequence_len = records[i].seq_len;
      uchar *cur_sequence = new uchar[cur_sequence_len];
      for (uint32 j = 0; j < cur_sequence_len; ++j)
      {
        // Symbols
        uint32 bit;
        sb_b_bit_stream.GetBits(bit, Huffman_sym->GetMinLen());
        int32 h_tmp = Huffman_sym->DecodeFast(bit);
        while (h_tmp < 0)
        {
          sb_b_bit_stream.GetBit(bit);
          h_tmp = Huffman_sym->Decode(bit);
        };

        cur_sequence[j] = symbols[h_tmp];
      } 
      //printf("else: %s", cur_sequence);
      //getchar();
      records[i].AppendData(cur_sequence, cur_sequence_len, 'D');
      records[i].AppendChar('\n','D');
      delete[] cur_sequence;
    }
  }

  if (Huffman_sym)
    delete Huffman_sym;

  sb_b_bit_stream.FlushInputWordBuffer();
}


// --------------------------------------------------------------------------------------------
void FetchQuality(BitStream &sb_b_bit_stream, std::vector<uchar> &qualities, uint32 no_qualities, 
  Record *&records, int64 no_records, uint32 fastq_flags, std::vector<uint32> &no_ambiguity,
  uint32 max_quality_length, int64 &prev_qua_len)
{
  uint32 tmp;
  HuffmanEncoder* huf;
  std::vector<HuffmanEncoder*> Huffman_qua;

  // Quality symbols
  qualities.resize(no_qualities);
  for (uint32 i = 0; i < no_qualities; ++i)
  {
    sb_b_bit_stream.GetByte(tmp);
    my_assert(tmp < 256);
    qualities[i] = (uchar) tmp;
  }

  // Quality Huffman codes
  for (uint32 i = 0; i <= max_quality_length; ++i)
  {
    Huffman_qua.push_back(huf = new HuffmanEncoder((uint32) qualities.size()));
    HuffmanEncoder::LoadTree(sb_b_bit_stream, *huf);
  }

  sb_b_bit_stream.FlushInputWordBuffer();

  // Quality data
  no_ambiguity.resize(no_records);
  for (uint32 i = 0; i < no_records; ++i)
  {
    no_ambiguity[i] = 0;
    uint32 cur_qua_len = records[i].qua_len;
    uchar *cur_quality = new uchar[cur_qua_len];

    for (uint32 j = 0; j < cur_qua_len; ++j)
    {
      uint32 bit;
      int32 h_tmp;

      sb_b_bit_stream.GetBits(bit, Huffman_qua[j+1]->GetMinLen());
      h_tmp = Huffman_qua[j+1]->DecodeFast(bit);

      while (h_tmp < 0)
      {
        sb_b_bit_stream.GetBit(bit);
        h_tmp = Huffman_qua[j+1]->Decode(bit);
      };

      if ((cur_quality[j] = qualities[h_tmp]) >= 128)
      {
        no_ambiguity[i]++;
      }
    }

    // This for loop only works if you have QUALITY_PLAIN_TRUNC
    for (uint32 j = cur_qua_len; j < cur_qua_len; ++j)
    {
      cur_quality[j] = '#';
    }

    records[i].AppendData(cur_quality, cur_qua_len, 'Q');
    records[i].AppendChar('\n','Q');
    delete[] cur_quality;

    records[i].prev_seq_qua_len = prev_qua_len;
    prev_qua_len += ((records[i].qua_len + 1) * 2) + 2;
  }
  sb_b_bit_stream.FlushInputWordBuffer();


  if (huf)
    delete huf;
}

// --------------------------------------------------------------------------------------------
void MakeFooter(BitStream &footer_bit_stream, int32 g_size, int64 FASTQ_size, int32 n_blocks, 
  int32 n_subblocks, std::multimap<double, int32> &blocks_order, uint32 *lb_sizes, uchar *b_data_buffer, 
  int32 b_data_size)  
{
  // Bits Encoding Processes Size
  int32 BEPS = floor(log2(g_size)) + 1;
  // Bits Encoding File Size
  int32 BEFS = floor(log2(FASTQ_size)) + 1;           
  // Bits Encoding Block Size
  int32 BEBS = floor(log2(n_blocks)) + 1;    
  // Bits Encoding SubBlock Size
  int32 BESS = floor(log2(n_subblocks)) + 1;     
  // Bits Encoding Last Blocks Size
  int32 BELB = floor(log2(*std::max_element(lb_sizes, lb_sizes+g_size))) + 1;
  // Are Last Blocks of Equal Size?
  const bool LBES = *std::max_element(lb_sizes, lb_sizes+g_size) != *std::min_element(lb_sizes, lb_sizes+g_size) ? 0 : 1; 
  // Create Footer buffer
  footer_bit_stream.Create(0);

  footer_bit_stream.PutBits(BEPS,4);
  footer_bit_stream.PutBits(BEFS,6);
  footer_bit_stream.PutBits(BEBS,5);
  footer_bit_stream.PutBits(BESS,5);
  footer_bit_stream.PutBits(BELB,5);
  footer_bit_stream.PutBit(LBES);

  footer_bit_stream.PutBits(g_size, BEPS);

  // When dealing with files bigger than 4Gb
  if (BEFS > 32)
  {
    footer_bit_stream.PutBits(FASTQ_size >> 32, BEFS-32);
    footer_bit_stream.PutBits(FASTQ_size & 0xFFFFFFFF, 32);
  } 
  else
  {
    footer_bit_stream.PutBits(FASTQ_size, BEFS);
  }

  footer_bit_stream.PutBits(n_blocks, BEBS);
  footer_bit_stream.PutBits(n_subblocks, BESS);

  // Recalculate to know the minimun amount of bits to encode the Processes ranks
  // from P0 to P(p-1)
  BEPS = ceil(log2(g_size));
  // Bitpack block order
  for (auto& bo: blocks_order)
  {
    footer_bit_stream.PutBits(bo.second, BEPS);
  }
  // Bitpack size of last blocks
  if (LBES == 0)
  {
    for (int32 i = 0; i < g_size; ++i)
    {
      footer_bit_stream.PutBits(lb_sizes[i], BELB);
    }
  }

  footer_bit_stream.FlushPartialWordBuffer();
  // d_data_buffer is already bitpacked, so add to the footer
  footer_bit_stream.PutBytes(b_data_buffer, b_data_size);
  // Round up to the next complete byte (< than 8 bits wasted)
  footer_bit_stream.FlushPartialWordBuffer();
  // Add footer size at the end (2bytes)
  footer_bit_stream.Put2Bytes(footer_bit_stream.GetIO_Buffer_Pos());
}

// // --------------------------------------------------------------------------------------------
// void MakeHeader(BitStream &header_bit_stream, BlockHeader &p_block_header)  
// {
//   // Create header buffer
//   header_bit_stream.Create(0);
//   // Number of subblocks in this block
//   p_block_header.NOSB = p_block_header.SBOL.size();

//   header_bit_stream.PutBits(p_block_header.WRID, p_block_header.BEWR);
//   header_bit_stream.PutBits(p_block_header.BHS, 12);
//   header_bit_stream.PutBits(p_block_header.NOSB, 6);
//   header_bit_stream.PutBits(p_block_header.BESO, 5);
//   header_bit_stream.PutBits(p_block_header.BCSS, 2);

//   // printf("WRID=%03d::BHS=%d::NOSB=%d::BESO=%d::BCSS=%d\n", p_block_header.WRID, p_block_header.BHS, p_block_header.NOSB, p_block_header.BESO, p_block_header.BCSS);

//   // Bitpack subblock's offsets
//   for (int32 i = 0, j = 0; i < p_block_header.NOSB; ++i)
//   {
//     header_bit_stream.PutBits(p_block_header.SBOL.at(i), p_block_header.BESO);
//     // printf("(%03d)SBOL[%d]=%d\n", p_block_header.WRID, i, p_block_header.SBOL.at(i));
//     if (i == 0 && p_block_header.BCSS > 1)
//       continue;

//     // printf("(%03d)QSPS[%d]=%d\n", p_block_header.WRID, j, p_block_header.QSPS.at(j));
//     header_bit_stream.PutBits(p_block_header.QSPS.at(j++), p_block_header.BESO);
//   }
//   header_bit_stream.FlushPartialWordBuffer();
// }

// --------------------------------------------------------------------------------------------
void MakeHeader(BitStream &header_bit_stream, std::vector<BlockHeader> &p_blocks_data)
{
  // Create header buffer
  header_bit_stream.Create(0);

  for (std::vector<BlockHeader>::iterator b = p_blocks_data.begin(); b != p_blocks_data.end(); ++b)
  {
    BlockHeader sb = *b;
    int32 no_sb    = sb.SBOL.size();
    int32 BESO     = 22; // Make this into a Global Variable

    // Put the number of subblocks in the bitstream
    // Encode with 6 bits (up to 64 sb's)
    header_bit_stream.PutBits(no_sb, 6);
    // Block contains splitted subblocks?
    header_bit_stream.PutBits(sb.BCSS, 2);

    // Bitpack subblock's offsets
    for (int32 i = 0; i < no_sb; ++i)
      header_bit_stream.PutBits(sb.SBOL.at(i), BESO);
  }
  header_bit_stream.FlushPartialWordBuffer();
}

// --------------------------------------------------------------------------------------------
void ReadFooter(BitStream &footer_bit_stream, Footer &footer, std::vector<SubBlock> &p_subblocks, 
  MPI_File &input_NGSC, int32 g_size, int32 p_rank)  
{
  uint32 bits;

  //Read BEPS
  footer_bit_stream.GetBits(bits, 4);
  footer.BEPS = bits;
  //Read BEFS
  footer_bit_stream.GetBits(bits, 6);
  footer.BEFS = bits;
  //Read BEBS
  footer_bit_stream.GetBits(bits, 5);
  footer.BEBS = bits;
  //Read BESS
  footer_bit_stream.GetBits(bits, 5);
  footer.BESS = bits;
  //Read BELB
  footer_bit_stream.GetBits(bits, 5);
  footer.BELB = bits;
  //Read LBES
  footer_bit_stream.GetBits(bits, 1);
  footer.LBES = bits;

  // Read PS: # of P's used to compressed
  footer_bit_stream.GetBits(bits, footer.BEPS);
  footer.PS = bits;

  // Read FS: FASTQ file size
  if (footer.BEFS > 32)
  {
    footer_bit_stream.GetBits(bits, footer.BEFS - 32);
    footer.FS = bits;
    footer_bit_stream.GetBits(bits, 32);
    footer.FS = (footer.FS << 32) + bits;
  } 
  else
  {
    footer_bit_stream.GetBits(bits, footer.BEFS);
    footer.FS = bits;
  }

  // Read BS: # of Compressed Blocks in NGSC file
  footer_bit_stream.GetBits(bits, footer.BEBS);
  footer.BS = bits;
  // Read SS: # of Subblocks in NGSC file
  footer_bit_stream.GetBits(bits, footer.BESS);
  footer.SS = bits;

  // Minimun bits used to encode each CBO
  int32 CBO_bits = ceil(log2(footer.PS));
  // Resize block's counter (counts blocks per BS)
  footer.BC.resize(footer.PS, 0);
  // Auxiliary variable for indexes of CBO
  std::vector<std::vector<int32> > indexCBO(footer.BS);

  // Get array of blocks order
  for (int32 i = 0; i < footer.BS; ++i)
  {
    footer_bit_stream.GetBits(bits, CBO_bits);
    footer.CBO.push_back(bits);
    footer.BC[bits]++;
    // indexCBO.at(bits).push_back(i);
  }

  // Get last blocks sizes if footer.LBES == 0
  if (footer.LBES == 0)
  {
    for (int32 i = 0; i < footer.PS; ++i)
    {
      footer_bit_stream.GetBits(bits, footer.BELB);
      footer.LBS.push_back(bits); 
    }
  }

  footer.ABS.resize(footer.BS);
  uint64 pos  = 0;
  char *buff  = new char[2];
  std::vector<int32> tmp_BC(footer.BC);

  // Check in case clocks were not synch
  for (int32 i = 0; i < footer.BS; ++i)
  {
    MPI_File_read_at(input_NGSC, pos, buff, 2, MPI_CHAR, MPI_STATUS_IGNORE);

    int32 WRID = buff[0];
    WRID = (WRID << 8) + buff[1];

    footer.ABS.at(i) = pos;

    // Correct block's order
    if (WRID != footer.CBO[i])      
      footer.CBO[i] = WRID;

    // Get correct index
    indexCBO.at(footer.CBO[i]).push_back(i);

    // Correct all block sizes
    if (tmp_BC[footer.CBO[i]] > 1)
      pos += 8 << 20;
    else 
      pos += footer.LBS[footer.CBO[i]];

    tmp_BC[footer.CBO[i]]--;
  }

  delete[] buff;

  // Devide subblocks equally among Ps
  int32 p_no_sb  = footer.SS / g_size;
  int32 sb_count = 0;
  // This part will make sure that the remainder sb
  // get distributed equally
  int32 remain = footer.SS % g_size;
  int32 plus_s = 0, plus_e = 0;
  if (remain && p_rank >= (g_size - remain) && p_no_sb)
  {
    plus_s = p_rank - (g_size - remain);
    plus_e = 1;
  }

  if (!p_no_sb)
    p_no_sb = 1;

  int32 p_sb_start = p_rank * p_no_sb + plus_s;
  int32 p_sb_end   = p_sb_start + p_no_sb - 1 + plus_e;

  // printf("R::%d S::%d E::%d SS(%d)\n", p_rank, p_sb_start, p_sb_end, footer.SS );

  footer_bit_stream.FlushInputWordBuffer();

  for (int32 p = 0; p < footer.PS; ++p)
  {
    for (int b = 0; b < footer.BC[p]; ++b)
    {
      footer_bit_stream.GetBits(bits, 6);
      int32 no_sb = bits;
      footer_bit_stream.GetBits(bits, 2);
      int32 BCSS  = bits;

      int32 idx         = indexCBO.at(p).at(b);
      uint64 global_pos = footer.ABS.at(idx) + 2;

      for (int sb = 0; sb < no_sb; ++sb)
      {
        footer_bit_stream.GetBits(bits, 22);
    
        if (sb_count >= p_sb_start && sb_count <= p_sb_end)
        {
          if (sb == 0 && (BCSS & 2))
          {
            p_subblocks.back().sb_cont_pos    = global_pos;
            p_subblocks.back().sb_cont_length = bits;
            p_subblocks.back().is_splitted    = 1;
          } 
          else
          {
            SubBlock curr_sb;
            curr_sb.wr_id        = p;
            curr_sb.sb_start_pos = global_pos; 
            curr_sb.sb_length    = bits;
            curr_sb.is_splitted  = 0;

            p_subblocks.push_back(curr_sb);
          }
        }

        global_pos += bits; 

        // If readding a splitted subblock skip incrementing counter
        if (sb == no_sb-1 && (BCSS & 1))
          continue;

        ++sb_count;

        if (sb_count > p_sb_end)
          return;
      }
    }
    footer_bit_stream.FlushInputWordBuffer();
  }
}