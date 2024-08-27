/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#ifndef _COMPRESSION_H
#define _COMPRESSION_H

#include "defs.h"

// --------------------------------------------------------------------------------------------
void CompressData(const char * in_file, const char * out_file, int32 g_size, int32 p_rank, int32 no_threads);

#endif