/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
  Incompresso Implementation Author: Cole Koryto
*/

#ifndef _INCOMPRESSO_H
#define _INCOMPRESSO_H

#include "defs.h"

// --------------------------------------------------------------------------------------------
bool FindFirst(const char * in_file, int32 g_size, int32 p_rank, char* pat, bool to_print, int32 no_threads);

// --------------------------------------------------------------------------------------------
void FindAll(const char * in_file, int32 g_size, int32 p_rank, char* pat, bool to_print, int32 no_threads);

// --------------------------------------------------------------------------------------------
int KMPSearch(char* pat, char* txt, int* lps);

// --------------------------------------------------------------------------------------------
void ComputeLPSArray(char* pat, int M, int* lps);

// --------------------------------------------------------------------------------------------
void ToFASTA (const char * in_file, const char * out_file, int32 g_size, int32 p_rank, int32 no_threads);


#endif