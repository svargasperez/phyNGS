/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
  Debug Author: Juniper Pasternak
*/

#ifndef _DEV_MODES_H
#define _DEV_MODES_H

// --------------------------------------------------------------------------------------------
bool get_debug();

// --------------------------------------------------------------------------------------------
void debug_on();

// --------------------------------------------------------------------------------------------
void debug_off();

// --------------------------------------------------------------------------------------------
void debug_print(const char *format, ...);


#endif