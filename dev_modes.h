/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
  Development Modes Author: Juniper Pasternak
*/

#ifndef _DEV_MODES_H
#define _DEV_MODES_H

#include "defs.h"

enum class Mode { DEBUG = 0, TEST };

// --------------------------------------------------------------------------------------------
bool get_mode_status(Mode mode);

// --------------------------------------------------------------------------------------------
void set_mode_status(Mode mode, bool status);

// --------------------------------------------------------------------------------------------
void mode_print(Mode mode, const char * format, ...);


#endif