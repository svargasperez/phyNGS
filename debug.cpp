/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
  Debug Author: Juniper Pasternak
*/

// #include "defs.h"
#include <stdarg.h>
#include <stdio.h>

bool debug_status = false;

bool get_debug()
{
    return debug_status;
}

void debug_on()
{
    debug_status = true;
}

void debug_off()
{
    debug_status = false;
}

void debug_print(const char *format, ...)
{
    // Return if debug mode is active
    if (debug_status == false)
        return;

    // Handle a variable number of arguments and print
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}
