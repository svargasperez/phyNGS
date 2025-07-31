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
#include <sys/stat.h>
#include <string>
#include <fstream>
#include "defs.h"

using std::string;


bool debug_status = false;
const string PERFORMANCE_PREFIX = "performanceOutput-";


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
    // Return if debug mode is inactive
    if (debug_status == false)
        return;

    // Handle a variable number of arguments and print
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}

void write_performance_file(double runtime, const string mode, int32 g_size, int32 no_threads, const string pattern)
{
    string folder_name = PERFORMANCE_PREFIX + mode;

    // Creates output directory if needed
    if (mkdir(folder_name.c_str(), 0777) != -1)
    {
        fflush(stdout);
        printf("Performance output directory made.\n");
    }

    // Construct file path using folder, process count, thread count, and possibly pattern
    string file_path = "./" + folder_name + "/" + std::to_string(g_size) + "," + std::to_string(no_threads);

    if (pattern != "") // Use pattern if string is set
        file_path += "," + pattern;

    file_path += ".txt";

    std::ofstream test_file;
    test_file.open(file_path);
    test_file << runtime;
}