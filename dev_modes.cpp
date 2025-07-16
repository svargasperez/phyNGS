// #include "defs.h"
#include <stdarg.h>
#include <stdio.h>

// Create enum to represent different modes
enum class Mode { DEBUG = 0, TEST };

// Array represents whether modes are enabled, indexed by (int)Mode
bool mode_statuses[2] = {false, false};


bool get_mode_status(Mode mode)
{
    return mode_statuses[(int)mode];
}

void set_mode_status(Mode mode, bool status)
{
    mode_statuses[(int)mode] = status;
}

void mode_print(Mode mode, const char * format, ...)
{
    // Return if mode is not active
    if (mode_statuses[(int)mode] == false)
        return;

    // Handle a variable number of arguments and print
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}
