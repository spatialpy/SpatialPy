#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#if DEBUG_LEVEL > 0
    // http://support.raisonance.com/content/how-remove-file-path
    #define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
    #define SSA_LOG(level, fmt, ...) \
        do { if (DEBUG_LEVEL >= level) \
            printf("%s %s - ssa-sdpd.%s - %s - %i - " fmt, __DATE__, __TIME__, __FILENAME__, ssa_level(level), __LINE__, ##__VA_ARGS__); } while(0)
#else
    #define SSA_LOG(level, fmt, ...)
#endif

const char *ssa_level(int level);
