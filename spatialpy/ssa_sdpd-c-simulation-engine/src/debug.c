#include "debug.h"

const char *ssa_level(int level){
    switch (level){
        case 1:
            return "INFO";
        case 2:
            return "DEBUG";
        default:
            return "UNKNOWN LEVEL";
    }
}
