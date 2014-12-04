#ifndef PARAM_H_INCLUDED
#define PARAM_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESSION LEVELS FOR HEADERS(-H), DNA (-D) and SCORES (-S). 
// USAGE: DNA (-D) CONTEXT DEPTHS MUST BE ORDERED [ASCENDING].
//
#define LEVEL_0                                           \
  "-H 1 "                                                 \
  "-D 12 "                                                \
  "-S 1"
#define LEVEL_1                                           \
  "-H 1 "                                                 \
  "-D 4 -D 13 "                                           \
  "-S 1"
#define LEVEL_2                                           \
  "-H 1 "                                                 \
  "-D 4 -D 14 "                                           \
  "-S 2"        
#define LEVEL_3                                           \
  "-H 2 "                                                 \
  "-D 4 -D 9 -D 13 "                                      \
  "-S 2"  
#define LEVEL_4                                           \
  "-H 2 "                                                 \
  "-D 4 -D 7 -D 11 -D 14 "                                \
  "-S 3"
#define LEVEL_5                                           \
  "-H 2 "                                                 \
  "-D 4 -D 16 "                                           \
  "-S 3"
#define LEVEL_6                                           \
  "-H 3 "                                                 \
  "-D 4 -D 8 -D 12 -D 16 "                                \
  "-S 3"
#define LEVEL_7                                           \
  "-H 3 "                                                 \
  "-D 4 -D 9 -D 13 -D 18 "                                \
  "-S 3" 
#define LEVEL_8                                           \
  "-H 3 "                                                 \
  "-D 4 -D 9 -D 13 "                                      \
  "-S 3"
#define LEVEL_9                                           \
  "-H 3 "                                                 \
  "-D 3 -D 5 -D 7 -D 9 -D 11 -D 13 -D 18 "                \
  "-S 4" 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// BASIC INPUT PARAMETERS
//
typedef struct{
  uint8_t verbose;     // GIVES MORE INFORMATION
  uint8_t force;       // FORCE REWRITE OVER FILES
  uint8_t level;       // COMPRESSION LEVEL
  uint8_t fLow;        // LOWER CONTEXT ORDER FOR FILTER 
  uint8_t fHigh;       // HIGH CONTEXT ORDER FOR FILTER
  uint8_t filter;      // USE FILTER
  uint8_t inverse;     // USE INVERSIONS
  uint8_t adjust;      // ADJUST MODELS TO DATA [AFTER FAST ANALYSIS]
  uint8_t mode;        // USE HASH OR BLOOM FOR DEEP CONTEXT MODELS
  }
PARAM;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void    Msg        (PARAM *, char *);
char    *GetParam  (uint8_t);
void    PrintParam (void);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

