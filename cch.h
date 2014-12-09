#ifndef CCH_H_INCLUDED
#define CCH_H_INCLUDED

#include "defs.h"
#include "hash.h"

typedef uint8_t C_CCH;    // CCH COUNTERS PRECISION
#define CCH_SIZE 67108879 // 1<<27 //1193711 //16777259 
#define MAX_CCH   (((uint64_t)1<<(sizeof(C_CCH)*8))-1)   // MAX CCH COUNTER

typedef struct{
  C_CCH    *array;      // ARRAY COUNTERS FOR CCH, LENGTH = SIZE * NSYM
  uint32_t size;        // SIZE OF THE ARRAY FOR CCH COUNTERS
  uint32_t nSym;   
  HFAM     *H;          // HASH FAMILY AGREGGATED TO CCH
  }
CCH;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CCH       *CreateCCH   (uint32_t, uint32_t);
C_CCH     *SearchCCH   (CCH *, uint64_t);
void      DeleteCCH    (CCH *);
void      UpdateCCH    (CCH *, uint64_t, uint8_t);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

