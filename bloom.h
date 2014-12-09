#ifndef BLOOM_H_INCLUDED
#define BLOOM_H_INCLUDED

#include "defs.h"
#include "hash.h"

typedef uint32_t    BCC;   // BLOOM COUNTERS PRECISION
#define N_HASH_FUNC 3
#define BLOOM_SIZE  1<<25 //1193711 //16777259 // NEXT PRIME AFTER 24 BITS
#define MAX_BCC     (((uint64_t)1<<(sizeof(BCC)*8))-1)   // MAX BLOOM COUNTER

typedef struct{
  BCC      *array;      // ARRAY COUNTERS FOR BLOOM
  uint32_t size;        // SIZE OF THE ARRAY FOR BLOOM COUNTERS
  uint32_t nSym;        // NUMBER OF SYMBOLS
  uint32_t *freqs;      // BLOOM FREQS
  HFAM     *H;          // HASH FAMILY AGREGGATED TO BLOOM
  }
BLOOM;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

BLOOM     *CreateBloom    (uint32_t, uint32_t, uint32_t);
void      SearchBloom     (BLOOM *, uint64_t);
void      DeleteBloom     (BLOOM *);
void      UpdateBloom     (BLOOM *, uint64_t);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

