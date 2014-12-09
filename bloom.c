#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "bloom.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

BLOOM *CreateBloom(uint32_t k, uint32_t s, uint32_t n){
  BLOOM *B = (BLOOM *) Calloc(1, sizeof(BLOOM));
  B->size  = s;
  B->nSym  = n;
  B->array = (BCC *) Calloc(s, sizeof(BCC));
  B->freqs = (uint32_t *) Calloc(B->nSym, sizeof(uint32_t)); 
  B->H     = CreateHFamily(k, 2147483647); // 2147483647 => LARGE PRIME (2^30)-1
  return B;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void DeleteBloom(BLOOM *B){
  DeleteHFamily(B->H);
  Free(B->array, B->size * sizeof(BCC));
  Free(B, sizeof(BLOOM));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void SearchBloom(BLOOM *B, uint64_t i){
  uint32_t s, n, x, min;
  i = i>>2<<2; // i &= 0xfffffffc [FIXME: THIS ONLY WORKS FOR 4 SYMBOLS]
  for(s = 0 ; s < B->nSym ; ++s){
    x = Hash(B->H, i, 0) % B->size;
    min = B->array[x];
    for(n = 1 ; n < B->H->k ; ++n){
      x = Hash(B->H, i, n) % B->size;
      if(B->array[x] < min)
        min = B->array[x];
      }
    B->freqs[s] = min;
    ++i;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateBloom(BLOOM *B, uint64_t i){
  uint32_t n, k, x;
  for(n = 0 ; n < B->H->k ; ++n){
    x = Hash(B->H, i, n) % B->size;
    if(B->array[x]++ == MAX_BCC){
      //fprintf(stderr, "Renormalizing Bloom counters\n");
      for(k = B->size ; k-- ; )
        B->array[k] >>= 1;
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

