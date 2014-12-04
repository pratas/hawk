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
  B->array = (BCC *) Calloc(B->size * B->nSym, sizeof(BCC));
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

BCC *SearchBloom(BLOOM *B, uint64_t i){
  uint32_t s, n, k = 0, min = MAX_BCC, sum[B->H->k];
  BCC *bc;
  for(n = 0 ; n < B->H->k ; ++n){
    bc = &B->array[(Hash(B->H, i, n) % B->size) * B->nSym];
    sum[n] = 0;
    for(s = 0 ; s < B->nSym ; ++s) 
      sum[n] += bc[s];
    if(sum[n] < min){
      min = sum[n];
      k = n;
      }
    }
  return &B->array[(Hash(B->H, i, k) % B->size) * B->nSym]; // TODO: computed...
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateBloom(BLOOM *B, uint64_t i, uint8_t s){
  uint32_t n, k;
  for(n = 0 ; n < B->H->k ; ++n){
    BCC *bc = &B->array[(Hash(B->H, i, n) % B->size) * B->nSym];
    if(++bc[s] == MAX_BCC){
      //fprintf(stderr, "Renormalizing Bloom counters\n");
      for(k = B->nSym ; k-- ; )
        bc[k] >>= 1;
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

