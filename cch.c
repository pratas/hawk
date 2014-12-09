#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cch.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CCH *CreateCCH(uint32_t s, uint32_t n){
  CCH   *C = (CCH *) Calloc(1, sizeof(CCH));
  C->size  = s;
  C->nSym  = n;
  C->array = (C_CCH *) Calloc(C->size * C->nSym, sizeof(C_CCH));
  C->H     = CreateHFamily(1, 2147483647); // 2147483647 => LARGE PRIME (2^30)-1
  return C;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void DeleteCCH(CCH *C){
  DeleteHFamily(C->H);
  Free(C->array, C->size * sizeof(C_CCH));
  Free(C, sizeof(CCH));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

C_CCH *SearchCCH(CCH *C, uint64_t i){
  return &C->array[(Hash(C->H, i, 0) % C->size) * C->nSym]; 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateCCH(CCH *C, uint64_t i, uint8_t s){
  uint32_t k;
  C_CCH *x = &C->array[(Hash(C->H, i, 0) % C->size) * C->nSym];
  if(++x[s] == MAX_CCH){
    for(k = C->nSym ; k-- ; )
      x[k] >>= 1;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

