#include <stdio.h>
#include <stdlib.h>
#include "sfcm.h"
#include "arith_aux.h"
#include "arith.h"
#include "bitio.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SFCM *CreateSFCM(int size, int nSym){
  SFCM *M  = (SFCM *) Calloc(1, sizeof(SFCM));
  M->array = (AP *)   Calloc(size<<6, sizeof(AP));
  M->freqs = (int *)  Calloc(nSym, sizeof(int));
  M->size  = size;
  M->n     = nSym;
  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void EncodeSFCM(FILE *F, SFCM *M, int last, int q){
  int n;
  uint32_t total = 0;

  AP *a = &M->array[last<<6];

  for(n = 0 ; n < M->n ; ++n){
    M->freqs[n] = 8 * a[n] + 1;    //INITIALIZE TO 1 !! YOU WONT NEED TO DO THIS
    total      += M->freqs[n];
    }

  AESym(q, (int *) M->freqs, total, F);

  if(++a[q] == MAXAP){
    for(n = 0 ; n < M->n ; ++n)
      a[n] -= a[n]>>1;    //INITIALIZE TO 1 !! YOU WONT NEED TO DO THIS
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void DeleteSFCM(SFCM *M){
  Free(M->array, M->size * sizeof(AP));
  Free(M->freqs, M->n    * sizeof(int));
  Free(M, sizeof(SFCM));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

