#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "gun.h"
#include "misc.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INITIALIZE SHOTGUN
//
SHOTGUN *CreateShotgun(uint32_t m, uint32_t l, uint32_t s){
  uint32_t a, b;
  SHOTGUN *G = (SHOTGUN *) Calloc(1, sizeof(SHOTGUN));
  G->sym = (uint32_t *) Calloc(l, sizeof(uint32_t));
  G->bits = (uint32_t *) Calloc(m, sizeof(uint32_t));
  G->freqs = (uint32_t ***) Calloc(m, sizeof(uint32_t **));
  for(a=0 ; a<m ; ++a){
    G->freqs[a] = (uint32_t **) Calloc(l, sizeof(uint32_t *));
    for(b=0 ; b<l ; ++b)
      G->freqs[a][b] = (uint32_t *) Calloc(s+1, sizeof(uint32_t));
    }
  return G;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// DELETE SHOTGUN
//
void DeleteShotgun(SHOTGUN *G, uint32_t m, uint32_t l, uint32_t s){
  uint32_t a, b;
  Free(G->sym, l * sizeof(uint32_t));
  Free(G->bits, m * sizeof(uint32_t));
  for(a=0 ; a<m ; ++a){
    for(b=0 ; b<l ; ++b)
      Free(G->freqs[a][b], (s+1) * sizeof(uint32_t));
    Free(G->freqs[a], l * sizeof(uint32_t *));
    }
  Free(G->freqs, m * sizeof(uint32_t **));
  Free(G, sizeof(SHOTGUN));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CALCULATE THE BEST IN GUN
//
uint32_t BestInGun(uint32_t *b, uint32_t n){
  uint32_t x, min = b[0], i = 0;
  for(x = 1 ; x < n ; ++x){
    if(min > b[x]){
      min = b[x];
      i = x;
      }
    b[x] = 0;
    }
  return i;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPUTE SHOTGUN PROBABILITIES LOG
//
uint32_t CompGunProbs(uint32_t *f, uint32_t s){
  return Log(f[4]/f[s]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

