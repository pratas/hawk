#ifndef GUN_H_INCLUDED
#define GUN_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// SHOTGUN SCHEME [STORE COMPETING FEATURES]

typedef struct{
  uint32_t ***freqs;
  uint32_t *sym;
  uint32_t *bits;
  }
SHOTGUN;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SHOTGUN   *CreateShotgun  (uint32_t, uint32_t, uint32_t);
void      DeleteShotgun   (SHOTGUN *, uint32_t, uint32_t, uint32_t);
uint32_t  BestInGun       (uint32_t *, uint32_t);
uint32_t  CompGunProbs    (uint32_t *, uint32_t);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

