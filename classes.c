#include <stdio.h>
#include <stdlib.h>
#include "classes.h"
#include "param.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INITIALIZE AND FREE CLASSES
//
CLASSES *InitClasses(void){
  CLASSES *C = (CLASSES *) Calloc(1, sizeof(CLASSES));
  C->nReads  = 0;
  C->length  = 0;
  return C;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint32_t *CreateStates(uint32_t n, uint32_t m){
  uint32_t k;
  uint32_t *states = (uint32_t *) Calloc(m+GUARD, sizeof(uint32_t));
  for(k = 0 ; k < m ; ++k)
    states[k] = k * n / m;  // SET STATES MASK    
  return states;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void SetValues(CLASSES *C, PARAM *P){
  C->H.extra      = 1; // 1=>EXISTS EXTRA HEADER (ALGORITHM FASTER W/ INVERSION)
  C->H.maxLine    = 1;
  C->H.nStates    = P->hNStates;
  C->S.nStates    = P->sNStates;
  C->S.maxLine    = 1;
  C->S.dynamic    = 0; // 1=>IS DYNAMIC
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CreateAuxStates(CLASSES *C){
  C->H.states = CreateStates(C->H.nStates, C->H.maxLine);
  C->S.states = CreateStates(C->S.nStates, C->S.maxLine);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void DeleteAuxStates(CLASSES *C){
  Free(C->H.states, (C->H.maxLine+GUARD) * sizeof(uint32_t));
  Free(C->S.states, (C->S.maxLine+GUARD) * sizeof(uint32_t));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FreeClasses(CLASSES *C){
  Free(C, sizeof(CLASSES));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

