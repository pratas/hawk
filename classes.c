#include <stdio.h>
#include <stdlib.h>
#include "classes.h"
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

void FreeClasses(CLASSES *C){
  Free(C, sizeof(CLASSES));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

