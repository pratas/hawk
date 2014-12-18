#ifndef SFCM_H_INCLUDED
#define SFCM_H_INCLUDED

#include "defs.h"
#include "mem.h"

typedef uint16_t AP;
#define MAXAP    (((uint64_t)1<<(sizeof(AP)*8))-1)

typedef struct{
  AP  *array;
  int *freqs;
  int size;
  int n;
  }
SFCM;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SFCM   *CreateSFCM  (int, int);
void   EncodeSFCM   (FILE *, SFCM *, int, int);
void   DeleteSFCM   (SFCM *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
