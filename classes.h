#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#include "defs.h"
#include "models.h"

typedef struct{
  HEADERS  H;         // HEADERS POINTER
  DNA      D;         // DNA POINTER
  SCORES   S;         // SCORES POINTER
  char     *f;        // FILE NAME FOR COMPRESS/UNCOMPRESS
  uint64_t nReads;    // NUMBER OF FASTQ READS
  uint64_t length;    // NUMBER OF CHARS IN THE FILE
  }
CLASSES;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CLASSES  *InitClasses  (void);
void     FreeClasses   (CLASSES *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

