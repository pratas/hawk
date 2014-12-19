#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#include "defs.h"
#include "param.h"
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

CLASSES  *InitClasses    (void);
void     InitAlphabets   (CLASSES *);
void     FreeAlphabets   (CLASSES *);
void     PrintStreamInfo (CLASSES *);
void     ParseFile       (CLASSES *, FILE *, PARAM *);
void     CreateModels    (CLASSES *, char *, FILE *, PARAM *);
void     FreeModels      (CLASSES *, PARAM *);
void     SetValues       (CLASSES *, PARAM *);  
void     CreateAuxStates (CLASSES *);
void     DeleteAuxStates (CLASSES *);
void     FreeClasses     (CLASSES *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

