#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

#include "defs.h"

#define MAX_STR   4096
#define MPW       1072632447

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char     *ReplaceSubStr (char *, char *, char *);
uint32_t Log            (uint64_t);
double   PW             (double, double);
FILE     *Fopen         (const char *, const char *);
void     CheckFile      (uint8_t, char *);
uint64_t FNBytes        (FILE *);
void     Progress       (uint64_t, uint64_t);
void     ReverseStr     (uint8_t *, uint8_t *, uint32_t);
uint8_t  Pack8bits      (uint8_t *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

