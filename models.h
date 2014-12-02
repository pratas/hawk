#ifndef MODELS_H_INCLUDED
#define MODELS_H_INCLUDED

#include <string.h>
#include "defs.h"
#include "param.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// MODELS SCHEME
//
typedef uint16_t ACC;        // SIZE OF COUNTERS FOR COUNTER TABLE [8|16|32]
typedef uint8_t  HCC;        // SIZE OF COUNTERS FOR HASH TABLE
typedef uint16_t ENTMAX;     // ENTRY SIZE (NKEYS FOR EACH HINDEX)
typedef HCC      HCCs[4];    // DEFAULT COUNTERS
#define MAXACC_C (((uint64_t)1<<(sizeof(ACC)*8))-1)
#define MAXHCC_C (((uint64_t)1<<(sizeof(HCC)*8))-1)
#define MAXHCC_H 3
#define MAXHSIZE (((uint64_t)1<<(sizeof(ENTMAX)*8))-1)

typedef struct{       // ENTRY FOR 4 SYMBOLS
  uint16_t key;       // THE KEY (INDEX / HASHSIZE) STORED IN THIS ENTRY
  HCC      cnts;      // SMALL COUNTERS [2 BITS P/ SYMBOL]
  }
ENTRY;

typedef struct{       // ENTRY FOR MULTIPLE FILES
  uint32_t key;       // THE KEY (INDEX / HASHSIZE) STORED IN THIS ENTRY
  HCC      *cnts;     // COUNTERS
  }
GENTRY;

typedef struct{       // HASH FOR 4 SYMBOLS
  ENTMAX   *size;     // NUMBER OF KEYS FOR EACH ENTRY
  ENTRY    **ent;     // THE HEADS OF THE HASH TABLE LISTS
  HCCs     **cnts;    // THE CONTEXT COUNTERS
  }
HASH;

typedef struct{       // HASH FOR MULTIPLE SYMBOLS
  ENTMAX   *size;     // NUMBER OF KEYS FOR EACH ENTRY
  GENTRY   **ent;     // THE HEADS OF THE HASH TABLE LISTS
  HCCs     **cnts;    // THE CONTEXT COUNTERS
  uint64_t sum;       // TOTAL NUMBER OF CNTS
  }
GHASH;

typedef struct{
  ACC      *cnts;     // TABLE COUNTERS
  }
ARRAY;

typedef struct{
  uint64_t nPMod;     // MAXIMUM NUMBER OF PROBABILITY MODELS
  uint64_t mult;      // MULTIPLICATOR TO CALCULATE INDEX
  uint64_t idx;       // CURRENT CONTEXT INDEX
  uint64_t idxRev;    // CURRENT INVERTED REPEAT CONTEXT INDEX
  uint32_t aDen;      // ESTIMATOR DENOMINATOR
  uint32_t ctx;       // CONTEXT ORDER FOR THE FCM
  uint32_t *freqs;    // FCM SYMBOL PROBABILITIES
  uint8_t  rev;       // INVERSIONS USAGE
  uint8_t  nSym;      // FCM NUMBER OF SYMBOLS
  uint8_t  mode;      // USING HASH-TABLES OR NOT [COUNTER=0]
  ARRAY    A;         // COUNTER-TABLE ARRAY
  HASH     H;         // HASH-TABLE LINK
  }
FCM;

typedef struct{
  uint64_t nPMod;     // MAXIMUM NUMBER OF PROBABILITY MODELS
  uint64_t mult;      // MULTIPLICATOR TO CALCULATE INDEX
  uint64_t idx;       // CURRENT CONTEXT INDEX
  uint32_t aDen;      // ESTIMATOR DENOMINATOR
  uint32_t ctx;       // CONTEXT ORDER FOR THE FCM
  uint32_t *freqs;    // FCM SYMBOL PROBABILITIES
  uint8_t  nSym;      // FCM NUMBER OF SYMBOLS
  uint8_t  mode;      // USING HASH-TABLES OR NOT [COUNTER=0]
  ARRAY    A;         // COUNTER-TABLE ARRAY
  GHASH    H;         // GENERAL HASH-TABLE LINK
  }
GFCM;

typedef struct{
  uint64_t length;    // NUMBER OF SEQUENCE SYMBOLS
  uint8_t  nSym;      // NUMBER OF SYMBOLS [ALPHABET CARDINALITY]
  uint8_t  *symbolic; // SYMBOLIC ALPHABET ARRAY
  uint8_t  *numeric;  // NUMERIC ALPHABET ARRAY
  }
ALPHA;

typedef struct{
  GFCM     **M;       // HEADERS FCMs POINTER          
  ALPHA    A;         // ALPHABET POINTER  
  uint8_t  nFCM;      // NUMBER OF FCMs
  uint8_t  type;      // TYPE OF HEADERS
  uint8_t  extra;     // PRESENCE OR NOT OF THE EXTRA HEADER
  }
HEADERS;

typedef struct{
  FCM      **M;       // DNA FCMs POINTER
  ALPHA    A;         // ALPHABET POINTER  
  uint8_t  nFCM;      // NUMBER OF FCMs
  uint8_t  *bica;     // BINARY INFORMATION CONTENT ARRAY [LOW OR HIGH ENTROPY]
  uint8_t  lfb;       // LESS FREQUENT BASE
  }
DNA;

typedef struct{
  GFCM     **M;       // SCORES FCMs POINTER
  ALPHA    A;         // ALPHABET POINTER  
  uint8_t  nFCM;      // NUMBER OF FCMs
  uint64_t maxLine;   // MAXIMUM LINE NUMBER
  uint8_t  dynamic;   // IF NUMBER OF LINE IS DYNAMIC
  uint8_t  not0N;     // ANY SCORE ZERO MATCH WITH 'N' BASE? IF NOT => 1 
  }
SCORES;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t    CalcAlphaDen    (uint8_t, uint32_t);
uint32_t    CalcAlphaDenT   (uint8_t, uint32_t);
uint32_t    AdjustContext   (uint8_t, uint32_t);
void        Init4DnaHashTab (FCM *);
void        InitGHashTab    (GFCM *);
void        Reset4DnaModel  (FCM *);
void        Free4DnaModel   (FCM *);
void        FreeGModel      (GFCM *);
void        Update4DnaFCM   (FCM *, uint32_t, uint8_t);
void        UpdateGFCM      (GFCM *, uint32_t);
FCM         *Create4DnaFCM  (uint32_t, uint32_t, uint8_t, uint8_t, PARAM *);
GFCM        *CreateGFCM     (uint32_t, uint32_t, uint8_t, PARAM *);
inline void Compute4DnaFCM  (FCM *);
inline void ComputeGun      (FCM *, uint32_t *);
double      CompProbs       (FCM *, uint32_t);
void        ComputeGFCM     (GFCM *);
inline void GetIdx4DnaRev   (uint8_t *, FCM *);
inline void GetIdx4Dna      (uint8_t *, FCM *);
inline void GetIdx          (uint8_t *, GFCM *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif




