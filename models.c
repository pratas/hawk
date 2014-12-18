#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "models.h"
#include "misc.h"
#include "dna.h"

static uint64_t XHASH(uint64_t x){
  return (x * 12582917 /*786433*/ + 196613) % 68719476735;
     // 19m = 12582917 (but lower error)
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CALCULATE ALPHA DEN BASED ON CONTEXT AND ALPHABET CARDINALITY 
// 
uint32_t CalcAlphaDen(uint8_t n, uint32_t c){
  uint64_t x = (uint64_t) PW(n, c);
  return x<268435456?1:x<4294967297?1:x<274877906944?50:100;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CALCULATE ALPHA DEN FOR TRAINNING PHASE
//
uint32_t CalcAlphaDenT(uint8_t n, uint32_t c){
  return 1;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// ADJUST CONTEXT MODEL IN FUNCTION OF THE ALPHABET
//
uint32_t AdjustContext(uint8_t n, uint32_t c){
  if(n < 4 && c>14) return 14;
  switch(n){
    case 4:  if(c>13) return 13;
    case 5:  if(c>11) return 11;
    case 6:  if(c>10) return 10;
    case 7:  if(c>9)  return 9;
    default: if(c>4)  return 4;
    }
  return c;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static HCCs hzeroCnts = {0x00, 0x00, 0x00, 0x00};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INITIALIZE 4 SYMBOL ARRAY TABLE
//
static void Init4DnaArrayTab(FCM *M){
  M->A.cnts = (ACC *) Calloc(M->nPMod<<2, sizeof(ACC));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INITIALIZE GENERIC ARRAY TABLE
//
static void InitGArrayTab(GFCM *M){
  M->A.cnts = (GACC *) Calloc(M->nPMod*M->nSym, sizeof(GACC));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RESET MODEL
//
void Reset4DnaModel(FCM *M){
  switch(M->mode){
    case CCH_TABLE:
      DeleteCCH(M->B);
      M->B = CreateCCH(CCH_SIZE, M->nSym);
    break;
    case HASH_TABLE:
      ResetPHash(M->H);
    break;
    default:
      Free(M->A.cnts, (M->nPMod<<2) * sizeof(ACC)); 
      Init4DnaArrayTab(M); 
    break;
    }
  M->idx = 0; M->idxRev = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FREE MODEL
//
void Free4DnaModel(FCM *M){
  switch(M->mode){
    case CCH_TABLE:
      DeleteCCH(M->B);
    break;
    case HASH_TABLE:
      DeletePHash(M->H);
    break;
    default: 
      Free(M->A.cnts, (M->nPMod<<2) * sizeof(ACC));
    break;
    }
  Free(M->freqs, 5 * sizeof(uint32_t));
  Free(M, sizeof(FCM));  
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FREE GLOBAL MODEL
//
void FreeGModel(GFCM *M){
  switch(M->mode){
    case CCH_TABLE:
      DeleteCCH(M->B);
    break;
    case HASH_TABLE:
      fprintf(stderr, "[x] Error: currently not implemented!\n");
      exit(1);
    break;
    default:
      Free(M->A.cnts, M->nPMod * M->nSym * sizeof(GACC));
    break;
    }
  Free(M->freqs, (M->nSym+1) * sizeof(uint32_t));
  Free(M, sizeof(GFCM));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UPDATE FCM COUNTERS
//
void Update4DnaFCM(FCM *M, uint32_t c, uint8_t ir){
  ACC *ac;
  uint64_t idx = (ir == 0) ? M->idx : M->idxRev;
  if(M->mode == CCH_TABLE){
    UpdateCCH(M->B, idx, c);
    }
  else if(M->mode == HASH_TABLE){
    UpdatePHash(M->H, XHASH(idx), c);
    }
  else{
    #ifdef PREFETCHING
    _mm_prefetch((ACC *)&M->A.cnts[idx<<2], _MM_HINT_T1);
    #endif
    ac = &M->A.cnts[idx<<2];
    if(++ac[c] == MAXACC_C){ ac[0]>>=1; ac[1]>>=1; ac[2]>>=1; ac[3]>>=1; }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UPDATE GENERAL FCM
//
void UpdateGFCM(GFCM *M, uint32_t c){
  GACC *ac;
  if(M->mode == CCH_TABLE){
    UpdateCCH(M->B, M->idx, c);
    }
  else{
    uint32_t n;
    ac = &M->A.cnts[M->idx*M->nSym];
    if(++ac[c] == MAXGACC_C){ 
      for(n = 0 ; n < M->nSym ; ++n)
        ac[n]>>=1; 
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATES 4 SYMBOL FCM
//
FCM *Create4DnaFCM(uint32_t c, uint32_t a, uint8_t i, uint8_t n, PARAM *A){
  FCM    *M = (FCM *) Calloc(1, sizeof(FCM));
  M->nSym   = n;
  M->ctx    = c;
  M->rev    = i;
  M->aDen   = a;
  M->mult   = CalcMult(c, n);
  M->idx    = 0;
  M->idxRev = 0;
  M->freqs  = (uint32_t *) Calloc(M->nSym+1, sizeof(uint32_t));
  M->nPMod  = (uint64_t) pow(n,c);

  if(M->nPMod < DEEP_CTX){
    Init4DnaArrayTab(M);
    M->mode = 0;
    }
  else{
    if(A->mode == 0){
      M->B = CreateCCH(CCH_SIZE, M->nSym);
      M->mode = 1;
      }
    else{
      M->H = CreatePHash();
      M->mode = 2;
      }
    }

  if(A->verbose == 1){
    fprintf(stderr, "  [+] mode .................. %s\n", M->mode == 0 ? 
    "Counter-table":(A->mode==0?"Counter Coast Hash-table":"Hash-table"));
    fprintf(stderr, "  [+] context ............... %u\n", c);
    fprintf(stderr, "  [+] alpha ................. %.3g\n", 1.0/a);
    if(i == 1) fprintf(stderr, "  [+] using inversions ...... yes\n");
    }
  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATES GENERAL FCM
//
GFCM *CreateGFCM(uint32_t c, uint32_t a, uint8_t n, PARAM *A){
  GFCM   *M = (GFCM *) Calloc(1, sizeof(GFCM));
  M->nSym   = n;
  M->ctx    = c;
  M->aDen   = a;
  M->mult   = CalcMult(c, n);
  M->idx    = 0;
  M->freqs  = (uint32_t *) Calloc(M->nSym+1, sizeof(uint32_t));
  M->nPMod  = (uint64_t) pow(n,c);

  if(M->nPMod < DEEP_CTX){
    InitGArrayTab(M); 
    M->mode = 0;
    }
  else{
    M->B = CreateCCH(CCH_SIZE, M->nSym);
    M->mode = 1;
    }

  if(A->verbose == 1){
    fprintf(stderr, "  [+] mode .................. %s\n", M->mode == 0 ?
    "Counter-table" : "Counter Coast Hash-table");
    fprintf(stderr, "  [+] context ............... %u\n", c);
    fprintf(stderr, "  [+] alpha ................. %.3g\n", 1.0/a);
    }
  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPUTE SPECIFIC 4 SYMBOL FCM PROBABILITIES
//
inline void Compute4DnaFCM(FCM *M){
  HCC *h;
  ACC *a;
  if(M->mode == HASH_TABLE){
    if(!(h = GetPHashCounters(M->H, XHASH(M->idx)))) h = hzeroCnts;
    M->freqs[0] = 1+M->aDen*h[0]; M->freqs[1] = 1+M->aDen*h[1];
    M->freqs[2] = 1+M->aDen*h[2]; M->freqs[3] = 1+M->aDen*h[3];
    }
  else{
    a = &M->A.cnts[M->idx<<2];
    M->freqs[0] = 1+M->aDen*a[0]; M->freqs[1] = 1+M->aDen*a[1];
    M->freqs[2] = 1+M->aDen*a[2]; M->freqs[3] = 1+M->aDen*a[3];
    }
  M->freqs[4] = M->freqs[0]+M->freqs[1]+M->freqs[2]+M->freqs[3];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPUTE SPECIFIC 4 SYMBOL FCM PROBABILITIES
//
inline void ComputeGun(FCM *M, uint32_t *f){
  HCC   *h;
  C_CCH *c;
  ACC   *a;
  switch(M->mode){
    case CCH_TABLE:
      c = SearchCCH(M->B, M->idx);
      f[0] = 1+M->aDen*c[0]; f[1] = 1+M->aDen*c[1];
      f[2] = 1+M->aDen*c[2]; f[3] = 1+M->aDen*c[3];
    break;
    case HASH_TABLE:
      if(!(h = GetPHashCounters(M->H, XHASH(M->idx)))) h = hzeroCnts;
      f[0] = 1+M->aDen*h[0]; f[1] = 1+M->aDen*h[1];
      f[2] = 1+M->aDen*h[2]; f[3] = 1+M->aDen*h[3];
    break;
    default:
      a = &M->A.cnts[M->idx<<2];
      f[0] = 1+a[0]; f[1] = 1+a[1]; f[2] = 1+a[2]; f[3] = 1+a[3];
    break;
    }
  f[4] = f[0]+f[1]+f[2]+f[3];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPUTE GENERAL FCM PROBABILITIES
//
void ComputeGFCM(GFCM *M){
  GACC *a;
  uint32_t x;
  if(M->mode == HASH_TABLE){
    fprintf(stderr, "[x] Error: currently not implemented!\n");
    exit(1);
    }
  else{
    a = &M->A.cnts[M->idx*M->nSym];
    for(x = 0 ; x < M->nSym ; ++x)
      M->freqs[x] = 1 + M->aDen * a[x]; 
    }
  M->freqs[M->nSym] = M->freqs[0];
  for(x = 1 ; x < M->nSym ; ++x)
    M->freqs[M->nSym] += M->freqs[x];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPUTE PROBABILITIES
//
uint32_t CompProbs(FCM *M, uint32_t s){
  return Log(M->freqs[M->nSym]/M->freqs[s]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REVERSE COMPLEMENT INDEX BASED ON PAST SYMBOLS FOR FCM WITH 4 SYMBOLS
//
inline void GetIdx4DnaRev(uint8_t *p, FCM *M){
  M->idxRev = (M->idxRev>>2)+Comp(*p)*M->mult;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INDEX CALC BASED ON PAST SYMBOLS FOR FCM WITH 4 SYMBOLS
//
inline void GetIdx4Dna(uint8_t *p, FCM *M){
  M->idx = ((M->idx-*(p-M->ctx)*M->mult)<<2)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INDEX CALC BASED ON PAST SYMBOLS FOR GENERAL FCMs 
//
inline void GetIdx(uint8_t *p, GFCM *M){
  M->idx = ((M->idx-*(p-M->ctx)*M->mult)*M->nSym)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CALC INDEX AND RETURN
//
inline uint64_t GetIdxA(uint8_t *p, GFCM *M){
  return (M->idx = ((M->idx-*(p-M->ctx)*M->mult)*M->nSym)+*p);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

