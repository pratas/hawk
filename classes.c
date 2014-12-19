#include <stdio.h>
#include <stdlib.h>
#include "classes.h"
#include "param.h"
#include "args.h"
#include "dna.h"
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
// INITIALIZE AND ALPHABETS
//
void InitAlphabets(CLASSES *C){
  C->H.A.symbolic = (uint8_t *) Calloc(MAX_AL, sizeof(uint8_t));
  C->D.A.symbolic = (uint8_t *) Calloc(MAX_AL, sizeof(uint8_t));
  C->S.A.symbolic = (uint8_t *) Calloc(MAX_AL, sizeof(uint8_t));
  C->H.A.numeric  = (uint8_t *) Calloc(MAX_AL, sizeof(uint8_t));
  C->D.A.numeric  = (uint8_t *) Calloc(MAX_AL, sizeof(uint8_t));
  C->S.A.numeric  = (uint8_t *) Calloc(MAX_AL, sizeof(uint8_t));
  C->H.A.length   = 0;
  C->D.A.length   = 0;
  C->S.A.length   = 0;
  C->H.A.nSym     = 0;
  C->D.A.nSym     = 0;
  C->S.A.nSym     = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FREE ALPHABETS
//
void FreeAlphabets(CLASSES *C){
  Free(C->H.A.symbolic, MAX_AL * sizeof(uint8_t));
  Free(C->D.A.symbolic, MAX_AL * sizeof(uint8_t));
  Free(C->S.A.symbolic, MAX_AL * sizeof(uint8_t));
  Free(C->H.A.numeric,  MAX_AL * sizeof(uint8_t));
  Free(C->D.A.numeric,  MAX_AL * sizeof(uint8_t));
  Free(C->S.A.numeric,  MAX_AL * sizeof(uint8_t));
  C->H.A.length   = 0;
  C->D.A.length   = 0;
  C->S.A.length   = 0;
  C->H.A.nSym     = 0;
  C->D.A.nSym     = 0;
  C->S.A.nSym     = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PRINT FILE INFORMATION [PROPERTIES & ALPHABET]
//
void PrintStreamInfo(CLASSES *C){
  uint32_t k;
  char y[]="yes", n[]="no";
  fprintf(stderr, "==[ FILE INFORMATION ]=======\n");
  fprintf(stderr, "File has got ................ "); PrintHRBytes(C->length);
  fprintf(stderr, " (%"PRIu64" Bytes)\n", C->length);
  fprintf(stderr, "Number of reads ............. %"PRIu64"\n", C->nReads);
  fprintf(stderr, "Headers information:\n");
  fprintf(stderr, "  [+] Length ................ %"PRIu64"\n", C->H.A.length);
  fprintf(stderr, "  [+] Cardinality ........... %u\n", C->H.A.nSym);
  fprintf(stderr, "  [+] Alphabet .............. ");
  for(k = 0 ; k < C->H.A.nSym ; ++k)
    C->H.A.symbolic[k] == '\n' ? fprintf(stderr, "\\n") :
    fprintf(stderr, "%c", C->H.A.symbolic[k]);
  fprintf(stderr, "\n  [+] Maximum line size ..... %"PRIu64"\n", C->H.maxLine);
  fprintf(stderr, "  [+] Number of states ...... %u\n", C->H.nStates);
  fprintf(stderr, "  [+] Extra header .......... %s\n", C->H.extra?y:n);
  fprintf(stderr, "DNA information:\n");
  fprintf(stderr, "  [+] Length ................ %"PRIu64"\n", C->D.A.length);
  fprintf(stderr, "  [+] Cardinality ........... %u\n", C->D.A.nSym);
  fprintf(stderr, "  [+] Alphabet .............. ");
  for(k = 0 ; k < C->D.A.nSym ; ++k)
    C->D.A.symbolic[k] == '\n' ? fprintf(stderr, "\\n") :
    fprintf(stderr, "%c", C->D.A.symbolic[k]);
  fprintf(stderr, "\nScores information:\n");
  fprintf(stderr, "  [+] Length ................ %"PRIu64"\n", C->S.A.length);
  fprintf(stderr, "  [+] Cardinality ........... %u\n", C->S.A.nSym);
  fprintf(stderr, "  [+] Alphabet .............. ");
  for(k = 0 ; k < C->S.A.nSym ; ++k)
    C->S.A.symbolic[k] == '\n' ? fprintf(stderr, "\\n") :
    fprintf(stderr, "%c", C->S.A.symbolic[k]);
  fprintf(stderr, "\n  [+] Maximum line size ..... %"PRIu64"\n", C->S.maxLine);
  fprintf(stderr, "  [+] Dynamic line size ..... %s\n", C->S.dynamic?y:n);
  fprintf(stderr, "  [+] Lower score match N ... %s\n", C->S.not0N?n:y);
  fprintf(stderr, "  [+] Number of states ...... %u\n", C->S.nStates);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PARSE FILE & BUILD ALPHABET
//
void ParseFile(CLASSES *C, FILE *F, PARAM *A){
  int64_t k;
  uint64_t i, pos = 0, size = BIN_LINE;
  uint32_t line = 0, e = 0, off = 0, low = 500;
  uint8_t s, binH[MAX_AL], binD[MAX_AL], binS[MAX_AL], *buf, *lD, *bases;
  for(k = 0 ; k < MAX_AL ; ++k){ binH[k] = 0; binD[k] = 0; binS[k] = 0; }

  fprintf(stderr, "[>] Running small file analysis ...\n");
  InitAlphabets(C);

  bases = (uint8_t *) Calloc(NSYM,     sizeof(uint8_t));
  buf   = (uint8_t *) Calloc(BUF_SIZE, sizeof(uint8_t));
  lD    = (uint8_t *) Calloc(BIN_LINE, sizeof(uint8_t));
  while((k = fread(buf, 1, BUF_SIZE, F)))
    for(i = 0 ; i < k ; ++i){
      s = *(buf+i);
      switch(line){
        case 0:
          if(C->H.maxLine < ++pos) C->H.maxLine = pos; 
          if(s == '\n'){ line = 1; pos = 0; break; }
          C->H.A.length++; binH[s] = 1; 
        break;
        case 1: 
          switch(s){
            case '\n': line = 2; break;
            case 'N' : case '.': binD[s] = 1; lD[pos] = 1; C->D.A.length++; 
            break;
            default  : binD[s] = 1; lD[pos] = 0; AssignLowerBase(bases, s);
            C->D.A.length++;
            }
          if(++pos == size)
            lD = (uint8_t *) Realloc(lD, (size += BIN_LINE) * sizeof(uint8_t), 
                 BIN_LINE * sizeof(uint8_t));
        break;
        case 2:
          ++e;
          if(s == '\n'){ if(e == 3) C->H.extra = 0; e = 0; line = 3; pos = 0; }
        break;
        case 3:
          switch(off){ // VERIFY IF: LOWER SCORE MATCH 'N'. If not: S->not0N = 1
            case 0: 
              if(lD[pos] == 1){ if(binS[s] == 1){ C->S.not0N = 1; off = 2; }
                else{ low = s; off = 1; }}
            break;           
            case 1:          
              if((s == low && lD[pos] == 0) || (lD[pos] == 1 && s != low)){ 
                C->S.not0N = 1; off = 2; }
            break;          
            case 2: break;
            }
          if(C->S.maxLine < ++pos) C->S.maxLine = pos;
          binS[s] = 1;
          if(s == '\n'){
            C->nReads++; line = 0; if(pos != C->S.maxLine) C->S.dynamic = 1;
            pos = 0;
            break;
            }
          C->S.A.length++;
        break;
        }
     }
  fprintf(stderr, "[>] Done!                     \n");
 
  for(k = 0 ; k < MAX_AL ; ++k){
    if(binH[k] == 1){ C->H.A.symbolic[C->H.A.nSym]=k; 
      C->H.A.numeric[k]=C->H.A.nSym++; } 
    else C->H.A.numeric[k] = (uint8_t) INVALID_S;
    if(binD[k] == 1){ C->D.A.symbolic[C->D.A.nSym]=k; 
      C->D.A.numeric[k]=C->D.A.nSym++; }
    else C->D.A.numeric[k] = (uint8_t) INVALID_S;
    if(binS[k] == 1){ C->S.A.symbolic[C->S.A.nSym]=k; 
      C->S.A.numeric[k]=C->S.A.nSym++; }
    else C->S.A.numeric[k] = (uint8_t) INVALID_S;
    }
    
  rewind(F); if(A->verbose) PrintStreamInfo(C); C->D.lfb = CalcLessFreq(bases);
  // ADJUST DNA TO 4 SYM &( USE EXTRA STREAM OR USE GENERAL FCM)
  if(binD['A']==1 && binD['C']==1 && binD['G']==1 && binD['T']==1){
    if(C->D.A.nSym>4){
      for(k = 0 ; k < MAX_AL ; ++k) binD[k] = 0;
      binD['A']=1; binD['C']=1; binD['G']=1; binD['T']=1; C->D.A.nSym=0;
      for(k = 0 ; k < MAX_AL ; ++k){
        if(binD[k] == 1){ C->D.A.symbolic[C->D.A.nSym] = k;
          C->D.A.numeric[k] = C->D.A.nSym++; }
        }
      }
    }  
  else if(A->filter == 1){
    Msg(A, "[!] Found non-regular DNA data  [filtering usage: off].\n");
    Msg(A, "[!] General purpose FCM will be used [inversions: off].\n");
    A->filter = 0; A->inverse = 0; 
    }

  // REMOVE '\n' IF QUALITY STREAMS HAVE FIXED SIZES
  if(C->S.dynamic == 0){
    binS['\n']  = 0;
    C->S.A.nSym = 0;
    for(k = 0 ; k < MAX_AL ; ++k){
        if(binS[k] == 1){ C->S.A.symbolic[C->S.A.nSym] = k;
          C->S.A.numeric[k] = C->S.A.nSym++; }
        }
    }

  Free(bases,   NSYM * sizeof(uint8_t));
  Free(lD,      size * sizeof(uint8_t));
  Free(buf, BUF_SIZE * sizeof(uint8_t));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATE CLASSES WITH MODELS
// 
void CreateModels(CLASSES *C, char *p, FILE *F, PARAM *A){
  char **xv;
  int32_t n, state, nD = 0, ctx, xc = StrToArgv(p, &xv);

  if(A->verbose) fprintf(stderr, "==[ MODELS ]=================\n");
  C->H.M = (GFCM **) Calloc(C->H.nStates, sizeof(GFCM *));
  C->H.nFCM = C->H.nStates;
  for(n = 0 ; n < xc ; ++n)
    if(!strcmp("-H", xv[n])){
      if((ctx = atoi(xv[n+1])) > HMAX_CTX || ctx < HMIN_CTX){
        fprintf(stderr, "Error: invalid header context (-H)!\n");
        exit(1);
        }
      for(state = 0 ; state < C->H.nStates ; ++state){
        if(A->verbose == 1)
          fprintf(stderr, "Creating Headers model [state %u of %u]:\n", state+1,
          C->H.nStates);
        C->H.M[state] = CreateGFCM(ctx, CalcAlphaDen(C->H.A.nSym, ctx),
        C->H.A.nSym, A);
        }
      break;
      }
  C->H.nFCM = C->H.nStates;

  for(n = 0 ; n < xc ; ++n) if(!strcmp("-D", xv[n])) ++nD;
  C->D.M = (FCM **) Calloc(nD, sizeof(FCM *));
  for(n = 0 ; n < xc ; ++n)
  if(!strcmp("-D", xv[n])){
    if((ctx = atoi(xv[n+1])) > DMAX_CTX || ctx < DMIN_CTX){
      fprintf(stderr, "Error: invalid DNA context (-D)!\n");
      exit(1);
      }
    Msg(A, "Creating DNA model:\n");
    if(A->adjust) ctx = AdjustContext(C->D.A.nSym, ctx);
    C->D.M[C->D.nFCM++] = Create4DnaFCM(ctx, CalcAlphaDen(C->D.A.nSym, ctx),
    0, C->D.A.nSym, A);
    }

  C->S.M = (GFCM **) Calloc(C->S.nStates, sizeof(GFCM *));
  C->S.nFCM = C->S.nStates;
  for(n = 0 ; n < xc ; ++n)
    if(!strcmp("-S", xv[n])){
      if((ctx = atoi(xv[n+1])) > SMAX_CTX || ctx < SMIN_CTX){
        fprintf(stderr, "Error: invalid Scores context (-S)!\n");
        exit(1);
        }
      for(state = 0 ; state < C->S.nStates ; ++state){
        if(A->verbose == 1)
          fprintf(stderr, "Creating Scores model [state %u of %u]:\n", state+1, 
          C->S.nStates);
        C->S.M[state] = CreateGFCM(ctx, CalcAlphaDen(C->S.A.nSym, ctx),
        C->S.A.nSym, A);
        } 
      break;
      }

  Free(*xv, xc * sizeof(char *)); 
  Free(xv, sizeof(char));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FREE ALL REMAINING MODELS
// 
void FreeModels(CLASSES *C, PARAM *A){
  uint32_t n;
  if(A->filter == 1) Free(C->D.bica, C->nReads * sizeof(uint8_t));
  for(n = 0 ; n < C->H.nStates ; ++n) FreeGModel(C->H.M[n]);
  for(n = 0 ; n < C->D.nFCM ; ++n) Free4DnaModel(C->D.M[n]);
  for(n = 0 ; n < C->S.nStates ; ++n) FreeGModel(C->S.M[n]);
  Free(C->H.M, C->H.nStates * sizeof(GFCM *));
  Free(C->D.M, C->D.nFCM    * sizeof(FCM  *));
  Free(C->S.M, C->S.nStates * sizeof(GFCM *));
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

