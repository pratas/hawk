////////////////////////////////////////////////////////////////////////////////
//............................................................................//
//...........................,~=~~:?$7?+II+=:I=~..............................//
//.......................,:+=II+=++I~:===~?=,I7=~?~...........................//
//.....................,=?+:++?+=?I?+:~,:=+:=7Z7I?Z?~.........................//
//...................,?=++?==~:+==+=~,,=~~?:.,=I+=7?II,.......................//
//.................,=:~~,~+=?I=~?II~=..::.,...~I=~:+~=:~,.....................//
//................,~7+~I$I$~?O?8+ZN?88OZI7$I~+:::~+?=7?..?....................//
//..............,~??I+?OZI?+IOMNNZ$8=7$DOOD8DOO$77I:~7?~:+?,..................//
//.............,,~:I=:~??+D8MNMOMMMON$IZZOIOZ$O7=?$$?,,??I~$7.................//
//............,+~,:=I~7$I?7I?MNI$D88NI?7OI?7?+$==I?$7+.=~,+7?:................//
//...........,==:~~=II$7ZZ=?+=NNNMNZ$777Z7++=?I++==7II~:+=I=IZO,..............//
//...........7~~==N8?$I++I+=7D7OOOO87I$?$IZ??=+++=?+??Z++?II$Z~Z..............//
//..........$:~:=~+?OI====I==?+=77I,+II~~?~:+=I7?+==+$$7~:7I?=$~?.............//
//.........M$=I?===+::~,~:~=?$$??I:?+~=~~=~:=:~+=++??~+IIIZZ$?~$ZI............//
//.........Z777Z8?=+++=?$=78ND8NN+78Z$$O$I=??7Z77I++==7I+==:~=~Z7$$...........//
//.........O+Z?7OI+$I7I$Z$+?ON?7DOZ8OD8$?+~~~,~78Z7+.==~==?~+:.:+~.~..........//
//.........O77MI8OND+~,+8DD8DDZO88II=$NDM888=OO877+==++=~?+~:==:?8Z=$.........//
//.........II7N......:....:O8D8ODNNDDODO$D8ZO7Z$7++?,:=I=~~:.+II?+~Z$7........//
//.........Z=7..............Z8OOOZ8OZ$$7II7$:,~~=~=~=Z?+=II=:,..:77:.~=.......//
//..........?$...............8NZ8OZDO+=~,~:::::~+7$:~++I?=?7$+~7?~ZO,=$8,.....//
//...........8..............+ZD8OZZNI?~=:,?~:::+?I=......,.::=+IO8?OOO+8O=....//
//.........................,OZ8Z$8O?$?=~7~~+::,~,,...........:III+==+I,.~8O~..//
//.........................$ZZO8OZOO7I7::~~~,,,,:,..::~,......?7II~7+7,:88+O$~//
//........................,$OZZZZZDZI?OI::$~::?=,=+~:,I=:..,:=+++77ZOO8O$8OO$Z//
//........................ZZZZZOO8Z$IZ7+=+?~=~~I~:I~,:~I?~,,.~$$OO$=:~O?++IO$Z//
//.......................,$8$ZZOOOZZOZO++?=?++++=:~~,:,++=:.......,7$$7ZI=+=+?//
//.......................$$DDZZOOZZZZZZ$77++~++?~~~:::,,:?~,..  ....7I7I=,:..,//
//.......................$ZDD$$$$ZOZO8OO$77++==?=~::::::::~:.........::,,~,:..//
//......................,N$O8Z$$$ZO888OOZ7?+=~=?~+~~:::::=7::,.,~..,..:~+~$:=://
////////////////////////////////////////////////////////////////////////////////
//............................................................................//
//.........@@@.....@@@.......@@@.....@@@.............@@@..@@@....@@@..........//
//.........@@@.....@@@......@@@@@....@@@.....@@@.....@@@..@@@...@@@...........//
//.........@@@@@@@@@@@.....@@@.@@@....@@@...@@@@@...@@@...@@@@@@@@............//
//.........@@@@@@@@@@@....@@@@@@@@@....@@@@@@@.@@@@@@@....@@@@@@@@............//
//.........@@@.....@@@...@@@.....@@@....@@@@@...@@@@@.....@@@...@@@...........//
//.........@@@.....@@@..@@@.......@@@....@@@.....@@@......@@@....@@@..........//
//............................................................................//
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     HAWK: A COMPRESSOR FOR FASTQ FILES                     //
//                                                                            //
//                     IEETA @ UNIVERSITY OF AVEIRO, 2015                     //
//                                                                            //
//                     TO COMPILE : make                                      //
//                     TO RUN HELP: ./hawk -h                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <malloc.h>
#include <unistd.h>

#include "defs.h"
#include "mem.h"
#include "misc.h"
#include "args.h"
#include "param.h"
#include "reads.h"
#include "dna.h"
#include "hash.h"
#include "cch.h"
#include "classes.h"
#include "models.h"
#include "gun.h"
#include "info.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INITIALIZE AND FREE ALPHABETS
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
  int32_t n, nH = 0, nD = 0, nS = 0, ctx, xc = StrToArgv(p, &xv);

  if(A->verbose) fprintf(stderr, "==[ MODELS ]=================\n");
  for(n = 0 ; n < xc ; ++n) if(!strcmp("-H", xv[n])) ++nH;
  C->H.M = (GFCM **) Calloc(nH, sizeof(GFCM *));
  for(n = 0 ; n < xc ; ++n)
    if(!strcmp("-H", xv[n])){
      if((ctx = atoi(xv[n+1])) > HMAX_CTX || ctx < HMIN_CTX){
        fprintf(stderr, "Error: invalid header context (-H)!\n");
        exit(1);
        }
      Msg(A, "Creating Header model:\n");
      C->H.M[C->H.nFCM++] = CreateGFCM(ctx, CalcAlphaDen(C->H.A.nSym, ctx), 
      C->H.A.nSym, A);
      }

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

  //for(n = 0 ; n < xc ; ++n) if(!strcmp("-S", xv[n])) ++nS;
  C->S.M = (GFCM **) Calloc(C->S.nStates, sizeof(GFCM *));
  for(n = 0 ; n < xc ; ++n)
    if(!strcmp("-S", xv[n])){
      if((ctx = atoi(xv[n+1])) > SMAX_CTX || ctx < SMIN_CTX){
        fprintf(stderr, "Error: invalid Scores context (-S)!\n");
        exit(1);
        }
      uint32_t state;
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
  for(n = 0 ; n < C->H.nFCM ; ++n) FreeGModel(C->H.M[n]);
  for(n = 0 ; n < C->D.nFCM ; ++n) Free4DnaModel(C->D.M[n]);
  for(n = 0 ; n < C->S.nFCM ; ++n) FreeGModel(C->S.M[n]);
  Free(C->H.M, C->H.nFCM * sizeof(GFCM *));
  Free(C->D.M, C->D.nFCM * sizeof(FCM  *));
  Free(C->S.M, C->S.nFCM * sizeof(GFCM *));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// DNA SEQUENCE TO DISK
// 
void DNA2Disk(CLASSES *C, PARAM *A, FILE *FI, FILE *DSK, FILE *DSKR){
  fprintf(stderr, "[>] DNA sequence to disk ...\n");
  uint64_t k, i, r = 0, x = 0;
  int64_t size = C->D.A.length+C->nReads;
  uint32_t line = 0;
  uint8_t *buf, *out, *rev, s;
  buf = (uint8_t *) Calloc(BUF_SIZE, sizeof(uint8_t));
  out = (uint8_t *) Calloc(BUF_SIZE, sizeof(uint8_t));
  rev = (uint8_t *) Calloc(BUF_SIZE, sizeof(uint8_t));
  while((k=fread(buf, 1, BUF_SIZE, FI)))
    for(i=0 ; i<k ; ++i){
      s=*(buf+i);
      switch(line){
        case 0: if(s == '\n'){ line = 1; } break;
        case 1:
          if(s == '\n'){ line = 2; } out[x]=s; 
          if(++x==BUF_SIZE){ fwrite(out, 1, x, DSK); ReverseStr(out, rev, x);
            fseek(DSKR, size -= x, SEEK_SET); fwrite(rev, 1, x, DSKR); x=0; } 
        break;
        case 2: if(s == '\n')  line = 3; break;
        case 3: if(s == '\n'){ line = 0; Progress(C->nReads, ++r); } break;
        }
      }
  if(x){ fwrite(out, 1, x, DSK); ReverseStr(out, rev, x);
    fseek(DSKR, size -= x, SEEK_SET); fwrite(rev, 1, x, DSKR); }

  Free(buf, BUF_SIZE * sizeof(uint8_t)); 
  Free(rev, BUF_SIZE * sizeof(uint8_t)); 
  Free(out, BUF_SIZE * sizeof(uint8_t)); 
  rewind(FI); rewind(DSK); rewind(DSKR);
  fprintf(stderr, "[>] Done!                         \n"); // SPACES ARE VALID 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// EVALUATE LEFT TO RIGHT "ACGT" DNA SEQUENCE [4 SYMBOL]
// 
void EvaluateLRDna(CLASSES *C, PARAM *A, FILE *F, FCM *F1, FCM *F2){
  fprintf(stderr, "[>] Evaluating LR: M1[ctx:%u] vs "
  "M2[ctx:%u] ...\n", F1->ctx, F2->ctx);
  double ca = 0, cb = 0;
  int64_t x = 0;
  uint64_t k, i, r = 0, z = 0;
  uint8_t *buf, *seq, s, n;

  buf  = (uint8_t *) Calloc(BUF_SIZE, sizeof(uint8_t));
  seq  = (uint8_t *) Calloc(BUF_SIZE+GUARD, sizeof(uint8_t));
  seq += GUARD;
  while((k=fread(buf, 1, BUF_SIZE, F)))
    for(i=0 ; i<k ; ++i){
      s=*(buf+i);
      if(s=='\n'){ 
        if(C->D.bica[r]==0 && ca>cb)
          C->D.bica[r]=1; 
        ca=cb=z=0; 
        Progress(C->nReads, ++r); continue; 
        } 
      else if(s=='A'||s=='C'||s=='G'||s=='T') n=C->D.A.numeric[s];
      else n = C->D.lfb;
      seq[x]=n;
      GetIdx4Dna(seq+x-1, F1); GetIdx4Dna(seq+x-1, F2);
      if(A->inverse == 1) GetIdx4DnaRev(seq+x, F2);
      if(C->D.bica[r]==0){
        Compute4DnaFCM(F1); Compute4DnaFCM(F2);
        ca+=CompProbs(F1, n); cb+=CompProbs(F2, n);
        }
      Update4DnaFCM(F1, n, 0); Update4DnaFCM(F2, n, 0);
      if(A->inverse == 1) Update4DnaFCM(F2, 3-(seq[x-F2->ctx]), 1);
      if(++x==BUF_SIZE){ memcpy(seq-GUARD, seq+x-GUARD, GUARD); x=0; } 
      ++z;
      }
  
  Reset4DnaModel(F1); Reset4DnaModel(F2);
  Free(seq-GUARD, (BUF_SIZE+GUARD) * sizeof(uint8_t)); 
  Free(buf, BUF_SIZE * sizeof(uint8_t)); 
  rewind(F);
  fprintf(stderr, "[>] Done!                         \n"); // SPACES ARE VALID 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// EVALUATE RIGHT TO LEFT "ACGT" DNA SEQUENCE [4 SYMBOL]
// 
void EvaluateRLDna(CLASSES *C, PARAM *A, FILE *F, FCM *F1, FCM *F2){
  fprintf(stderr, "[>] Evaluating RL: M1[ctx:%u] vs "
  "M2[ctx:%u] ...\n", F1->ctx, F2->ctx);
  double ca = 0, cb = 0;
  int64_t x = 0;
  uint64_t k, i, r = 0, bk = C->nReads-1, z = 0;
  uint8_t *buf, *seq, s, n;

  buf  = (uint8_t *) Calloc(BUF_SIZE, sizeof(uint8_t));
  seq  = (uint8_t *) Calloc(BUF_SIZE+GUARD, sizeof(uint8_t));
  seq += GUARD;
  while((k=fread(buf, 1, BUF_SIZE, F)))
    for(i=0 ; i<k ; ++i){
      s=*(buf+i);
      if(s=='\n'){
        if(i==0 && r==0) continue;
        if(C->D.bica[bk]==0 && ca>cb) C->D.bica[bk]=1;
        ca=cb=z=0;
        Progress(C->nReads, ++r); 
        --bk; 
        continue;
        }
      else if(s=='A'||s=='C'||s=='G'||s=='T') n=C->D.A.numeric[s];
      else n = C->D.lfb;
      seq[x]=n;
      GetIdx4Dna(seq+x-1, F1); GetIdx4Dna(seq+x-1, F2);
      if(A->inverse == 1) GetIdx4DnaRev(seq+x, F2);
      if(C->D.bica[bk]==0){
        Compute4DnaFCM(F1); Compute4DnaFCM(F2);
        ca+=CompProbs(F1, n); cb+=CompProbs(F2, n);
        }
      Update4DnaFCM(F1, n, 0); Update4DnaFCM(F2, n, 0);
      if(A->inverse == 1) Update4DnaFCM(F2, 3-(seq[x-F2->ctx]), 1);
      if(++x==BUF_SIZE){ memcpy(seq-GUARD, seq+x-GUARD, GUARD); x=0; }
      ++z;
      }

  Free(seq-GUARD, (BUF_SIZE+GUARD) * sizeof(uint8_t)); 
  Free(buf, BUF_SIZE * sizeof(uint8_t)); 
  rewind(F);
  fprintf(stderr, "[>] Done!                         \n"); // SPACES ARE VALID 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESS STREAM
//
void CompressStream(GFCM *M, FILE *F, CBUF *B, uint8_t sym, uint8_t nSym){
  B->buf[B->idx] = sym;
  GetIdx(B->buf+B->idx-1, M);
  ComputeGFCM(M);
  AESym(sym, (int *)M->freqs, (int)M->freqs[nSym], F);
  UpdateGFCM(M, sym);
  UpdateCBuffer(B);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESS HEADER
//
void CompressHeader(FILE *W, CLASSES *C, Read *R, CBUF *B){
  uint32_t x, s = strlen((char *) R->header1[1]);
  uint8_t sym;
  for(x = 0 ; x < s ; ++x){
    B->buf[B->idx] = sym = C->H.A.numeric[R->header1[1][x]];

    GetIdx(B->buf+B->idx-1, C->H.M[0]);
    ComputeGFCM(C->H.M[0]);
    AESym(sym, (int *)C->H.M[0]->freqs, (int)C->H.M[0]->freqs[C->H.A.nSym], W);
    UpdateGFCM(C->H.M[0], sym);

    UpdateCBuffer(B);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESS SCORES WITH CONSTANT LENGTH
//
void CompressCSScores(FILE *W, CLASSES *C, Read *R, CBUF *B){
  uint32_t x, s = strlen((char *)R->scores)-1, state;
  uint8_t sym;

  // REMOVE "KILLER BEES" [ONLY FROM CONSTANT SIZE READS]
  while(s > 0 && R->scores[s] == '#') --s;

  for(x = 0 ; x < s ; ++x){
    B->buf[B->idx] = sym = C->S.A.numeric[R->scores[x]];
    state = C->S.states[x];

    GetIdx(B->buf+B->idx-1, C->S.M[state]);
    ComputeGFCM(C->S.M[state]);
    AESym(sym, (int *) C->S.M[state]->freqs, (int)
    C->S.M[state]->freqs[C->S.A.nSym], W);
    UpdateGFCM(C->S.M[state], sym);

    UpdateCBuffer(B);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESS DNA BASES
//
void CompressBases(FILE *W, CLASSES *C, Read *R, CBUF *B, CBUF *BE, CBUF *BM,
PARAM *A, SHOTGUN *Gun, GFCM *ME, GFCM *MM){

  uint32_t z, x, m, best = 0, size = strlen((char *) R->bases)-1;
  uint8_t s, n, update;

  for(x = 0 ; x < size ; ++x){

    s = R->bases[x]; 
    n = (s == 'N' ? C->D.lfb : C->D.A.numeric[s]); // BREAK ON C->D.lfb
    B->buf[B->idx] = n;
    update = (A->filter == 1 ? C->D.bica[C->D.idx] : 1);
    Gun->sym[x] = n;

    for(m = 0 ; m < C->D.nFCM ; ++m){

      GetIdx4Dna(B->buf+B->idx-1, C->D.M[m]);
      if(A->inverse == 1) GetIdx4DnaRev(B->buf+B->idx, C->D.M[m]);
      #ifdef REVERSE 
      if(A->reverse == 1) GetIdx4DnaRev(B->buf+B->idx, C->D.M[m]);
      #endif

      ComputeGun(C->D.M[m], Gun->freqs[m][x]);
      Gun->bits[m] += CompGunProbs(Gun->freqs[m][x], n);

      if(update == 1 || C->D.M[m]->ctx < HIGH_CTXBG)
        Update4DnaFCM(C->D.M[m], n, 0);
      if(A->inverse == 1)
        Update4DnaFCM(C->D.M[m], 3-B->buf[B->idx-C->D.M[m]->ctx], 1);
      #ifdef REVERSE
      if(A->reverse == 1)
        Update4DnaFCM(C->D.M[m], B->buf[B->idx-C->D.M[m]->ctx], 1);
      #endif

      #ifdef MEMORY
      if(C->D.M[m]->mode == HASH_TABLE && PeakMem() > A->memory){
        fprintf(stderr, "Reseting DNA Hash model ...\n");
        Reset4DnaModel(C->D.M[m]);
        RestartPeakAndRS();
        fprintf(stderr, "Done!\n");
        }
      #endif
      }
    UpdateCBuffer(B);
    }

  if(A->filter == 1)
    CompressStream(ME, W, BE, C->D.bica[C->D.idx++], 2);

  if(C->D.nFCM != 1){
    best = BestInGun(Gun->bits, C->D.nFCM);
    CompressStream(MM, W, BM, best, C->D.nFCM);
    for(z = 0 ; z < size ; ++z) AESym(Gun->sym[z], (int *)
      Gun->freqs[best][z], (int) Gun->freqs[best][z][4], W);
    }
  else{
    for(z = 0 ; z < size ; ++z) AESym(Gun->sym[z], (int *)
    Gun->freqs[0][z], (int) Gun->freqs[0][z][4], W);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESSION
// 
void Compress(CLASSES *C, PARAM *A, FILE *F, char *fn, char *cn){
  FILE *W = NULL;
  SHOTGUN *Gun = CreateShotgun(C->D.nFCM, C->S.maxLine, 4);
  Read *Read = CreateRead(C->H.maxLine+GUARD, C->S.maxLine+GUARD);
  CBUF *BH   = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD); // HEADERS
  CBUF *BS   = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD); // SCORES
  CBUF *BD   = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD); // DNA
  CBUF *BM   = NULL, *BE = NULL;
  GFCM *ME   = NULL, *MM = NULL;
  uint32_t readIdx = 0;
  uint8_t  *tmp;
  
  if(A->filter == 1){
    if(A->verbose == 1)
      fprintf(stderr, "Creating side information model [Filter]:\n");
    ME = CreateGFCM(DEF_ENT_CTX, DEF_ENT_DEN, 2, A);
    BE = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD); // ENTROPY
    }

  if(C->D.nFCM > 1){
    if(A->verbose == 1)
      fprintf(stderr, "Creating side information model [Models]:\n");
    MM = CreateGFCM(DEF_MOD_CTX, DEF_MOD_DEN, C->D.nFCM, A);
    BM = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD); // MODELS
    }

  Msg(A, "==[ COMPRESSION ]============\n");
  fprintf(stderr, "[>] Compressing ...\n");
  W = Fopen(cn, "w");
  startoutputtingbits();
  start_encode();
  EncodeParameters(C, A, W);

  while(GetRead(F, Read)){

    CompressHeader   (W, C, Read, BH);
    CompressCSScores (W, C, Read, BS);
    CompressBases    (W, C, Read, BD, BE, BM, A, Gun, ME, MM);

    tmp = Read->header1[1];
    Read->header1[1] = Read->header1[0];
    Read->header1[0] = tmp;
    Progress(C->nReads, ++readIdx);
    }

  finish_encode(W);
  doneoutputtingbits(W);
  fclose(W);
  DeleteShotgun(Gun, C->D.nFCM, C->S.maxLine, 4);
  FreeRead(Read);
  RemoveCBuffer(BH);
  RemoveCBuffer(BS);
  RemoveCBuffer(BD);
  if(A->filter == 1){
    RemoveCBuffer(BE);
    FreeGModel(ME);
    }
  if(C->D.nFCM > 1){
    RemoveCBuffer(BM);
    FreeGModel(MM);
    }

  fprintf(stderr, "[>] Done!                           \n"); // SPACES ARE VALID 
  fprintf(stderr, "[i] Compressed %"PRIu64" bases using %llu bytes (%.3g)\n", 
  C->D.A.length, _bytes_output, (double) _bytes_output / C->D.A.length);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UNCOMPRESSION
// 
void Uncompress(CLASSES *C, PARAM *A){
  fprintf(stderr, "[>] Uncompressing ...\n");


  fprintf(stderr, "[>] Done!                         \n"); // SPACES ARE VALID 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// ACTION FOR HEADERS, DNA & SCORES COMPRESSION
//
void ActionC(PARAM *A, char *fn){
  FILE    *F  = Fopen(fn, "r");
  uint64_t b  = FNBytes(F); 
  char    *p  = GetParam(A->level), y[]="yes", n[]="no";
  char    *cn = (char *) Calloc(MFILENM, sizeof(char));
  sprintf(cn, "%s.hawk", fn);
  CheckFile(A->force, cn);

  if(A->verbose){
    fprintf(stderr, "==[ SETTINGS ]===============\n");
    fprintf(stderr, "Compression level ........... %u\n", A->level);
    fprintf(stderr, "Parameters .................. %s\n", p); 
    fprintf(stderr, "Force rewrite ............... %s\n", A->force?   y : n);
    fprintf(stderr, "Use filter .................. %s\n", A->filter?  y : n);
    fprintf(stderr, "Use inversions .............. %s\n", A->inverse? y : n); 
    fprintf(stderr, "Adjust models to data ....... %s\n", A->adjust?  y : n); 
    }
  CLASSES *C = InitClasses();
  C->length = b;
  SetValues(C, A);
  ParseFile(C, F, A);

  // IF nFCM = 1 IT USES THE FILTER, HOWEVER ONLY IN CREATEMODELS THE nFCM IS
  // LOADED. PERHAPS, FIX THIS WITH PRE-LOAD PARAMETERS.
  // TODO: IF nFCM(FROM PARAMETERS) < 2 THEN A->filter = 0 

  // DNA DEEP EVALUATION: AN ALGORITHMIC INFORMATION CONTENT FILTER, FROM PAPER: 
  // D. PRATAS AND A. J. PINHO. EXPLORING DEEP MARKOV MODELS IN GENOMIC DATA 
  // COMPRESSION USING SEQUENCE PRE-ANALYSIS. SIGNAL PROCESSING CONFERENCE 
  // (EUSIPCO), PROCEEDINGS OF THE 22ND EUROPEAN. PP.2395,2399. IEEE, 2014.
  if(A->filter == 1){
    char *nmn = (char *) Calloc(MFILENM, sizeof(char));
    char *nmr = (char *) Calloc(MFILENM, sizeof(char));
    uint32_t seed = time(NULL);
    sprintf(nmn, "Diskn.%ux.hawk", seed);
    sprintf(nmr, "Diskr.%ux.hawk", seed);
    FILE *DK = Fopen(nmn, "w"), *DKR = Fopen(nmr, "w");
    DNA2Disk(C, A, F, DK, DKR);
    fclose(DK); fclose(DKR); 
    DK = Fopen(nmn, "r"); DKR = Fopen(nmr, "r");
    Msg(A, "==[ EVALUATION ]=============\n");
    FCM *F1=Create4DnaFCM(A->fLow,  CalcAlphaDenT(4, A->fLow),  0, 4, A);
    FCM *F2=Create4DnaFCM(A->fHigh, CalcAlphaDenT(4, A->fHigh), 0, 4, A);
    C->D.bica = (uint8_t *) Calloc(C->nReads, sizeof(uint8_t));
    C->D.idx  = 0;
    EvaluateLRDna(C, A, DK,  F1, F2);
    EvaluateRLDna(C, A, DKR, F1, F2);
    Free4DnaModel(F1);
    Free4DnaModel(F2);
    fclose(DK); fclose(DKR);
    unlink(nmn); unlink(nmr); 
    Free(nmn, MFILENM * sizeof(char)); 
    Free(nmr, MFILENM * sizeof(char));
    }

  RestartPeak();      // RESTART MEMORY PEAK [IT WILL IGNORE THE FILTERING PART]
  CreateAuxStates(C);
  CreateModels(C, p, F, A);      // CREATE ALL FCM MODELS [FOR COMPRESSION MODE] 

  // DO NOT USE FILTERING IF THE NUMBER OF MODELS IS NOT AT LEAST 2 AND REPORT A
  // ERROR IF EQUAL OR HIGHER THAN MAX_AL
  if(C->D.nFCM < 2) A->filter = 0;
  else if(C->D.nFCM >= MAX_AL){
    fprintf(stderr, "[x] Error: too many models!\n"); exit(1);
    }

  C->D.idx  = 0;
  // HEADERS, DNA SEQUENCE AND QUALITY-SCORES COMPRESSION
  Compress(C, A, F, fn, cn);

  // FREE MODELS & PRINT MEMORY INFORMATION FOR THE UNCOMPRESSION PROCESS
  FreeModels(C, A);
  FreeAlphabets(C);
  DeleteAuxStates(C);
  FreeClasses(C);
  Free(cn, MFILENM * sizeof(char));
  #ifdef MEMORY
  PrintRAM(A->memory);
  #endif
  fclose(F);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// ACTION FOR UNCOMPRESSION
//
void ActionD(PARAM *A, char *fn){
  FILE    *F = Fopen(fn, "r");
  char    *dn = ReplaceSubStr(fn, ".hawk", ".d");
  CheckFile(A->force, dn);

  CLASSES *C = InitClasses();
  Uncompress(C, A);
  //DecodeHeader(F);
  //CreateModels(C, p, F, A);

  FreeClasses(C); 
  Free(dn, MAX_STR * sizeof(char));
  fclose(F);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// MAIN $ COMMAND INTERFACE $
//
int main(int argc, char *argv[]){
  char **p = *&argv;
  PARAM *A = Calloc(1, sizeof(PARAM));

  if(argc < 2 || ArgBin(DEF_HELP, argv, argc, "-h")){
    fprintf(stderr, 
    "Usage: Hawk [OPTION]... [FILE]                                 \n" 
    "Compress or uncompress a FASTQ file (by default, compress).    \n"
    "                                                               \n"
    "Non-mandatory arguments:                                       \n"
    "                                                               \n"
    "  -h            give this help,                                \n"
    "  -v            verbose mode (more information),               \n" 
    "  -V            display version number,                        \n"
    "  -P            display internal parameters from levels,       \n"
    "  -f            force overwrite of output,                     \n"
    "  -fl <ORDER>   lower context order for filter,                \n"
    "  -fh <ORDER>   higher context order for filter,               \n"
    "  -fn           do NOT use filter,                             \n" 
    "  -sh <STATES>  number of multi-states for headers,            \n" 
    "  -ss <STATES>  number of multi-states for quality scores,     \n" 
    "  -l  <LEVEL>   compression level [1,...,9],                   \n" 
    "  -a            adjust context to data,                        \n" 
    "  -i            use inversions (DNA sequence only),            \n" 
    #ifdef REVERSE
    "  -r            use reversions (DNA sequence only),            \n" 
    #endif
    #ifdef MEMORY
    "  -m  <MEMORY>  maximum hash memory for deepest model (in MB). \n" 
    #endif
    "  -c            use CCH instead of hash (deepest contexts).    \n" 
    "                                                               \n"
    "Mandatory compress arguments:                                  \n"
    "                                                               \n"
    "  <FILE>       file to compress (last argument).               \n"
    "                                                               \n"
    "Mandatory uncompress arguments:                                \n"
    "                                                               \n"
    "  -d           uncompress mode,                                \n" 
    "  <FILE>       file to uncompress (last argument).             \n"
    "                                                               \n"
    "Report bugs to <{pratas,ap}@ua.pt>.                            \n");
    return EXIT_SUCCESS;
    }

  if(ArgBin(DEF_VERSION, p, argc, "-V")){
    fprintf(stderr, "Hawk %u.%u\n"
    "Copyright (C) 2015 University of Aveiro.\nThis is Free software. \nYou "
    "may redistribute copies of it under the terms of the GNU General \n"
    "Public License v2 <http://www.gnu.org/licenses/gpl.html>.\nThere is NO "
    "WARRANTY, to the extent permitted by law.\nWritten by Diogo Pratas and "
    "Armando J. Pinho.\n", RELEASE, VERSION);
    return EXIT_SUCCESS;
    }

  if(ArgBin(DEF_PARAM, p, argc, "-P")){
    PrintParam();
    return EXIT_SUCCESS;
    }

  A->verbose  = ArgBin(DEF_VERBOSE, p, argc, "-v");
  A->force    = ArgBin(DEF_FORCE,   p, argc, "-f");
  A->level    = ArgNum(DEF_LEVEL,   p, argc, "-l",  MIN_LEVEL, MAX_LEVEL);
  A->fLow     = ArgNum(DEF_FL_CTX,  p, argc, "-fl", MIN_F_CTX, MAX_F_CTX);
  A->fHigh    = ArgNum(DEF_FH_CTX,  p, argc, "-fh", MIN_F_CTX, MAX_F_CTX);
  A->filter   = ArgBin(DEF_FILTER,  p, argc, "-fn");
  A->hNStates = ArgNum(DEF_STATES,  p, argc, "-sh", MIN_STATES, MAX_STATES);
  A->sNStates = ArgNum(DEF_STATES,  p, argc, "-ss", MIN_STATES, MAX_STATES);
  A->inverse  = ArgBin(DEF_INVERSE, p, argc, "-i");
  #ifdef REVERSE
  A->reverse  = ArgBin(DEF_REVERSE, p, argc, "-r");
  #endif
  A->adjust   = ArgBin(DEF_ADJUST,  p, argc, "-a");
  A->mode     = ArgBin(DEF_MODE,    p, argc, "-c");
  #ifdef MEMORY
  A->memory   = (uint64_t) ArgNum(DEF_MEM,  p, argc, "-m", MIN_MEM, MAX_MEM) *
               1048576;
  #endif

  if(ArgBin(DEF_ACTION, p, argc, "-d")) ActionD(A, argv[argc-1]);
  else{ ActionC(A, argv[argc-1]); }

  Free(A, sizeof(PARAM));
  PrintCurrMem();
  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
