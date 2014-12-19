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
#include <assert.h>

#include "defs.h"
#include "mem.h"
#include "misc.h"
#include "args.h"
#include "param.h"
#include "reads.h"
#include "dna.h"
#include "hash.h"
#include "phash.h"
#include "cch.h"
#include "sfcm.h"
#include "classes.h"
#include "models.h"
#include "gun.h"
#include "info.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

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
  uint32_t x, s = strlen((char *) R->header1[1]), state;
  uint8_t sym;
  uint64_t idx = C->H.M[0]->idx;
  
  for(x = 0 ; x < s ; ++x){
    B->buf[B->idx] = sym = C->H.A.numeric[R->header1[1][x]];
    state = (uint32_t) C->H.states[x];

    idx = GetIdxA(B->buf+B->idx-1, C->H.M[state]);
    ComputeGFCM(C->H.M[state]);
    AESym(sym, (int *) C->H.M[state]->freqs, (int)
    C->H.M[state]->freqs[C->H.A.nSym], W);
    UpdateGFCM(C->H.M[state], sym);
    UpdateCBuffer(B);
    if(x < s-1) C->H.M[(uint32_t) C->H.states[x+1]]->idx = idx;
    }
  C->H.M[0]->idx = idx;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESS SCORES WITH CONSTANT LENGTH
//
void CompressScores1(FILE *W, CLASSES *C, Read *R, CBUF *B){
  uint32_t x, s = strlen((char *)R->scores), state;
  uint8_t sym;
  uint64_t idx = C->S.M[0]->idx;

  while(s > 0 && R->scores[s-1] == '#') --s; // REMOVE "KILLER BEES"

  for(x = 0 ; x < s ; ++x){
    B->buf[B->idx] = sym = C->S.A.numeric[R->scores[x]];
    state = (uint32_t) C->S.states[x];
    idx = GetIdxA(B->buf+B->idx-1, C->S.M[state]);
    ComputeGFCM(C->S.M[state]);
    AESym(sym, (int *) C->S.M[state]->freqs, (int)
    C->S.M[state]->freqs[C->S.A.nSym], W);
    UpdateGFCM(C->S.M[state], sym);
    UpdateCBuffer(B);
    if(x < s-1) C->S.M[(uint32_t) C->S.states[x+1]]->idx = idx;
    }
  C->S.M[0]->idx = idx;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESS SCORES
//
// IT HAS BEEN (SLIGHTLY) ADAPTED FROM JAMES BONFIELD PROGRAM: FQZ_COMP (SCORES 
// MODEL IS THE SAME).
//
// BONFIELD, JAMES K., AND MATTEW V. MAHONEY. "COMPRESSION OF FASTQ AND SAM 
// FORMAT SEQUENCING DATA." PLOS ONE 8.3 (2013): e59190.
// 
// https://sourceforge.net/projects/fqzcomp
//
#define QMAX  64 // 128 // KEEP AS POWER OF 2
#define QBITS 12 // SIZE OF THE CONTEXT 2x QUALS
#define QSIZE (1<<QBITS)

void CompressScores2(FILE *W, SFCM *M, char *qual, int len){
  unsigned int last = 0;
  int delta = 5, i, len2 = len, q1 = 0, q2 = 0;

  while(len2 > 0 && qual[len2-1] == '#') --len2; // REMOVE "KILLER BEES"

  for(i = 0 ; i < len2 ; ++i){
    uint8_t q = (qual[i]-'!') & (QMAX-1);
    EncodeSFCM(W, M, last, q);

    last   = ((MAX(q1, q2)<<6)+q)&((1<<QBITS)-1);
    last  += (q1==q2)<<QBITS;
    delta += (q1>q)*(q1-q);
    last  += (MIN(7*8, delta)&0xf8)<<(QBITS-2);
    last  += (MIN(i+15, 127)&(15<<3))<<(QBITS+1);

    q2 = q1; 
    q1 = q;
    }

  if(len != len2) EncodeSFCM(W, M, last, QMAX-1);
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
    n = (s == 'N' ? C->D.lfb : C->D.A.numeric[s]); // TODO: BREAK ON C->D.lfb
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
  SFCM *S    = CreateSFCM(1048576, QMAX);
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
    CompressHeader  (W, C, Read, BH);
//  CompressScores1 (W, C, Read, BS);
    CompressScores2 (W, S, (char*) Read->scores, strlen((char*)Read->scores)-1);
    CompressBases   (W, C, Read, BD, BE, BM, A, Gun, ME, MM);

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
  DeleteSFCM(S);
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
    "  -m  <MEMORY>  maximum hash memory for deepest model (in MB), \n" 
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
