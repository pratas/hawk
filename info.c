#include <stdio.h>
#include <stdlib.h>
#include "info.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// ENCODE PARAMETERS
//
void EncodeParameters(CLASSES *C, PARAM *A, FILE *W){
  uint32_t x;
  // BASIC HAWK INFORMATION
  //
  WriteNBits(HAWKFILE,             21, W);
  WriteNBits(VERSION,               8, W);
  WriteNBits(RELEASE,               8, W);
  WriteNBits(YEAR,                 12, W);
  //
  // BASIC READS INFORMATION
  WriteNBits(C->length,            64, W);
  WriteNBits(C->nReads,            64, W);
  //
  // HEADER:
  WriteNBits(C->H.nFCM,             8, W);
  WriteNBits(C->H.type,             8, W);
  WriteNBits(C->H.extra,            8, W);
  WriteNBits(C->H.A.length,        64, W);
  WriteNBits(C->H.A.nSym,           8, W);
  for(x = 0 ; x < C->H.A.nSym ; ++x)
    WriteNBits(C->H.A.symbolic[x],  8, W);
  for(x = 0 ; x < C->H.nFCM ; ++x){
    WriteNBits(C->H.M[x]->ctx,     32, W);
    WriteNBits(C->H.M[x]->aDen,    32, W);
    WriteNBits(C->H.M[x]->nSym,     8, W);
    }
  //
  // DNA SEQUENCE
  WriteNBits(A->filter,             1, W);
  if(A->filter == 1){
    WriteNBits(DEF_ENT_CTX,         8, W);
    WriteNBits(DEF_ENT_DEN,        32, W);
    }
  WriteNBits(C->D.nFCM,             8, W);
  if(C->D.nFCM>1){
    WriteNBits(DEF_MOD_CTX,         8, W);
    WriteNBits(DEF_MOD_DEN,        32, W);
    }
  WriteNBits(C->D.lfb,              8, W);
  WriteNBits(C->D.A.length,        64, W);
  WriteNBits(C->D.A.nSym,           8, W);
  for(x = 0 ; x < C->D.A.nSym ; ++x)
    WriteNBits(C->D.A.symbolic[x],  8, W);
  for(x = 0 ; x < C->D.nFCM ; ++x){
    WriteNBits(C->D.M[x]->ctx,     32, W);
    WriteNBits(C->D.M[x]->aDen,    32, W);
    WriteNBits(C->D.M[x]->nSym,     8, W);
    WriteNBits(C->D.M[x]->rev,      8, W);
    }
  //
  // QUALITY SCORES:
  WriteNBits(C->S.nFCM,             8, W);
  WriteNBits(C->S.maxLine,         64, W);
  WriteNBits(C->S.dynamic,          8, W);
  WriteNBits(C->S.not0N,            8, W);
  WriteNBits(C->S.A.nSym,           8, W);
  for(x = 0 ; x < C->S.A.nSym ; ++x)
    WriteNBits(C->S.A.symbolic[x],  8, W);
  for(x = 0 ; x < C->S.nFCM ; ++x){
    WriteNBits(C->S.M[x]->ctx,     32, W);
    WriteNBits(C->S.M[x]->aDen,    32, W);
    WriteNBits(C->S.M[x]->nSym,     8, W);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

