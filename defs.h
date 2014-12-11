#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>

#define MAX(a,b) (((a)>(b))?(a):(b))

uint64_t garbage;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// HAWK INFORMATION
//
#define HAWKFILE    1337411
#define RELEASE     1
#define VERSION     1
#define YEAR        2015

// SYSTEM DEFENITIONS
//
#define MAX_BUF     1000000
#define ONE_GB      1073741824
#define SCACHE      32
#define GUARD       30
#define NSYM        4
#define MAXC        65535
#define BUF_SIZE    65535
#define MFILENM     256
#define BIN_LINE    512
#define DEF_HELP    0
#define DEF_ACTION  0    // DEFAULT = COMPRESSION
#define DEF_VERBOSE 0
#define DEF_PARAM   0
#define DEF_FILTER  1
#define DEF_FORCE   0
#define DEF_LEVEL   5
#define DEF_MODE    1     // DEFAULT = HASH MODE
#define DEF_INVERSE 0
#define DEF_REVERSE 0
#define DEF_ADJUST  0
#define DEF_VERSION 0
#define DEF_ENT_CTX 5
#define DEF_ENT_DEN 1
#define DEF_MOD_CTX 5
#define DEF_MOD_DEN 1
#define DEF_STATES  3
#define MIN_STATES  1
#define MAX_STATES  200
#define MAX_AL      256
#define INVALID_S   256
#define HMIN_CTX    1
#define HMAX_CTX    5
#define EDMIN_CTX   1
#define EDMAX_CTX   20
#define DMIN_CTX    1
#define DMAX_CTX    20
#define SMIN_CTX    1
#define SMAX_CTX    5
#define MIN_F_CTX   1
#define MAX_F_CTX   18
#define DEF_FL_CTX  4
#define DEF_FH_CTX  13
#define MIN_MEM     1024
#define MAX_MEM     131072
#define DEF_MEM     5120
#define HIGH_CTXBG  13 
#define SCACHE      32
#define MIN_LEVEL   0
#define MAX_LEVEL   9
#define HASH_TABLE  2
#define CCH_TABLE   1
#define ARRAY_TABLE 0
#define DEEP_CTX    268435457  // (4^14)+1

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

