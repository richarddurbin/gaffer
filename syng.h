/*  File: syng.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 29 18:34 2024 (rd109)
 * Created: Mon May 29 08:19:18 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "hash.h"
#include "kmerhash.h"
#include "ONElib.h"

typedef struct {
  int w, k, seed ;
} Params ;

typedef struct {
  int n ;		// number of copies of the node
} Node ;

typedef struct {
  I64 sa, sb ;          // syncs, encode as -x-1 if -ve orientation, made canonical as below
  I64 n ;		// number of copies of the link
  I64 overlap ;
} Link ;

Params params ;
static int PARAMS_K_DEFAULT = 16 ;
static int PARAMS_W_DEFAULT = 1023 ;
static int PARAMS_SEED_DEFAULT = 7 ;
Array    nodes ;
Array    links ;

static char *syngSchemaText =
  "1 3 def 1 0               schema for syng\n"
  ".\n"
  "P 3 seq                   SEQUENCE\n"
  "S 4 sync                  syncmer sequence\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  ".\n"
  "O S 1 3 DNA               sequence of the syncmer\n"
  "D T 1 3 INT               KmerHash table index for read/write of hash table\n"
  ".\n"
  "O N 1 3 INT               length of graph node\n"
  "D K 1 3 INT               kmer count\n"
  "D E 4 4 CHAR 3 INT 3 INT 3 INT   +/-, following node (- if reversed), offset, count\n"
  "D B 2 4 CHAR 8 INT_LIST          +/-, GBWT list of node indices from E lines\n"
  "D C 2 4 CHAR 8 INT_LIST          +/-, GBWT list of run-length counts\n"
  ".\n"
  "P 7 syncseq               SYNCMER SEQUENCE\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  "O S 1 8 INT_LIST          sequence: list of syncmer ids\n"
  "D P 1 8 INT_LIST          positions of the syncmers\n"
  "D D 1 6 STRING            orientations of the syncmers\n" // better compression with directions
  "D R 2 3 INT 3 INT         origin of sequence - reference file number and position in file \n"
  ".\n"
  "P 5 khash                 KMER HASH\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  "D t 3 3 INT 3 INT 3 INT   max, len, dim for KmerHash table\n"
  "O S 1 3 DNA               packed sequences aligned to 64-bit boundaries\n" 
  "D L 1 8 INT_LIST          locations in the table\n"
  "D C 1 8 INT_LIST          kmer counts\n"
  ".\n"
  ;

/****************** end of file ********************/
