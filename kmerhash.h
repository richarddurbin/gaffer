/*  File: kmerhash.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: fixed length DNA string hash set package (e.g. syncmers)
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 25 14:38 2024 (rd109)
 * Created: Tue Sep  3 19:39:02 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"  // includes utils.h

typedef struct {
  int    len ;    // length of dna sequences stored
  int    dim ;    // dimension of table (size is 2^dim)
  U64   *table ;  // the main table indexed by the hash
  U64    mask ;   // mask limited to tsize bits
  U64    max ;	  // current number of entries
  int    plen ;   // packed sequence length = (len+31) >> 5
  U64    psize ;  // max number of elements in pack before doubling
  U64   *pack ;   // the packed sequences
  char  *seqbuf ; // size len+1, for printing out sequences
  U64    finds ;  // stats: number of finds
  U64    deltas ; // stats: number of deltas (not in remapping)
} KmerHash ;

// idea is to make a custom library based on dict.[ch]
// for now don't support removal
// use SeqPack from seqio.h
// flips to canonical and reports orientation of match - if palindromic always reports false

KmerHash *kmerHashCreate (U64 initialSize, int len) ;
void     kmerHashDestroy (KmerHash *kh) ;
bool     kmerHashAdd (KmerHash *kh, char *dna, U64 *index, bool *isRC) ; // true if added, always fill index
bool     kmerHashFind (KmerHash *kh, char *dna, U64 *index, bool *isRC) ; // true if found
char*    kmerHashSeq (KmerHash *kh, U64 i) ;
#define  kmerHashMax(kh)  ((kh)->max)

/*********** end of file ***********/



