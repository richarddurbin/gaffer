/*  File: kmerhash.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: fixed length DNA string hash set package (e.g. syncmers)
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 29 13:07 2024 (rd109)
 * Created: Tue Sep  3 19:39:02 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"  // includes utils.h

// fixed length kmer hasing, based on dict.[ch]
// does not support removal
// pack sequences to 2bit encoding as in ONElib, so can read/write ONEfiles natively
// kmers are stored in canonical orientation: kmer < reverseComplement(kmer) 
// new kmers are added from 1 .. kmerHashMax() inclusive: so 1-based numbering
// this allows to return matches in negative orientation as -index
// kmerHashAdd() and kmerHashFind() are not threadsafe, but can use kmerHasFindThreadSafe() in threads

typedef struct {
  int    len ;    // length of dna sequences stored
  int    dim ;    // dimension of table (size is 2^dim)
  I64   *table ;  // the main table indexed by the hash
  I64    mask ;   // mask limited to tsize bits
  I64    max ;	  // current number of entries
  int    plen ;   // packed sequence length = (len+31) >> 5
  U64    psize ;  // max number of elements in pack before doubling
  U64   *pack ;   // the packed sequences
  char  *seqbuf ; // size len+1, for printing out sequences
  U64    finds ;  // stats: number of finds
  U64    deltas ; // stats: number of deltas (not in remapping)
} KmerHash ;

KmerHash *kmerHashCreate (U64 initialSize, int len) ;
void      kmerHashDestroy (KmerHash *kh) ;
bool      kmerHashAdd (KmerHash *kh, char *dna, I64 *index) ;  // true if addded - always fill *index
bool      kmerHashFind (KmerHash *kh, char *dna, I64 *index) ; // true if found
bool      kmerHashFindThreadSafe (KmerHash *kh, char *dna, I64 *index, U64 *buf) ;
// buf must point to user memory of size kh->plen or larger
char*     kmerHashSeq (KmerHash *kh, U64 i) ; // retrieve the i'th sequence
#define   kmerHashMax(kh)  ((kh)->max)        // count of stored kmers

#include "ONElib.h"

bool      kmerHashWriteOneFile (KmerHash *kh, OneFile *of) ;
KmerHash *kmerHashReadOneFile (OneFile *of) ;

/*********** end of file ***********/



