/*  File: sbwt.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 16 13:04 2024 (rd109)
 * Created: Tue Apr 16 08:28:54 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"

#define CHECK // this does some QC checks should not be necessary when all is working

typedef unsigned int Sync ;
typedef unsigned char U8 ;

// we use 0 as a terminator for lists of syncs, so can't have a real sync with value 0
// we will use Sync == 1 for terminators $, and Sync == 2 for repeats

#define TERMINATOR 1
#define REPEAT 2

typedef struct { // one of these per Sync
  Sync *syncs ;        // 0-terminated list of the possible preceding symbols
  int  *occ ;          // for each of these, how many syncs < this followed it
  U8   *bwt ;          // 0-terminated array of pairs (n,k), where k is an index in sym
} SBWT ;
// strictly, bwt does not need to be 0-terminated, and we don't need to check *b in backwards
// we do this only when CHECK is set

typedef struct {
  int   nSync ;
  SBWT *sbwt ;         // length nSync+1 (0 is unused)
  int   nSeq ;
  Sync *last ;         // length nSeq: last symbol in sequence[i]
  int  *occ ;          // length nSeq: number of last[i] occurrences at last position for j < i
  int  *rev ;          // index of i in lexicographic ordering of the reads
} SBWTset ;

bool backwardsMatch1 (SBWT *sb, Sync sx, int *i)
{ // returns true and updates *i to position in sbwt[sx] if sx is found
  // returns false if sx is not found
  
  // first find kx = the index of sx in sb->syncs
  int   kx = 0 ;
  Sync *s = sb->syncs ;
  while (*s && *s != sx) { ++kx ; ++s ; }
  if (!*s) return false ;  // symbol sx not found

  // next find r = how many kx are in the bwt before i
  int r = sb->occ[kx] ; // will be return value
  U8 *b = sb->bwt ;
  int n = 0 ;              // position in bwt
#ifdef CHECK
  if (i1 < 0) die ("i1 %d < 0 in backwards1", i1) ;
  while (*b && (n + *b < *i)) { n += *b ; if (b[1] == kx) r += *b ;  b += 2 ; }
  if (!*b) die ("i %d > n %d in backwards1", *i, n) ; // this should not happen
#else
  while (n + *b < *i) { n += *b ; if (b[1] == kx) r += *b ;  b += 2 ; }
#endif
  if (b[1] == kx) r += *i - n ;
  *i = r ;
  return true ;
}

int backwardsMatch2 (SBWT *b, Sync sx, int *i1, int *i2) // returns numMatches and updates i1, i2
{ // similar to backwards1(), but updates two pointers, as required for string matching

  int   kx = 0 ;
  Sync *s = sb->syncs ;
  while (*s && *s != sx) { ++kx ; ++s ; }
  if (!*s) return 0 ;  // symbol sx not found

  int r = sb->occ[kx] ; // will be return value
  U8 *b = sb->bwt ;
  int n = 0 ;              // position in bwt
#ifdef CHECK
  if (i1 < 0) die ("i1 %d < 0 in backwards2", i1) ;
  while (*b && (n + *b < *i1)) { n += *b ; if (b[1] == kx) r += *b ;  b += 2 ; }
  if (!*b) die ("i1 %d > n %d in backwards2", *i1, n) ; // this should not happen
#else
  while (n + *b < *i1) { n += *b ; if (b[1] == kx) r += *b ;  b += 2 ; }
#endif
  if (b[1] == kx) *i1 = r + *i1 - n ; else *i1 = r ;
  
#ifdef CHECK
  if (i2 < i1) die ("i2 %d < i1 %d in backwards2", i2, i1) ;
  while (*b && (n + *b < *i2)) { n += *b ; if (b[1] == kx) r += *b ;  b += 2 ; }
  if (!*b) die ("i2 %d > n %d in backwards2", *i2, n) ;
#else
  while (n + *b < *i2) { n += *b ; if (b[1] == kx) r += *b ;  b += 2 ; }
#endif
  if (b[1] == kx) *i2 = r + *i2 - n ; else *i2 = r ;
  
  return *i2 - *i1 ;
}

Sync backwards (SBWT *sb, int *i) // returns the Sync at i and updates i
{ 
  U8 *b = sb->sbwt ;
  int n = 0 ;
#ifdef CHECK
  while (*b && (n + *b < *i)) { n += *b ; b += 2 ; }
  if (!*b) die ("i %d > n %d in backwards1", *i, n) ; // this should not happen
#else
  while (n + *b < *i) { n += *b ; b += 2 ; }
#endif
  int kx = b[1] ;
  b = sb->sbwt ;
  n = 0 ;
  int r = sb->occ[kx] ;
  while (n + *b < *i) { n += *b ; if (b[1] == kx) r += *b ; b += 2 ; }
  *i = r + *i - n ;
  return sb->syncs[kx] ;
}

void sbwtSize (SBWT *sb)
{
  int n = 0 ;
  U8 *b = sb->bwt ;
  while (*b) { n += *b ; b += 2 ; }
  return n ;
}

int countMatches (Sync *string, int len, SBWTset *set)
{
  if (len <= 0) return 0 ;
  int i1 = 0, i2 = sbwtSize (set->sbwt + string[--len]) ;
  while (len)
    { if (!backwardsMatch2 (set->sbwt + string[len], string[len-1], &i1, &i2)) return 0 ;
      --len ;
    }
  return i2-i1 ;
}

Sync *extractSequence (SBWTset *set, int k) // extracts a 0-terminated Sync string for the sequence
{
  static Array a = 0 ; if (!a) arrayCreate (256, Sync) ; else a->max = 0 ;
  Sync sx = set->last[k] ;
  int   i = set->occ[k] ;
  while (sx != TERMINATOR)
    { array(a,arrayMax(a),Sync) = sx ;
      sx = backwards (set->sbwt + sx, &i) ;
    }
#ifdef CHECK
  if (set->rev[i] != k) die ("some problem in extractRead: rev[i %d] != k %d", i, k) ;
#endif

  //now reverse the array and zero-terminate, then return it
  i = 0 ;
  int j = arrayMax(a) - 1 ;
  for (i = 0 ; i < j ; ++i ; --j)
    { sx = arr(a,i,Sync) ; arr(a,i,Sync) = arr(a,j,Sync) ; arr(a,j,Sync) = sx ; }
  array(a,arrayMax(a),Sync) = 0 ;
  return arrp(a,0,Sync) ; // 0-terminated array of Syncs
}

int findRead (SBWTset *set, SBWT *sb, int i, int *pos) // returns the index of the read at i in sb
{
  *pos = 0 ;
  Sync sx = backwards (sb, &i) ;
  while (sx != TERMINATOR)
    { ++*pos ;
      sx = backwards (set->sbwt + sx, &i) ;
    }
  return set->rev[i] ;
}

SBWTset *sbwtSetCreate (Sync **reads, int nReads)
{
  
}



 

  
