/*  File: kmerhash.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: fixed length DNA string hash set package (e.g. syncmers)
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 27 22:51 2024 (rd109)
 * Created: Tue Sep  3 19:38:07 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "kmerhash.h"

static inline void Compress_DNA (int len, char *s, U64 *u) ;    // from s to u
static inline void Compress_DNA_RC (int len, char *s, U64 *u) ; // from s to u
static inline void Uncompress_DNA (int len, U64 *u, char *t) ;  // from u to t

KmerHash *kmerHashCreate (U64 initialSize, int len)
{
  if (!(len & 0x01)) die ("kmerHash len %d must be odd", len) ;
  KmerHash *kh = new0(1, KmerHash) ;
  kh->len = len ;
  kh->dim = 20 ;
  U64 size ;
  for (size = 1 << kh->dim ; size < initialSize ; ++kh->dim, size <<= 1) ;
  kh->table = new0(size, U64) ;
  kh->mask = size - 1 ;
  kh->plen = (len+31) >> 5 ;
  kh->psize = size * 0.3 ;
  kh->pack = new0(kh->plen*kh->psize, U64) ;
  kh->seqbuf = new0(len+1,char) ;
  return kh ; 
}

void kmerHashDestroy (KmerHash *kh)
{
  free (kh->table) ;  totalAllocated -= ((U64)1 << kh->dim) * sizeof(U64) ;
  free (kh->pack) ;   totalAllocated -= (kh->psize*kh->plen) * sizeof(U64) ;
  free (kh->seqbuf) ; totalAllocated -= kh->len + 1 ;
  free (kh) ;         totalAllocated -= sizeof (KmerHash) ;
}

static U8 comp[] = {   /* sends N (indeed any non-CGT) to A, except 0,1,2,3 are maintained */
   3,   2, 1,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0,   0, 0,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0,   0, 0,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0,   0, 0,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0, 'T', 0, 'G',   0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0,
   0,   0, 0,   0, 'A', 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0,
   0, 't', 0, 'g',   0, 0, 0, 'c', 0, 0, 0, 0, 0, 0, 'n', 0,
   0,   0, 0,   0, 'a', 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0
} ;

static inline void packDNA (char *dna, U64 *u, int len, bool *isRC)
{
  int x = -1, y = len ;
  while (dna[++x] == comp[dna[--y]]) ; // first find orientation
  if (dna[x] < comp[dna[y]])
    { if (isRC) *isRC = false ;
      Compress_DNA (len, dna, u) ;
    }
  else
    { if (isRC) *isRC = true ;
      Compress_DNA_RC (len, dna, u) ;
    }
}

static inline bool isMatch (U64 *u, U64 *v, int n)
{
  while (n--) if (*u++ != *v++) return false ;
  // { printf (" mismatch %d %llx != %llx", n, *--u, *--v) ; return false ; }
  return true ;
}

static inline U64 hashDelta (U64 *u, int plen, int dim)
{ // very simple hash is enough here I think
  U64 v = *u ;
  while (--plen) v = v ^ u[plen] ;
  plen = 64 / dim ;
  while (--plen) v = v ^ (v >> dim) ;
  v |= 1 ; // ensure that delta is odd so coprime with 2^dim
  return v ;
}

#define packseq(kh,i) ((kh)->pack + (i)*(kh)->plen)

static inline bool find (KmerHash *kh, U64 *u, U64 *index, bool *isRC, U64 *ploc)
{
  ++kh->finds ;
  U64 loc = *u & kh->mask ; // initially try the simplest thing
  U64 delta = 0 ;
  U64 x = kh->table[loc] ;
  //  printf ("find: loc %0llx x %llu", loc, x) ;
  while (x)
    { if (isMatch (u, packseq(kh,x-1), kh->plen)) // found
	{ if (index) *index = x-1 ; return true ; }
	//	{ if (index) *index = x-1 ; putchar ('\n') ; return true ; }
      if (!delta) delta = hashDelta (u, kh->plen, kh->dim) ;
      loc = (loc + delta) & kh->mask ;
      ++kh->deltas ;
      x = kh->table[loc] ;
      //      printf (" loc %0llx x %llu", loc, x) ;
    }
  *ploc = loc ;
  //  putchar ('\n') ;
  return false ; // not found
}

bool kmerHashFind (KmerHash *kh, char *dna, U64 *index, bool *isRC)
{
  U64 *u = packseq(kh,kh->max) ; // pack into the next free slot in pack[]
  packDNA (dna, u, kh->len, isRC) ;
  U64 loc ;
  return find (kh, u,  index, isRC, &loc) ;
}

bool kmerHashFindThreadSafe (KmerHash *kh, char *dna, U64 *index, bool *isRC, U64 *u)
{
  packDNA (dna, u, kh->len, isRC) ; // pack into user-provided space
  U64 loc ;
  return find (kh, u,  index, isRC, &loc) ;
}

bool kmerHashAdd (KmerHash *kh, char *dna, U64 *index, bool *isRC)
{
  U64 *u = packseq(kh,kh->max) ; // pack into the next free slot in pack[]
  packDNA (dna, u, kh->len, isRC) ;
  U64 loc ;
  //  printf ("add %llu\n", kh->max) ;
  bool isFound = find (kh, u, index, isRC, &loc) ;
  if (isFound) return false ;
  
  if (index) *index = kh->max++ ;
  kh->table[loc] = kh->max ; // NB this is *index + 1, as required
  //  printf ("setting kh->table[%0llx] = %llu\n", loc, kh->max) ;
  
  if (kh->max == kh->psize)  // need to double
    { ++kh->dim ;
      U64 *newTable = new0 (1 << kh->dim, U64) ;
      kh->mask = (1 << kh->dim) - 1 ;
      U64  i ;
      for (i = 0 ; i < kh->max ; ++i) // remap all the packed sequences into newTable
	{ loc = *packseq(kh,i) & kh->mask ;
	  U64 delta = 0 ;
	  while (newTable[loc])
	    { if (!delta) delta = hashDelta(packseq(kh,i), kh->plen, kh->dim) ;
	      loc = (loc + delta) & kh->mask ;
	    }
	  newTable[loc] = i+1 ;
	}
      free (kh->table) ; totalAllocated -= ((U64)1 << (kh->dim-1)) * sizeof(U64) ;
      kh->table = newTable ;
      resize (kh->pack, kh->plen*kh->psize, 2*kh->plen*kh->psize, U64) ;
      kh->psize *= 2 ;
      // printf ("doubled at max = %llu to %llu\n", kh->max, kh->psize) ;
    }

  // printf ("add: loc %llx index %lld max %lld\n", loc, *index, kh->max) ;
  return true ;
}

char* kmerHashSeq (KmerHash *kh, U64 i)
{
  if (i >= kh->max) die ("out of range in kmerHashSeq: %lld >= %lld", i, kh->max) ;
  Uncompress_DNA (kh->len, packseq(kh,i), kh->seqbuf) ;
  return kh->seqbuf ;
}

#ifdef TEST

#include <stdlib.h>

int len ;

static inline void seqSim (char *s)
{ int i ;
  for (i = 0 ; i < len ; ++i) *s++ = random() & 0x3 ;
}

int main (int argc, char *argv[])
{
  U64 i, count ;
  int target ;
  timeUpdate (stdout) ;
  if (argc != 4) die ("Usage: dnahash <len> <number> <target>") ;
  len = atoi(argv[1]) ; count = atoi(argv[2]) ; target = atoi(argv[3]) ;

  KmerHash *kh = kmerHashCreate (0, len) ;

  /* test the compression/decompression code
  char *dna = "cagtatcttagatcctatgggagtacgagcgagcactagcggagggaggagcatcagcacgagcgagggcacaatcta" ;
  char buf[50] ;
  U64 u[2] ;
  printf ("    %s\n", dna) ;
  for (i = 1 ; i < 40 ; ++i)
    { Compress_DNA (i, dna, u) ; Uncompress_DNA (i, u, buf) ; printf ("%2d  %-40s", (int)i, buf) ;
      Compress_DNA_RC (i, dna, u) ; Uncompress_DNA (i, u, buf) ; printf ("  %40s\n", buf) ;
    }
  printf ("    %s\n", dna) ;
  exit (0) ;
  */
  
  char *data = new(count*len, char) ;
  char *s = data ;
  for (i = 0 ; i < count ; ++i, s += len) seqSim(s) ;
  printf ("simulated %llu sequences of length %d\n", count, len) ;
  timeUpdate (stdout) ;
  
  bool isRC ;
  U64 oriTotal[2] ; oriTotal[0] = 0 ; oriTotal[1] = 0 ;
  U64 index ;
  s = data ;
  for (i = 0 ; i < count ; ++i, s += len)
    { kmerHashAdd (kh, s, &index, &isRC) ;
      if (i < 4) printf ("orientation[%llu] is %c\n", i, isRC ? '-' : '+') ;
      ++oriTotal[isRC] ;
    }

  char *acgt = "acgt" ;
  s = data + target*len ;
  printf ("seq %d is  ", target) ; for (i = 0 ; i < len ; ++i) putchar (acgt[s[i]]) ; putchar ('\n') ;
  printf ("hash %d is %s\n", target, kmerHashSeq (kh,target)) ;
  
  printf ("loaded: max %lld finds %lld deltas %lld + %llu - %llu\n",
	  kh->max, kh->finds, kh->deltas, oriTotal[0], oriTotal[1]) ;
  index = 0 ; kmerHashFind (kh, data + target*len, &index, &isRC) ;
  printf ("found target %d at %llu orientation %d\n", target, index, isRC) ;
  timeUpdate (stdout) ;
  kmerHashDestroy (kh) ;
  free (data) ; totalAllocated -= count*len ;
}

#endif

/*********** following ONElib.c, so representation is binary compatible ************/
/*********** NB this means bases are naturally ordered in a little-endian machine **/

static U8 pack[] = {    /* sends N (indeed any non-CGT) to A, except 0,1,2,3 are maintained */
   0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
} ;

static inline void Compress_DNA (int len, char *s, U64 *u)
{ int i ;
  U8 *s0, *s1, *s2, *s3 ;
  U8 *t = (U8*) u ;

  s0 = (U8 *) s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  len -= 3;
  for (i = 0; i < len; i += 4)
    *t++ = pack[s0[i]] | (pack[s1[i]] << 2) | (pack[s2[i]] << 4) | (pack[s3[i]] << 6) ;
  switch (i-len)
  { case 0: *t++ = pack[s0[i]] | (pack[s1[i]] << 2) | (pack[s2[i]] << 4) ; break;
    case 1: *t++ = pack[s0[i]] | (pack[s1[i]] << 2) ; break;
    case 2: *t++ = pack[s0[i]] ; break;
    default: break;
  }
}

static U8 packRC[] = {   /* sends N (indeed any non-CGT) to A, except 0,1,2,3 are maintained */
   3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
} ;

static inline void Compress_DNA_RC (int len, char *s, U64 *u)
{ int    i, j;
  U8 *s0, *s1, *s2, *s3;
  U8 *t = (U8*) u ;

  s0 = (U8 *) s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  for (i = len-4, j = 0; i >= 0; i -= 4)
    t[j++] = packRC[s3[i]] | (packRC[s2[i]] << 2) | (packRC[s1[i]] << 4) | (packRC[s0[i]] << 6) ;
  switch (-i)
  { case 1: t[j++] = packRC[s3[i]] | (packRC[s2[i]] << 2) | (packRC[s1[i]] << 4) ; break;
    case 2: t[j++] = packRC[s3[i]] | (packRC[s2[i]] << 2) ; break;
    case 3: t[j++] = packRC[s3[i]] ; break;
    default: break;
  }
}

static char Base[4] = { 'a', 'c', 'g', 't' };

static inline void Uncompress_DNA (int len, U64 *u, char *t)
{ int   i, tlen, byte;
  char *t0, *t1, *t2, *t3;
  U8 *s = (U8*) u ;

  t0 = t;
  t1 = t0+1;
  t2 = t1+1;
  t3 = t2+1;

  tlen = len-3;
  for (i = 0; i < tlen; i += 4)
    { byte = *s++;
      t0[i] = Base[byte & 0x3];
      t1[i] = Base[(byte >> 2) & 0x3];
      t2[i] = Base[(byte >> 4) & 0x3];
      t3[i] = Base[(byte >> 6) & 0x3];
    }

  switch (i-tlen)
  { case 0:
      byte = *s++;
      t0[i] = Base[byte & 0x3];
      t1[i] = Base[(byte >> 2) & 0x3];
      t2[i] = Base[(byte >> 4) & 0x3];
      break;
    case 1:
      byte = *s++;
      t0[i] = Base[byte & 0x3];
      t1[i] = Base[(byte >> 2) & 0x3];
      break;
    case 2:
      byte = *s++;
      t0[i] = Base[byte & 0x3];
      break;
    default:
      break;
  }
  t[len] = 0 ; // 0 terminate
}

/*********** end of file ***********/
