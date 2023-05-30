/*  File: test.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2022
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 21 00:12 2022 (rd109)
 * Created: Mon Dec 19 20:29:29 2022 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "seqio.h"

int test1 (int argc, char *argv[])
{
  int i, j ;
  char *a = "accgttaaatcgagg" ;
  char *b = "cgttaatcgacgt" ;
  U8 ua[8], ub[8] ;
  SeqPack *sp = seqPackCreate ('a') ;
  seqPack (sp, a, ua, strlen(a)) ;
  seqPack (sp, b, ub, strlen(b)) ;
  { char *s1 = seqUnpack (sp, ub, 0, 0, 6), *s2 = seqUnpack (sp, ub, 0, 1, 6) ;
    printf ("b at 0 %s at 1 %s\n", s1, s2) ;
    free (s1) ; free (s2) ;
  }
  for (i = 0 ; i < 5 ; ++i)
    for (j = 0 ; j < 5 ; ++j)
      { int d = seqMatchPacked (ua, i, ub, j, 6) ;
	printf ("%s %d %s %d - %d\n", a, i, b, j, d) ;
      }
  return 0 ;
}

int test2 (int argc, char *argv[])
{
  char x[16] ;
  U8 y[8] ;
  U16 *z16 = (U16*)y ;
  U32 *z32 = (U32*)y ;
  U64 *z64 = (U64*)y ;
  int i, j ;

  for (i = 0 ; i < 16 ; ++i) x[i] = i ;
  for (j = 0 ; j < 8 ; ++j) y[j] = (x[2*j+1] << 4) | x[2*j] ;
  for (i = 0 ; i < 16 ; ++i) printf (" %1x", x[i]) ; putchar ('\n') ;
  for (i = 0 ; i < 8 ; ++i) printf (" %02x", y[i]) ; putchar ('\n') ;
  for (i = 0 ; i < 4 ; ++i) printf (" %4x", z16[i]) ; putchar ('\n') ;
  for (i = 0 ; i < 2 ; ++i) printf (" %8x", z32[i]) ; putchar ('\n') ;
  for (i = 0 ; i < 1 ; ++i) printf (" %16llx", z64[i]) ; putchar ('\n') ;
  return 0 ;
}

int test3 (int argc, char *argv[])
{
  char *fname = "-" ;
  --argc ; ++argv ;
  if (argc) fname = *argv ;
  SeqIO *si = seqIOopenRead (fname, dna2textConv, false) ;
  SeqPack *sp = seqPackCreate ('a') ;
  while (seqIOread (si))
    { int len = si->seqLen ;
      printf ("S  %.*s\n", len, sqioSeq(si)) ;
      char *s = seqRevComp (sqioSeq(si), len) ;
      printf ("R  %.*s\n", len, s) ;
      U8 *u = seqPack (sp, sqioSeq(si), 0, len) ;
      printf ("SP %.*s\n", len, seqUnpack (sp, u, 0, 0, len)) ;
      U8 *v = seqRevCompPacked (u, 0, len) ;
      printf ("RP %.*s\n\n", len, seqUnpack (sp, v, 0, 0, len)) ;
    }
  return 0 ;
}

int test4 (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (argc < 2) die ("need 2 files to compare") ;
  SeqIO *si1 = seqIOopenRead (argv[0], dna2textConv, false) ;
  SeqIO *si2 = seqIOopenRead (argv[1], dna2textConv, false) ;
  SeqPack *sp = seqPackCreate ('a') ;
  int x = 0, y ;
  while (seqIOread (si1) && seqIOread (si2))
    { if ((y = seqMatchPacked (seqPack(sp, sqioSeq(si1), 0, si1->seqLen), 0,
			       seqPack(sp, sqioSeq(si2), 0, si2->seqLen), 0, si1->seqLen)) != 0)
	printf ("mismatch at %d in seq %d : %.20s %.20s\n", y, x, sqioSeq(si1), sqioSeq(si2)) ;
      ++x ;
    }
  return 0 ;
}

int main (int argc, char *argv[]) { return test4 (argc, argv) ; }
