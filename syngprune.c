/*  File: syngprune.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 29 20:15 2023 (rd109)
 * Created: Mon May 29 08:22:38 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "syng.h"
#include "ONElib.h"

typedef struct {
  int seg ;
  int count ;
  int pos ;
  char ori ;
} Unit ;

static inline char flip (char c) { return (c == '+') ? '-' : '+' ; }

void view (int object, Array units, int min, int max)
{
  int i, j, n = arrayMax(units) ;
  
  for (i = 0 ; i < n ; ++i)
    { Unit *u = arrp(units,i,Unit) ;
      if (u->count >= min && u->count <= max)
	{ if (u->ori == '+')
	    { printf ("%6d %4d ", u->seg, u->count) ;
	      if (i) printf ("L %6d %4d %c", u[-1].seg, u[-1].count, u[-1].ori) ;
	      else printf ("L      .    . .") ;
	      if (i < n-1) printf (" R %6d %4d %c |", u[1].seg, u[1].count, u[1].ori) ;
	      else printf (" R      .    . . |") ;
	      printf (" OBJ %6d |", object) ;
	      for (j = i+1 ; j < n ; ++j) printf (" %d", arrp(units,j,Unit)->count) ;
	      printf ("\n") ;
	    }
	  else
	    { printf ("%6d %4d ", u->seg, u->count) ;
	      if (i < n-1) printf ("L %6d %4d %c", u[1].seg, u[1].count, flip(u[1].ori)) ;
	      else printf ("L      .    . .") ;
	      if (i) printf (" R %6d %4d %c |", u[-1].seg, u[-1].count, flip(u[-1].ori)) ;
	      else printf (" R      .    . . |") ;
	      printf (" OBJ %6d |", object) ;
	      for (j = i-1 ; j >= 0 ; --j) printf (" %d", arrp(units,j,Unit)->count) ;
	      printf ("\n") ;
	    }
	}
    }
}

OneFile *readsIn, *readsOut ;

typedef struct {
  OneFile *readsyn ;
  OneFile *readsIn ;
  OneFile *readsOut ; // share the kCounts Array as a read-only global
} Thread_info ;

static size_t fieldSize[128] ;

static inline void transferLine (OneFile *vfIn, OneFile *vfOut)
{ memcpy (vfOut->field, vfIn->field, fieldSize[(int)vfIn->lineType]) ;
  oneWriteLine (vfOut, vfIn->lineType, oneLen(vfIn), oneString(vfIn)) ;
  char *s = oneReadComment (vfIn) ; if (s) oneWriteComment (vfOut, "%s", s) ;
}

bool filter (Array units, int min, int max)
{
  static U8 *dna = 0 ;
  static I64 dnaLen = 0, dnaMax = 0 ;
  int i, n = arrayMax(units) ;
  bool isPrint = true ;
  Unit *u = arrp(units,0,Unit) ;
  for (i = 0 ; i < n ; ++i,++u) if (u->count >= min && u->count <= max) isPrint = false ; 
  if (isPrint)
    { if (dna) oneWriteLineDNA2bit (readsOut, 'S', dnaLen, dna) ;
      while (oneReadLine(readsIn) && readsIn->lineType != 'S')
	transferLine (readsIn, readsOut) ;
    }
  else
    while (oneReadLine(readsIn) && readsIn->lineType != 'S') ;
  if (readsIn->lineType == 'S')
    { dnaLen = oneLen(readsIn) ;
      if (dnaLen > dnaMax) { if (dna) free(dna) ; dnaMax = 2*dnaLen ; dna = new(dnaMax,U8) ; }
      memcpy (dna, oneDNA2bit(readsIn), dnaLen) ;
    }
  return !isPrint ;
}

static char usage[] =
  "Usage: syngprune <options> <onecode file prefix> <read sequence file>\n"
  "  -min <min>    : [0]\n"
  "  -max <max>    : [1]\n"
  "  -v            : view\n"
  "  -T <nthreads> : [1]\n" ;

int main (int argc, char **argv)
{
  OneFile *seg, *readsyn ;
  int min = 0, max = 1 ;
  int i, object ;
  bool isView = false ;
  int nFiltered = 0 ;
  int nThreads = 1 ;

  timeUpdate (0) ;

  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;
  if (!argc) { printf ("%s",usage) ; exit (0) ; }
  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-min") && argc > 1) { min = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-max") && argc > 1) { max = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-T") && argc > 1) { nThreads = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-v")) { isView = true ; argc -= 1 ; argv += 1 ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  if (!argc || argc > 2) die ("bad number of arguments %d\n%s", argc, usage) ;

  if (isView && nThreads > 1) die ("can only use 1 thread in view mode\n") ;
  if (nThreads < 1 || nThreads > 16) die ("nThreads %d must be between 1 and 16", nThreads) ;
  
  if (!(seg = oneFileOpenRead (fnameTag (argv[0],"1seg"), 0, "seg", 1)))
    die ("failed to open %s.1seg", argv[0]) ;
  Array kCounts = arrayCreate(seg->info['K']->given.count, int) ;
  while (oneReadLine (seg))
    if (seg->lineType == 'K') array(kCounts, arrayMax(kCounts), int) = oneInt(seg,0) ;
  oneFileClose (seg) ;
  fprintf (stderr, "read counts for %d segments\n", arrayMax(kCounts)) ;
  
  if (!isView)
    { if (argc != 2) die ("need two arguments when pruning, not %d\n%s", argc, usage) ;
      if (!(readsIn = oneFileOpenRead (argv[1], 0, "seq", nThreads)))
	die ("failed to open reads file %s", argv[1]) ;
      for (i = 0 ; i < 128 ; ++i)
	if (readsIn->info[i]) fieldSize[i] = readsIn->info[i]->nField*sizeof(OneField) ;
      if (!(readsOut = oneFileOpenWriteFrom ("-", readsIn, true, nThreads)))
	die ("failed to open stdout to write to", argv[1]) ;
      oneAddProvenance (readsOut, "syngprune", "0.0", getCommandLine()) ;
    }

  Array units = arrayCreate (64,Unit) ;
  if (!(readsyn = oneFileOpenRead (fnameTag (argv[0],"1readsyn"), 0, "seqsyn", nThreads)))
    die ("failed to open %s.1readsyn to read", argv[0]) ;
  while (oneReadLine (readsyn))
    if (readsyn->lineType == 'S')
      { I64 *x = oneIntList(readsyn) ;
	if (isView) view (object, units, min, max) ; else if (filter (units, min, max)) ++nFiltered ;
	arrayMax(units) = 0 ;
	object = readsyn->object ;
	for (i = 0 ; i < oneLen(readsyn) ; ++i) arrayp(units,i,Unit)->seg = x[i] ;
 	for (i = 0 ; i < oneLen(readsyn) ; ++i) arrayp(units,i,Unit)->count = arr(kCounts,(int)x[i],int) ;
      }
    else if (readsyn->lineType == 'O')
      { char*s = oneString(readsyn) ;
	for (i = 0 ; i < oneLen(readsyn) ; ++i) arrayp(units,i,Unit)->ori = s[i] ;
      }
    else if (readsyn->lineType == 'P')
      { I64 *p = oneIntList(readsyn) ;
	for (i = 0 ; i < oneLen(readsyn) ; ++i) arrayp(units,i,Unit)->pos = p[i] ;
      }
  if (isView)
    view (object, units, min, max) ;
  else
    { if (filter (units, min, max)) ++nFiltered ;
      fprintf (stderr, "filtered %d of %d objects\n", nFiltered, (int)readsyn->object+1) ;
      oneFileClose (readsIn) ;
      oneFileClose (readsOut) ;
    }
  oneFileClose (readsyn) ;
  
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}
