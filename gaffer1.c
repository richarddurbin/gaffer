/*  File: gaffer.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2022
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 11 13:02 2022 (rd109)
 * Created: Thu Mar 24 01:02:39 2022 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "seqio.h"
#include "ONElib.h"

typedef struct Seqstruct {
  I64 len ;
  char* s ;
  int rd ; // read depth from input unitig file - do we want this?
} Seq ;

typedef struct {
  int s1, s2 ;
  bool d1, d2 ;
  int nBad ;
  int overlap ;
} Link ;

typedef struct {
} Walk ;

typedef struct {
  DICT  *seqName ;
  Array seq ;			// of Seq
  Array link ;                  // of Link
  Array walk ;                  // of Array of Walk
  bool isPerfect ;
} Gfa ;

static char comp[128] ;
static char *gfaSchemaText ;

void gfaParseSL (Gfa *gf, char *filename)
{
  FILE *f = fopen (filename, "r") ;

  if (!f) die ("failed to open GFA file %s", filename) ;
  int line = 1 ;
  while (!feof (f))
    { char *word = fgetword (f) ;
      if (feof(f)) break ;
      if (strlen(word) > 1) die ("line %d starts with %s not a single char", line, word) ;
      U64 index ;
      switch (*word)
	{
	case 'S':
	  word = fgetword (f) ;
	  if (!dictAdd (gf->seqName, word, &index)) die ("duplicate S line %d", line) ;
	  Seq *seq = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	  word = fgetword (f) ;
	  if (*word != '*') seq->s = strdup (word) ;
	  word = fgetword (f) ;
	  if (strncmp (word, "LN:i:", 5)) die ("bad LN field line %d", line) ;
	  seq->len = atoi (&word[5]) ;
	  if (seq->s && seq->len != strlen(seq->s)) die ("length error line %d", line) ;
	  break ;
	case 'L':
	  word = fgetword (f) ;
	  if (!dictFind (gf->seqName, word, &index)) die ("bad s1 %d line %d", line) ;
	  Link *link = arrayp(gf->link, arrayMax(gf->link), Link) ;
	  link->s1 = index ;
	  word = fgetword (f) ;
	  if (*word == '+') link->d1 = true ;
	  else if (*word != '-') die ("bad d1 line %d", line) ;
	  word = fgetword (f) ;
	  if (!dictFind (gf->seqName, word, &index)) die ("bad s2 %d line %d", line) ;
	  link->s2 = index ;
	  word = fgetword (f) ;
	  if (*word == '+') link->d2 = true ;
	  else if (*word != '-') die ("bad d2 line %d", line) ;
	  word = fgetword (f) ; if (!*word) die ("empty match line %d", line) ;
	  char *s = word ; while (*s) ++s ; --s ;
	  if (*s != 'M') die ("not an expected match line %d", line) ;
	  *s = 0 ;
	  link->overlap = atoi (word) ;
	  break ;
	case 'A':
	  break ;
	default:
	  die ("bad gfa line type %c (ASCII %d) line %d", *word, *word, line) ;
	}
      char c ;
      while ((c = fgetc(f)) && !feof(f) && c != '\n') ; // flush end of line
      ++line ;
    }
  fclose (f) ;
  printf ("read %d S lines and %d L lines from GFA file %s\n",
	  arrayMax(gf->seq), arrayMax(gf->link), filename) ;
}

void readSeqFile (Gfa *gf, char *filename)
{
  SeqIO *si = seqIOopenRead (filename, 0, false) ;
  if (!si) die ("failed to open sequence file %s to read", filename) ;
  U64 index ;
  U64 total = 0 ;
  while (seqIOread (si))
    if (dictFind (gf->seqName, sqioId(si), &index))
      { Seq *seq = arrp (gf->seq, index, Seq) ;
	if (si->seqLen != seq->len)
	  die ("length mismatch for %s seq %lld != gfa %lld",
	       sqioId(si), si->seqLen, seq->len) ;
	seq->s = new (si->seqLen, char) ;
	memcpy (seq->s, sqioSeq(si), si->seqLen) ;
	// seq->s = strdup (sqioSeq(si)) ;
	total += si->seqLen ;
      }
    else
      die ("unknown sequence %s", sqioId(si)) ;
  printf ("read %lld sequences total %lld from %s file %s\n",
	  si->nSeq, total, seqIOtypeName[si->type], filename) ;
  seqIOclose (si) ;
}

static void linkReport (Array a, int i, int n)
{ while (i < n)
    { Link *l = arrp(a, i, Link) ;
      printf ("link %d: l1 %c%d %c%d %d\n", i,
	      l->d1?'+':'-', l->s1, l->d2?'+':'-', l->s2, l->overlap) ;
      ++i ;
    }
}

int checkOrder1 (const void *a, const void *b)
{ Link *la = (Link*)a, *lb = (Link*)b ;
  if (la->s1 < lb->s1) return -1 ; else if (la->s1 > lb->s1) return 1 ;
  if (la->s2 < lb->s2) return -1 ; else if (la->s2 > lb->s2) return 1 ;
  return 0 ;
}

int checkOrder2 (const void *a, const void *b)
{ Link *la = (Link*)a, *lb = (Link*)b ;
  if (la->s2 < lb->s2) return -1 ; else if (la->s2 > lb->s2) return 1 ;
  if (la->s1 < lb->s1) return -1 ; else if (la->s1 > lb->s1) return 1 ;
  return 0 ;
}

void fixOverlaps (Gfa *gf)
{
  int i, j ;

  // first check all the links to see if they are correct
  int nPerfect = 0, nImperfect = 0 ;
  for (i = 0 ; i < arrayMax (gf->link) ; ++i)
    { Link *link = arrp (gf->link, i, Link) ;
      Seq *s1 = arrp(gf->seq,link->s1,Seq) ;
      Seq *s2 = arrp(gf->seq,link->s2,Seq) ;
      char *c1 = link->d1 ? &(s1->s[s1->len - link->overlap]) : &(s1->s[link->overlap]) ;
      char *c2 = link->d2 ? s2->s : &(s2->s[s2->len]) ;
      link->nBad = 0 ;
      if (link->d1 && link->d2)
	for (j = 0 ; j < link->overlap ; ++j) { if (*c1++ != *c2++) ++link->nBad ; }
      else if (link->d1 && !link->d2)
	for (j = 0 ; j < link->overlap ; ++j) { if (*c1++ != comp[*--c2]) ++link->nBad ; }
      else if (!link->d1 && link->d2)
	for (j = 0 ; j < link->overlap ; ++j) { if (comp[*--c1] != *c2++) ++link->nBad ; }
      else if (!link->d1 && !link->d2)
	for (j = 0 ; j < link->overlap ; ++j)
	  { if (*--c1 != *--c2) ++link->nBad ; } // could comp[] them both, but unnecessary
      if (link->nBad) ++nImperfect ; else ++nPerfect ;
    }

  // next check perfect pairs are reciprocal
  // in some cases (surprisingly) one is perfect and the other is not: fix these
  int nFixed = 0 ;
  arraySort (gf->link, checkOrder1) ;
  { Array a2 = arrayCopy (gf->link) ;
    arraySort (a2, checkOrder2) ;
    for (i = 0 ; i < arrayMax(a2) ; ++i)
      { Link *l1 = arrp(gf->link, i, Link), *l2 = arrp(a2,i,Link) ;
	if (l1->s1 != l2->s2 || l1->s2 != l2->s1)
	  die ("link mismatch at %d: l1 %c%d %c%d %d, l2 %c%d %c%d %d", i,
	       l1->d1?'+':'-', l1->s1, l1->d2?'+':'-', l1->s2, l1->overlap,
	       l2->d1?'+':'-', l2->s1, l2->d2?'+':'-', l2->s2, l2->overlap) ;
	if (l1->nBad && !l2->nBad) // strange case where it is perfect one way but not the other
	  { l1->nBad = 0 ; l1->overlap = l2->overlap ; }
	if (!l1->nBad && l2->nBad) // the reciprocal strange case - requires more work
	  { int j ;
	    for (j = 0 ; j < arrayMax(gf->link) ; ++j) // linear search is bad, but I hope rare
	      { Link *lj = arrp(gf->link, j, Link) ;
		if (lj->s1 != l2->s1 || lj->s2 != l2->s2) continue ;
		lj->nBad = 0 ; lj->overlap = l1->overlap ;
		--nImperfect ; ++nPerfect ; ++nFixed ;
		break ;
	      }
	  }
      }
    arrayDestroy (a2) ;
  }
  
  // finally only keep the perfect links, for now
  Array l2 = arrayCreate (arrayMax(gf->link), Link) ;
  for (i = 0 ; i < arrayMax(gf->link) ; ++i)
    if (!arrp(gf->link,i,Link)->nBad) array(l2,arrayMax(l2),Link) = arr(gf->link,i,Link) ;
  arrayDestroy (gf->link) ; gf->link = l2 ;
  
  printf ("%d imperfect overlaps removed, %d fixed, %d remain\n", nImperfect, nFixed, nPerfect) ;
}

static int intSort (const void *a, const void *b) { return *(int*)a - *(int*)b ; }

Gfa *bluntify (Gfa *gf1)
{
  Gfa *gf2 = new0 (1, Gfa) ;
  gf2->seq = arrayCreate (arrayMax(gf1->seq)+arrayMax(gf1->link), Seq) ;
  gf2->link = arrayCreate (2*arrayMax(gf1->link), Link) ;
  gf2->walk = arrayCreate (arrayMax(gf1->seq), Array) ;

  // the new segs in gf2 will be blunt-ended subsegments of gf1 segs - call them frags here
  // the old segs in gf1 will become paths in gf2

  // strategy is to represent each frag by its occurrence in the least s containing it
  // use overlapping frags from s for s < s1, 
  // then partition the remaining sequence from this seg using overlaps for s > s1
  arraySort (gf1->link, checkOrder1) ; // must be sorted on s1, s2
  Array aCut = arrayCreate (256, int) ; // just for sorting cutpoints
  Link *l = arrp(gf1->link, 0, Link) ;
  Link *lEnd = arrp(gf1->link, arrayMax(gf1->link), Link) ;
  while (l < lEnd) 
    { int s1 = l->s1 ;
      int s1len = arrp(gf1->seq, s1, Seq)->len ;
      int maxLeft = 0, minRight = s1len ; // s1 runs from left to right
      Link *lLeft = 0, *lRight = 0 ; // real ones can't be 0
      while (l->s1 == s1 && l->s2 < s1)
	{ if (l->d1) // forward overlap, at end of s
	    { int cut = s1len - l->overlap ;
	      if (cut < minRight)
		{ minRight = (cut < maxLeft) ? maxLeft : cut ;
		  lRight = l ;
		}
	    }
	  else // reverse overlap, at start of s
	    { int cut = l->overlap ;
	      if (cut > maxLeft)
		{ maxLeft = (cut > minRight) ? minRight : cut ;
		  lLeft = l ;
		}
	    }
	  ++l ;
	}
      
      // now sort the cutpoints for this seg between maxLeft and minRight ;
      arrayMax(aCut) = 0 ; // clear aCut
      array(aCut,0,int) = maxLeft ; array(aCut,1,int) = minRight ;
      while (l < lEnd && l->s1 == s1)
	{ int cut = l->d1 ? s1len - l->overlap : l->overlap ;
	  if (cut > maxLeft && cut < minRight) array(aCut,arrayMax(aCut),int) = cut ;
	  ++l ;
	}
      if (arrayMax(aCut) > 2) arraySort (aCut, intSort) ; // sort if anything added
      
      // now assemble the walk in gfa2 corresponding to s1
      Array *aw1 = arrayp(gf2->walk, s1, Array) ;
      if (lLeft) // add frags from lLeft->s2
	{ Array *aw2 = arrp(gf2->walk, l->s2, Array) ;
	  
	}
    }
  
  return gf2 ;
}

void usage (void)
{
  fprintf (stderr, "gaffer <commands>\n") ;
  fprintf (stderr, "       -sl <gfa file with S and L lines>\n") ;
  fprintf (stderr, "       -seq <sequence file matching S line names>\n") ;
  fprintf (stderr, "       -fixOverlaps : (for now) remove imperfect overlaps\n") ;
  fprintf (stderr, "       -blunt  : makes new non-blunt graph\n") ;
  fprintf (stderr, "       -depth <gfa file with A lines>\n") ;
  fprintf (stderr, "       -extend : adds match blocks for shared incoming edges - needs seq\n") ;
  exit (0) ;
}

int main (int argc, char *argv[])
{
  Gfa *gf = new0 (1, Gfa) ;
  gf->seqName = dictCreate (10000) ;
  gf->seq = arrayCreate (10000,Seq) ;
  gf->link = arrayCreate (10000,Link) ;

  comp['a'] = 't' ; comp['c'] = 'g' ; comp['g'] = 'c' ; comp['t'] = 'a' ;
  
  timeUpdate (0) ;

  argc-- ; ++argv ; // swallow the program name
  if (!argc) usage() ;
  while (argc)
    { if (**argv != '-') die ("argument %s is expected to start with '-'", *argv) ;
      if (!strcmp (*argv, "-sl") && argc >= 2)
	{ gfaParseSL (gf, argv[1]) ;
	  fprintf (stderr, "%s: ", *argv+1) ; argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-seq") && argc >= 2)
	{ readSeqFile (gf, argv[1]) ;
	  fprintf (stderr, "%s: ", *argv+1) ; argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-fixOverlaps"))
	{ fixOverlaps (gf) ;
	  fprintf (stderr, "%s: ", *argv+1) ; argc-- ; argv++ ;
	}
      else if (!strcmp (*argv, "-blunt"))
	{ Gfa *gf2 = bluntify (gf) ;
	  fprintf (stderr, "%s: ", *argv+1) ; argc-- ; argv++ ;
	}
      else die ("unrecognised option %s - run without args for help", argv[1]) ;
      timeUpdate (stderr) ;
    }

  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
}


static char *gfaSchemaText =
  "1 3 def 1 0   schema for GFA\n"
  ".\n"
  ". segments don't need to have explicit sequences - so we make them their own class\n"
  ". if segment sequences are available use the sgs subtype of seq\n" 
  ".\n"
  "P 3 seg                 SEGMENT\n"
  "O S 2 6 STRING 3 INT    id, length - required fields\n"
  "D R 1 3 INT             RC read count\n"
  "D F 1 3 INT             FC fragment count\n"
  "D K 1 3 INT             KC k-mer count\n"
  "D H 1 6 STRING          SH checksum in hex\n"
  "D U 1 6 STRING          UR URL for sequence\n"
  ".\n"
  "P 3 seq                 SEQUENCE\n"
  "S 3 sgs                 segment sequences - objects are 1:1 with those in seg file\n"
  "S 3 srd                 read sequences\n"
  "O S 1 3 DNA             sequence: the DNA string\n"
  "D I 1 6 STRING          id - sequence identifier; only use for reads, not segments\n"
  "D Q 1 6 STRING          base quality - Q values (ascii string = q+33)\n"
  "D P 1 6 STRING          kmer ploidy estimates - 0 for error, 1,2... for ploidy, R for repeat\n"
  ".\n"
  "P 3 lnk                             LINK\n"
  "O L 4 3 INT 4 CHAR 3 INT 4 CHAR     s1 dir1 s2 dir2 - s1,2 are indices in seg file, dir=+|-\n"
  "D G 1 6 STRING                      cigar string - else presume abut\n"
  "D Q 1 3 INT                         MQ mapping quality\n"
  "D M 1 3 INT                         NM number of mismatches\n"
  "D R 1 3 INT                         RC read count\n"
  "D F 1 3 INT                         FC fragment count\n"
  "D K 1 3 INT                         KC k-mer count\n"
  "D I 1 6 STRING                      ID edge identifier (deprecated)\n"
  ".\n"
  "P 3 ctn                                  CONTAINMENT - contained and container are both segs\n"
  "O L 5 3 INT 4 CHAR 3 INT 4 CHAR 3 INT    s1 dir1 s2 dir2 pos - s1,2 in seg file dir=+|- start position\n"
  "D G 1 6 STRING                           cigar string - else presume abut\n"
  "D I 1 6 STRING                           ID edge identifier (deprecated)\n"
  "D R 1 3 INT                              RC read count\n"
  "D M 1 3 INT                              NM number of mismatches\n"
  ".\n"
  "P 3 pth                PATH\n"
  "O P 1 8 INT_LIST       list of segments\n"
  "D D 1 9 CHAR_LIST      list of directions of each segment (+ or -) - required\n"
  "D G 1 11 STRING_LIST   list of cigar strings for overlaps - optional/deprecated\n"
  "D I 1 6 STRING         path identifier (required by GFA)\n"
  ".\n"
  "P 3 wlk                WALK - in GFA 1.1 only - requires non-overlapping segments\n"
  "O P 1 8 INT_LIST       list of segments\n"
  "D D 1 9 CHAR_LIST      list of directions of each segment (+ or -) - required\n"
  "D I 1 6 STRING         identifier (required by GFA)\n"
  "D J 1 6 STRING         sample identifier (required by GFA)\n"
  "D H 1 3 INT            haplotype : 1..<ploidy> - missing if 0 in GFA = haploid or unknown\n"
  "D S 1 3 INT            start - defaults to 0 = start of first segment if missing0\n"
  "D E 1 3 INT            end (open, i.e. 1 + end coord so end-start = length) - defaults to end of last segment if missing\n"
  ".\n"
  "P 3 aln                  ALIGNMENT - of sequences from a sequence file, e.g. reads\n"
  "O A 2 3 INT 8 INT_LIST   index of seq to align, list of segments as in WALK\n"
  "D D 1 9 CHAR_LIST        list of directions of each segment (+ or -) - required\n"
  "D S 1 3 INT              start - as in WALK\n"
  "D E 1 3 INT              end - as in WALK\n"
  "D G 1 6 STRING           cigar string - for mapping\n"
  "D Q 1 3 INT              mapping quality\n"
  "D M 1 3 INT              number of mismatches\n"
  ".\n" ;
  
/******************* end of file *******************/
