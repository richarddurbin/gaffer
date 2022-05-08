/*  File: gaffer.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2022
 *  License: MIT license (see file LICENSE)
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May  8 08:15 2022 (rd109)
 * Created: Thu Mar 24 01:02:39 2022 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "seqio.h"
#include "ONElib.h"

typedef struct Seqstruct {
  int len ;
  char* s ;
} Seq ;

typedef struct {
  int s1, s2 ;
  int overlap ;
} Link ;

typedef struct {
  Array as ;                    // of int = Seq index
  int start, end ;              // start in as[0], end in as[max-1]
  int len ;
} Walk ;

typedef struct {
  DICT  *seqName ;
  Array seq ;			// of Seq
  Array link ;                  // of Link
  Array walk ;                  // of Array of Walk
  bool isPerfect ;              // all links match exactly
  bool isOwnSeq ;               // does this object own its sequences (used in destroy)
} Gfa ;

// have every sequence forwards and reverse complemented adjacent in gf->seq
#define RC(si) ((si)%2 ? (si)-1 : (si)+1)

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
      int index ;
      switch (*word)
	{
	case 'S':
	  word = fgetword (f) ; // name
	  if (!dictAdd (gf->seqName, word, &index)) die ("duplicate S line %d", line) ;
	  Seq *seq = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	  word = fgetword (f) ; // sequence
	  if (*word != '*')
	    { seq->len = strlen (word) ;
	      seq->s = new (seq->len, char) ;
	      memcpy (seq->s, word, seq->len) ;
	    }
	  word = fgetword (f) ; // length
	  if (strncmp (word, "LN:i:", 5)) die ("bad LN field line %d", line) ;
	  if (seq->s)
	    { if (atoi(&word[5]) != seq->len) die ("length error line %d", line) ; }
	  else
	    seq->len = atoi (&word[5]) ;
	  // now make the reverse complement as the (n+1)th Seq
	  Seq *seq2 = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	  seq2->len = seq->len ;
	  if (seq->s) seq2->s = sqioRevComp (seq->s, seq->len) ;
	  break ;
	case 'L':
	  word = fgetword (f) ; // seq1
	  if (!dictFind (gf->seqName, word, &index)) die ("bad s1 %d line %d", line) ;
	  Link *link = arrayp(gf->link, arrayMax(gf->link), Link) ;
	  link->s1 = 2*index ;
	  word = fgetword (f) ; // dir1
	  if (*word == '-') ++link->s1 ;
	  else if (*word != '+') die ("bad d1 line %d", line) ;
	  word = fgetword (f) ; // seq2
	  if (!dictFind (gf->seqName, word, &index)) die ("bad s2 %d line %d", line) ;
	  link->s2 = 2*index ;
	  word = fgetword (f) ; // dir2
	  if (*word == '-') ++link->s2 ;
	  else if (*word != '+') die ("bad d2 line %d", line) ;
	  word = fgetword (f) ; // cigar - this code only parses exact match
	  if (!*word) die ("empty match line %d", line) ;
	  char *s = word ; while (*s) ++s ; --s ;
	  if (*s != 'M') die ("not an expected match line %d", line) ;
	  *s = 0 ;
	  link->overlap = atoi (word) ;
	  // now add the reverse complement link
	  Link *l2 = arrayp(gf->link, arrayMax(gf->link), Link) ;
	  l2->s1 = RC(link->s2) ; l2->s2 = RC(link->s1) ; l2->overlap = link->overlap ;
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

int linkOrder1 (const void *a, const void *b) // by s1, then s2
{ Link *la = (Link*)a, *lb = (Link*)b ;
  if (la->s1 < lb->s1) return -1 ; else if (la->s1 > lb->s1) return 1 ;
  if (la->s2 < lb->s2) return -1 ; else if (la->s2 > lb->s2) return 1 ;
  return 0 ;
}

void linkRemoveDuplicates (Gfa *gf)
{
  arraySort (gf->link, linkOrder1) ;
  int i ;
  Link *lu = arrp (gf->link, 0, Link) ;
  for (i = 1 ; i < arrayMax(gf->link) ; ++i)
    { Link *l = arrp (gf->link, i, Link) ;
      if (lu->s1 != l->s1 || lu->s2 != l->s2 || lu->overlap != l->overlap) *++lu = *l ;
    }
  int newMax = lu - arrp(gf->link, 0, Link) + 1 ;
  printf ("removed %d duplicate links, leaving %d\n", arrayMax(gf->link) - newMax, newMax) ;
  arrayMax (gf->link) = newMax ;
}

void readSeqFile (Gfa *gf, char *filename)
{
  SeqIO *si = seqIOopenRead (filename, 0, false) ;
  if (!si) die ("failed to open sequence file %s to read", filename) ;
  int index ;
  U64 total = 0 ;
  while (seqIOread (si))
    if (dictFind (gf->seqName, sqioId(si), &index))
      { index *= 2 ;
	Seq *seq = arrp (gf->seq, index, Seq) ;
	if (si->seqLen != seq->len)
	  die ("length mismatch for %s seq %lld != gfa %lld", sqioId(si), si->seqLen, seq->len) ;
	seq->s = new (si->seqLen, char) ;
	memcpy (seq->s, sqioSeq(si), si->seqLen) ;
	++seq ; seq->s = sqioRevComp (sqioSeq(si), si->seqLen) ;
	total += si->seqLen ;
      }
    else
      die ("unknown sequence %s", sqioId(si)) ;
  printf ("read %lld sequences total %lld from %s file %s\n",
	  si->nSeq, total, seqIOtypeName[si->type], filename) ;
  seqIOclose (si) ;
}

static void linkReport (Array a, int i, int n) // previously used for debugging
{ while (i < n)
    { Link *l = arrp(a, i, Link) ;
      printf ("link %d: l1 %d %d %d\n", i, l->s1, l->s2, l->overlap) ;
      ++i ;
    }
}

void linkRemoveBad (Gfa *gf)
{
  int i, j ;

  int nPerfect = 0, nImperfect = 0 ;
  Link *ul = arrp (gf->link, 0, Link) ;
  for (i = 0 ; i < arrayMax (gf->link) ; ++i)
    { Link *l = arrp (gf->link, i, Link) ;
      Seq *seq1 = arrp(gf->seq,l->s1,Seq) ;
      char *s1 = &(seq1->s[seq1->len - l->overlap]) ;
      char *s2 = arrp(gf->seq,l->s2,Seq)->s ;
      for (j = 0 ; j < l->overlap ; ++j) { if (*s1++ != *s2++) break ; }
      if (j == l->overlap) *ul++ = *l ; // pefect match
    }

  int newMax = ul - arrp(gf->link, 0, Link) ;
  printf ("%d imperfect overlaps removed, %d remain\n", arrayMax(gf->link) - newMax, newMax) ;
  arrayMax(gf->link) = newMax ;

  gf->isPerfect = true ;
}

static int intOrder (const void *a, const void *b) { return *(int*)a - *(int*)b ; }

int linkOrder2 (const void *a, const void *b) // by s2 then s1
{ Link *la = (Link*)a, *lb = (Link*)b ;
  if (la->s2 < lb->s2) return -1 ; else if (la->s2 > lb->s2) return 1 ;
  if (la->s1 < lb->s1) return -1 ; else if (la->s1 > lb->s1) return 1 ;
  return 0 ;
}

static int linkOrderSize (const void *a, const void *b) // by overlap size small to large
{ Link *la = (Link*)a, *lb = (Link*)b ;
  if (la->overlap < lb->overlap) return -1 ; else if (la->overlap > lb->overlap) return 1 ;
  return 0 ;
}

typedef struct {
  int x ;
  int s ;	// index of lowest seq in which the cut ending at x is found
  int sx ;      // position in s
} Cut ;

static int cutOrder (const void *a, const void *b) // by x, else s
{ Cut *ca = (Cut*)a, *cb = (Cut*)b ;
  if (ca->x < cb->x) return -1 ; else if (ca->x > cb->x) return 1 ;
  return (ca->s - cb->s) ;
}

Gfa *bluntify (Gfa *gf1) // returns a list of int arrays of cutpoints per seq
{
  int nSeq = arrayMax(gf1->seq) ; // for convenience
  Cut *c, *c1, *c2, *cEnd ; // utility variables
  int X = 218, Y = 4 ; // for debugging

  Array *cut = new(nSeq, Array) ; // working space to store cutpoints in each seq
  int is ;
  for (is = 0 ; is < nSeq ; ++is) // initialise it
    { cut[is] = arrayCreate (8, Cut) ; // start with endpoint of each seq
      c = arrayp(cut[is], 0, Cut) ; c->x = c->sx = 0 ; c->s = is ;
      c = arrayp(cut[is], 1, Cut) ; c->x = c->sx = arrp(gf1->seq, is, Seq)->len ; c->s = is ;
      if (is == X || is == Y) printf ("seq %d len %d\n", is, arrp(gf1->seq, is, Seq)->len) ;
    }

  // next add the direct cutpoints from the links into cut[]
  arraySort (gf1->link, linkOrder1) ;
  Link *l, *lEnd = arrp(gf1->link, arrayMax(gf1->link), Link) ;
  for (l = arrp(gf1->link, 0, Link) ; l < lEnd ; ++l)
    { c1 = arrayp(cut[l->s1], arrayMax(cut[l->s1]), Cut) ; // at start of overlap in s1
      c2 = arrayp(cut[l->s2], arrayMax(cut[l->s2]), Cut) ; // at end of overlap in s2
      int c1len = arrp(gf1->seq,l->s1,Seq)->len ;
      c1->x = c1len - l->overlap ;
      c2->x = l->overlap ;
      if (l->s1 < l->s2) { c1->s = c2->s = l->s1 ; c1->sx = c1->x ; c2->sx = c1len ; }
      else { c1->s = l->s1 ; c2->s = l->s2 ; c1->sx = c1->x ; c2->sx = c2->x ; }
      if (l->s1 == X || l->s2 == X)
	{ printf ("link %d to %d overlap %d\n", l->s1,l->s2,l->overlap) ;
	  printf ("  c1 = (%d,%d = %d,%d)\n", l->s1, c1->x, c1->s, c1->sx) ;
	  printf ("  c2 = (%d,%d = %d,%d)\n", l->s2, c2->x, c2->s, c2->sx) ;
	}
    }

  // now iterate until all cuts match across all links
  int nTouch = nSeq ;
  int *cutMax = new0(nSeq,int) ; // arrayMax(cut[i]) prior to updating
  while (nTouch)
    { bool *isTouched = new0(nSeq,bool) ;
      // first sort and compactify the cut lists - keep Cut at each x with lowest ->s
      for (is = 0 ; is < nSeq ; ++is)
	{ arraySort(cut[is], cutOrder) ;
	  if (cutMax[is] < arrayMax(cut[is])) isTouched[is] = true ;
	  c1 = c2 = arrp(cut[is], 0, Cut) ;
	  cEnd = arrp(cut[is], arrayMax(cut[is]), Cut) ;
	  while (++c2 < cEnd)
	    if (c2->x != c1->x) { ++c1 ; c1->x = c2->x ; c1->s = c2->s ; c1->sx = c2->sx ; }
	  cutMax[is] = arrayMax(cut[is]) = ++c1 - arrp(cut[is], 0, Cut) ;
	}
      // next check whether cut sites on l->s1 and l->s2 match up: if not add them
      nTouch = 0 ;
      for (l = arrp(gf1->link, 0, Link) ; l < lEnd ; ++l)
	if (isTouched[l->s1] || isTouched[l->s2]) // no need to check if they didn't change
	  { Array cut1 = cut[l->s1], cut2 = cut[l->s2] ;
	    int off = arrp(gf1->seq,l->s1,Seq)->len - l->overlap ;
	    int i1 = 0, i2 = 0 ;
	    while (arrp(cut1,i1,Cut)->x < off) ++i1 ;
	    assert (arrp(cut1,i1,Cut)->x == off) ; // off must be in cut1 - we added it above
	    while (i1 < cutMax[l->s1]) // only up to end of previous list
	      { c1 = arrp(cut1, i1, Cut) ; c2 = arrp(cut2, i2, Cut) ;
		if (l->s1 == X || l->s2 == X)
		  printf ("comparing c1 (%d,%d = %d,%d) to c2 (%d,%d = %d,%d)\n",
			  l->s1,c1->x, c1->s,c1->sx, l->s2,c2->x, c2->s,c2->sx) ;
		if (c1->x - off > c2->x) // add c2 into cut1
		  { c1 = arrayp(cut1, arrayMax(cut1), Cut) ; ++nTouch ;
		    c1->x = c2->x + off ; c1->s = c2->s ; c1->sx = c2->sx ;
		    if (l->s1 < c2->s) // need to update ->s, ->sx for new c1 and c2
		      { c1->s = l->s1 ; c1->sx = c1->x ;
			c2 = arrayp(cut2, arrayMax(cut2), Cut) ; ++nTouch ;
			c2->x = c1->x - off ; c2->s = c1->s ; c2->sx = c1->sx ;
			if (l->s1 == X || l->s2 == X)
			  printf ("  adding c1 (%d,%d = %d,%d) and c2 (%d,%d = %d,%d)\n",
				  l->s1,c1->x, c1->s,c1->sx, l->s2,c2->x, c2->s,c2->sx) ;
		      }
		    else if (l->s1 == X || l->s2 == X)
		      printf ("  adding c1 (%d,%d = %d,%d)\n", l->s1,c1->x, c1->s,c1->sx) ;
		    ++i2 ;
		  }
		else if (c1->x - off < c2->x) // add c1 into cut2
		  { c2 = arrayp(cut2, arrayMax(cut2), Cut) ; ++nTouch ;
		    c2->x = c1->x - off ; c2->s = c1->s ; c2->sx = c1->sx ;
		    if (l->s2 < c1->s) // need to update ->s, ->sx for new c2 and c1
		      { c2->s = l->s2 ; c2->sx = c2->x ;
			c1 = arrayp(cut1, arrayMax(cut1), Cut) ; ++nTouch ;
			c1->x = c2->x + off ; c1->s = c2->s ; c1->sx = c2->sx ;
			if (l->s1 == X || l->s2 == X)
			  printf ("  adding c1 (%d,%d = %d,%d) and c2 (%d,%d = %d,%d)\n",
				  l->s1,c1->x, c1->s,c1->sx, l->s2,c2->x, c2->s,c2->sx) ;
		      }
		    else if (l->s1 == X || l->s2 == X)
		      printf ("  adding c2 (%d,%d = %d,%d)\n", l->s2,c2->x, c2->s,c2->sx) ;
		    ++i1 ;
		  }
		else
		  { if (c1->s < c2->s) // check if we need to update ->s for either c1 or c2
		      { c2 = arrayp(cut2, arrayMax(cut2), Cut) ; ++nTouch ;
			c2->x = c1->x - off ; c2->s = c1->s ; c2->sx = c1->sx ;
			if (l->s1 == X || l->s2 == X)
			  printf ("  adding c2 (%d,%d = %d,%d)\n", l->s2,c2->x, c2->s,c2->sx) ;
		      }
		    else if (c1->s > c2->s && c2->sx > 0)
		      { c1 = arrayp(cut1, arrayMax(cut1), Cut) ; ++nTouch ;
			c1->x = c2->x + off ; c1->s = c2->s ; c1->sx = c2->sx ;
			if (l->s1 == X || l->s2 == X)
			  printf ("  adding c1 (%d,%d = %d,%d)\n", l->s1,c1->x, c1->s,c1->sx) ;
		      }
		    ++i1 ; ++i2 ;
		  }
	      }
	  }
      free (isTouched) ;
    }
  free (cutMax) ;

  // next create gfa2: seqs from the segments between cuts, and links and walks from the seqs
  Gfa *gf2 = new0 (1, Gfa) ;
  gf2->seq = arrayCreate (arrayMax(gf1->seq)+arrayMax(gf1->link), Seq) ;
  gf2->link = arrayCreate (8*arrayMax(gf1->link), Link) ;
  gf2->walk = arrayCreate (arrayMax(gf1->seq), Walk) ;
  
  Hash newSeqHash = hashCreate (2*arrayMax(gf1->link)) ;
  for (is = 0 ; is < nSeq ; ++is)
    { int i ;
      int s1 = -1, s2 ;
      int index, lastIndex = -1 ;
      Walk *wi = arrayp(gf2->walk, is, Walk) ;
      wi->len = arrp(gf1->seq,is,Seq)->len ;
      wi->start = 0 ;
      wi->as = arrayCreate (8, int) ;
      c1 = arrp(cut[is],0,Cut) ;
      if (c1->s == nSeq) c1->s = is ;
      if (is == X || is == RC(X)) // debug section
	{ printf ("seq %d len %d when building gf2\n  cut", is, wi->len) ;
	  for (c = arrp(cut[is],0,Cut) ; c < arrp(cut[is],arrayMax(cut[is]),Cut) ; ++c)
	    printf (" (%d,%d,%d)", c->x, c->s, c->sx) ;
	  putchar ('\n') ;
	}
      for (i = 1 ; i < arrayMax(cut[is]) ; ++i)
	{ c2 = arrp(cut[is],i,Cut) ;
	  Seq *sfs = arrp(gf1->seq,c2->s,Seq) ;
	  if (!hashFind (newSeqHash, HASH_INT2(c2->s,c2->sx), &index))
	    { hashAdd (newSeqHash, HASH_INT2(c2->s,c2->sx), &index) ;
	      hashAdd (newSeqHash, HASH_INT2(RC(c2->s),sfs->len - c2->sx), 0) ; // add RC seq
	    }
	  Seq *sNew = arrayp(gf2->seq,index,Seq) ;
	  if (!sNew->len)
	    { sNew->len = c2->x - c1->x ;
	      if (sfs->s) sNew->s = &(sfs->s[c2->sx - sNew->len]) ;
	    }
	  array(wi->as,arrayMax(wi->as),int) = index ;
	  if (lastIndex >= 0)
	    { Link *lNew = arrayp(gf2->link, arrayMax(gf2->link), Link) ;
	      lNew->s1 = lastIndex ; lNew->s2 = index ; lNew->overlap = 0 ;
	    }
	  lastIndex = index ;
	  c1 = c2 ;
	}
    }
  hashDestroy (newSeqHash) ;

  // finally sort and compress the links
  arraySort (gf2->link, linkOrder1) ;
  arrayCompress (gf2->link) ;

  printf ("made blunt gfa with %d seqs, %d links and %d walks\n",
	  arrayMax(gf2->seq), arrayMax(gf2->link), arrayMax(gf2->walk)) ;

  // now check the walks
  for (is = 0 ; is < nSeq ; ++is)
    { char *s = arrp(gf1->seq, is, Seq)->s ;
      if (!s) continue ;
      Walk *w = arrp(gf2->walk, is, Walk) ;
      int i, j, n = 0 ;
      for (i = 0 ; i < arrayMax(w->as) ; ++i)
	{ Seq *fs = arrp(gf2->seq, arr(w->as,i,int), Seq) ;
	  char *t = fs->s ;
	  //	  printf ("seq %d pos %d starting cut %d len %d\n", is, n, i, fs->len) ;
	  for (j = 0 ; j < fs->len ; ++j, ++n)
	    if (*s++ != *t++)
	      die ("walk mismatch: seq %d len %d pos %d frag %d base %d",
		   is, arrp(gf1->seq, is, Seq)->len, n, i, j) ;
	}
    }
  printf ("all walks check out\n") ;

  return gf2 ;
}

void usage (void)
{
  fprintf (stderr, "gaffer <commands>\n") ;
  fprintf (stderr, "       -sl <gfa file with S and L lines>\n") ;
  fprintf (stderr, "       -seq <sequence file matching S line names>\n") ;
  fprintf (stderr, "       -removeBadLinks : (for now) remove imperfect overlaps\n") ;
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
  
  timeUpdate (0) ;

  argc-- ; ++argv ; // swallow the program name
  if (!argc) usage() ;
  while (argc)
    { if (**argv != '-') die ("argument %s is expected to start with '-'", *argv) ;
      if (!strcmp (*argv, "-sl") && argc >= 2)
	{ gfaParseSL (gf, argv[1]) ;
	  linkRemoveDuplicates (gf) ;
	  fprintf (stderr, "%s: ", *argv+1) ; argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-seq") && argc >= 2)
	{ readSeqFile (gf, argv[1]) ;
	  fprintf (stderr, "%s: ", *argv+1) ; argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-removeBadLinks"))
	{ linkRemoveBad (gf) ;
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
