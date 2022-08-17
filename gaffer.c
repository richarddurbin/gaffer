/*  File: gaffer.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2022
 *  License: MIT license (see file LICENSE)
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 17 13:25 2022 (rd109)
 * Created: Thu Mar 24 01:02:39 2022 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "seqio.h"
#include "ONElib.h"

typedef struct Seqstruct {
  int len ;
  char* dna ;
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
  OneFile *vf ;                 // just to carry header info
  DICT  *seqName ;
  Array seq ;			// of Seq
  Array link ;                  // of Link
  Array walk ;                  // of Array of Walk
  bool isPerfect ;              // all links match exactly
  bool isOwnSeq ;               // does this object own its sequences (used in destroy)
  bool hasDNA ;                 // true if the dna fields of seq[i] are filled
} Gfa ;

// have every sequence forwards and reverse complemented adjacent in gf->seq
#define RC(si) ((si)%2 ? (si)-1 : (si)+1)

static char *gfaSchemaText ;

Gfa *gfaCreate (int nSeq, int nLink)
{
  Gfa *gf = new0 (1, Gfa) ;
  gf->seqName = dictCreate (nSeq) ;
  gf->seq = arrayCreate (nSeq, Seq) ;
  gf->link = arrayCreate (nLink, Link) ;
  return gf ;
}

void gfaDestroy (Gfa *gf)
{ int i ;
  dictDestroy (gf->seqName) ;
  if (gf->hasDNA) for (i = arrayMax(gf->seq) ; i-- ;) free (arrp(gf->seq,i,Seq)->dna) ;
  arrayDestroy (gf->seq) ;
  arrayDestroy (gf->link) ;
  if (gf->walk)
    { for (i = arrayMax(gf->walk) ; i-- ;) arrayDestroy(arrp(gf->walk,i,Walk)->as) ;
      arrayDestroy (gf->walk) ;
    }
  if (gf->vf) oneFileClose (gf->vf) ; // will just free data as we set isWrite == false
  free (gf) ;
}

Gfa *readOneFiles (char *stem)
{
  OneSchema *schema = oneSchemaCreateFromText (gfaSchemaText) ;
  int stemLen = strlen(stem) ;
  char *fileName = new (stemLen+6, char) ; strcpy (fileName, stem) ;

  // general strategy is to check the schema just for the lines that we will read

  // have to start by opening both the seg and links files to find numbers of each to make gfa
  
  strcpy (fileName+stemLen, ".1seg") ;
  OneFile *vfs = oneFileOpenRead (fileName, schema, "seg", 1) ;
  if (!vfs) die ("can't open %s to read", fileName) ;
  if (!oneFileCheckSchema (vfs, "P 3 seg\nO S 1 3 INT\n"))
    die ("schema mismatch %s", fileName) ;
  int nSeg = vfs->info[vfs->objectType]->given.count ;

  strcpy (fileName+stemLen, ".1lnk") ;
  OneFile *vfl = oneFileOpenRead (fileName, schema, "lnk", 1) ;
  if (!vfl) die ("can't open %s to read", fileName) ;
  if (!oneFileCheckSchema (vfl, "P 3 lnk\nO L 4 3 INT 4 CHAR 3 INT 4 CHAR\nD O 1 3 INT\n"))
    die ("schema mismatch %s", fileName) ;
  int nLink = vfl->info[vfl->objectType]->given.count ;

  Gfa *gf = gfaCreate (nSeg, nLink) ;
  
  while (oneReadLine (vfs))
    if (vfs->lineType == 'S')
      { Seq *s = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	s->len = oneLen(vfs) ;
	s = arrayp(gf->seq, arrayMax(gf->seq), Seq) ; // and make one for the reverse complement
	s->len = oneLen(vfs) ;
      }
    else if (vfs->lineType == 'I') // assumes one per 'S' line - should be true 
      { if (!dictAdd (gf->seqName, oneString(vfs), 0)) die ("duplicate S line %d", vfs->line) ;
      }
  // create shell OneFile to maintain header information
  gf->vf = oneFileOpenWriteFrom ("-", vfs, false, 1) ; gf->vf->isWrite = false ;
  printf ("read file %s with %d segs\n", fileName, (int)vfs->object) ;
  oneFileClose (vfs) ;

  Link *l ;
  while (oneReadLine (vfl))
    if (vfl->lineType == 'L')
      { l = arrayp(gf->link, arrayMax(gf->link), Link) ;
	l->s1 = 2 * oneInt(vfl,0) ;
	if (oneChar(vfl,1) == '<') ++l->s1 ;
	l->s2 = 2 * oneInt(vfl,2) ;
	if (oneChar(vfl,3) == '<') ++l->s2 ;
      }
    else if (vfl->lineType == 'O')
      l->overlap = oneInt(vfl,0) ;
  printf ("read file %s with %d links\n", fileName, (int)vfl->object) ;
  oneFileClose (vfl) ;

  strcpy (fileName+stemLen, ".1sgs") ;
  OneFile *vfq = oneFileOpenRead (fileName, schema, "seq", 1) ;
  if (vfq)
    { if (!oneFileCheckSchema (vfq, "P 3 seq\nO S 1 3 DNA\n"))
	die ("schema mismatch %s", fileName) ;
      while (oneReadLine (vfq))
	if (vfq->lineType == 'S')
	  { Seq *s = arrp(gf->seq, 2*(vfq->object-1), Seq) ;
	    assert (oneLen(vfq) == s->len) ;
	    s->dna = new(s->len, char) ; memcpy (s->dna, oneString(vfq), s->len) ;
	    s[1].dna = new(s->len, char) ; s[1].dna = sqioRevComp(s->dna, s->len) ;
	  }
      printf ("read file %s with %d sequences\n", fileName, (int)vfq->object) ;
      oneFileClose (vfq) ;
      gf->hasDNA = true ;
    }

  strcpy (fileName+stemLen, ".1wlk") ;
  OneFile *vfw = oneFileOpenRead (fileName, schema, "wlk", 1) ;
  if (vfw)
    { if (!oneFileCheckSchema (vfw,"P 3 wlk\nO W 2 3 INT 1 8 INT_LIST\n"
			       "D D 1 6 STRING\nD S 1 3 INT\nD E 1 3 INT\n"))
	die ("schema mismatch %s", fileName) ;
      Walk *w ;
      while (oneReadLine (vfw))
	switch (vfw->lineType)
	  {
	  case 'W':
	    w = arrayp(gf->walk, arrayMax(gf->walk), Walk) ;
	    w->len = oneInt(vfw,0) ;
	    int i, n = oneLen(vfw) ;
	    I64 *i64 = oneIntList(vfw) ;
	    w->as = arrayCreate (n,int) ;
	    for (i = 0 ; i < n ; ++i) arr(w->as,i,int) = 2 * i64[i] ;
	    arrayMax(w->as) = n ;
	    break ;
	  case 'D':
	    n = oneLen(vfw) ; // inherit declaration of i, n from above
	    char *c = oneString(vfw) ;
	    for (i = 0 ; i < n ; ++i) if (c[i] == '<') ++arr(w->as,i,int) ;
	    break ;
	  case 'S': w->start = oneInt(vfw,0) ; break ;
	  case 'E': w->end = oneInt(vfw,0) ; break ;
	  }
      oneFileClose (vfw) ;
    }

  free (fileName) ;
  oneSchemaDestroy (schema) ;
  return gf ;
}

void writeOneFiles (Gfa *gf, char *stem)
{
  int i ;
  int stemLen = strlen(stem) ;
  char *fileName = new (stemLen+5, char) ; strcpy (fileName, stem) ;
  OneSchema *schema = oneSchemaCreateFromText (gfaSchemaText) ;

  strcpy (fileName+stemLen, ".1seg") ;
  OneFile *vfseg = oneFileOpenWriteFrom (fileName, gf->vf, true, 1) ; // use existing header
  if (!vfseg) die ("can't open % to write", fileName) ;
  oneWriteHeader (vfseg) ;
  for (i = 0 ; i < arrayMax(gf->seq) ; i += 2) /* only write sequences in one direction */
    { Seq *s = arrp(gf->seq, i, Seq) ;
      oneInt(vfseg,0) = s->len ;
      oneWriteLine (vfseg, 'S', 0, 0) ;
      if (dictMax(gf->seqName) > i/2)
	{ char *name = dictName(gf->seqName, i/2) ;
	  oneWriteLine (vfseg, 'I', strlen(name), name) ;
	}	  
    }
  fprintf (stderr, "wrote %d objects to %s\n", (int)vfseg->object, fileName) ;
  oneFileClose (vfseg) ;

  if (gf->hasDNA)
    { strcpy (fileName+stemLen, ".1sgs") ;
      OneFile *vfseq = oneFileOpenWriteNew (fileName, schema, "sgs", true, 1) ;
      if (!vfseq) die ("can't open % to write", fileName) ;
      oneWriteHeader (vfseq) ;
      for (i = 0 ; i < arrayMax(gf->seq) ; i += 2) /* only write sequences in one direction */
	{ Seq *s = arrp(gf->seq, i, Seq) ;
	  oneWriteLine (vfseq, 'S', s->len, s->dna) ;
	}
      fprintf (stderr, "wrote %d objects to %s\n", (int)vfseq->object, fileName) ;
      oneFileClose (vfseq) ;
    }

  strcpy (fileName+stemLen, ".1lnk") ;
  OneFile *vfl = oneFileOpenWriteNew (fileName, schema, "lnk", true, 1) ;
  if (!vfl) die ("can't open %s", fileName) ;
  oneWriteHeader (vfl) ;
  for (i = 0 ; i < arrayMax(gf->link) ; ++i)
    { Link *l = arrp(gf->link, i, Link) ;
      oneInt(vfl,0) = l->s1 >> 1 ;
      oneChar(vfl,1) = (l->s1 & 1) ? '<' : '>' ;
      oneInt(vfl,2) = l->s2 >> 1 ;
      oneChar(vfl,3) = (l->s2 & 1) ? '<' : '>' ;
      oneWriteLine (vfl, 'L', 0, 0) ;
      if (l->overlap) { oneInt(vfl,0) = l->overlap ; oneWriteLine (vfl, 'O', 0, 0) ; }
    }
  fprintf (stderr, "wrote %d objects to %s\n", (int)vfl->object, fileName) ;
  oneFileClose (vfl) ;

  if (gf->walk)
    { strcpy (fileName+stemLen, ".1wlk") ;
      OneFile *vfw = oneFileOpenWriteNew (fileName, schema, "wlk", true, 1) ;
      if (!vfw) die ("can't open %s", fileName) ;
      oneWriteHeader (vfw) ;
      Array aI64 = arrayCreate (64, I64) ;
      Array aDir = arrayCreate (64, char) ;
      int j ;
      for (i = 0 ; i < arrayMax(gf->walk) ; ++i)
	{ Walk *w = arrp(gf->walk, i, Walk) ;
	  int j, n = arrayMax(w->as) ;
	  array(aI64,n-1,I64) = 0 ; // ensure space so can use arr inside loop
	  array(aDir,n-1,char) = 0 ; // ensure space so can use arr inside loop
	  if (i < 10) fprintf (stderr, "walk %d length %d -", i, n) ;
	  int *asj = arrp(w->as, 0, int) ;
	  for (j = 0 ; j < n ; ++j, ++asj)
	    { arr(aI64,j,I64) = *asj >> 1 ; // can use arr() not array() since space confirmed
	      arr(aDir,j,char) = (*asj & 0x1) ? '<' : '>' ;
	      if (i < 10) fprintf (stderr, " %d", *asj) ;
	    }
	  if (i < 10) fprintf (stderr, "\n") ;
	  arrayMax(aI64) = n ; arrayMax(aDir) = n ; // must set by hand because not set by arr
	  oneInt(vfw,0) = w->len ; oneWriteLine (vfw, 'W', n, arrp(aI64,0,I64)) ;
	  oneWriteLine (vfw, 'D', n, arrp(aDir,0,char)) ;
	  if (w->start) { oneInt(vfw,0) = w->start ; oneWriteLine (vfw, 'S', 0, 0) ; }
	  if (w->end) { oneInt(vfw,0) = w->end ; oneWriteLine (vfw, 'E', 0, 0) ; }
	}
      fprintf (stderr, "wrote %d objects to %s\n", (int)vfw->object, fileName) ;
      oneFileClose (vfw) ;
    }
  
  free (fileName) ;
  oneSchemaDestroy (schema) ;
}

Gfa *gfaParseSL (char *filename)
{
  FILE *f = fzopen (filename, "r") ;
  if (!f) die ("failed to open GFA file %s", filename) ;

  Gfa *gf = gfaCreate (10000, 10000) ;
  // make a stub OneFile, into which provenance and comments can be written
  OneSchema *schema = oneSchemaCreateFromText (gfaSchemaText) ;
  gf->vf = oneFileOpenWriteNew ("-", schema, "seg", false, 1) ; gf->vf->isWrite = false ;
  oneSchemaDestroy (schema) ;
  
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
	      seq->dna = new (seq->len, char) ;
	      memcpy (seq->dna, word, seq->len) ;
	      gf->hasDNA = true ; // assume every seq will have DNA...
	    }
	  word = fgetword (f) ; // length
	  if (strncmp (word, "LN:i:", 5)) die ("bad LN field line %d", line) ;
	  if (seq->dna)
	    { if (atoi(&word[5]) != seq->len) die ("length error line %d", line) ; }
	  else
	    seq->len = atoi (&word[5]) ;
	  // now make the reverse complement as the (n+1)th Seq
	  Seq *seq2 = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	  seq2->len = seq->len ;
	  if (seq->dna) seq2->dna = sqioRevComp (seq->dna, seq->len) ;
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
	  if (link->s1 == link->s2)
	    { warn ("** self link from/to %s %d overlap %d - ignoring",
		    dictName (gf->seqName, index), link->s1, link->overlap) ;
	      --arrayMax(gf->link) ;
	    }
	  else	  // add the reverse complement link
	    { Link *l2 = arrayp(gf->link, arrayMax(gf->link), Link) ;
	      l2->s1 = RC(link->s2) ; l2->s2 = RC(link->s1) ; l2->overlap = link->overlap ;
	    }
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
  return gf ;
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

void readSeqFile (Gfa *gf, char *filename) // to be used in conjunction with GFA files
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
	seq->dna = new (si->seqLen, char) ;
	memcpy (seq->dna, sqioSeq(si), si->seqLen) ;
	++seq ; seq->dna = sqioRevComp (sqioSeq(si), si->seqLen) ;
	total += si->seqLen ;
      }
    else
      die ("unknown sequence %s", sqioId(si)) ;
  printf ("read %lld sequences total %lld from %s file %s\n",
	  si->nSeq, total, seqIOtypeName[si->type], filename) ;
  seqIOclose (si) ;
  gf->hasDNA = true ;
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

  if (!gf->hasDNA)
    { fprintf (stderr, "no DNA sequence loaded yet - aborting removeBadLinks\n") ;
      return ;
    }
  
  int nPerfect = 0, nImperfect = 0 ;
  Link *ul = arrp (gf->link, 0, Link) ;
  for (i = 0 ; i < arrayMax (gf->link) ; ++i)
    { Link *l = arrp (gf->link, i, Link) ;
      Seq *seq1 = arrp(gf->seq,l->s1,Seq) ;
      char *s1 = arrp(gf->seq,l->s1,Seq)->dna + (seq1->len - l->overlap) ;
      char *s2 = arrp(gf->seq,l->s2,Seq)->dna ;
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
  int x ;	// position of the cut
  int s ;	// index of lowest seq in which the cut ending at x is found
  int sx ;      // position in s
  bool isLeft ; // true if s, sx apply to the left, false if to the right
} Cut ;

static int cutOrder (const void *a, const void *b) // by x, else isLeft, else s
{ Cut *ca = (Cut*)a, *cb = (Cut*)b ;
  if (ca->x < cb->x) return -1 ; else if (ca->x > cb->x) return 1 ;
  if (ca->isLeft && !cb->isLeft) return -1 ; else if (!ca->isLeft && cb->isLeft) return 1 ;
  return (ca->s - cb->s) ;
}

static int X = 1 ; /* debugging */

void newCut (Array *cut, int i, int x, int s, int sx, bool isLeft)
{
  Cut *c = arrayp(cut[i], arrayMax(cut[i]), Cut) ;
  c->x = x ; c->s = s ; c->s = s ; c->sx = sx, c->isLeft = isLeft ;
  if (i == X) // debug
    printf ("  new cut (%d,%d = %d,%d)%c\n", i, x, s, sx, isLeft?'L':'R') ;
}

Gfa *bluntify (Gfa *gf1) // returns a list of int arrays of cutpoints per seq
{
  int nSeq = arrayMax(gf1->seq) ; // for convenience
  Cut *c, *c1, *c2, *cEnd ; // utility variables

  Array *cut = new(nSeq, Array) ; // working space to store cutpoints in each seq
  int is ;
  for (is = 0 ; is < nSeq ; ++is) // initialise it
    { cut[is] = arrayCreate (8, Cut) ; // start with endpoints of each seq
      int len = arrp(gf1->seq, is, Seq)->len ;
      if (is == X) printf ("seq %d len %d\n", is, len) ;
      newCut (cut, is, 0, is, 0, false) ;
      newCut (cut, is, len,  is, len, true) ;
    }

  // next add the direct cutpoints from the links into cut[]
  arraySort (gf1->link, linkOrder1) ;
  Link *l, *lEnd = arrp(gf1->link, arrayMax(gf1->link), Link) ;
  for (l = arrp(gf1->link, 0, Link) ; l < lEnd ; ++l)
    { if (l->s1 == X || l->s2 == X) printf ("link %d to %d overlap %d\n", l->s1,l->s2,l->overlap) ;
      int s1len = arrp(gf1->seq, l->s1,Seq)->len ;
      int s1x = s1len - l->overlap ;
      newCut (cut, l->s1, s1x, l->s1, s1x, true) ; // left of cut in s1
      newCut (cut, l->s2, l->overlap, l->s2, l->overlap, false) ; // right of cut in s2
      if (l->s1 < l->s2)
	{ newCut (cut, l->s1, s1x, l->s1, s1x, false) ; // right of cut in s1
	  newCut (cut, l->s2, l->overlap, l->s1, s1len, true) ; // left of cut in s2
	  newCut (cut, l->s2, 0, l->s1, s1x, false) ; // replace the start of s2
	}
      else
	{ newCut (cut, l->s1, s1x, l->s2, 0, false) ; // right of cut in s1
	  newCut (cut, l->s2, l->overlap, l->s2, l->overlap, true) ; // left of cut in s2
	  newCut (cut, l->s1, s1len, l->s2, l->overlap, true) ; // replace the end of s1
	}
    }

  // now iterate until all cuts match across all links
  int nTouch = nSeq ;
  int *cutMax = new0(nSeq,int) ; // arrayMax(cut[i]) prior to updating
  while (nTouch)
    { bool *isTouched = new0(nSeq,bool) ;
      // first sort and compactify the cut lists - keep Cut at each x and direction with lowest ->s
      for (is = 0 ; is < nSeq ; ++is)
	{ arraySort(cut[is], cutOrder) ;
	  if (cutMax[is] < arrayMax(cut[is])) isTouched[is] = true ;
	  c1 = c2 = arrp(cut[is], 0, Cut) ;
	  cEnd = arrp(cut[is], arrayMax(cut[is]), Cut) ;
	  while (++c2 < cEnd)
	    if (c2->x != c1->x || c2->isLeft != c1->isLeft)
	      { ++c1 ; c1->x = c2->x ; c1->s = c2->s ; c1->sx = c2->sx ; c1->isLeft = c2->isLeft ; }
	  cutMax[is] = arrayMax(cut[is]) = ++c1 - arrp(cut[is], 0, Cut) ;
	}
      // next for each link l check whether cut sites on l->s1 and l->s2 match up: if not add them
      nTouch = 0 ;
      for (l = arrp(gf1->link, 0, Link) ; l < lEnd ; ++l)
	if (isTouched[l->s1] || isTouched[l->s2]) // no need to check if they didn't change
	  { Array cut1 = cut[l->s1], cut2 = cut[l->s2] ;
	    int off = arrp(gf1->seq,l->s1,Seq)->len - l->overlap ;
	    int i1 = 0, i2 = 0 ;
	    if (l->s1 == X || l->s2 == X)
	      printf ("comparing link %d to %d off %d\n", l->s1, l->s2, off) ;
	    while (arrp(cut1,i1,Cut)->x < off || arrp(cut1,i1,Cut)->isLeft) ++i1 ;
	    assert (arrp(cut1,i1,Cut)->x == off) ; // off right must be in cut1 - we added it above
	    while (i1 < cutMax[l->s1]) // only up to end of previous list
	      { c1 = arrp(cut1, i1, Cut) ; c2 = arrp(cut2, i2, Cut) ;
		if (l->s1 == X || l->s2 == X)
		  printf ("  comparing cuts c1 %d (%d,%d = %d,%d)%c to c2 %d (%d,%d = %d,%d)%c\n",
			  i1, l->s1, c1->x, c1->s, c1->sx, c1->isLeft?'L':'R',
			  i2, l->s2, c2->x, c2->s, c2->sx, c2->isLeft?'L':'R') ;
		if (c1->x - off > c2->x) // add c2 into cut1
		  { if (l->s1 < l->s2)
		      { newCut (cut, l->s1, c2->x + off, l->s1, c2->x + off, c2->isLeft) ;
			newCut (cut, l->s2, c2->x, l->s1, c2->x + off, c2->isLeft) ;
		      }
		    else
		      newCut (cut, l->s1, c2->x + off, l->s2, c2->x, c2->isLeft) ;
		    ++i2 ; ++nTouch ;
		  }
		else if (c1->x - off < c2->x) // add c1 into cut2
		  { if (l->s2 < l->s1)
		      { newCut (cut, l->s2, c1->x - off, l->s2, c1->x - off, c1->isLeft) ;
			newCut (cut, l->s1, c1->x, l->s2, c1->x - off, c1->isLeft) ;
		      }
		    else
		      newCut (cut, l->s2, c1->x - off, l->s1, c1->x, c1->isLeft) ;
		    ++i1 ; ++nTouch ;
		  }
		else
		  { assert (c1->isLeft == c2->isLeft) ;
		    if (c1->s < c2->s) // check if we need to update ->s for either c1 or c2
		      { c2 = arrayp(cut2, arrayMax(cut2), Cut) ; ++nTouch ;
			c2->x = c1->x - off ; c2->s = c1->s ; c2->sx = c1->sx ; c2->isLeft = c1->isLeft ;
			if (l->s1 == X || l->s2 == X)
			  printf ("  adding c2 (%d,%d = %d,%d)\n", l->s2,c2->x, c2->s,c2->sx) ;
		      }
		    else if (c1->s > c2->s)
		      { c1 = arrayp(cut1, arrayMax(cut1), Cut) ; ++nTouch ;
			c1->x = c2->x + off ; c1->s = c2->s ; c1->sx = c2->sx ; c1->isLeft = c2->isLeft ;
			if (l->s1 == X || l->s2 == X)
			  printf ("  adding c1 (%d,%d = %d,%d)\n", l->s1,c1->x, c1->s,c1->sx) ;
		      }
		    ++i1 ; ++i2 ; ++nTouch ;
		  }
	      }
	  }
      free (isTouched) ;
    }
  free (cutMax) ;

  // next create gf2: seqs from the segments between cuts, and links and walks from the seqs
  Gfa *gf2 = gfaCreate (arrayMax(gf1->seq)+arrayMax(gf1->link), 8*arrayMax(gf1->link)) ;
  gf2->walk = arrayCreate (arrayMax(gf1->seq), Walk) ;
  
  Hash newSeqHash = hashCreate (2*arrayMax(gf1->link)) ;
  for (is = 0 ; is < nSeq ; ++is)
    { Walk *wi = arrayp(gf2->walk, is, Walk) ;
      wi->len = arrp(gf1->seq,is,Seq)->len ;
      wi->start = 0 ;
      wi->as = arrayCreate (arrayMax(cut[is])/2, int) ;
      if (is == X || is == RC(X)) // debug section
	{ printf ("seq %d len %d when building gf2\n  cut", is, wi->len) ;
	  for (c = arrp(cut[is],0,Cut) ; c < arrp(cut[is],arrayMax(cut[is]),Cut) ; ++c)
	    printf (" (%d,%d,%d)%c", c->x, c->s, c->sx, c->isLeft?'L':'R') ;
	  putchar ('\n') ;
	}
      int i, index, lastIndex = -1 ;
      for (i = 0 ; i < arrayMax(cut[is]) ; i += 2)
	{ c1 = arrp(cut[is],i,Cut) ;
	  c2 = arrp(cut[is],i+1,Cut) ;
	  assert (!c1->isLeft && c2->isLeft) ;
	  assert (c1->s == c2->s) ; 
	  Seq *sfs = arrp(gf1->seq,c1->s,Seq) ;
	  if (!hashFind (newSeqHash, HASH_INT2(c1->s,c1->sx), &index))
	    { hashAdd (newSeqHash, HASH_INT2(c1->s,c1->sx), &index) ;
	      hashAdd (newSeqHash, HASH_INT2(RC(c1->s),sfs->len - c1->sx), 0) ; // add RC seq
	    }
	  Seq *sNew = arrayp(gf2->seq,index,Seq) ;
	  if (sNew->len)
	    assert (sNew->len == c2->x - c1->x) ;
	  else
	    { sNew->len = c2->x - c1->x ;
	      sNew->dna = &(sfs->dna[c1->sx]) ;
	    }
	  array(wi->as,arrayMax(wi->as),int) = index ;
	  if (lastIndex >= 0)
	    { Link *lNew = arrayp(gf2->link, arrayMax(gf2->link), Link) ;
	      lNew->s1 = lastIndex ; lNew->s2 = index ; lNew->overlap = 0 ;
	    }
	  lastIndex = index ;
	}
    }
  hashDestroy (newSeqHash) ;

  // finally sort and compress the links
  arraySort (gf2->link, linkOrder1) ;
  printf ("%d initial new links before compression\n", arrayMax(gf2->link)) ;
  arrayCompress (gf2->link) ;

  printf ("made blunt gfa with %d seqs, %d links and %d walks\n",
	  arrayMax(gf2->seq), arrayMax(gf2->link), arrayMax(gf2->walk)) ;

  // now check the walks
  for (is = 0 ; is < nSeq ; ++is)
    { char *s = arrp(gf1->seq, is, Seq)->dna ;
      if (!s) continue ;
      Walk *w = arrp(gf2->walk, is, Walk) ;
      int i, j, n = 0 ;
      for (i = 0 ; i < arrayMax(w->as) ; ++i)
	{ Seq *fs = arrp(gf2->seq, arr(w->as,i,int), Seq) ;
	  char *t = fs->dna ;
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
  fprintf (stderr, "Usage: gaffer <commands>\n") ;
  fprintf (stderr, "   -read  <stem>   : read ONE files\n") ;
  fprintf (stderr, "   -write <stem>   : write ONE files\n") ;
  fprintf (stderr, "   -schema         : print ONE schema\n") ;
  fprintf (stderr, "   -readGfa <gfa file> : only reads S and L lines for now\n") ;
  fprintf (stderr, "   -dna <sequence file matching S line names>\n") ;
  fprintf (stderr, "   -removeBadLinks : (for now) remove imperfect overlaps\n") ;
  fprintf (stderr, "   -blunt          : makes new non-overlapping graph\n") ;
  fprintf (stderr, "   -depth <gfa file with A lines>\n") ;
  fprintf (stderr, "   -extend         : adds match blocks for shared incoming edges - needs seq\n") ;
  exit (0) ;
}

int main (int argc, char *argv[])
{
  Gfa *gf = 0 ;
  
  timeUpdate (0) ;

  argc-- ; ++argv ; // swallow the program name
  if (!argc) usage() ;
  while (argc)
    { if (**argv != '-') die ("argument %s is expected to start with '-'", *argv) ;
      char *command = *argv+1 ;
      if (gf)
	{
	}
      if (!strcmp (*argv, "-readGfa") && argc >= 2)
	{ if (gf) { fprintf (stderr, "removing existing gf\n") ; gfaDestroy (gf) ; }
	  gf = gfaParseSL (argv[1]) ;
	  linkRemoveDuplicates (gf) ;
	  argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-dna") && argc >= 2)
	{ if (gf) readSeqFile (gf, argv[1]) ;
	  else fprintf (stderr, "can't read sequences without a graph\n") ;
	  argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-read") && argc >= 2)
	{ if (gf) { fprintf (stderr, "removing existing gf\n") ; gfaDestroy (gf) ; }
	  gf = readOneFiles (argv[1]) ;
	  argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-write") && argc >= 2)
	{ if (gf) writeOneFiles (gf, argv[1]) ;
	  else fprintf (stderr, "can't write without a graph\n") ;
	  argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-schema"))
	{ printf ("%s", gfaSchemaText) ;
	  argc-- ; argv++ ; 
	}
      else if (!strcmp (*argv, "-removeBadLinks"))
	{ if (gf) linkRemoveBad (gf) ;
	  else fprintf (stderr, "can't remove bad links without a graph\n") ;
	  argc-- ; argv++ ;
	}
      else if (!strcmp (*argv, "-blunt"))
	{ if (gf)
	    { Gfa *gf2 = bluntify (gf) ;
	      gf2->vf = gf->vf ; gf->vf = 0 ; 
	      gfaDestroy (gf) ;
	      gf = gf2 ;
	    }
	  else fprintf (stderr, "can't bluntify without a graph\n") ;
	  argc-- ; argv++ ;
	}
      else die ("unrecognised option %s - run without args for help", *argv) ;
      fprintf (stderr, "%s: ", command) ; timeUpdate (stderr) ;
    }

  if (gf) gfaDestroy (gf) ;
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
}

// here is the working schema derived from the GFA spec at

static char *gfaSchemaText =
  "1 3 def 1 0   schema for GFA\n"
  ".\n"
  ". segments don't need to have explicit sequences - so we make them their own class\n"
  ". if segment sequences are available use the sgs subtype of seq\n" 
  ".\n"
  "P 3 seg                 SEGMENT\n"
  "O S 1 3 INT             length\n"
  "D I 1 6 STRING          identifier - required by GFA\n"
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
  "D I 1 6 STRING          id - sequence identifier; unnecessary for segments\n"
  "D Q 1 6 STRING          base quality - Q values (ascii string = q+33)\n"
  "D P 1 8 INT_LIST        kmer ploidy estimates - 0 for error, 1,2... for ploidy, R for repeat\n"
  ".\n"
  "P 3 lnk                             LINK\n"
  "O L 4 3 INT 4 CHAR 3 INT 4 CHAR     s1 dir1 s2 dir2 - s1,2 are indices in seg file, dir=+|-\n"
  "D O 1 3 INT                         overlap - else presume abut\n"
  "D G 1 6 STRING                      cigar string - else presume exact (only meaningful if O)\n"
  "D Q 1 3 INT                         MQ mapping quality\n"
  "D M 1 3 INT                         NM number of mismatches\n"
  "D R 1 3 INT                         RC read count\n"
  "D F 1 3 INT                         FC fragment count\n"
  "D K 1 3 INT                         KC k-mer count\n"
  "D I 1 6 STRING                      ID edge identifier (deprecated)\n"
  ".\n"
  "P 3 wlk                 WALK - in GFA 1.1 only - requires non-overlapping segments\n"
  "O W 2 3 INT 8 INT_LIST  length, list of segments\n"
  "D D 1 6 STRING          list of directions of each segment (+ or -) - required\n"
  "D I 1 6 STRING          identifier (required by GFA)\n"
  "D J 1 6 STRING          sample identifier (required by GFA)\n"
  "D H 1 3 INT             haplotype : 1..<ploidy> - missing if 0 in GFA = haploid or unknown\n"
  "D S 1 3 INT             start in first segment - defaults to 0 if missing\n"
  "D E 1 3 INT             end in last segment - defaults to end of segment if missing\n"
  ".\n"
  "P 3 ctn                                  CONTAINMENT - contained and container are both segs\n"
  "O L 5 3 INT 4 CHAR 3 INT 4 CHAR 3 INT    s1 dir1 s2 dir2 pos - s1,2 in seg file dir=+|- start position\n"
  "D G 1 6 STRING                           cigar string - else presume exact match\n"
  "D I 1 6 STRING                           ID edge identifier (deprecated)\n"
  "D R 1 3 INT                              RC read count\n"
  "D M 1 3 INT                              NM number of mismatches\n"
  ".\n"
  "P 3 pth                PATH\n"
  "O P 1 8 INT_LIST       list of segments\n"
  "D D 1 6 STRING         list of directions of each segment (+ or -) - required\n"
  "D G 1 11 STRING_LIST   list of cigar strings for overlaps - optional/deprecated\n"
  "D I 1 6 STRING         path identifier (required by GFA)\n"
  ".\n"
  "P 3 aln                  ALIGNMENT - of sequences from a sequence file, e.g. reads\n"
  "O A 2 3 INT 8 INT_LIST   index of seq to align, list of segments as in WALK\n"
  "D D 1 6 STRING           list of directions of each segment (+ or -) - required\n"
  "D S 1 3 INT              start - as in WALK\n"
  "D E 1 3 INT              end - as in WALK\n"
  "D G 1 6 STRING           cigar string - for mapping\n"
  "D Q 1 3 INT              mapping quality\n"
  "D M 1 3 INT              number of mismatches\n"
  ".\n" ;
  
/******************* end of file *******************/
