/*  File: gaffer.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2022
 *  License: MIT license (see file LICENSE)
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 15 15:57 2023 (rd109)
 * Created: Thu Mar 24 01:02:39 2022 (rd109)
 *-------------------------------------------------------------------
 */

#define VERSION "1.1"

#include "utils.h"
#include "seqio.h"
#include "ONElib.h"

typedef struct Seqstruct {
  int len ;
  U8* dna ;
} Seq ;

typedef struct {
  int s1, s2 ;
  int overlap ;
} Link ;

typedef struct {
  Array as ;                    // of int = Seq index
  int start, end ;              // start in as[0], end in as[max-1]
} Path ;

typedef struct {
  OneFile *vf ;                 // just to carry header info
  Array seq ;			// of Seq
  Array link ;                  // of Link
  Array path ;                  // of Array of Path
  DICT *seqName ;
  DICT *pathName ;
  bool  isPerfect ;             // all links match exactly
  char *buf ;			// for temporary operations,large enough for longest DNA sequence
  int   bufSize ;		// size of buf
} Gfa ;

// have every sequence forwards and reverse complemented adjacent in gf->seq
#define RC(si) ((si)%2 ? (si)-1 : (si)+1)

static SeqPack* SP ;		// create in main(), unpacking to acgt
static char *gfaSchemaText ;	// at end of file

Gfa *gfaCreate (int nSeq, int nLink)
{
  Gfa *gf = new0 (1, Gfa) ;
  gf->seq = arrayCreate (nSeq, Seq) ;
  gf->link = arrayCreate (nLink, Link) ;
  return gf ;
}

void gfaDestroy (Gfa *gf)
{ int i ;
  for (i = arrayMax(gf->seq) ; i-- ;)
    if (arrp(gf->seq,i,Seq)->dna) free (arrp(gf->seq,i,Seq)->dna) ;
  arrayDestroy (gf->seq) ;
  arrayDestroy (gf->link) ;
  if (gf->seqName) dictDestroy (gf->seqName) ;
  if (gf->pathName) dictDestroy (gf->pathName) ;
  if (gf->path)
    { for (i = arrayMax(gf->path) ; i-- ;) arrayDestroy(arrp(gf->path,i,Path)->as) ;
      arrayDestroy (gf->path) ;
    }
  if (gf->vf) oneFileClose (gf->vf) ; // will just free data as we set isWrite == false
  if (gf->buf) free (gf->buf) ;
  free (gf) ;
}

/*************************************************************************/

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
  int nSeg = vfs->info['S']->given.count ;

  strcpy (fileName+stemLen, ".1lnk") ;
  OneFile *vfl = oneFileOpenRead (fileName, schema, "lnk", 1) ;
  if (!vfl) die ("can't open %s to read", fileName) ;
  if (!oneFileCheckSchema (vfl, "P 3 lnk\nO L 4 3 INT 4 CHAR 3 INT 4 CHAR\nD O 1 3 INT\n"))
    die ("schema mismatch %s", fileName) ;
  int nLink = vfl->info['I']->given.count ;

  Gfa *gf = gfaCreate (nSeg, nLink) ;

  if (vfs->info['I'] && vfs->info['I']->given.count)
    { if (vfs->info['I']->given.count != nSeg)
	die ("name number %d does not match seq number %d", vfs->info['I']->given.count, nSeg) ;
      gf->seqName = dictCreate (nSeg) ;
    }
  
  while (oneReadLine (vfs))
    if (vfs->lineType == 'S')
      { Seq *s = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	s->len = oneLen(vfs) ;
	s = arrayp(gf->seq, arrayMax(gf->seq), Seq) ; // and make one for the reverse complement
	s->len = oneLen(vfs) ;
	if (s->len > gf->bufSize)
	  { gf->bufSize = s->len ;
	    if (gf->buf) free (gf->buf) ;
	    gf->buf = new(gf->bufSize,char) ;
	  }
      }
    else if (vfs->lineType == 'I') // assumes one per 'S' line - should be true 
      if (!dictAdd (gf->seqName, oneString(vfs), 0))
	die ("duplicate seg name line %d", vfs->line) ;
  // create shell OneFile to maintain header information
  gf->vf = oneFileOpenWriteFrom ("-", vfs, false, 1) ; gf->vf->isWrite = false ;
  printf ("read file %s with %d segs\n", fileName, (int)vfs->object) ;
  oneFileClose (vfs) ;

  Link *l ;
  while (oneReadLine (vfl))
    if (vfl->lineType == 'L')
      { l = arrayp(gf->link, arrayMax(gf->link), Link) ;
	l->s1 = 2 * (oneInt(vfl,0) - 1) ;
	if (oneChar(vfl,1) == '<') ++l->s1 ;
	l->s2 = 2 * (oneInt(vfl,2) - 1) ;
	if (oneChar(vfl,3) == '<') ++l->s2 ;
      }
    else if (vfl->lineType == 'O')
      l->overlap = oneInt(vfl,0) ;
  printf ("read file %s with %d links\n", fileName, (int)vfl->object) ;
  oneFileClose (vfl) ;

  strcpy (fileName+stemLen, ".1sgs") ;
  OneFile *vfd = oneFileOpenRead (fileName, schema, "seq", 1) ;
  if (vfd)
    { if (!oneFileCheckSchema (vfd, "P 3 seq\nO S 1 3 DNA\n"))
	die ("schema mismatch %s", fileName) ;
      while (oneReadLine (vfd))
	if (vfd->lineType == 'S')
	  { Seq *s = arrp(gf->seq, 2*(vfd->object-1), Seq) ;
	    assert (oneLen(vfd) == s->len) ;
	    s->dna = new((s->len+3)/4, U8) ; memcpy (s->dna, oneDNA2bit(vfd), (s->len+3)/4) ;
	    s[1].dna = seqRevCompPacked (s->dna, 0, s->len) ;
	  }
      printf ("read file %s with %d sequences\n", fileName, (int)vfd->object) ;
      oneFileClose (vfd) ;
    }

  strcpy (fileName+stemLen, ".1pth") ;
  OneFile *vfp = oneFileOpenRead (fileName, schema, "pth", 1) ;
  if (vfp)
    { if (!oneFileCheckSchema (vfp,"P 3 wlk\nO P 1 8 INT_LIST\nD D 1 6 STRING\n"))
	die ("schema mismatch %s", fileName) ;
      gf->path = arrayCreate (vfp->info['P']->given.count, Path) ;
      if (vfp->info['I'] && vfp->info['I']->given.count)
	gf->pathName = dictCreate (vfp->info['I']->given.count) ;
      Path *p = 0 ;
      while (oneReadLine (vfp))
	switch (vfp->lineType)
	  {
	  case 'P':
	    p = arrayp(gf->path, arrayMax(gf->path), Path) ;
	    int i, n = oneLen(vfp) ;
	    I64 *i64 = oneIntList(vfp) ;
	    p->as = arrayCreate (n,int) ;
	    for (i = 0 ; i < n ; ++i) arr(p->as,i,int) = 2 * (i64[i] - 1) ;
	    arrayMax(p->as) = n ;
	    break ;
	  case 'D':
	    n = oneLen(vfp) ; // inherit declaration of i, n from above
	    char *c = oneString(vfp) ;
	    for (i = 0 ; i < n ; ++i) if (c[i] == '<') ++arr(p->as,i,int) ;
	    break ;
	  case 'I':
	    if (!dictAdd (gf->pathName, oneString(vfs), 0))
	      die ("duplicate S line %d", vfs->line) ;
	    break ;
	  case 'S': p->start = oneInt(vfp,0) ; break ;
	  case 'E': p->end = oneInt(vfp,0) ; break ;
	  }
      printf ("read file %s with %d paths\n", fileName, (int)vfp->object) ;
      oneFileClose (vfp) ;
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
  for (i = 0 ; i < arrayMax(gf->seq) ; i += 2) /* only write sequences in one direction */
    { Seq *s = arrp(gf->seq, i, Seq) ;
      oneInt(vfseg,0) = s->len ;
      oneWriteLine (vfseg, 'S', 0, 0) ;
      if (gf->seqName)
	{ char *name = dictName(gf->seqName, i/2) ;
	  oneWriteLine (vfseg, 'I', strlen(name), name) ;
	}	  
    }
  fprintf (stderr, "wrote %d objects to %s\n", (int)vfseg->object, fileName) ;
  oneFileClose (vfseg) ;

  if (arrayMax(gf->seq) && arrp(gf->seq, 0, Seq)->dna)
    { strcpy (fileName+stemLen, ".1sgs") ;
      OneFile *vfseq = oneFileOpenWriteNew (fileName, schema, "sgs", true, 1) ;
      if (!vfseq) die ("can't open % to write", fileName) ;
      for (i = 0 ; i < arrayMax(gf->seq) ; i += 2) /* only write sequences in one direction */
	{ Seq *s = arrp(gf->seq, i, Seq) ;
	  oneWriteLineDNA2bit (vfseq, 'S', s->len, s->dna) ;
	}
      fprintf (stderr, "wrote %d objects to %s\n", (int)vfseq->object, fileName) ;
      oneFileClose (vfseq) ;
    }

  strcpy (fileName+stemLen, ".1lnk") ;
  OneFile *vfl = oneFileOpenWriteNew (fileName, schema, "lnk", true, 1) ;
  if (!vfl) die ("can't open %s", fileName) ;
  for (i = 0 ; i < arrayMax(gf->link) ; ++i)
    { Link *l = arrp(gf->link, i, Link) ;
      oneInt(vfl,0) = 1 + (l->s1 >> 1) ;
      oneChar(vfl,1) = (l->s1 & 1) ? '<' : '>' ;
      oneInt(vfl,2) = 1 + (l->s2 >> 1) ;
      oneChar(vfl,3) = (l->s2 & 1) ? '<' : '>' ;
      oneWriteLine (vfl, 'L', 0, 0) ;
      if (l->overlap) { oneInt(vfl,0) = l->overlap ; oneWriteLine (vfl, 'O', 0, 0) ; }
    }
  fprintf (stderr, "wrote %d objects to %s\n", (int)vfl->object, fileName) ;
  oneFileClose (vfl) ;

  if (gf->path)
    { strcpy (fileName+stemLen, ".1pth") ;
      OneFile *vfp = oneFileOpenWriteNew (fileName, schema, "pth", true, 1) ;
      if (!vfp) die ("can't open %s", fileName) ;
      Array aI64 = arrayCreate (64, I64) ;
      Array aDir = arrayCreate (64, char) ;
      int j ;
      for (i = 0 ; i < arrayMax(gf->path) ; ++i)
	{ Path *p = arrp(gf->path, i, Path) ;
	  int j, n = arrayMax(p->as) ;
	  array(aI64,n-1,I64) = 0 ; // ensure space so can use arr inside loop
	  array(aDir,n-1,char) = 0 ; // ensure space so can use arr inside loop
	  int *asj = arrp(p->as, 0, int) ;
	  for (j = 0 ; j < n ; ++j, ++asj)
	    { arr(aI64,j,I64) = 1 + (*asj >> 1) ; // can use arr() since space confirmed
	      arr(aDir,j,char) = (*asj & 0x1) ? '<' : '>' ;
	    }
	  arrayMax(aI64) = n ; arrayMax(aDir) = n ; // must set by hand because not set by arr
	  oneWriteLine (vfp, 'P', n, arrp(aI64,0,I64)) ;
	  oneWriteLine (vfp, 'D', n, arrp(aDir,0,char)) ;
	  if (gf->pathName)
	    { char *name = dictName (gf->pathName, i) ;
	      oneWriteLine (vfp, 'I', strlen(name), name) ;
	    }
	  if (p->start) { oneInt(vfp,0) = p->start ; oneWriteLine (vfp, 'S', 0, 0) ; }
	  if (p->end) { oneInt(vfp,0) = p->end ; oneWriteLine (vfp, 'E', 0, 0) ; }
	}
      fprintf (stderr, "wrote %d objects to %s\n", (int)vfp->object, fileName) ;
      oneFileClose (vfp) ;
    }
  
  free (fileName) ;
  oneSchemaDestroy (schema) ;
}

/*************************************************************************/

Gfa *gfaParseSL (char *filename)
{
  FILE *f = fzopen (filename, "r") ;
  if (!f) die ("failed to open GFA file %s", filename) ;

  Gfa *gf = gfaCreate (10000, 10000) ;
  gf->seqName = dictCreate (10000) ; // gfa parsing always requires sequence names
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
	case '#': // comment line
	  break ; 
	case 'H': // header line
	  word = fgetword (f) ;
	  if (line > 1) die ("H header line not at start of GFA file - line %d", line) ;
	  if (strlen (word) >= 6 && !strncmp (word, "VN:Z:", 5) && word[5] != '1')
	    die ("GFA file version %d - can only parse version 1 for now", word[5]) ;
	  break ;
	case 'S':
	  word = fgetword (f) ; // name
	  if (!dictAdd (gf->seqName, word, &index)) die ("duplicate S id line %d", line) ;
	  Seq *seq = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	  word = fgetword (f) ; // sequence
	  if (*word != '*')
	    { seq->len = strlen (word) ;
	      seq->dna = seqPack (SP, word, 0, seq->len) ;
	    }
	  word = fgetword (f) ; // length
	  if (strncmp (word, "LN:i:", 5)) die ("bad LN field line %d", line) ;
	  if (seq->dna)
	    { if (atoi(&word[5]) != seq->len) die ("length error line %d", line) ; }
	  else
	    seq->len = atoi (&word[5]) ;
	  if (seq->len > gf->bufSize)
	    { gf->bufSize = seq->len ;
	      if (gf->buf) free (gf->buf) ;
	      gf->buf = new(gf->bufSize,char) ;
	    }
	  // now make the reverse complement as the (n+1)th Seq
	  Seq *seq2 = arrayp(gf->seq, arrayMax(gf->seq), Seq) ;
	  seq2->len = seq->len ;
	  if (seq->dna) seq2->dna = seqRevCompPacked (seq->dna, 0, seq->len) ;
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
	  if (!*word) die ("empty cigar field line %d", line) ;
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
	case 'P':
	  if (!gf->path) gf->path = arrayCreate (10000, Path) ;
	  if (!gf->pathName) gf->pathName = dictCreate (10000) ;
	  word = fgetword (f) ; // name
	  if (!dictAdd (gf->pathName, word, &index)) die ("duplicate P id in line %d", line) ;
	  Path *p = arrayp(gf->path, arrayMax(gf->path), Path) ;
	  word = fgetword (f) ; // path name
	  int n = 1, i ;
	  for (s = word ; *s ; ++s) if (*s == ',') { ++n ; *s = 0 ; }
	  p->as = arrayCreate (n, int) ; arrayMax(p->as) = n ; // will use arr() below
	  for (i = 0, s = word ; i < n ; ++i)
	    { int k = strlen (s) - 1 ;
	      arr(p->as,i,int) = (s[k] == '+') ? 0 : 1 ;
	      s[k] = 0 ;
	      if (!dictFind (gf->seqName, s, &index)) die ("unknown seq %s in line %d", s, line) ;
	      arr(p->as,i,int) += 2*index ;
	    }
	  fgetword (f) ; // ignore the cigar field
	  while ((word = fgetword (f)) && *word) // optional fields
	    if (!strncmp (word, "ST:i:", 5)) p->start = atoi (&word[5]) ;
	    else if (!strncmp (word, "EN:i:", 5)) p->end = atoi (&word[5]) ;
	    else die ("unknown additional field %s in line %d", word, line) ;
	  break ;
	case 'A': // HiFiAsm makes these - I will ignore them for now
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

// need to write names and directions of segments and paths into GFA file
// write 1-based integers if they have no names

static inline char  seqDir (int i) { return (i&1) ? '-' : '+' ; }
static inline char* seqName (Gfa *gf, int i)
{ static char buf[64] ;
  if (gf->seqName) return dictName(gf->seqName, i/2) ;
  else { sprintf (buf, "%d", 1 + i/2) ; return buf ; }
}
static inline char* pathName (Gfa *gf, int i)
{ static char buf[64] ;
  if (gf->pathName) return dictName(gf->pathName, i) ;
  else { sprintf (buf, "%d", arrayMax(gf->seq)/2 + 1 + i) ; return buf ; }
}

void gfaWrite (Gfa *gf, char *file)
/* example from https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
H	VN:Z:1.0
S	11	ACCTT
S	12	TCAAGG
S	13	CTTGATT
L	11	+	12	-	4M
L	12	-	13	+	5M
L	11	+	13	+	3M
P	14	11+,12-,13+	4M,5M 
// in the P lines it is possible to have * for cigar, and get overlaps from links - I write this

// GFA 1.1 has path 'W' lines, which are nicer in some ways, but there is less support for these
// W walks need to have 0M end join (non-overlapping) segments
H	VN:Z:1.1
S	s11	ACCTT
S	s12	TC
S	s13	GATT
L	s11	+	s12	-	0M
L	s12	-	s13	+	0M
L	s11	+	s13	+	0M
W	NA12878	1	chr1	0	11	>s11<s12>s13
*/
{
  int i, j ;
  FILE *f = fzopen (file, "w") ;
  if (!f) die ("failed to open %s to write gfa file", file) ;

  for (i = 0 ; i < arrayMax(gf->seq) ; i += 2) /* only write sequences in one direction */
    { Seq *s = arrp(gf->seq, i, Seq) ;
      fprintf (f, "S\t%s", seqName (gf, i)) ;
      if (s->dna)
	fprintf (f, "\t%*s\n", s->len, seqUnpack (SP, s->dna, gf->buf, 0, s->len)) ;
      else fprintf (f, "\t*\tLN:i:%d\n", s->len) ;
    }

  for (i = 0 ; i < arrayMax(gf->link) ; ++i)
    { Link *l = arrp(gf->link, i, Link) ;
      fprintf (f, "L\t%s\t%c\t%s\t%c\t%dM\n",
	       seqName(gf,l->s1), seqDir(l->s1), seqName(gf,l->s2), seqDir(l->s2), l->overlap) ;
    }

  if (gf->path)
    for (i = 0 ; i < arrayMax(gf->path) ; ++i)
      { Path *p = arrp(gf->path, i, Path) ;
	fprintf (f, "P\t%s", pathName(gf,i)) ;
	assert (p->as && arrayMax(p->as) > 0) ;
	fprintf (f, "\t%s%c", seqName(gf,arr(p->as,0,int)), seqDir(arr(p->as,0,int))) ;
	for (j = 1 ; j < arrayMax(p->as) ; ++j)
	  fprintf (f, ",%s%c", seqName(gf,arr(p->as,j,int)), seqDir(arr(p->as,j,int))) ;
	fprintf (f, "\t*") ;
	if (p->start) fprintf (f, "\tST:i:%d", p->start) ;
	if (p->end)   fprintf (f, "\tEN:i:%d", p->end) ;
	fputc ('\n', f) ;
      }

  fclose (f) ;
}

/*************************************************************************/

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
	seq->dna = seqPack (SP, sqioSeq(si), 0, si->seqLen) ;
	++seq ; seq->dna = seqRevCompPacked (seq->dna, 0, si->seqLen) ;
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
      if (!seqMatchPacked (arrp(gf->seq,l->s1,Seq)->dna, seq1->len - l->overlap,
			   arrp(gf->seq,l->s2,Seq)->dna, 0, l->overlap))
	*ul++ = *l ; // perfect match
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

static int X = -1 ; /* debugging */

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

  // next create gf2: seqs from the segments between cuts, and links and paths from the seqs
  Gfa *gf2 = gfaCreate (arrayMax(gf1->seq)+arrayMax(gf1->link), 8*arrayMax(gf1->link)) ;
  gf2->path = arrayCreate (arrayMax(gf1->seq), Path) ;
  
  Hash newSeqHash = hashCreate (2*arrayMax(gf1->link)) ;
  for (is = 0 ; is < nSeq ; ++is)
    { Path *wi = arrayp(gf2->path, is, Path) ;
      wi->start = 0 ;
      wi->as = arrayCreate (arrayMax(cut[is])/2, int) ;
      if (is == X || is == RC(X)) // debug section
	{ printf ("seq %d len %d when building gf2\n  cut", is, arrp(gf1->seq,is,Seq)->len) ;
	  for (c = arrp(cut[is],0,Cut) ; c < arrp(cut[is],arrayMax(cut[is]),Cut) ; ++c)
	    printf (" (%d,%d,%d)%c", c->x, c->s, c->sx, c->isLeft?'L':'R') ;
	  putchar ('\n') ;
	}
      int i, index, indexRC, lastIndex = -1 ;
      for (i = 0 ; i < arrayMax(cut[is]) ; i += 2)
	{ c1 = arrp(cut[is],i,Cut) ;
	  c2 = arrp(cut[is],i+1,Cut) ;
	  assert (!c1->isLeft && c2->isLeft) ;
	  assert (c1->s == c2->s) ; 
	  Seq *sfs = arrp(gf1->seq,c1->s,Seq) ;
	  if (!hashFind (newSeqHash, HASH_INT2(c1->s,c1->sx), &index))
	    { hashAdd (newSeqHash, HASH_INT2(c1->s,c1->sx), &index) ;
	      hashAdd (newSeqHash, HASH_INT2(RC(c1->s),sfs->len - c2->sx), &indexRC) ;
	      if (is == X || is == RC(X))
		printf ("adding seq2 %d from %d-%d in %d and RC %d from %d to %d in %d\n",
			index, c1->x, c2->x, c1->s,
			indexRC, sfs->len - c2->sx, sfs->len - c1->x, RC(c1->s)) ;
	    }
	  Seq *sNew = arrayp(gf2->seq,index,Seq) ;
	  if (sNew->len)
	    assert (sNew->len == c2->x - c1->x) ;
	  else
	    { sNew->len = c2->x - c1->x ;
	      if (sNew->len > gf2->bufSize)
		{ gf2->bufSize = sNew->len ;
		  if (gf2->buf) free (gf2->buf) ;
		  gf2->buf = new(gf2->bufSize,char) ;
		}
	      char* tmp = seqUnpack (SP, sfs->dna, gf1->buf, c1->sx, sNew->len) ;
	      sNew->dna = seqPack (SP, tmp, 0, sNew->len) ;
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

  printf ("made blunt gfa with %d seqs, %d links and %d paths\n",
	  arrayMax(gf2->seq), arrayMax(gf2->link), arrayMax(gf2->path)) ;

  // now check the paths
  for (is = 0 ; is < nSeq ; ++is)
    { Seq *ps = arrp(gf1->seq, is, Seq) ; if (!ps->dna) continue ;
      Path *p = arrp(gf2->path, is, Path) ;
      if (is == X || is == RC(X)) // debug section
	printf ("is %d path %d segs %d start %d end\n",
		is, arrayMax(p->as), p->start, p->end) ;
      assert (p->start == 0 && p->end == 0) ; // code below currently requires this
      int i, j, n = 0 ;
#ifdef OLD
      char *s = seqUnpack (SP, ps->dna, gf1->buf, 0, ps->len) ;
      for (i = 0 ; i < arrayMax(p->as) ; ++i)
	{ Seq *fs = arrp(gf2->seq, arr(p->as,i,int), Seq) ;
	  if (is == X || is == RC(X)) // debug section
	    printf ("  seg %d len %d\n", arr(p->as,i,int), fs->len) ;
	  char *t = seqUnpack (SP, fs->dna, gf2->buf, 0, fs->len) ;
	  //	  printf ("seq %d pos %d starting cut %d len %d\n", is, n, i, fs->len) ;
	  for (j = 0 ; j < fs->len ; ++j, ++n)
	    if (*s++ != *t++)
	      die ("path mismatch: seq %d len %d pos %d frag %d base %d",
		   is, arrp(gf1->seq, is, Seq)->len, n, i, j) ;
	}
#else
      for (i = 0 ; i < arrayMax(p->as) ; ++i)
	{ Seq *fs = arrp(gf2->seq, arr(p->as,i,int), Seq) ;
	  if (is == X || is == RC(X)) // debug section
	    printf ("  seg %d len %d\n", arr(p->as,i,int), fs->len) ;
	  if ((j = seqMatchPacked (ps->dna, n, fs->dna, 0, fs->len)))
	    die ("path mismatch: seq %d len %d pos %d frag %d base %d",
		 is, arrp(gf1->seq, is, Seq)->len, n+j, i, j) ;
	  n += fs->len ;
	}
#endif
    }
  
  printf ("all paths check out\n") ;

  return gf2 ;
}

Gfa *chain (Gfa *gf) { return gf ; } // a stub for now

void usage (void)
{
  fprintf (stderr, "Usage: gaffer <commands>\n") ;
  fprintf (stderr, "   -read  <stem>   : read ONE files\n") ;
  fprintf (stderr, "   -write <stem>   : write ONE files\n") ;
  fprintf (stderr, "   -schema         : print ONE schema\n") ;
  fprintf (stderr, "   -readGfa <gfa file> : only reads S and L lines for now\n") ;
  fprintf (stderr, "   -writeGfa <gfa file>\n") ;
  fprintf (stderr, "   -readDna <sequence file matching S line names>\n") ;
  fprintf (stderr, " // options below change the graph\n") ;
  fprintf (stderr, "   -removeBadLinks : (for now) remove imperfect overlaps\n") ;
  fprintf (stderr, "   -blunt          : makes new non-overlapping graph\n") ;
  fprintf (stderr, "   -chain          : chain 1-1 links into unitigs\n") ;
  fprintf (stderr, "   -extend         : adds match blocks for shared incoming edges - needs seq\n") ;
  exit (0) ;
}

int main (int argc, char *argv[])
{
  Gfa *gf = 0 ;
  
  timeUpdate (0) ;
  SP = seqPackCreate ('a') ;

  argc-- ; ++argv ; // swallow the program name
  if (!argc) usage() ;
  while (argc)
    { if (**argv != '-') die ("argument %s is expected to start with '-'", *argv) ;
      char *command = *argv+1 ;
      if (!strcmp (*argv, "-read") && argc >= 2)
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
      else if (!strcmp (*argv, "-readGfa") && argc >= 2)
	{ if (gf) { fprintf (stderr, "removing existing gf\n") ; gfaDestroy (gf) ; }
	  gf = gfaParseSL (argv[1]) ;
	  linkRemoveDuplicates (gf) ;
	  oneAddProvenance (gf->vf, "gaffer", VERSION, "readGfa %s // nS %d nL %d nP %d", argv[1],
			    arrayMax(gf->seq), arrayMax(gf->link),
			    gf->path ? arrayMax(gf->path) : 0) ;
	  argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-writeGfa") && argc >= 2)
	{ gfaWrite (gf, argv[1]) ;
	  oneAddProvenance (gf->vf, "gaffer", VERSION, "writeGfa %s", argv[1]) ;
	  argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-readDna") && argc >= 2)
	{ if (gf)
	    { readSeqFile (gf, argv[1]) ;
	      oneAddProvenance (gf->vf, "gaffer", VERSION, "readDna %s", argv[1]) ;
	    }
	  else fprintf (stderr, "can't read sequences without a graph\n") ;
	  argc -= 2 ; argv += 2 ; 
	}
      else if (!strcmp (*argv, "-removeBadLinks"))
	{ if (gf)
	    { linkRemoveBad (gf) ;
	      oneAddProvenance (gf->vf, "gaffer", VERSION, "removeBadLinks // nL %d",
				arrayMax (gf->link)) ;
	    }
	  else fprintf (stderr, "can't remove bad links without a graph\n") ;
	  argc-- ; argv++ ;
	}
      else if (!strcmp (*argv, "-blunt"))
	{ if (gf)
	    { Gfa *gf2 = bluntify (gf) ;
	      oneAddProvenance (gf->vf, "gaffer", VERSION, "bluntify // nS %d nL %d nP %d",
				arrayMax(gf->seq), arrayMax(gf->link),
				gf->path ? arrayMax(gf->path) : 0) ;
	      gf2->vf = gf->vf ; gf->vf = 0 ; 
	      gfaDestroy (gf) ;
	      gf = gf2 ;
	    }
	  else fprintf (stderr, "can't bluntify without a graph\n") ;
	  argc-- ; argv++ ;
	}
      else if (!strcmp (*argv, "-chain"))
	{ if (gf)
	    { Gfa *gf2 = chain (gf) ;
	      oneAddProvenance (gf->vf, "gaffer", VERSION, "chain // nS %d nL %d",
				arrayMax(gf->seq), arrayMax(gf->link)) ;
	      gf2->vf = gf->vf ; gf->vf = 0 ; 
	      gfaDestroy (gf) ;
	      gf = gf2 ;
	    }
	  else fprintf (stderr, "can't chain without a graph\n") ;
	  argc-- ; argv++ ;
	}
      else die ("unrecognised option %s - run without args for help", *argv) ;
      fprintf (stderr, "%s: ", command) ; timeUpdate (stderr) ;
    }

  fprintf (stderr, "cleaning up\n") ;
  seqPackDestroy (SP) ;
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
  "P 3 pth                 PATH\n"
  "O P 1 8 INT_LIST        list of segments\n"
  "D D 1 6 STRING          list of directions of each segment (+ or -) - required\n"
  "D I 1 6 STRING          path identifier (required by GFA) - optional here\n"
  "D G 1 11 STRING_LIST    list of cigar strings for overlaps - optional/deprecated\n"
  "D S 1 3 INT             start in first segment - defaults to 0 if missing\n"
  "D E 1 3 INT             end distance from end of last segment - defaults to 0 if missing\n"
  ".\n"
  "P 3 wlk                 WALK - in GFA 1.1 only - requires non-overlapping segments\n"
  "O W 2 3 INT 8 INT_LIST  length, list of segments\n"
  "D D 1 6 STRING          list of directions of each segment (+ or -) - required\n"
  "D I 1 6 STRING          identifier (required by GFA)\n"
  "D J 1 6 STRING          sample identifier (required by GFA)\n"
  "D H 1 3 INT             haplotype : 1..<ploidy> - missing if 0 in GFA = haploid or unknown\n"
  "D S 1 3 INT             start in first segment - defaults to 0 if missing\n"
  "D E 1 3 INT             end distance from end of last segment - defaults to 0 if missing\n"
  ".\n"
  "P 3 ctn                                  CONTAINMENT - contained and container are both segs\n"
  "O L 5 3 INT 4 CHAR 3 INT 4 CHAR 3 INT    s1 dir1 s2 dir2 pos - s1,2 in seg file dir=+|- start position\n"
  "D G 1 6 STRING                           cigar string - else presume exact match\n"
  "D I 1 6 STRING                           ID edge identifier (deprecated)\n"
  "D R 1 3 INT                              RC read count\n"
  "D M 1 3 INT                              NM number of mismatches\n"
  ".\n"
  "P 3 aln                 ALIGNMENT - of sequences from a sequence file, e.g. reads\n"
  "O A 2 3 INT 8 INT_LIST  index of seq to align, list of segments as in WALK\n"
  "D D 1 6 STRING          list of directions of each segment (+ or -) - required\n"
  "D S 1 3 INT             start - as in WALK\n"
  "D E 1 3 INT             end - as in WALK\n"
  "D G 1 6 STRING          cigar string - for mapping\n"
  "D Q 1 3 INT             mapping quality\n"
  "D M 1 3 INT             number of mismatches\n"
  ".\n" ;
  
/******************* end of file *******************/
