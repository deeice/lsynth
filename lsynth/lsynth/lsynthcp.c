//---------------------------------------------------------------------------

#pragma hdrstop

/*
 * Program Name: LSynth
 * Program Description:
 *   LDraw CAD system compatible flexible part synthesizer.  This program
 *   reads your LDraw file, searching for unofficial META commands (structured
 *   comments), that specify things you want synthesized.
 *
 *   LSynth has two primary forms of synthesis.  The first kind of synthesis
 *   creates hose like things including:
 *     rubbed, rubber, pneumatic, and flex system hoses, as well as electric
 *     and fiber optic cables, flex system cables, flexible axles, string, and
 *     minifig chain.
 *
 *   It also creates things that travel around circular lego parts.  Things
 *   like:
 *     rubber band, rubber belt, technic chain, technic plastic tread, technic
 *     rubber tread.
 *
 *   The files tube.c, tube.h, curve.c and curve.h perform hose synthesis.
 *   The files band.c and band.h perform band synthesis.
 *
 *   This file (main.c) contains the main entry/exit points for the program.
 *   It opens and scans the LDraw file provided, identifies synthesis
 *   synthesis specifications and hands them off to the appropriate synthesis
 *   methodology.
 */

#pragma hdrstop

#include "lsynthcp.h"
#include "hose.h"
#include "band.h"
#include "ctype.h"
#include "string.h"

//---------------------------------------------------------------------------
void strupper(char *s) {
  char *p;
  for (p = s; *p; p++) {
    *p = toupper(*p);
  }
}

//---------------------------------------------------------------------------
/* If this code works, it was written by Lars C. Hassing. */
/* If not, I don't know who wrote it.                     */

/* Like fgets, except that 1) any line ending is accepted (\n (unix),
\r\n (DOS/Windows), \r (Mac (OS9)) and 2) Str is ALWAYS zero terminated
(even if no line ending was found) */
//---------------------------------------------------------------------------
char *L3fgets(char *Str, int n, FILE *fp)
{
   register int   c;
   int            nextc;
   register char *s = Str;

   while (--n > 0)
   {
      if ((c = getc(fp)) == EOF)
         break;
      if (c == '\032')
         continue;              /* Skip CTRL+Z                               */
      if (c == '\r' || c == '\n')
      {
         *s++ = '\n';
         /* We got CR or LF, eat next character if LF or CR respectively */
         if ((nextc = getc(fp)) == EOF)
            break;
         if (nextc == c || (nextc != '\r' && nextc != '\n'))
            ungetc(nextc, fp);  /* CR-CR or LF-LF or ordinary character      */
         break;
      }
      *s++ = c;
   }
   *s = 0;

   /* if (ferror(fp)) return NULL; if (s == Str) return NULL; */
   if (s == Str)
      return NULL;

   return Str;
}

//---------------------------------------------------------------------------
char *
fgetline(
  char *line,
  int   len,
  FILE *file)
{
  char *rc;
  while (rc = L3fgets(line,len,file)) {
    char caps[256];
    char *nonwhite;

    strcpy(caps,line);
    strupper(caps);

    nonwhite = caps + strspn(caps," \t");

    if (strncmp(nonwhite,"0 ROTATION C",strlen("0 ROTATION C")) == 0 ||
        strncmp(nonwhite,"0 COLOR",strlen("0 COLOR")) == 0) {
      continue;
    }

    nonwhite = line + strspn(line," \t");
    if (strncmp(nonwhite,"0 WRITE ",strlen("0 WRITE ")) == 0) {
      strcpy(nonwhite + 2, nonwhite + strlen("0 WRITE "));
    }
    break;
  }
  return rc;
}

//---------------------------------------------------------------------------
void
strclean(char *str)
{
  if (strncmp(str,"0 WRITE ",strlen("0 WRITE ")) == 0) {
    strcpy(str + 2, str + strlen("0 WRITE "));
  }
}

//---------------------------------------------------------------------------

/*
 * Skip over MLCad ROTATION and COLOR statements
 */

static int skip_rot(char *line, int n_line, FILE *dat, FILE *temp)
{
  char *nonwhite;
  char  caps[512];

  strcpy(caps,line);

  nonwhite = caps + strspn(caps," \t");

  while (strncmp(nonwhite,"0 ROTATION C",strlen("0 ROTATION C")) == 0 ||
         strncmp(nonwhite,"0 COLOR",strlen("0 COLOR")) == 0) {
    fputs(line,temp);
    L3fgets(line,n_line,dat);  /* FIXME: check L3fgets rc */
    strcpy(caps,line); strupper(caps); nonwhite = caps + strspn(caps," \t");
  }

  return 0;
}


/*****************************************************************************
 *
 * Read in and parse up synthesis descriptions and constraints.
 *
 ****************************************************************************/

int
parse_descr(char *fullpath_progname)
{
  char filename[256];
  FILE *mpd;
  char line[256];

  strcpy(filename,fullpath_progname);
  {
    char *l, *p;

    for (l = p = filename; *p; *p++) {
      if ((*p == '\\') || (*p == '/')) {
        l = p+1;
      }
    }
    *l = '\0';
  }
  strcat(filename,"lsynth.mpd");

  mpd = fopen(filename,"r");

  if (mpd == NULL) {
    printf("Failed to open lsynth.mpd for reading\n");
    return -1;
  }

  while(fgetline(line,sizeof(line),mpd)) {
    char stretch[64];
    char type[64];
    int  d,st,i;
    PRECISION s,t;

    if (sscanf(line,"0 SYNTH BEGIN DEFINE %s HOSE %s %d %d %f\n",
        type,stretch,&d,&st,&t) == 5) {
      if (strcmp(stretch,"STRETCH") == 0) {
        hose_types[n_hose_types].fill = STRETCH;
      } else if (strcmp(stretch,"FIXED") == 0) {
        hose_types[n_hose_types].fill = FIXED;
      } else {
        printf("Error: Unrecognized fill type %s for hose type %s.  Aborting\n",
          stretch,type);
        fclose(mpd);
        return -1;
      }

      strcpy(hose_types[n_hose_types].type,type);
      hose_types[n_hose_types].diameter = d;
      hose_types[n_hose_types].stiffness = st;
      hose_types[n_hose_types].twist = t;

      for (i = 0; i < 3; i++) {
        part_t *part;

        if (i == 0) {
          part = &hose_types[n_hose_types].start;
        } else if (i == 1) {
          part = &hose_types[n_hose_types].mid;
        } else {
          part = &hose_types[n_hose_types].end;
        }

        if (fgetline(line,sizeof(line),mpd)) {

          int n;

          n = sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
                &part->attrib,
                &part->offset[0],    &part->offset[1],    &part->offset[2],
                &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
                &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
                &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
                 part->type);

          if (n != 14) {
            printf("Error: Failed to parse type one line.  Got this instead:\n");
            printf(line);
            fclose(mpd);
            return -1;
          }
        } else {
          printf("Error: Unexpected end of file\n");
          fclose(mpd);
          return -1;
        }
      }

      if (fgetline(line,sizeof(line),mpd)) {
        if (strcmp(line,"0 SYNTH END\n") != 0) {
          printf("Error: Expected SYNTH END, got this instead\n");
          printf(line);
          fclose(mpd);
          return -1;
        }
      } else {
        printf("Error: Unexepcted end of file\n");
        fclose(mpd);
        return -1;
      }
      n_hose_types++;
    } else if (strcmp(line,"0 SYNTH BEGIN DEFINE HOSE CONSTRAINTS\n") == 0) {
      while(fgetline(line,sizeof(line),mpd)) {
        part_t *part = &hose_constraints[n_hose_constraints];
        if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          &part->attrib,
          &part->offset[0],    &part->offset[1],    &part->offset[2],
          &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
          &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
          &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
           part->type) == 14) {
          n_hose_constraints++;
        } else if (strcmp(line,"0 SYNTH END\n") == 0) {
          break;
        } else {
          printf("Error: Failed to parse type one line.  Got this instead:\n");
          printf(line);
          fclose(mpd);
          return -1;
        }
      }
    } else if (sscanf(line,"0 SYNTH BEGIN DEFINE %s BAND %s %f %f\n",
                 type,stretch,&s,&t) == 4) {
      int n;

      if (strcmp(stretch,"STRETCH") == 0) {
        band_types[n_band_types].fill = STRETCH;
        n = 2;
      } else if (strcmp(stretch,"FIXED") == 0) {
        band_types[n_band_types].fill = FIXED;
        n = 2;
      } else if (strcmp(stretch,"FIXED3") == 0) {
        band_types[n_band_types].fill = FIXED3;
        n = 4;
      } else {
        printf("Error: Unrecognized fill type %s for hose type %s.  Aborting\n",
          stretch,type);
        fclose(mpd);
        return -1;
      }

      strcpy(band_types[n_band_types].type,type);
      band_types[n_band_types].scale = s;
      band_types[n_band_types].thresh = t;

      for (i = 0; i < n; i++) {
        part_t *part;

        if (i == 0) {
          part = &band_types[n_band_types].tangent;
        } else if (i == 1) {
          part = &band_types[n_band_types].arc;
        } else if (i == 2) {
          part = &band_types[n_band_types].start_trans;
        } else {
          part = &band_types[n_band_types].end_trans;
        }

        if (fgetline(line,sizeof(line),mpd)) {
          if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
            &part->attrib,
            &part->offset[0],    &part->offset[1],    &part->offset[2],
            &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
            &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
            &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
             part->type) != 14) {
            printf("Error: Failed to parse type one line.  Got this instead:\n");
            printf(line);
            fclose(mpd);
            return -1;
          }
        } else {
          printf("Error: Unexpected end of file\n");
          fclose(mpd);
          return -1;
        }
      }

      if (L3fgets(line,sizeof(line),mpd)) {
        if (strcmp(line,"0 SYNTH END\n") != 0) {
          printf("Error: Expected SYNTH END, got this instead\n");
          printf(line);
          fclose(mpd);
          return -1;
        }
      } else {
        printf("Error: Unexepcted end of file\n");
        fclose(mpd);
        return -1;
      }
      n_band_types++;

    } else if (strcmp(line,"0 SYNTH BEGIN DEFINE BAND CONSTRAINTS\n") == 0) {
      while(fgetline(line,sizeof(line),mpd)) {
        part_t *part = &band_constraints[n_band_constraints];
        if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          &part->attrib,
          &part->offset[0],    &part->offset[1],    &part->offset[2],
          &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
          &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
          &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
           part->type) == 14) {
          n_band_constraints++;
        } else if (strcmp(line,"0 SYNTH END\n") == 0) {
          break;
        } else {
          printf("Error: Failed to parse type one line.  Got this instead:\n");
          printf(line);
          fclose(mpd);
          return -1;
        }
      }
    }
  }

  fclose(mpd);
  return 0;
}

/*
 * Skip over results rom previous syntesis efforts.
 */

int skip_synthesized(FILE *dat, char *line, int sizeof_line)
{
  int rc;
  char *nonwhite;

  if (strcmp(line,"0 SYNTHESIZED BEGIN\n") == 0) {
    int rc;

    while (L3fgets(line,sizeof(line),dat)) {
      strupper(line); nonwhite = line + strspn(line," \t");
      if (strcmp(nonwhite,"0 SYNTHESIZED END\n") == 0) {
        return 0;
      }
    }
    return -1;
  } else {
    return 0;
  }
}

/*
 * Gather constraints from hose synthesis
 */

static
int synth_hose_class(
  char *nonwhite,
  FILE *dat,
  FILE *temp)
{
  char   type[256];
  char   line[512];
  char   caps[512];
  part_t constraints[128];
  int    constraint_n = 0;
  int    color;
  char   start_type[64];
  float  a,b,c, d,e,f, g,h,i, j,k,l;
  int    rc = 0;
  int    hide = 1;
  int    hose_color;
  int    ghost = 0;

  memset(constraints, 0, 128*sizeof(part_t));

  /* parse up remainder of line and convert units to LDU */

  if (sscanf(nonwhite,"%s %d\n",type,&hose_color) != 2) {
    return -1;
  }

  /* gather up the constraints */

  while (L3fgets(line,sizeof(line), dat)) {
    strcpy(caps,line); strupper(caps); nonwhite = caps + strspn(caps," \t");

    skip_rot(line,sizeof(line),dat,temp);

    strcpy(caps,line); strupper(caps); nonwhite = caps + strspn(caps," \t");

    if (strncmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      nonwhite += strspn(nonwhite," \t");
    }

    strclean(nonwhite);

    if (sscanf(nonwhite,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s",
        &color, &a,&b,&c, &d,&e,&f, &g,&h,&i, &j,&k,&l, start_type) == 14) {

      strupper(start_type);

      if (ishoseconstraint(start_type)) {

        if (hide) {
          fputs("0 ",temp);
        }
        fputs(line,temp);

        constraints[constraint_n].offset[0] = a;
        constraints[constraint_n].offset[1] = b;
        constraints[constraint_n].offset[2] = c;

        constraints[constraint_n].orient[0][0] = d;
        constraints[constraint_n].orient[0][1] = e;
        constraints[constraint_n].orient[0][2] = f;
        constraints[constraint_n].orient[1][0] = g;
        constraints[constraint_n].orient[1][1] = h;
        constraints[constraint_n].orient[1][2] = i;
        constraints[constraint_n].orient[2][0] = j;
        constraints[constraint_n].orient[2][1] = k;
        constraints[constraint_n].orient[2][2] = l;

        strcpy(constraints[constraint_n].type,start_type);
        constraint_n++;
      }
    } else if (strcmp(nonwhite,"0 SYNTH HIDE\n") == 0) {
      fputs(line,temp);
      hide = 1;
    } else if (strcmp(nonwhite,"0 SYNTH SHOW\n") == 0) {
      fputs(line,temp);
      hide = 0;
    } else {
      break;
    }
  }

  rc = skip_synthesized(dat, line, sizeof(line));

  if (rc == 0) {

    strcpy(caps,line); strupper(caps); nonwhite = caps + strspn(caps, " \t");
    strclean(nonwhite);

    if (strcmp(nonwhite,"0 SYNTH END\n") == 0) {
      rc = synth_hose(type,constraint_n,constraints,ghost,hose_color,temp);
      fputs(line,temp);
    }
  }
  return rc;
}

/*
 * Gather up rubber band constraints
 */

static
int synth_band_class(
  char *nonwhite,
  FILE *dat,
  FILE *temp)
{
  char type[256];
  char line[512];
  char caps[512];
  LSL_band_constraint constraints[128];
  int  constraint_n = 0;
  int  color = 0;
  int  rc = 0;
  int  hide = 0;
  int  ghost = 0;

  memset(constraints, 0, sizeof(LSL_band_constraint)*128);

  if (sscanf(nonwhite,"%s %d\n",type,&color) != 2) {
    return -1;
  }

  /* gather up the constraints */

  while (L3fgets(line,sizeof(line), dat)) {
    float a,b,c, d,e,f, g,h,i, j,k,l;
    char start_type[64];
    int t;

    strcpy(caps,line); strupper(caps); nonwhite = caps + strspn(caps, " \t");

    skip_rot(line,sizeof(line),dat,temp);

    strcpy(caps,line); strupper(caps); nonwhite = caps + strspn(caps, " \t");

    if (strncmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      nonwhite += strspn(nonwhite," \t");
    }

    if (sscanf(nonwhite,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s",
        &t, &a,&b,&c, &d,&e,&f, &g,&h,&i, &j,&k,&l, start_type) == 14) {

      if (isbandconstraint(start_type)) {
        if (hide) {
          fputs("0 ",temp);
        }
        fputs(line,temp);

        constraints[constraint_n].part.attrib = color;
        strcpy(constraints[constraint_n].part.type,start_type);

        constraints[constraint_n].part.offset[0] = a;
        constraints[constraint_n].part.offset[1] = b;
        constraints[constraint_n].part.offset[2] = c;

        constraints[constraint_n].part.orient[0][0] = d;
        constraints[constraint_n].part.orient[0][1] = e;
        constraints[constraint_n].part.orient[0][2] = f;
        constraints[constraint_n].part.orient[1][0] = g;
        constraints[constraint_n].part.orient[1][1] = h;
        constraints[constraint_n].part.orient[1][2] = i;
        constraints[constraint_n].part.orient[2][0] = j;
        constraints[constraint_n].part.orient[2][1] = k;
        constraints[constraint_n].part.orient[2][2] = l;

        strcpy(constraints[constraint_n].part.type,start_type);
        constraint_n++;
      }
    } else {

      strclean(nonwhite);

      if (strcmp(nonwhite,"0 SYNTH INSIDE\n") == 0) {
        fputs(line,temp);
        strcpy(constraints[constraint_n].part.type,"INSIDE");
        constraint_n++;
      } else if (strcmp(nonwhite,"0 SYNTH OUTSIDE\n") == 0) {
        fputs(line,temp);
        strcpy(constraints[constraint_n].part.type,"OUTSIDE");
        constraint_n++;
      } else if (strcmp(nonwhite,"0 SYNTH CROSS\n") == 0) {
        fputs(line,temp);
        strcpy(constraints[constraint_n].part.type,"CROSS");
        constraint_n++;
      } else if (strcmp(nonwhite,"0 SYNTH HIDE\n") == 0) {
        fputs(line,temp);
        hide = 1;
      } else if (strcmp(nonwhite,"0 SYNTH SHOW\n") == 0) {
        fputs(line,temp);
        hide = 0;
      } else {
        break;
      }
    }
  }

  rc = skip_synthesized(dat, line, sizeof(line));
  strcpy(caps,line); strupper(caps); strclean(caps);

  if (rc == 0) {

    nonwhite = caps + strspn(caps, " \t");

    if (strcmp(nonwhite,"0 SYNTH END\n") == 0 ||
        strcmp(nonwhite,"0 WRITE SYNTH END\n") == 0) {

      rc = synth_band(type,constraint_n,constraints,color,temp,ghost);
      fputs(line,temp);
    }
  }
  return 0;
}

PRECISION hose_res_angle = 0.01;
PRECISION band_res = 2;

//---------------------------------------------------------------------------

char * stripquotes(char *s)
{
  char *p;
  int i;

  // Strip away leading whitespace (spaces and tabs).
  s += strspn(s, " \t");

  // Remove leading quotes
  if (*s == '\"')
    s++;

  // Allocate memory so we can modify the end of the string.
  s = strdup(s);

  // Eliminate trailing whitespace

  for (p = s + (strlen(s)-1); p >= s; p--) {
    if ((*p == ' ') || (*p == '\t')) {
      *p = 0;
    } else {
      break;
    }
  }

  // Remove trailing quotes.
  if ((p = strrchr(s, '\"')) != NULL) {
    *p = 0;
  }

  return(s);
}

#pragma argsused
int main(int argc, char* argv[])
{
  char *dat_name = argv[1];
  char *dst_name = argv[2];
  char *synth_name = NULL;
  char  filename[512];
  FILE *outfile;
  FILE *synthfile;
  FILE *dat;
  char line[512];
  char caps[512];
  char *nonwhite;
  int   synthcount = 0;
  int   subfiles = 0;

  /*
   * Read in the descriptions and constraints for synthesis
   */

  if (parse_descr(argv[0])) {
    return 1;
  }

  /*
   * Output MLCad.ini for LSynth
   */

  if (argc == 2 && strcmp(argv[1],"-m") == 0) {
    extern void hose_ini(void);
    extern void band_ini(void);
    char path[256];
    char *l,*p;
    int i;

    strcpy(path,argv[0]);

    for (i = 0; i < 2; i++) {
      for (p = path; *p; p++) {
        if (*p == '\\') {
          l = p;
        }
      }
      *l = '\0';
    }

    printf("[LSYNTH]\n");
    printf("%%PATH = \"%s\"\n",path);
    hose_ini();
    band_ini();
    printf("Tangent Statement: INSIDE = SYNTH INSIDE\n");
    printf("Tangent Statement: OUTSIDE = SYNTH OUTSIDE\n");
    printf("Tangent Statement: CROSS = SYNTH CROSS\n");
    printf("Visibility Statement: SHOW = SYNTH SHOW\n");
    printf("Visibility Statement: HIDE = SYNTH HIDE\n");
    return 1;
  }

  printf("LSynth version 2.2 by Kevin Clague, kevin_clague@yahoo.com\n");
  printf("                   and Don Heyse\n");

  if (argc == 2 && strcmp(argv[1],"-v") == 0) {
    return 1;
  }

  /*
   * Print out help display
   */

  if (argc == 2 && strcmp(argv[1],"-h") == 0) {
    extern void list_hose_types(void);
    extern void list_band_types(void);
    extern void list_band_constraints(void);

    printf("LSynth is an LDraw compatible flexible part synthesizer\n");
    printf("  usage: lsynthcp [-v] [-h] [-m] <src> <dst>\n");
    printf("    -v - prints lsynthcp version\n");
    printf("    -h - prints this help message\n");
    printf("    -m - prints out the LSynth portion of the MLcad.ini for using\n");
    printf("         this program\n");
    printf("The easiest way to use LSynth is from within MLcad.  You need to\n");
    printf("make additions to MLCad.ini.  Please see Willy Tscager's tutorial\n");
    printf("page http://www.holly-wood.it/mlcad/lsynth-en.html.\n");
    printf("\n");
    printf("To create a flexible part, you put specifications for the part\n");
    printf("directly into your LDraw file, where the part is needed.\n");

    list_hose_types();
    list_hose_constraints();
    list_band_types();
    list_band_constraints();
    return 1;
  }

  if (argc < 3) {
    printf("usage: lsynth <input_file> <output_file>\n");
    return 1;
  }

  dat = fopen(dat_name,"r");

  if (dat == NULL) {
    printf("%s: Failed to open file %s for reading\n",argv[0],dat_name);
    return -1;
  }

  outfile = fopen(dst_name,"w");

  if (outfile == NULL) {
    printf("%s: Failed to open file %s for writing\n",argv[0],dst_name);
    return -1;
  }

  /*
   * Scan the input file looking for synthesis specifications
   */

  while (L3fgets(line,sizeof(line), dat)) {

    fputs(line,outfile);

    strcpy(caps,line); strupper(caps); nonwhite = caps + strspn(caps, " \t");
    strclean(nonwhite);

    if (strncmp(nonwhite,"0 SYNTH BEGIN ", strlen("0 SYNTH BEGIN ")) == 0) {
      char *option;

      nonwhite += strlen("0 SYNTH BEGIN ");

      synthfile = outfile; // By default synth data goes to main outfile.
      synthcount++; // Count the synth parts.

      // We may be writing synth data to separate subfiles.
      if (subfiles) {
        // Check if the subfile is named in the SYNTH BEGIN line.
        option = strstr(nonwhite, "NAME=");
        if (option) {
          option += strlen("NAME=");
          // Allocate memory for name, strip quotes and trailing spaces.
          synth_name = stripquotes(option);
        } else {
         // Otherwise just number the subfile.  Use lsynthN.ldr for stdout?
          char *p;

          // Remove the extension from dat_name if not already gone.
          if ((p = strrchr(dst_name, '.')) != NULL)
            *p = 0;
          // Build the new subfilename.
          sprintf(filename, "%s%d.ldr", dst_name, synthcount);
          synth_name = strdup(filename);
        }

        synthfile = fopen(synth_name,"w");
        if (synthfile == NULL) {
          printf("%s: Failed to open file %s for writing\n",argv[0],synth_name);
          return -1;
        }

        fputs(line, synthfile); // Warning, this line has been uppercased
      }
      option = strstr(nonwhite, "LENGTH=");
      if (option) {
        option += strlen("LENGTH=");
      }
      option = strstr(nonwhite, "UNITS=");
      if (option) {
        option += strlen("UNITS=");
      }

      /* check to see if it is a known synth command */

      if (ishosetype(nonwhite)) {
        synth_hose_class(nonwhite,dat,synthfile);

      } else if (isbandtype(nonwhite)) {
        synth_band_class(nonwhite,dat,synthfile);

      } else {
        printf("Unknown synthesis type %s\n",nonwhite);
      }

      // Close subfile and cleanup

      if (synth_name) {
        fclose(synthfile);
        free(synth_name);
        synth_name = NULL;
      }
    } else {
      float foo;

      if (sscanf(nonwhite,"0 SYNTH HOSE_RES %f",&foo) == 1) {
        hose_res_angle = foo;
      } else if (sscanf(nonwhite,"0 WRITE SYNTH HOSE_RES %f",&foo) == 1) {
        hose_res_angle = foo;
      } else if (sscanf(nonwhite,"0 SYNTH BAND_RES %f",&foo) == 1) {
        band_res = foo;
      } else if (sscanf(nonwhite,"0 WRITE SYNTH BAND_RES %s",&foo) == 1) {
        band_res = foo;
      }
    }
  }
  fclose(dat);
  fclose(outfile);

  printf("LSynth complete\n");
  return 0;
}

