//---------------------------------------------------------------------------

#pragma hdrstop

#include "lsynth.h"
#include "tube.h"
#include "band.h"
#include "ctype.h"

void strupper(char *s) {
  char *p;
  for (p = s; *p; p++) {
    *p = toupper(*p);
  }
}
//---------------------------------------------------------------------------

static int skip_rot(char *line, int n_line, FILE *dat, FILE *temp)
{
  char *nonwhite;

  nonwhite = line + strspn(line, " \t");

  while (strncmp(nonwhite,"0 ROTATION C",strlen("0 ROTATION C")) == 0) {
    fputs(line,temp);
    fgets(line,n_line,dat);  /* FIXME: check fgets rc */
    strupper(line);
    nonwhite = line + strspn(line, " \t");
  }

  return 0;
}

int skip_synthesized(FILE *dat, char *line, int sizeof_line)
{
  int rc;
  char *nonwhite;

  nonwhite = line + strspn(line, " \t");

  if (strcmp(nonwhite,"0 SYNTHESIZED BEGIN\n") == 0) {
    int rc;

    while (fgets(line,sizeof(line),dat)) {
      strupper(line);
      nonwhite = line + strspn(line, " \t");
      if (strcmp(nonwhite,"0 SYNTHESIZED END\n") == 0) {
        return 0;
      }
    }
    return -1;
  } else {
    return 0;
  }
}

static char *tube_types[] = {
  ALL_HOSES,
  "",
};

static
int synth_tube_like(
  char *type,
  char *resid,
  FILE *dat,
  FILE *temp)
{
  char line[512];
  char *nonwhite;
  LSL_part_usage constraints[128];
  int constraint_n = 0;
  int color;
  char start_type[64];
  float a,b,c, d,e,f, g,h,i, j,k,l;
  int rc = 0;
  float length;
  char units[64];
  int  hide = 1;
  int hose_color;
  int ghost = 0;

  /* parse up remainder of line and convert units to LDU */

  if (sscanf(resid,"%d\n",&hose_color) != 1) {
    return -1;
  }

  length = 0;
  units[0] = '\0';
  *resid = '\0';
  /* gather up the constraints */

  while (fgets(line,sizeof(line), dat)) {
    strupper(line);

    nonwhite = line + strspn(line, " \t");

    skip_rot(line,sizeof(line),dat,temp);

    nonwhite = line + strspn(line, " \t");

    if (strncmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      nonwhite += strspn(nonwhite, " \t");
    }

    if (sscanf(nonwhite,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s",
        &color, &a,&b,&c, &d,&e,&f, &g,&h,&i, &j,&k,&l, start_type) == 14) {
      strupper(start_type);
      if (strcmp(start_type,"LS00.DAT") == 0) {

        if (hide) {
          fputs("0 ",temp);
        }
        fputs(line,temp);

        constraints[constraint_n].loc.x = a;
        constraints[constraint_n].loc.y = b;
        constraints[constraint_n].loc.z = c;

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
      } else {
        fputs(line,temp);
      }
    } else if (strcmp(nonwhite,"0 SYNTH HIDE\n") == 0 ||
               strcmp(nonwhite,"0 WRITE SYNTH HIDE\n") == 0) {
      fputs(line,temp);
      hide = 1;
    } else if (strcmp(nonwhite,"0 SYNTH SHOW\n") == 0 ||
               strcmp(nonwhite,"0 WRITE SYNTH SHOW\n") == 0) {
      fputs(line,temp);
      hide = 0;
    } else {
      break;
    }
  }

  rc = skip_synthesized(dat, line, sizeof(line));

  if (rc == 0) {

    nonwhite = line + strspn(line, " \t");

    if (strcmp(nonwhite,"0 SYNTH END\n") == 0 ||
        strcmp(nonwhite,"0 WRITE SYNTH END\n") == 0) {

      rc = synth_tube(type,constraint_n,constraints,length,units,hose_color,temp,ghost);
      fputs(line,temp);
    }
  }
  return rc;
}

static char *band_types[] = {
  ALL_BANDS,
  "",
};

static
int synth_band_like(
  char *type,
  char *resid,
  FILE *dat,
  FILE *temp)
{
  char line[512];
  char *nonwhite;
  LSL_band_constraint constraints[128];
  int constraint_n = 0;
  int color = 0;
  int t;
  char start_type[64];
  float a,b,c, d,e,f, g,h,i, j,k,l;
  int rc = 0;
  float length;
  char units[64];
  int  hide = 0;
  int  ghost = 0;
  /* parse up remainder of line and convert units to LDU */

  if (sscanf(resid,"%d\n",&color) != 1) {
    return -1;
  }
  *resid = '\0';

  /* gather up the constraints */

  while (fgets(line,sizeof(line), dat)) {
    strupper(line);

    nonwhite = line + strspn(line, " \t");

    skip_rot(line,sizeof(line),dat,temp);

    nonwhite = line + strspn(line, " \t");

    if (strncmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      nonwhite += strspn(nonwhite, " \t");
    }

    if (sscanf(nonwhite,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s",
        &t, &a,&b,&c, &d,&e,&f, &g,&h,&i, &j,&k,&l, start_type) == 14) {

      if (hide) {
        fputs("0 ",temp);
      }
      fputs(line,temp);

      constraints[constraint_n].part.color = color;
      strcpy(constraints[constraint_n].part.type,start_type);

      constraints[constraint_n].part.loc.x = a;
      constraints[constraint_n].part.loc.y = b;
      constraints[constraint_n].part.loc.z = c;

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
    } else if (strcmp(nonwhite,"0 SYNTH INSIDE\n") == 0 ||
               strcmp(nonwhite,"0 WRITE SYNTH INSIDE\n") == 0) {
      fputs(line,temp);
      strcpy(constraints[constraint_n].part.type,"INSIDE");
      constraint_n++;
    } else if (strcmp(nonwhite,"0 SYNTH OUTSIDE\n") == 0 ||
               strcmp(nonwhite,"0 WRITE SYNTH OUTSIDE\n") == 0) {
      fputs(line,temp);
      strcpy(constraints[constraint_n].part.type,"OUTSIDE");
      constraint_n++;
    } else if (strcmp(nonwhite,"0 SYNTH CROSS\n") == 0 ||
               strcmp(nonwhite,"0 WRITE SYNTH CROSS\n") == 0) {
      fputs(line,temp);
      strcpy(constraints[constraint_n].part.type,"CROSS");
      constraint_n++;
    } else if (strcmp(nonwhite,"0 SYNTH HIDE\n") == 0 ||
               strcmp(nonwhite,"0 WRITE SYNTH HIDE\n") == 0) {
      fputs(line,temp);
      hide = 1;
    } else if (strcmp(nonwhite,"0 SYNTH SHOW\n") == 0 ||
               strcmp(nonwhite,"0 WRITE SYNTH SHOW\n") == 0) {
      fputs(line,temp);
      hide = 0;
    } else {
      break;
    }
  }

  rc = skip_synthesized(dat, line, sizeof(line));

  if (rc == 0) {

    nonwhite = line + strspn(line, " \t");

    if (strcmp(nonwhite,"0 SYNTH END\n") == 0 ||
        strcmp(nonwhite,"0 WRITE SYNTH END\n") == 0) {

      rc = synth_band(type,constraint_n,constraints,color,temp,ghost);
      fputs(line,temp);
    }
  }
  return 0;
}


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
#if 0
#include <libgen.h>
  if (i = strrspn(s, " \t"))
    s[i] = '\0';
#else
  for (p = s + (strlen(s)-1); p >= s; p--)
  {
    if ((*p == ' ') || (*p == '\t'))
      *p = 0;
    else
      break;
  }
#endif
  
  // Remove trailing quotes.
  if ((p = strrchr(s, '\"')) != NULL)
    *p = 0;

  return(s);
}

//---------------------------------------------------------------------------
#pragma argsused
int main(int argc, char* argv[])
{
  char *dat_name = argv[1];
  char *dst_name = argv[2];
  char *synth_name = NULL;
  char filename[512];
  FILE *outfile;
  FILE *synthfile;
  FILE *dat;
  char line[512];
  char *nonwhite;
  int   rc = 0;
  int  synthcount = 0;
  int  synthonly = 0;
  int  subfiles = 0;

  printf("LSynth version 2.1 by Kevin Clague, kevin_clague@yahoo.com\n");

  if (argc == 2 && strcmp(argv[1],"-v") == 0) {
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

  // Consider using L3fgets() to avoid OS dependent line termination issues.
  while (fgets(line,sizeof(line), dat)) {

    int format;

    if (!synthonly)
      fputs(line,outfile); // Copy line to output file before uppercasing it.

    strupper(line);

    // Strip away leading whitespace (spaces and tabs).
    nonwhite = line + strspn(line, " \t");

    format = (strncmp(nonwhite,
                     "0 SYNTH BEGIN ",
                     strlen("0 SYNTH BEGIN ")) == 0) * 2;

    format += strncmp(nonwhite,
                       "0 WRITE SYNTH BEGIN ",
                       strlen("0 WRITE SYNTH BEGIN ")) == 0;

    if (format) {

      int type_index;
      char *option;

      if (format == 1) {
        nonwhite += strlen("0 WRITE SYNTH BEGIN ");
      } else {
        nonwhite += strlen("0 SYNTH BEGIN ");
      }

      synthfile = outfile; // By default synth data goes to main outfile.
      synthcount++; // Count the synth parts.

      // We may be writing synth data to separate subfiles.
      if (subfiles)
      {
	// Check if the subfile is named in the SYNTH BEGIN line.
	if (option = strstr(nonwhite, "NAME="))
	{
	  option += strlen("NAME=");
	  // Allocate memory for name, strip quotes and trailing spaces.
	  synth_name = stripquotes(option); 
	}
	// Otherwise just number the subfile.  Use lsynthN.ldr for stdout?
	else
	{
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
      if (option = strstr(nonwhite, "LENGTH="))
      {
	option += strlen("LENGTH=");
      }
      if (option = strstr(nonwhite, "UNITS="))
      {
	option += strlen("UNITS=");
      }

      /* check to see if it is a known synth command */

      for (type_index = 0; tube_types[type_index][0]; type_index++) {
        if (strncmp(nonwhite,
              tube_types[type_index],
              strlen(tube_types[type_index])) == 0) {

          rc = synth_tube_like(
                 nonwhite,
                 nonwhite + strlen(tube_types[type_index]),
                 dat,synthfile);
          break;
        }
      }

      /* check to see if it is a known synth command */

      for (type_index = 0; band_types[type_index][0]; type_index++) {
        if (strncmp(nonwhite,
              band_types[type_index],
              strlen(band_types[type_index])) == 0) {

          rc = synth_band_like(
                 nonwhite,
                 nonwhite + strlen(band_types[type_index]),
                 dat,synthfile);
          break;
	}
      }

      // Close subfile and clean up.
      if (synth_name)
      {
	fclose(synthfile);
	free(synth_name);
	synth_name = NULL;
      }
    }
  }
  fclose(dat);
  fclose(outfile);

  printf("LSynth complete\n");
  return 1;
}
//---------------------------------------------------------------------------
