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

  for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
    if (*nonwhite != ' ') {
      break;
    }
  }

  while (strncmp(nonwhite,"0 ROTATION C",strlen("0 ROTATION C")) == 0) {
    fputs(line,temp);
    fgets(line,n_line,dat);  /* FIXME: check fgets rc */
    strupper(line);
    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }
  }

  return 0;
}

int skip_synthesized(FILE *dat, char *line, int sizeof_line)
{
  int rc;
  char *nonwhite;

  if (strcmp(nonwhite,"0 SYNTHESIZED BEGIN\n") == 0) {
    int rc;

    while (fgets(line,sizeof(line),dat)) {
      strupper(line);
      for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
        if (*nonwhite != ' ') {
          break;
        }
      }
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

    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }

    skip_rot(line,sizeof(line),dat,temp);

    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }

    if (strncmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      for ( ; *nonwhite != '\0'; nonwhite++) {
        if (*nonwhite != ' ') {
          break;
        }
      }
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

    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }

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

    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }

    skip_rot(line,sizeof(line),dat,temp);

    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }

    if (strncmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      for ( ; *nonwhite != '\0'; nonwhite++) {
        if (*nonwhite != ' ') {
          break;
        }
      }
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

    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }

    if (strcmp(nonwhite,"0 SYNTH END\n") == 0 ||
        strcmp(nonwhite,"0 WRITE SYNTH END\n") == 0) {

      rc = synth_band(type,constraint_n,constraints,color,temp,ghost);
      fputs(line,temp);
    }
  }
  return 0;
}


#pragma argsused
int main(int argc, char* argv[])
{
  char *dat_name = argv[1];
  char *dst_name = argv[2];
  FILE *temp;
  FILE *dat;
  char line[512];
  char *nonwhite;
  int   rc = 0;

  printf("LSynth version 2.0 by Kevin Clague, kevin_clague@yahoo.com\n");

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

  temp = fopen(dst_name,"w");

  if (temp == NULL) {
    printf("%s: Failed to open file %s for writing\n",argv[0],dst_name);
    return -1;
  }

  while (fgets(line,sizeof(line), dat)) {

    int format;

    strupper(line);

    for (nonwhite = line; *nonwhite != '\0'; nonwhite++) {
      if (*nonwhite != ' ') {
        break;
      }
    }

    fputs(line,temp);

    format = (strncmp(nonwhite,
                     "0 SYNTH BEGIN ",
                     strlen("0 SYNTH BEGIN ")) == 0) * 2;

    format += strncmp(nonwhite,
                       "0 WRITE SYNTH BEGIN ",
                       strlen("0 WRITE SYNTH BEGIN ")) == 0;;

    if (format) {

      int type_index;

      if (format == 1) {
        nonwhite += strlen("0 WRITE SYNTH BEGIN ");
      } else {
        nonwhite += strlen("0 SYNTH BEGIN ");
      }

      /* check to see if it is a known synth command */

      for (type_index = 0; tube_types[type_index][0]; type_index++) {
        if (strncmp(nonwhite,
              tube_types[type_index],
              strlen(tube_types[type_index])) == 0) {

          rc = synth_tube_like(
                 nonwhite,
                 nonwhite + strlen(tube_types[type_index]),
                 dat,temp);
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
                 dat,temp);
          break;
        }
      }
    }
  }
  fclose(dat);
  fclose(temp);

  printf("LSynth complete\n");
  return 1;
}
//---------------------------------------------------------------------------
