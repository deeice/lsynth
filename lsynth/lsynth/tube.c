/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynth.h"
#include "tube.h"
#include "curve.h"

/*
 * a 1x1 brick is 20 LDU wide and 24 LDU high
 *
 * hose length = 14 brick widths long = 280 LDU
 * number of ribs = 45
 * 6.2 LDU per rib
 *
 */

int tube_res = 8;

#define MAX_SEGMENTS 1024*8

static PRECISION constr_len(
  LSL_part_usage *start,
  LSL_part_usage *end,
  LSL_part_usage *segments,
  int             n_segments,
  LSL_tube_attributes *attrib,
  FILE *output)
{
  int i;
  PRECISION x,y,z;
  PRECISION length;

  synth_curve(start,end,segments,n_segments,attrib,output);

  length = 0;
  for (i = 0; i < n_segments - 1; i++) {
    x = segments[i + 1].loc.x - segments[i].loc.x;
    y = segments[i + 1].loc.y - segments[i].loc.y;
    z = segments[i + 1].loc.z - segments[i].loc.z;

    length += sqrt(x*x+y*y+z*z);
  }

  return length;
}

void
rotate( PRECISION res[3], PRECISION loc[3], PRECISION rot[3][3])
{
  res[0] = rot[0][0]*loc[0] + rot[0][1]*loc[1] + rot[0][2]*loc[2];
  res[1] = rot[1][0]*loc[0] + rot[1][1]*loc[1] + rot[1][2]*loc[2];
  res[2] = rot[2][0]*loc[0] + rot[2][1]*loc[1] + rot[2][2]*loc[2];
}

void
mat_mult( PRECISION res[3][3], PRECISION lft[3][3], PRECISION rht[3][3])
{
  int i,j,k;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      res[i][j] = 0.0;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        res[i][j] += lft[i][k] * rht[k][j];
      }
    }
  }
}

LSL_part_usage segments[MAX_SEGMENTS];

int
synth_tube(
  char                 *type,
  int                   n_constraints,
  LSL_part_usage       *constraints,
  int                   length, /* in ribs */
  char                 *units,
  int                   color,
  FILE *output,
  int                   ghost)
{
  LSL_tube_attributes attrib;
  int c, i, n_segments, first, last;
  PRECISION real_length;
  int segment;
  PRECISION pi = 2*atan2(1,0);
  PRECISION reverse[3][3];
  PRECISION t[3][3];
  PRECISION reversed[3][3];

  if (strcmp(type,RIBBED_HOSE)) {

    attrib.exit_start = 50;
    attrib.exit_stop  = 50;
  } else {

    attrib.exit_start = 40;
    attrib.exit_stop  = 40;
  }

  reverse[0][0] = 1;
  reverse[0][1] = 0;
  reverse[0][2] = 0;
  reverse[1][0] = 0;
  reverse[1][1] = cos(pi);
  reverse[1][2] = -sin(pi);
  reverse[2][0] = 0;
  reverse[2][1] = sin(pi);
  reverse[2][2] = cos(pi);

  // determine the total minimum length

  segment = 0;

  fprintf(output,"0 SYNTH SYNTHESIZED BEGIN\n");

  for (c = 0; c < n_constraints - 1; c++) {

    real_length = constr_len(&constraints[c],&constraints[c+1], segments, MAX_SEGMENTS, &attrib, output);

    if (n_segments > MAX_SEGMENTS) {
      n_segments = MAX_SEGMENTS;
    }

    if (strncmp(type,RIBBED_HOSE,strlen(RIBBED_HOSE)) == 0 ||
        strncmp(type,RIBBED_TUBE,strlen(RIBBED_TUBE)) == 0) {

      n_segments = real_length /= 6.6;
      n_segments *= 2;

      if (n_segments > MAX_SEGMENTS) {
        n_segments = MAX_SEGMENTS;
      }

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      if (c == 0) {
        first = 1;
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[0].loc.x,
          segments[0].loc.y,
          segments[0].loc.z,
          -segments[0].orient[0][0],
          -segments[0].orient[0][1],
          -segments[0].orient[0][2],
          -segments[0].orient[1][0],
          -segments[0].orient[1][1],
          -segments[0].orient[1][2],
          -segments[0].orient[2][0],
          -segments[0].orient[2][1],
          -segments[0].orient[2][2],
          "79.DAT");
        segment++;
      } else {
        first = 0;
      }

      if (c == n_constraints-2) {
        last = n_segments-2;
      } else {
        last = n_segments - 1;
      }

      for (i = first; i <  last; i++) {
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "80.DAT");
        segment++;
      }

      if (c == n_constraints-2) {
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "79.DAT");
        segment++;
      }
    } else if (strncmp(type,TEST_HOSE,strlen(TEST_HOSE)) == 0) {

      n_segments = real_length /= 17;

      if (n_segments > MAX_SEGMENTS) {
        n_segments = MAX_SEGMENTS;
      }

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      for (i = 0; i <  n_segments-1; i++) {

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "LS00.DAT");
        segment++;
      }
    } else if (strncmp(type,ELECTRIC_CABLE,strlen(ELECTRIC_CABLE)) == 0) {

      n_segments = real_length*2.0;

      if (n_segments > MAX_SEGMENTS) {
        n_segments = MAX_SEGMENTS;
      }

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      for (i = 0; i <  n_segments-1; i++) {

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "LS10.DAT");
        segment++;
      }

    } else if (strncmp(type,PNEUMATIC_HOSE,strlen(PNEUMATIC_HOSE)) == 0 ||
               strncmp(type,PNEUMATIC_TUBE,strlen(PNEUMATIC_TUBE)) == 0) {

      n_segments = real_length*2.0;

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      if (c == 0) {
        first = 1;
        mat_mult(reversed,segments[0].orient,reverse);
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[0].loc.x,
          segments[0].loc.y,
          segments[0].loc.z,
          reversed[0][0],
          reversed[0][1],
          reversed[0][2],
          reversed[1][0],
          reversed[1][1],
          reversed[1][2],
          reversed[2][0],
          reversed[2][1],
          reversed[2][2],
          "LS20.DAT");
        segment++;
      } else {
        first = 0;
      }

      if (c == n_constraints-2) {
        last = n_segments-2;
      } else {
        last = n_segments - 1;
      }

      for (i = first; i <  last; i++) {

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "LS21.DAT");
        segment++;
      }

      if (c == n_constraints-2) {
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "LS20.DAT");
        segment++;
      }
    } else if (strncmp(type,FLEXIBLE_HOSE,strlen(FLEXIBLE_HOSE)) == 0 ||
               strncmp(type,FLEXIBLE_TUBE,strlen(FLEXIBLE_TUBE)) == 0 ||
               strncmp(type,RIGID_HOSE,strlen(RIGID_HOSE)) == 0 ||
               strncmp(type,RIGID_TUBE,strlen(RIGID_TUBE)) == 0) {

      n_segments = real_length*1.4;

      if (n_segments > MAX_SEGMENTS) {
        n_segments = MAX_SEGMENTS;
      }

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      if (c == 0) {
        first = 1;
        mat_mult(reversed,segments[0].orient,reverse);
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[0].loc.x,
          segments[0].loc.y,
          segments[0].loc.z,
          reversed[0][0],
          reversed[0][1],
          reversed[0][2],
          reversed[1][0],
          reversed[1][1],
          reversed[1][2],
          reversed[2][0],
          reversed[2][1],
          reversed[2][2],
          "76.DAT");
        segment++;
      } else {
        first = 0;
      }

      if (c == n_constraints-2) {
        last = n_segments-2;
      } else {
        last = n_segments - 1;
      }

      for (i = first; i <  last; i++) {

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "77.DAT");
        segment++;
      }

      if (c == n_constraints-2) {
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "76.DAT");
        segment++;
      }

    } else if (strncmp(type,FLEX_CABLE,strlen(FLEX_CABLE)) == 0) {

      n_segments = real_length*3;

      if (n_segments > MAX_SEGMENTS) {
        n_segments = MAX_SEGMENTS;
      }

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      for (i = 0; i <  n_segments-1; i++) {

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "LS51.DAT");
        segment++;
      }

    } else if (strncmp(type,FLEXIBLE_AXLE,strlen(FLEXIBLE_AXLE)) == 0) {

      n_segments = real_length*1.5;

      if (n_segments > MAX_SEGMENTS) {
        n_segments = MAX_SEGMENTS;
      }

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      for (i = 0; i <  n_segments-1; i++) {

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "AXLE.DAT");
        segment++;
      }

    } else if (strncmp(type,FIBER_OPTIC_CABLE,strlen(FIBER_OPTIC_CABLE)) == 0 ||
               strncmp(type,FIBRE_OPTIC_CABLE,strlen(FIBRE_OPTIC_CABLE)) == 0 ||
               strncmp(type,RUBBER_BAND,strlen(RUBBER_BAND)) == 0) {

      PRECISION dx,dy,dz;

      synth_curve(&constraints[c],&constraints[c+1],segments,n_segments,&attrib,output);

      if (c == n_constraints-2) {
        last = n_segments-2;
      } else {
        last = n_segments - 1;
      }

      for (i = 0; i <  last; i++) {

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "LS31.DAT");
        segment++;
      }

      if (c == n_constraints-2) {
        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
          color,
          segments[i].loc.x,
          segments[i].loc.y,
          segments[i].loc.z,
          segments[i].orient[0][0],
          segments[i].orient[0][1],
          segments[i].orient[0][2],
          segments[i].orient[1][0],
          segments[i].orient[1][1],
          segments[i].orient[1][2],
          segments[i].orient[2][0],
          segments[i].orient[2][1],
          segments[i].orient[2][2],
          "LS32.DAT");
        segment++;
      }
    }
  }

  fprintf(output,"0 SYNTH SYNTHESIZED END\n");
  printf("Synthesized %s\n",type);
  return 0;
}

