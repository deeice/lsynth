/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague and Don Heyse
 */

#include "lsynthcp.h"
#include "hose.h"
#include "curve.h"
#include "math.h"
#include "mathlib.h"

#define PI 2*atan2(1,0)

/*
 * 0 SYNTH BEGIN DEFINE HOSE <fill> PNEUMATIC_HOSE "descr" <diameter> <stiffness> <twist>
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 0 SYNTH END
 */

int           n_hose_types = 0;
hose_attrib_t   hose_types[32];

#define N_HOSE_TYPES n_hose_types

/* In hoses, the attrib field in constraints, indicates that
 * LSynth should turn the final constraint around to get everything
 * to work correctly (e.g. flex-axle ends).
 */

/*
 * 0 SYNTH BEGIN DEFINE HOSE CONSTRAINTS
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 0 SYNTH END
 */

int    n_hose_constraints = 0;
part_t   hose_constraints[64];

#define N_HOSE_CONSTRAINTS n_hose_constraints

void
list_hose_types(void)
{
  int i;

  printf("\n\nHose like synthesizable parts\n");
  for (i = 0; i < N_HOSE_TYPES; i++) {
    printf("  %-20s %s\n",hose_types[i].type, hose_types[i].descr);
  }
}

void
list_hose_constraints(void)
{
  int i;

  printf("\n\nHose constraints\n");
  for (i = 0; i < N_HOSE_CONSTRAINTS; i++) {
    printf("    %11s\n",hose_constraints[i].type);
  }
}

void
hose_ini(void)
{
  int i;

  for (i = 0; i < N_HOSE_TYPES; i++) {
    printf("%-20s = SYNTH BEGIN %s 16\n",hose_types[i].type, hose_types[i].type);
  }
}

int
ishosetype(char *type)
{
  int i;

  for (i = 0; i < N_HOSE_TYPES; i++) {
    if (strncmp(hose_types[i].type,type,strlen(hose_types[i].type)) == 0) {
      return 1;
    }
  }
  return 0;
}

int
ishoseconstraint(char *type)
{
  int i;

  for (i = 0; i < N_HOSE_CONSTRAINTS;i++) {
    if (strcmp(hose_constraints[i].type,type) == 0) {
      return 1;
    }
  }
  return 0;
}

/*
 * a 1x1 brick is 20 LDU wide and 24 LDU high
 *
 * hose length = 14 brick widths long = 280 LDU
 * number of ribs = 45
 * 6.2 LDU per rib
 *
 */

void orient( int n_segments, part_t *segments);

PRECISION
line_angle(
  int a,
  int b,
  part_t *segments)
{
  PRECISION va[3]; /* line A */
  PRECISION vb[3]; /* line B */
  PRECISION denom;
  PRECISION theta;

  /* we get the line A as point[a] to point[a+1] */
  /* and line B as point[b] to point[b+1] */

  vectorsub3(va,segments[a+1].offset,segments[a].offset);
  vectorsub3(vb,segments[b+1].offset,segments[b].offset);

  denom = vectorlen(va)*vectorlen(vb);

  theta = (va[0]*vb[0]+va[1]*vb[1]+va[2]*vb[2])/denom;
  if (theta >= 1 || theta <= -1) {
    theta = 0;
  } else {
    theta = acos(theta);
  }
  return theta;
}

int
merge_segments_angular(
  part_t    *segments,
  int       *n_segments,
  PRECISION  max,
  FILE      *output)
{
  int a,b;
  int n;
  PRECISION theta;

  a = 0; b = 1;
  n = 0;

  while (b < *n_segments) {
    theta = line_angle(a,b,segments);
    if (theta < max) {
      b++;
    } else {
      segments[n++] = segments[a++];
      a = b++;
    }
  }
  if (n < 3) {
    segments[1] = segments[*n_segments-2];
    segments[2] = segments[*n_segments-1];
    n = 3;
  } else {
    segments[n++] = segments[b-1];
  }
  *n_segments = n;
  orient(n,segments);
  return 0;
}

int
merge_segments_length(
  part_t    *segments,
  int       *n_segments,
  PRECISION  max,
  FILE      *output)
{
  int a,b;
  int n;
  PRECISION d[3],l;

  a = 0; b = 1;
  n = 1;

  while (b < *n_segments) {

    vectorsub3(d,segments[a].offset,segments[b].offset);
    l = vectorlen(d);

    if (l + 0.5 < max) {
      b++;
    } else {
      a = b;
      segments[n++] = segments[b++];
    }
  }
  segments[n++] = segments[b-1];
  *n_segments = n;
  orient(n,segments);
  return 0;
}

#define MAX_SEGMENTS 1024*8

part_t segments[MAX_SEGMENTS];

void
output_line(
  FILE           *output,
  int             ghost,
  char           *group,
  int             color,
  PRECISION       a,
  PRECISION       b,
  PRECISION       c,
  PRECISION       d,
  PRECISION       e,
  PRECISION       f,
  PRECISION       g,
  PRECISION       h,
  PRECISION       i,
  PRECISION       j,
  PRECISION       k,
  PRECISION       l,
  char            *type)
{
  if (group) {
    fprintf(output,"0 MLCAD BTG %s\n",group);
  }
  fprintf(output,"%s1 %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %s\n",
      ghost ? "0 GHOST " : "",
      color,
      a,b,c,d,e,f,g,h,i,j,k,l,
      type);
  group_size++;
}

/*
 * Twist
 *    cos(t) 0 sin(t)
 *         0 1 0
 *   -sin(t) 0 cos(t)
 */

void
render_hose_segment(
  hose_attrib_t  *hose,
  int             ghost,
  char           *group,
  int             color,
  part_t         *segments,
  int             n_segments,
  PRECISION      *total_twist,
  int             first,
  int             last,
  FILE           *output)
{
  int i,j,k;
  PRECISION pi = 2*atan2(1,0);
  PRECISION scale[3][3];
  PRECISION scaled[3][3];
  PRECISION reversed[3][3];
  PRECISION twist[3][3];
  PRECISION twisted[3][3];
  char     *type;

  scale[0][0] = 1;
  scale[0][1] = 0;
  scale[0][2] = 0;
  scale[1][0] = 0;
  scale[1][2] = 0;
  scale[2][0] = 0;
  scale[2][1] = 0;
  scale[2][2] = 1;

  twist[0][0] = 1;
  twist[0][1] = 0;
  twist[0][2] = 0;
  twist[1][0] = 0;
  twist[1][1] = 1;
  twist[1][2] = 0;
  twist[2][0] = 0;
  twist[2][1] = 0;
  twist[2][2] = 1;

  for (i = 0; i <  n_segments-1; i++) {
    PRECISION tx,ty,tz,l,theta;

    if (hose->fill == FIXED) {
      l = 1;
      *total_twist += hose->twist;
    } else {
      PRECISION d[3];

      vectorsub3(d,segments[i+1].offset,segments[i].offset);

      l = vectorlen(d);

      if (n_segments > 2) {
        theta = line_angle(i+1,i,segments);
        theta = sin(theta)*(hose->diameter/2);

        l += theta;
      }
    }

    scale[1][1] = l;

    if (i == 0 && first) {
      type = hose->start.type;
      matrixmult3(scaled,hose->start.orient,scale);
    } else if (i == n_segments-2 && last) {
      PRECISION offset[3];
      type = hose->end.type;
      offset[0] = 0; offset[1] = l; offset[2] = 0;
      vectorrot(offset,segments[i].orient);
      vectoradd(segments[i].offset,offset);
      matrixmult3(scaled,hose->end.orient,scale);
    } else {
      type = hose->mid.type;
      matrixmult3(scaled,hose->mid.orient,scale);
    }

    if (hose->fill == FIXED) {
      PRECISION angle = *total_twist * pi / 360;
      twist[0][0] =   cos(angle);
      twist[0][2] =   sin(angle);
      twist[2][0] =  -sin(angle);
      twist[2][2] =   cos(angle);
      *total_twist += hose->twist;
    }
    matrixmult3(twisted,twist,scaled);
    matrixmult3(reversed,segments[i].orient,twisted);

    output_line(output,ghost,group,color,
      segments[i].offset[0], segments[i].offset[1], segments[i].offset[2],
      reversed[0][0], reversed[0][1], reversed[0][2],
      reversed[1][0], reversed[1][1], reversed[1][2],
      reversed[2][0], reversed[2][1], reversed[2][2],
      type);
  }
}

void
adjust_constraint(
  part_t *part,
  part_t *orig,
  int     last)
{
  int i;
  PRECISION m[3][3];

  *part = *orig;

  for (i = 0; i < N_HOSE_CONSTRAINTS; i++) {
    if (strcmp(part->type,hose_constraints[i].type) == 0) {

      vectorsub(part->offset,orig->offset);
      matrixinv(m,orig->orient);
      matrixmult(part->orient,m);

      vectoradd(part->offset,hose_constraints[i].offset);
      vectorrot(part->offset,hose_constraints[i].orient);

      matrixmult(part->orient,hose_constraints[i].orient);

      matrixmult(part->orient,orig->orient);

      vectoradd(part->offset,orig->offset);

      if (hose_constraints[i].attrib && last) {
        matrixcp(m,part->orient);
        matrixneg(part->orient,m);
      }
      break;
    }
  }
}

void
render_hose(
  hose_attrib_t  *hose,
  int             n_constraints,
  part_t         *constraints,
  PRECISION       angular_res,
  int             ghost,
  char           *group,
  int             color,
  FILE *output)
{
  int       c, n_segments;
  part_t    mid_constraint;
  PRECISION total_twist = 0;

  group_size = 0;

  fprintf(output,"0 SYNTH SYNTHESIZED BEGIN\n");

  mid_constraint = constraints[0];

  for (c = 0; c < n_constraints - 1; c++) {
    part_t first,second;

    // reorient imperfectly oriented or displaced constraint types

    adjust_constraint(&first,&mid_constraint,0);

    // reorient imperflectly oriented or displaced constraint types

    adjust_constraint(&second, &constraints[c+1],c == n_constraints-2);

    n_segments = MAX_SEGMENTS;

    // create an oversampled curve

    synth_curve(&first,&second,segments,n_segments,hose->stiffness,output);

    // Make sure final segment matches second constraint

    vectorcp(segments[n_segments-1].offset,second.offset);

    // reduce oversampled curve to fixed length chunks, or segments limit
    // by angular resolution

    if (hose->fill == FIXED) {
      merge_segments_length(segments,&n_segments,hose->mid.attrib,output);
    } else {
      merge_segments_angular(segments,&n_segments,angular_res, output);
    }

    mid_constraint = constraints[c+1];

    if (hose->fill == FIXED) {
      vectorcp(mid_constraint.offset,segments[n_segments-2].offset);
    }

    // output the result

    render_hose_segment(
      hose,
      ghost,
      group,
      color,
      segments,n_segments,
      &total_twist,
      c ==0,
      c == n_constraints-2,
      output);
  }

  if (group) {
    fprintf(output,"0 GROUP %d %s\n",group_size,group);
  }

  fprintf(output,"0 SYNTH SYNTHESIZED END\n");
  printf("Synthesized %s\n",hose->type);
}

int
synth_hose(
  char   *type,
  int     n_constraints,
  part_t *constraints,
  int     ghost,
  char   *group,
  int     color,
  FILE   *output)
{
  int i;

  for (i = 0; i < N_HOSE_TYPES; i++) {
    if (strcmp(hose_types[i].type,type) == 0) {
      render_hose(
        &hose_types[i],
        n_constraints,constraints,
        hose_res_angle,
        ghost,
        group,
        color,
        output);
      return 0;
    }
  }
  return 1;
}

