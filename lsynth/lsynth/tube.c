/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynth.h"
#include "tube.h"
#include "curve.h"
#include "math.h"

#define PI 2*atan2(1,0)

typedef struct {
  char     *part;          // LDraw part name
  PRECISION orient[3][3];
  PRECISION length;        // For FIXED synthesis
} part_t;

#define FIXED       0
#define STRETCH     1

typedef struct {
  char     *type;
  int       fill;         // FIXED or STRETCH
  PRECISION diameter;     // of cross section
  PRECISION exit_attrib;  // stiffness
  PRECISION twist;        // in degrees
  part_t    end;          // LDraw part for end of hose
  part_t    mid;          // LDraw part for mis section
} hose_attrib_t;

hose_attrib_t hose_type[] = {
  {
    "RIBBED_HOSE",
    FIXED,            // fixed length segments
    1,                // diameter
    90,               // enter and exit
    0,                // twist
    {
      "79.DAT",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "80.DAT",     // mid
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      },
      6.6,
    }
  },
  {
    "RUBBER_HOSE",
    FIXED,            // fixed length segments
    1,                // diameter
    90,               // enter and exit
    0,                // twist
    {
      "755.DAT",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "756.DAT",     // mid
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      },
      4.5,
    }
  },
  {
    "STRING",
    FIXED,            // fixed length segments
    2.0,              // diameter
    40,               // enter and exit
    10,               // twist
    {
      "STRING.DAT",   // end
      {
        {1, 0,0 },
        {0, 0,1 },
        {0,-1,0 }
      }
    },{
      "STRING.DAT",   // mid
      {
        { 1, 0,0 },
        { 0, 0,1 },
        { 0,-1,0 }
      },
      2.0
    }
  },
  {
    "MINIFIG_CHAIN",
    FIXED,            // fixed length segments
    2.0,              // diameter
    60,               // enter and exit
    90,               // twist
    {
      "209.DAT",   // end
      {
        {1, 0, 0 },
        {0, 1, 0 },
        {0, 0, 1 }
      },
    },{
      "209.DAT",   // mid
      {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 }
      },
      12.0,
    }
  },
  {
    "PNEUMATIC_HOSE",
    STRETCH,          // stretchable segments
    6.0,              // diameter
    40,               // enter and exit
    0,                // twist
    {
      "LS20.DAT",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "LS21.DAT",     // mid
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    }
  },
  {
    "ELECTRIC_CABLE",
    STRETCH,          // stretchable segments
    6.0,              // diameter
    40,               // enter and exit
    0,                // twist
    {
      "LS10.DAT",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "LS10.DAT",     // mid
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    }
  },
  {
    "FLEX_SYSTEM_HOSE_LD",
    STRETCH,          // stretchable segments
    6.0,              // diameter
    40,               // enter and exit
    0,                // twist
    {
      "76.DAT",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "4-4cyli.dat",// mid
      {
        { 4,0,0 },
        { 0,1,0 },
        { 0,0,4 }
      }
    }
  },
  {
    "FLEX_SYSTEM_HOSE",
    STRETCH,          // stretchable segments
    6.0,              // diameter
    40,               // enter and exit
    0,                // twist
    {
      "76.DAT",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "4-4cyli.dat",     // mid
      {
        { 4,0,0 },
        { 0,1,0 },
        { 0,0,4 }
      }
    }
  },
  {
    "FLEX_SYSTEM_CABLE",
    STRETCH,          // stretchable segments
    6.0,              // diameter
    40,               // enter and exit
    0,                // twist
    {
      "4-4cyli.dat",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "4-4cyli.dat",     // mid
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    }
  },
  {
    "FLEXIBLE_AXLE",
    STRETCH,          // stretchable segments
    6.0,              // diameter
    70,               // enter and exit
    0,                // twist
    {
      "AXLE_1.DAT",     // end
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    },{
      "AXLE_1.DAT",     // mid
      {
        { 1,0,0 },
        { 0,1,0 },
        { 0,0,1 }
      }
    }
  },
  {
    "FIBER_OPTIC_CABLE",
    STRETCH,          // stretchable segments
    6.0,              // diameter
    40,               // enter and exit
    0,                // twist
    {
      "LS32.DAT",     // end
      {
        {1,0,0 },
        {0,1,0 },
        {0,0,1 }
      }
    },{
      "4-4cyli.DAT",     // mid
      {
        { 2,0,0 },
        { 0,1,0 },
        { 0,0,2 }
      }
    }
  },
};

/*
 * a 1x1 brick is 20 LDU wide and 24 LDU high
 *
 * hose length = 14 brick widths long = 280 LDU
 * number of ribs = 45
 * 6.2 LDU per rib
 *
 * ToDo:
 *   Examine segments , gathering them until they exceed
 *   error angle.  Back up one segment, then treat gathered
 *   segments as one segment with cross section scaled for
 *   length.
 */

void orient( int n_segments, LSL_part_usage *segments);

PRECISION
line_angle(
  int a,
  int b,
  LSL_part_usage *segments)
{
  PRECISION c,d;
  PRECISION f1,g1,h1; /* line A */
  PRECISION f2,g2,h2; /* line B */
  PRECISION denom;
  PRECISION theta;

  /* we get the line A as point[a] to point[a+1] */
  /* and line B as point[b] to point[b+1] */

  f1 = segments[a+1].loc.x - segments[a].loc.x;
  g1 = segments[a+1].loc.y - segments[a].loc.y;
  h1 = segments[a+1].loc.z - segments[a].loc.z;

  f2 = segments[b+1].loc.x - segments[b].loc.x;
  g2 = segments[b+1].loc.y - segments[b].loc.y;
  h2 = segments[b+1].loc.z - segments[b].loc.z;

  denom = sqrt((f1*f1+g1*g1+h1*h1)*(f2*f2+g2*g2+h2*h2));
  theta = acos((f1*f2+g1*g2+h1*h2)/denom);
  return theta;
}

int
merge_segments_angular(
  LSL_part_usage       *segments,
  int                  *n_segments,
  PRECISION             max,
  FILE                 *output)
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
  segments[n++] = segments[b-1];
  *n_segments = n;
  orient(n,segments);
  return 0;
}

int
merge_segments_length(
  LSL_part_usage       *segments,
  int                  *n_segments,
  PRECISION             max,
  FILE                 *output)
{
  int a,b;
  int n;
  PRECISION dx,dy,dz,l;

  a = 0; b = 1;
  n = 1;

  while (b < *n_segments) {

    dx = segments[a].loc.x - segments[b].loc.x;
    dy = segments[a].loc.y - segments[b].loc.y;
    dz = segments[a].loc.z - segments[b].loc.z;
    l = sqrt(dx*dx+dy*dy+dz*dz);

    if (l < max) {
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

static PRECISION constr_len(
  LSL_part_usage *start,
  LSL_part_usage *end,
  LSL_part_usage *segments,
  int             n_segments,
  PRECISION       attrib,
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

void
output_line(
  FILE           *output,
  int             ghost,
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
  fprintf(output,"%s1 %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %s\n",
      ghost ? "0 GHOST " : "",
      color,
      a,b,c,d,e,f,g,h,i,j,k,l,
      type);
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
  int             color,
  LSL_part_usage *segments,
  int             n_segments,
  PRECISION      *total_twist,
  FILE           *output)
{
  int i,j,k;
  PRECISION pi = 2*atan2(1,0);
  PRECISION scale[3][3];
  PRECISION scaled[3][3];
  PRECISION reversed[3][3];
  PRECISION twist[3][3];
  PRECISION twisted[3][3];

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
    PRECISION tx,ty,tz,l;

    if (hose->fill == FIXED) {
      l = 1;
      *total_twist += hose->twist;
    } else {

      tx = segments[i+1].loc.x - segments[i].loc.x;
      ty = segments[i+1].loc.y - segments[i].loc.y;
      tz = segments[i+1].loc.z - segments[i].loc.z;

      l = sqrt(tx*tx+ty*ty+tz*tz);
    }

    scale[1][1] = l;

    mat_mult(scaled,hose->mid.orient,scale);

    if (hose->fill == FIXED) {
      PRECISION angle = *total_twist * pi / 360;
      twist[0][0] =   cos(angle);
      twist[0][2] =   sin(angle);
      twist[2][0] =  -sin(angle);
      twist[2][2] =   cos(angle);
      *total_twist += hose->twist;
    }
    mat_mult(twisted,twist,scaled);

    mat_mult(reversed,segments[i].orient,twisted);
    output_line(output,ghost,color,
      segments[i].loc.x, segments[i].loc.y, segments[i].loc.z,
      reversed[0][0], reversed[0][1], reversed[0][2],
      reversed[1][0], reversed[1][1], reversed[1][2],
      reversed[2][0], reversed[2][1], reversed[2][2],
      hose->mid.part);
  }
}

void
render_hose(
  hose_attrib_t  *hose,
  int             n_constraints,
  LSL_part_usage *constraints,
  PRECISION       angular_res,
  int             ghost,
  int             color,
  FILE *output)
{
  int c, n_segments;
  PRECISION pi = 2*atan2(1,0);
  PRECISION reverse[3][3];
  PRECISION reversed[3][3];
  LSL_part_usage mid_constraint;
  PRECISION total_twist = 0;

  reverse[0][0] = 1;
  reverse[0][1] = 0;
  reverse[0][2] = 0;
  reverse[1][0] = 0;
  reverse[1][1] = cos(pi);
  reverse[1][2] = -sin(pi);
  reverse[2][0] = 0;
  reverse[2][1] = sin(pi);
  reverse[2][2] = cos(pi);

  fprintf(output,"0 SYNTH SYNTHESIZED BEGIN\n");

  mid_constraint = constraints[0];

  for (c = 0; c < n_constraints - 1; c++) {

    n_segments = MAX_SEGMENTS;

    synth_curve(&mid_constraint,&constraints[c+1],
                segments,n_segments,hose->exit_attrib,output);

    segments[n_segments-1].loc = constraints[c+1].loc;

    if (hose->fill == FIXED) {
      merge_segments_length(segments,&n_segments,hose->mid.length,output);
    } else {
      merge_segments_angular(segments,&n_segments,angular_res, output);
    }
    mid_constraint = constraints[c+1];
    mid_constraint.loc = segments[n_segments-2].loc;

    if (c == 0) {
      mat_mult(reversed,segments[0].orient,reverse);
      output_line(
        output, ghost, color,
        segments[0].loc.x,segments[0].loc.y,segments[0].loc.z,
        reversed[0][0],reversed[0][1],reversed[0][2],
        reversed[1][0],reversed[1][1],reversed[1][2],
        reversed[2][0],reversed[2][1],reversed[2][2],
        hose->end.part);

      total_twist += hose->twist;

      render_hose_segment(hose,ghost,color,segments+1,n_segments-2,&total_twist,output);
    } else if (c == n_constraints-2) {
      int i = n_segments - 2;

      render_hose_segment(hose,ghost,color,segments,n_segments-2,&total_twist,output);

      output_line(output, ghost, color,
        segments[i].loc.x,segments[i].loc.y,segments[i].loc.z,
        segments[i].orient[0][0],
        segments[i].orient[0][1],
        segments[i].orient[0][2],
        segments[i].orient[1][0],
        segments[i].orient[1][1],
        segments[i].orient[1][2],
        segments[i].orient[2][0],
        segments[i].orient[2][1],
        segments[i].orient[2][2],
        hose->end.part);
    } else {
      render_hose_segment(hose,ghost,color,segments,n_segments-1,&total_twist,output);
    }
  }

  fprintf(output,"0 SYNTH SYNTHESIZED END\n");
  printf("Synthesized %s\n",hose->type);
}

int
synth_tube(
  char                 *type,
  int                   n_constraints,
  LSL_part_usage       *constraints,
  int                   ghost,
  int                   color,
  FILE *output)
{
  int i;

  for (i = 0; i < sizeof(hose_type)/sizeof(hose_attrib_t); i++) {
    if (strcmp(hose_type[i].type,type) == 0) {
      render_hose(
        &hose_type[i],
        n_constraints,constraints,
        hose_res_angle,
        ghost,
        color,
        output);
      return 0;
    }
  }
  return 1;
}

