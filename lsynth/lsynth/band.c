/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague and Don Heyse
 */

#include "lsynthcp.h"
#include "band.h"
#include "hose.h"
#include <float.h>

/*
 * 0 SYNTH BEGIN DEFINE BAND <fill> RUBBER_BAND "Descr" <scale> <thresh>
 * 1 <len>  a b c  d e f  g h i  j k l "name"
 * 1 <len>  a b c  d e f  g h i  j k l "name"
 * 0 SYNTH END
 *
 * 0 SYNTH BEGIN DEFINE BAND CONSTRAINTS
 * 1 <dia>  a b c  d e f  g h i  j k l  "name"
 * 0 SYNTH END
 */

#define STRETCH 0
#define FIXED   1
#define FIXED3  2

typedef struct {
  char     *type;        // name of thing being specified (e.g. RUBBER_BAND)
  char     *descr;
  int       fill;        // method for synthesizing
                         // FIXED - chain and tread are composed of fixed
                         //         length parts
                         // FIXED3 - special case for rubber tread.
                         // STRETCH - for rubber bands
  PRECISION scale;       // convert LDUs to number of parts
  PRECISION delta;       // used to compute number of steps on arc
  part_t    tangent;     // the part type used for tangent
  part_t    arc;         // the part type used for going around constraints
  part_t    start_trans; // for rubber treads, transition from arc to tangent
  part_t    end_trans;   // for rubber treads, transition from tangent to arc
} band_attrib_t;

/*
 * TYPE: <name> <scale> <delta> <mode>
 * 1 <disp> <rot> <type
 */
band_attrib_t band_types[] =
{
  {
    "RUBBER_BAND","Technic rubber band (with circular cross section)",
    STRETCH,
    2.0,
    1,
    {
      "4-4CYLI.DAT",
      {
        { 2,0,0 },
        { 0,1,0 },
        { 0,0,2 }
      }
    },
    {
      "4-4CYLI.DAT",
      {
        { 2,0,0 },
        { 0,1,0 },
        { 0,0,2 }
      }
    },
  },
  {
    "RUBBER_BELT","Technic rubber belt (with square cross section)",
    STRETCH,
    2.0,
    1,
    {
      "box4o8a.DAT",
      {
        { 1.414,0,1.414 },
        {     0,1,0 },
        {-1.414,0,1.414 }
      }
    },
    {
      "box4o8a.DAT",
      {
        { 1.414,0,1.414 },
        {     0,1,0 },
        {-1.414,0,1.414 }
      }
    },
  },
  {
    "CHAIN","Technic chain composed of 3711.dat chain links",
    FIXED,
    1.0/16,
    8,
    { // tangent
      "3711.DAT",
      {
        { 0, 1, 0 },
        { 0, 0,-1 },
        {-1, 0, 0 }
      },
      {
        0, 0, 16
      }
    },
    { // arc
      "3711.DAT",
      {
        { 0, -1, 0 },
        { 0,  0, 1 },
        {-1,  0, 0 }
      },
      {
        0, 0, 0
      }
    }
  },
  {
    "PLASTIC_TREAD","Technic plastic tread composed of 3873.DAT tread links",
    FIXED,
    1.0/16,
    8,
    { // tangent
      "3873.DAT",
      {
        {  0, 1, 0 },
        {  0, 0,-1 },
        { -1, 0, 0 }
      },
      {
        0, 0, 16
      }
    },
    { // arc
      "3873.DAT",
      {
        { 0, -1, 0 },
        { 0,  0, 1 },
        { -1, 0, 0 }
      },
      {
        0, 0, 0
      }
    }
  },
  {
    "RUBBER_TREAD","Technic rubber tread composed of 680.DAT, 681.DAT and 682.DAT",
    FIXED3,
    1.0/20,
    1,
    {
      "680.DAT",
      {
        { 0, -1, 0 },
        { 1,  0, 0 },
        { 0,  0, 1 }
      },
      {
        0, 32, 0
      }
    },
    {
      "682.DAT",
      {
        { -0.342,  0.940, 0 },
        { -0.940, -0.342, 0 },
        {      0,      0, 1 }
      },
      {
        0, 32, 0
      }
    },
    {
      "681.DAT",
      {
        { -0.342,  0.940, 0 },
        { -0.940, -0.342, 0 },
        {      0,      0, 1 }
      },
      {
        0, 32, 0
      }
    },
    {
      "681.DAT",
      {
        {  0.342,   0.940, 0 },
        {  0.940,  -0.342, 0 },
        {      0,      0,  1 }
      },
      {
        0, 32, 0
      }
    },
  }
};

#define N_BAND_TYPES (sizeof(band_types)/sizeof(band_types[0]))

struct {
  char     *type;         // LDraw part name
  char     *descr;
  PRECISION radius;       // Radius of circular part
  PRECISION orient[3][3]; // How to orient it
  PRECISION offset[3];    // How much to offset it.
} constraint_type[] = {

  // bushings/pulleys

  {    "3736.DAT", "Technic Pulley Large",
                   45, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "4185.DAT", "Technic Wedge Belt Wheel",
                   29, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "6539.DAT", "Technic Transmission Driving Ring",
                   15, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "2983.DAT", "Electric Technic Micromotor Pulley",
                   13, {{1, 0, 0},{0,0, 1},{ 0,-1,0}}, { 0, 0, -2}},
  {   "4265A.DAT", "Technic Bush 1/2 Type I",
                    9, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "4265B.DAT", "Technic Bush 1/2 Type II",
                    9, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "4265C.DAT", "Technic Bush 1/2 Smooth",
                    9, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "3713.DAT", "Technic Bush",
                    9, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "2736.DAT", "Technic Axle Towball",
                    6, {{0, 0, 1},{0,1, 0},{-1,0,0}}, { 0, 0, 2}},
  {    "6628.DAT",  "Technic Friction Pin with Towball",
                    6, {{0, 0, 1},{0,1, 0},{-1,0,0}}, { 0, 0, 2}},
  {   "32007.DAT", "Technic Tread Sprocket Wheel",
                   32, {{1, 0, 0},{0,1, 0},{ 0,0,1}}, { 0, 32, 0}},
  {   "32089.DAT", "Technic Tread Sprocket Wheel Thin",
                   32, {{0,0,1},{0,1,0},{-1,0,0}}},
  // axles

  {    "3704.DAT",  "Technic Axle 2",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {   "32062.DAT",  "Technic Axle 2 Notched",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "4519.DAT",  "Technic Axle 3",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "6587.DAT",  "Technic Axle 3 with Stud",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3705.DAT",  "Technic Axle 4",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  { "3705C01.DAT",  "Technic Axle 4 Threaded",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {   "32073.DAT",  "Technic Axle 5",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3706.DAT",  "Technic Axle 6",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3707.DAT",  "Technic Axle 8",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3737.DAT",  "Technic Axle 10",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  { "3737C01.DAT",  "Technic Axle 10 Threaded",
                    8, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "3708.DAT",  "Technic Axle 12",
                    8, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  // axle joiner

  {   "6538A.DAT",  "Technic Axle Joiner",
                    9, {{0, 0, 1},{0,1, 0},{ 0,0,1}}},
  // gears

  {   "73071.DAT", "Technic Differential",
                   35, {{1, 0, 0},{0,1, 0},{ 0,0,1}}, { 0, 0,19}},
  {    "6573.DAT", "Technic Differential New",
                   30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}, { 0, 0, 30}},
  {    "4019.DAT", "Technic Gear 16 Tooth",
                   18, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "6542.DAT", "Technic Gear 16 Tooth with Clutch",
                   18, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "3648.DAT", "Technic Gear 24 Tooth",
                   30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "60C01.DAT", "Technic Gear 24 Tooth",
                   30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "3650A.DAT", "Technic Gear 24 Tooth Crown",
                   30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "3649.DAT", "Technic Gear 40 Tooth",
                   50, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "2855.DAT", "Technic Turntable Top",
                   71, {{1, 0, 0},{0,0,-1},{ 0,1,0}}, { 0, 0,-10}},
  {   "48451.DAT", "Technic Turntable New with Hole Top ",
                   71, {{1, 0, 0},{0,0,-1},{ 0,1,0}}, { 0, 0,19}},
  {    "LS00.DAT", "~LSynth Constraint",
                    1, {{1, 0, 0},{0,0,-1},{ 0,1,0}}},
};

#define N_BAND_CONSTRAINTS (sizeof(constraint_type)/sizeof(constraint_type[0]))

void band_ini(void)
{
  int i;

  for (i = 0; i < N_BAND_TYPES; i++) {
    printf("%-20s = SYNTH BEGIN %s 16\n",band_types[i].type, band_types[i].type);
  }
}

void
list_band_types(void)
{
  int i;

  printf("\n\nBand type synthesizable parts\n");
  for (i = 0; i < N_BAND_TYPES; i++) {
    printf("  %-20s %s\n",band_types[i].type, band_types[i].descr);
  }
}

int
isbandtype(char *type)
{
  int i;

  for (i = 0; i < N_BAND_TYPES; i++) {
    if (strncmp(band_types[i].type,type,strlen(band_types[i].type)) == 0) {
      return 1;
    }
  }
  return 0;
}

int
isbandconstraint(char *type)
{
  int i;

  for (i = 0; i < N_BAND_CONSTRAINTS; i++) {
    if (strcmp(constraint_type[i].type,type) == 0) {
      return 1;
    }
  }
  return 0;
}

void
list_band_constraints(void)
{
  int i;

  printf("\n\nBand type synthesis constraints\n");
  for (i = 0; i < sizeof(constraint_type)/sizeof(constraint_type[0]); i++) {
    printf("    %11s %s\n",constraint_type[i].type,constraint_type[i].descr);
  }
}

void
calc_crosses(
  LSL_band_constraint *k,
  LSL_band_constraint *m,
  int                 *layer,
  FILE *output)
{
  PRECISION xlk = k->end_line[0] - k->start_line[0];
  PRECISION ylk = k->end_line[1] - k->start_line[1];
  PRECISION xnm = m->end_line[0] - m->start_line[0];
  PRECISION ynm = m->end_line[1] - m->start_line[1];
  PRECISION xmk = m->start_line[0] - k->start_line[0];
  PRECISION ymk = m->start_line[1] - k->start_line[1];

  PRECISION det = xnm*ylk - ynm*xlk;
  if (fabs(det) < ACCY) {
    /* parallel lines */
  } else {
    PRECISION detinv = 1.0/det;
    PRECISION s = (xnm*ymk - ynm*xmk)*detinv;
    PRECISION t = (xlk*ymk - ylk*xmk)*detinv;
    if (s >= 0 && s <= 1.0 && t >= 0 && t <= 1.0) {
      PRECISION x = k->start_line[0] + xlk*s;
      PRECISION y = k->start_line[1] + ylk*s;
      if (k->n_crossings < 8) {
        k->crossings[k->n_crossings][0] = x;
        k->crossings[k->n_crossings][1] = y;
        k->crossings[k->n_crossings][2] = 0;
        k->n_crossings++;
      }
      if (m->n_crossings < 8) {
        m->crossings[m->n_crossings][0] = x;
        m->crossings[m->n_crossings][1] = y;
        m->crossings[m->n_crossings][2] = 0;
        m->n_crossings++;
      }
      if (k->layer == -1) {
        k->layer = *layer;
        *layer = *layer + 1;
      }
      if (m->layer == -1) {
        m->layer = *layer;
        *layer = *layer + 1;
      }
    }
  }
}

void
calc_angles(
  band_attrib_t *type,
  LSL_band_constraint *k,
  FILE *output)
{
  PRECISION first_x, first_y, last_x, last_y;
  PRECISION dx,dy;
  PRECISION angle,ta;
  int i;
  float n;
  PRECISION pi = 2*atan2(1,0);

  if (k->cross || ! k->inside) {
    first_x = k->end_angle[0];
    first_y = k->end_angle[1];
    last_x  = k->start_angle[0];
    last_y  = k->start_angle[1];
  } else {
    first_x = k->start_angle[0];
    first_y = k->start_angle[1];
    last_x  = k->end_angle[0];
    last_y  = k->end_angle[1];
  }

  dx = (first_x - k->part.offset[0]);
  dy = (first_y - k->part.offset[1]);
  dx /= k->radius;
  dy /= k->radius;

  if (dy > 0) {
    angle = acos(dx);
  } else {
    angle = -acos(dx);
  }

  k->s_angle = angle;

  if (type->fill == FIXED3) {
    PRECISION circ;
    if (angle < 0) {
      angle = - angle;
    }
    circ = angle*2*k->radius;
    k->n_steps = circ*type->scale+0.5;

  } else if (type->fill == FIXED) {

    n = type->scale * 2 * pi * k->radius + 0.5;

    // circumference
    for (i = 0; i < n; i++) {
      PRECISION f;
      f = i;
      f /= n;
      ta = angle - 2*pi*(1-f);
      dx = k->radius*cos(ta) + k->part.offset[0] - last_x;
      dy = k->radius*sin(ta) + k->part.offset[1] - last_y;
      if (sqrt(dx*dx+dy*dy) < type->delta) {
        break;
      }
    }
    k->n_steps = i+1;
  } else {

    n =  2 * pi * k->radius/band_res + 0.5;

    // circumference
    for (i = 0; i < n; i++) {
      PRECISION f;
      f = i;
      f /= n;
      ta = angle - 2*pi*(1-f);
      dx = k->radius*cos(ta) + k->part.offset[0] - last_x;
      dy = k->radius*sin(ta) + k->part.offset[1] - last_y;
      if (sqrt(dx*dx+dy*dy) < type->delta) {
        break;
      }
    }
    k->n_steps = i+1;

  }
}

int intersect_line_circle_2D(
  PRECISION xo,
  PRECISION yo,
  PRECISION f,
  PRECISION g,
  PRECISION xj,
  PRECISION yj,
  PRECISION rj,
  PRECISION *x,
  PRECISION *y)
{
  PRECISION fsq, gsq, fgsq;
  PRECISION xjo, yjo;
  PRECISION fygx;
  PRECISION fxgy;
  PRECISION root;
  PRECISION t;

  fsq = f * f;
  gsq = g * g;
  fgsq = fsq + gsq;

  if (fgsq < ACCY) {
    printf("line coefficients are corrupt\n");
  }

  xjo = xj - xo;
  yjo = yj - yo;
  fygx = f*yjo - g*xjo;
  root = rj*rj*fgsq - fygx*fygx;

  if (root < -ACCY) {
    printf("line does not intersect with circle\n");
  }

  fxgy = f*xjo + g*yjo;

  t = fxgy/fgsq;

  *x = xo + f*t;
  *y = yo + g*t;
  return 0;
}

int calc_tangent_line(
  LSL_band_constraint *k,
  LSL_band_constraint *l,
  FILE                *output)
{
  int inside1, inside2;
  PRECISION rl,rk,rlk;
  PRECISION xlk,xlksq,ylk,ylksq;
  PRECISION denom;
  PRECISION radius;
  PRECISION angle;
  PRECISION rx,ry;

  inside1 = k->inside;
  inside2 = l->inside;

  if (l->was_cross) {
    if (l->cross) {
      inside2 ^= 1;
    } else {
      inside1 ^= 1;
    }
  }

  rl = l->radius;
  rk = k->radius;

  switch ((inside1 << 1) | inside2) {
    case 3: /* inside to inside */
      /* end angle for i , start for j */
      /* start point line of ij line stored in j */
      /* endpoint of line of ij line stored in j */
    break;
    case 2: /* inside to outside */
      rl = - rl;
    break;
    case 1: /* outside to inside */
      rk = -rk;
    break;
    case 0: /* outside to outside */
      rk = -rk;
      rl = -rl;
    break;
  }

  rlk = rl - rk;

  /* rotate xk and yk up to -45 degrees */

  xlk = l->part.offset[0] - k->part.offset[0];
  ylk = l->part.offset[1] - k->part.offset[1];

  xlksq = xlk*xlk;
  ylksq = ylk*ylk;

  denom = xlksq + ylksq;

  if (denom < ACCY) {
    /* circles are coincident - badness */
  } else {
    PRECISION root;

    root = denom - rlk*rlk;
    if (root < -ACCY) {
      /* tangent doesn exist */
    } else {

      PRECISION a,b,c;
      PRECISION deninv,factor;
      PRECISION xo,yo;
      PRECISION f,g;

      if (root < 0) {
        root = 0;
      }
      root = sqrt(root);
      deninv = 1.0/denom;
      a = (-rlk*xlk - ylk*root)*deninv;
      b = (-rlk*ylk + xlk*root)*deninv;
      c = -(rk + a*k->part.offset[0] + b*k->part.offset[1]);

      /* we have the normalized form of the tangent line */

      /* now we map the line to parametric form */
      root = 1.0/(a * a + b * b);
      factor = -c*root;
      xo = a*factor;
      yo = b*factor;

      root = sqrt(root);

      f =  b*root;
      g = -a*root;

      /* now line is in x = xo + f*t, y = yo + g*t form */
      /* calculate endpoints of each line */

      intersect_line_circle_2D(
        xo,
        yo,
        f,
        g,
        k->part.offset[0],
        k->part.offset[1],
        k->radius,
       &k->start_line[0],
       &k->start_line[1]);
        k->start_line[2] = 0;
      vectorcp(k->start_angle,k->start_line);

      intersect_line_circle_2D(
        xo,
        yo,
        f,
        g,
        l->part.offset[0],
        l->part.offset[1],
        l->radius,
       &k->end_line[0],
       &k->end_line[1]);
        k->end_line[2] = 0;
      vectorcp(l->end_angle,k->end_line);

      // this means we need our previous neighbor's end line and our
      // start line to know the arc
    }
  }
  return 0;
}

int draw_arc_line(
  band_attrib_t       *type,
  LSL_band_constraint *constraint,
  int                  color,
  int                  draw_line,
  FILE                *output,
  int                  ghost,
  part_t      *absolute,
  LSL_band_constraint *f_constraint)
{
  int       i,j,n;
  PRECISION dx,dy,dz;
  PRECISION L1,L2;
  int steps;

  if (draw_line) {

    for (j = 0; j < constraint->n_crossings - 1; j++) {
      PRECISION orient[3][3];

      // determine the orientation of the part in the XY plane

      dx = constraint->crossings[j+1][0] - constraint->crossings[j][0];
      dy = constraint->crossings[j+1][1] - constraint->crossings[j][1];
      dz = constraint->crossings[j+1][2] - constraint->crossings[j][2];

      L1 = sqrt(dx*dx + dy*dy);
      L2 = sqrt(dx*dx + dy*dy + dz*dz);
      if (L1 == 0) {

        orient[0][0] = 1;
        orient[1][0] = 0;
        orient[2][0] = 0;
        orient[0][1] = 0;
        orient[1][1] = 1;
        orient[2][1] = 0;
        orient[0][2] = 0;
        orient[1][2] = 0;
        orient[2][2] = 1;
      } else {
        orient[0][0] =  dy/L1;  //  cos
        orient[1][0] = -dx/L1;  //  sin
        orient[2][0] =  0;
        orient[0][1] =  dx/L2;  // -sin
        orient[1][1] =  dy/L2;  //  cos
        orient[2][1] =  dz/L2;
        orient[0][2] = -dx*dz/(L1*L2);
        orient[1][2] = -dy*dz/(L1*L2);
        orient[2][2] =  L1/L2;
      }

      if (type->fill == STRETCH) {
        n = 1;
        steps = 0;
      } else if (type->fill == FIXED3) {
        n = L2*type->scale+0.5;
        steps = 1;
      } else {
        n = L2*type->scale+0.5;
        steps = 0;
      }

      for (i = steps; i < n; i++) {
        part_t part;
        PRECISION tm[3][3];
        PRECISION foffset[3];

        part = type->tangent;

        if (type->fill == STRETCH) {
          PRECISION scale[3][3];
          int i,j;

          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              scale[i][j] = 0;
            }
          }
          scale[0][0] = 1;
          scale[1][1] = L2;
          scale[2][2] = 1;
          matrixmult(part.orient,scale);
        }

        // We performed this to start:
        //   1.  move the assembly so the first constraint is at the origin
        //   2.  rotate offsetations about the inverse of the first constraint's
        //       orientation bringing everything into the X/Y plane

        // Now we're putting together the segments and putting them back into place
        //   1.  multiply tangent part orientation times the orient
        //       so we know how to orient the tangent part in the X/Y plane.
        //   2.  rotate the tangent offsets into the X/Y plane.

        matrixcp(tm,part.orient);
        matrixmult3(part.orient,orient,tm);

        vectorrot3(part.offset,type->tangent.offset,part.orient);

        //   3.  Calculate the tangent part's offsetation in the X/Y plane

        foffset[0] = constraint->crossings[j][0] + dx * i / n - part.offset[0];
        foffset[1] = constraint->crossings[j][1] + dy * i / n - part.offset[1];
        foffset[2] = constraint->crossings[j][2] + dz * i / n - part.offset[2];

        vectorcp(part.offset,foffset);

        //   4.  Reoffsetate the tangent part in the X/Y plane based on the first
        //       constraint's offset (for things like technic turntable top,
        //       where the gear plane does not go through the origin.

        vectoradd( part.offset,constraint_type[f_constraint->constraint_type_n].offset);

        //   5.  Orient tangent part based on the orientation of the first
        //       constraint's orientation (for things like technic turntable
        //       who's gear plane is perpendicular to the plane of say the 24T
        //       gears.

        vectorrot( part.offset,constraint_type[f_constraint->constraint_type_n].orient);

        //   6.  Orient the tangent part offsetation back to the absolute 3D
        //       offsetation.

        vectorrot( part.offset,absolute->orient);

        //   7.  Now add the absolute offsetation of the first constraint

        vectoradd( part.offset,absolute->offset);

        //   8.  Change the tangent part orientation based on the first
        //       constraint's orientation (again for things like the
        //       technic turntable top, where the gear plane is perpendicular
        //       to the standard 24T gear's gear plane

        matrixcp(tm,part.orient);
        matrixmult3(part.orient,
                    constraint_type[f_constraint->constraint_type_n].orient,
                    tm);

        //   9.  Now rotate the part back into its correct orientation in 3D
        //       space.

        matrixcp(tm,part.orient);
        matrixmult3(part.orient,absolute->orient,tm);

        output_line(
          output,
          ghost,
          color,
          part.offset[0],part.offset[1],part.offset[2],
          part.orient[0][0],
          part.orient[0][1],
          part.orient[0][2],
          part.orient[1][0],
          part.orient[1][1],
          part.orient[1][2],
          part.orient[2][0],
          part.orient[2][1],
          part.orient[2][2],
          type->tangent.type);
      }
    }
  }

  // Create the arc

  {
    PRECISION pi = 2*atan2(1,0);
    PRECISION f[3];

    /* now for the arc */

    if (type->fill == STRETCH) {
      n = 2*pi*constraint->radius/band_res;
    } else {
      n = 2*pi*constraint->radius*type->scale;
    }

    // vector for the first arc part

    f[0] = constraint->radius * cos(constraint->s_angle);
    f[1] = constraint->radius * sin(constraint->s_angle);
    f[2] = 0;

    if (type->fill == FIXED3) {
      steps = constraint->n_steps + 2;
    } else {
      steps = constraint->n_steps;
    }

    for (i = 1; i < steps; i++) {
      PRECISION orient[3][3];
      PRECISION foffset[3];
      PRECISION tm[3][3];
      PRECISION s[3];
      part_t    part;

      // for FIXED3 (e.g. rubber tread), the first and last segments are
      // treated as transition pieces, otherwise just use arc parts

      if (type->fill == FIXED3 && i == 1) {
        part = type->start_trans;
      } else if (type->fill == FIXED3 && i+1 == steps) {
        part = type->end_trans;
      } else {
        part = type->arc;
      }

      // rotate the arc part so it hugs the current constraint

      s[0] = constraint->radius * cos(constraint->s_angle + 2*pi*i/n);
      s[1] = constraint->radius * sin(constraint->s_angle + 2*pi*i/n);
      s[2] = 0;

      dx = s[0] - f[0];
      dy = s[1] - f[1];
      dz = s[2] - f[2];

      L1 = sqrt(dx*dx + dy*dy);
      L2 = sqrt(dx*dx + dy*dy + dz*dz);

      if (L1 == 0) {
        orient[0][0] =  1;
        orient[1][0] =  0;
        orient[2][0] =  0;
        orient[0][1] =  0;
        orient[1][1] =  1;
        orient[2][1] =  0;
        orient[0][2] =  0;
        orient[1][2] =  0;
        orient[2][2] =  1;
      } else {
        orient[0][0] =  dy/L1;
        orient[1][0] = -dx/L1;
        orient[2][0] =  0;
        orient[0][1] =  dx/L2;
        orient[1][1] =  dy/L2;
        orient[2][1] =  dz/L2;
        orient[0][2] = -dx*dz/(L1*L2);
        orient[1][2] = -dy*dz/(L1*L2);
        orient[2][2] =  L1/L2;
      }

      if (type->fill == STRETCH) {
        PRECISION scale[3][3];
        PRECISION angle = 2 * pi / n;
        PRECISION l = type->scale*sin(angle);
        int i,j;

        for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            scale[i][j] = 0;
          }
        }
        scale[0][0] = 1;
        scale[1][1] = L2+l;
        scale[2][2] = 1;
        matrixmult(part.orient,scale);
      }

      // Now we're putting together the segments and putting them back into place
      //   1.  multiply arc part orientation times the orient
      //       so we know how to orient the arc part in the X/Y plane.
      //   2.  rotate the arc offsets into the X/Y plane.

      matrixcp(tm,part.orient);
      matrixmult3(part.orient,orient,tm);

      vectorrot3(part.offset,type->arc.offset,part.orient);

      //   3.  Calculate the arc part's offsetation in the X/Y plane

      foffset[0] = f[0] + constraint->part.offset[0] - part.offset[0];
      foffset[1] = f[1] + constraint->part.offset[1] - part.offset[1];
      foffset[2] = f[2] + constraint->part.offset[2] - part.offset[2];

      vectorcp(part.offset,foffset);

      //   4.  Rotate the arc part in the X/Y plane based on the first
      //       constraint's offset (for things like technic turntable top,
      //       where the gear plane does not go through the origin.)

      vectoradd(part.offset,constraint_type[f_constraint->constraint_type_n].offset);
      vectoradd(part.offset,constraint_type[f_constraint->constraint_type_n].offset);

      //   5.  Orient arc part based on the orientation of the first
      //       constraint's orientation (for things like technic turntable
      //       who's gear plane is perpendicular to the plane of say the 24T
      //       gears.)

      vectorrot(part.offset,constraint_type[f_constraint->constraint_type_n].orient);

      //   6.  Orient the arc part offsetation back to the absolute 3D
      //       offsetation.

      vectorrot(part.offset,absolute->orient);

      //   7.  Now add the absolute offsetation of the first constraint

      vectoradd( part.offset,absolute->offset);

      //   8.  Change the arc part orientation based on the first
      //       constraint's orientation (again for things like the
      //       technic turntable top, where the gear plane is perpendicular
      //       to the standard 24T gear's gear plane)

      matrixcp(tm,part.orient);
      matrixmult3(part.orient,
                  constraint_type[f_constraint->constraint_type_n].orient,
                  tm);

      //   9.  Now rotate the part back into its correct orientation in 3D
      //       space.

      matrixcp(tm,part.orient);
      matrixmult3(part.orient,absolute->orient,tm);

      output_line(
        output,
        ghost,
        color,
        part.offset[0],part.offset[1],part.offset[2],
        part.orient[0][0],
        part.orient[0][1],
        part.orient[0][2],
        part.orient[1][0],
        part.orient[1][1],
        part.orient[1][2],
        part.orient[2][0],
        part.orient[2][1],
        part.orient[2][2],
        part.type);

      vectorcp(f,s);
    }
  }
  return 0;
}

void
showconstraints(
  FILE                *output,
  LSL_band_constraint *constraints,
  int                  n_constraints,
  int                  color)
{
#if 0
  int i;

  for (i = 0; i < n_constraints; i++) {
    output_line(
      output,
      0,
      color,
      constraints[i].part.offset[0],
      constraints[i].part.offset[1],
      constraints[i].part.offset[2],
      constraints[i].part.orient[0][0],
      constraints[i].part.orient[0][1],
      constraints[i].part.orient[0][2],
      constraints[i].part.orient[1][0],
      constraints[i].part.orient[1][1],
      constraints[i].part.orient[1][2],
      constraints[i].part.orient[2][0],
      constraints[i].part.orient[2][1],
      constraints[i].part.orient[2][2],
      constraints[i].part.type);
  }
#endif
}

static void
rotate_constraints(
  LSL_band_constraint *constraints,
  int                  n_constraints,
  PRECISION            m[3][3])
{
  int i;
  PRECISION t[3][3];

  for (i = 0; i < n_constraints; i++) {
    vectorrot(constraints[i].part.offset,m);
    matrixcp(t,constraints[i].part.orient);
    matrixmult3(constraints[i].part.orient,m,t);
  }
}

/*
 * This subroutine synthesizes planar rubber bands, chain, and treads.
 *
 * We do all the arc and tangent analysis in the X/Y plane.
 *
 * The synthesis plane is defined by the first constraint.  First we move the
 * first constaint to the origin. Then we calculate the inverse of the first
 * constraints orientation, and multiply all the constraints' offsetations by
 * the inverse of the first constraint's orientation.
 *
 * Some of the gears' are oriented in the X/Y plane, while other gears are
 * oriented in the X/Z plane.  The constraint_type array above, describes each of the
 * supported constraint types, and their orientation.  We multiply all the
 * constraints by the first constraint's constraint_type orientation.
 *
 * Also, in some cases, some of the constraint types described in constraint_type
 * need to be offset to get the place where the band should hit onto the
 * X/Y plane.
 */
int
synth_band(
  char *type,
  int n_constraints,
  LSL_band_constraint *constraints,
  int color,
  FILE *output,
  int ghost)
{
  int i;
  int cross = 0;
  int was_cross = 0;
  int inside = 1;
  int first,last;
  part_t absolute;
  PRECISION inv[3][3],trot[3][3];
  int layer = 0;
  band_attrib_t *band_type = NULL;

  /* Search for band type */

  for (i = 0; i < N_BAND_TYPES; i++) {
    if (strcmp(type,band_types[i].type) == 0) {
      band_type = &band_types[i];
      break;
    }
  }
  if (band_type == NULL) {
    return 0;
  }

  first = -1;

  for (i = 0; i < n_constraints; i++) {
    constraints[i].radius       = 0;
    constraints[i].inside       = inside;
    constraints[i].cross        = cross;
    constraints[i].was_cross    = was_cross;
    constraints[i].n_crossings  = 0;
    constraints[i].layer        = -1;
    was_cross = 0;
    if (strcmp(constraints[i].part.type,"INSIDE") == 0) {
      inside = 1;
    } else if (strcmp(constraints[i].part.type,"OUTSIDE") == 0) {
      inside = 0;
    } else if (strcmp(constraints[i].part.type,"CROSS") == 0) {
      inside ^= 1;

    } else {
      int k;

      // search the constraints table

      for (k = 0; k < sizeof(constraint_type)/sizeof(constraint_type[0]); k++) {
        if (strcmp(constraints[i].part.type,constraint_type[k].type) == 0) {
          constraints[i].constraint_type_n = k;

          constraints[i].radius = constraint_type[k].radius;
          break;
        }
      }
    }
    if (first == -1 && constraints[i].radius) {
      first = i;
    }
  }

  /* create an extra constraint that represents the final state of the
   * first pulley */

  memcpy(&constraints[i],&constraints[first],sizeof(constraints[i]));
  constraints[i].inside    = inside;
  constraints[i].cross     = cross;
  constraints[i].was_cross = was_cross;

  /* record the first constraint in its original form */

  absolute = constraints[first].part;

  /* 1. move the first constraint to the origin */

  for (i = 0; i <= n_constraints; i++) {
    if (constraints[i].radius) {
      vectorsub(constraints[i].part.offset,absolute.offset);
    }
  }

  showconstraints(output,constraints,n_constraints,14);

  /* 2. bring the entire assembly into the part's natural orientation */

  matrixinv(inv,absolute.orient);

  rotate_constraints(constraints,n_constraints+1,inv);

  showconstraints(output,constraints,n_constraints,4);

  /* 3. bring the assembly into the X/Y plane (necessary for first constraints
   *    who's gear plane is different that the default gear plane used by
   *    simple gears, like technic turntable).
   */

  rotate_constraints(constraints,n_constraints+1,
    constraint_type[constraints[first].constraint_type_n].orient);

  showconstraints(output,constraints,n_constraints,15);

  /* 4. offset constraint's who's gear plane does not intersect with the
   *    origin.
   */

  for (i = 0; i <= n_constraints; i++) {
    if (constraints[i].radius) {
      vectorsub(constraints[i].part.offset,
        constraint_type[constraints[i].constraint_type_n].offset);
    }
  }

  showconstraints(output,constraints,n_constraints,3);

  /* figure out the tangents' intersections with circles */

  first = -1;
  last  = -1;

  for (i = 0; i < n_constraints; ) {
    if (constraints[i].radius) {
      int j;
      if (first == -1) {
        first = i;
      }
      for (j = i+1; j <= n_constraints; j++) {
        if (constraints[j].radius) {
          calc_tangent_line(&constraints[i],&constraints[j],output);
          i = j;
          last = j;
          break;
        }
      }
      if (j == n_constraints+1) {
        i++;
      }
    } else {
      i++;
    }
  }
  vectorcp(constraints[first].end_angle,constraints[last].end_angle);

  /* calculate intersections between band straight line segments, so
   * we can make the line segments go around each other.
   */

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      vectorcp(constraints[i].crossings[0],constraints[i].start_line);
      constraints[i].n_crossings = 1;
    }
  }
  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      int j;
      for (j = i+1; j <= n_constraints; j++) {
        if (constraints[j].radius) {
          calc_crosses(&constraints[i],&constraints[j],&layer,output);
        }
      }
    }
  }

#define BAND_DIAM 4

  /* calculate the depth for each band at crossing */

  if (layer > 0) {
    {
      PRECISION layer_offset;
      layer_offset = (((layer) / 2) % 2) * BAND_DIAM/2;
      //layer_offset = 0;
      for (i = 0; i < n_constraints; i++) {
        if (constraints[i].radius) {
          int j;
          for (j = 1; j < constraints[i].n_crossings; j++) {
            constraints[i].crossings[j][2] =
              (constraints[i].layer-layer/2)*BAND_DIAM + layer_offset;
          }
        }
      }
    }
  }

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      vectorcp(constraints[i].crossings[constraints[i].n_crossings++],
               constraints[i].end_line);
    }
  }

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      calc_angles(band_type,&constraints[i],output);
    }
  }

  /*****************************************************
   * rotate everything back to whence it came
   * and put all parts back to original absolute
   * coordinates.
   *****************************************************/

  fprintf(output,"0 SYNTH SYNTHESIZED BEGIN\n");

  /* now draw out the rubber band in terms of lines and arcs */
  for (i = 0; i <= n_constraints; ) {
    if (constraints[i].radius != 0) {
      int j;
      for (j = i+1; j <= n_constraints; j++) {

        if (constraints[j].radius) {
          draw_arc_line(
            band_type,
            &constraints[i],
            color,
            n_constraints > 1,
            output,
            ghost,
            &absolute,
            &constraints[first]);
          i = j;
          break;
        }
      }
      if (j >= n_constraints) {
        i++;
      }
    } else {
      i++;
    }
  }

  fprintf(output,"0 SYNTH SYNTHESIZED END\n");
  printf("Synthesized %s\n",type);
  return 0;
}
