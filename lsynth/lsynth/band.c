/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynth.h"
#include "band.h"
#include "tube.h"
#include <float.h>

#define FIXED   0
#define FIXED3  1
#define STRETCH 2

typedef struct {
  char     *part;         // LDraw part type
  PRECISION orient[3][3]; // how to orient/scale it into the result
  PRECISION offset[3];    // how to displace it into the result
} part_t;

typedef struct {
  char     *type;   // name of thing being specified (e.g. RUBBER_BAND)
  int       fill;   // method for synthesizing
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
    "RUBBER_BAND ",
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
    "RUBBER_BELT ",
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
    "CHAIN ",
    FIXED,
    1.0/16,
    8,
    { // tangent
      "3711.DAT",
      {
        { 0, 1, 0 },
        { 0, 0,-1 },
        { -1,0, 0 }
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
        { -1, 0, 0 }
      },
      {
        0, 0, 0
      }
    }
  },
  {
    "PLASTIC_TREAD ",
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
    "RUBBER_TREAD ",
    FIXED3,
    1.0/20,
    4,
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
  },
};

struct {
  char     *type;         // LDraw part name
  PRECISION radius;       // Radius of circular part
  PRECISION orient[3][3]; // How to orient it
  PRECISION offset[3];    // How much to offset it.
} radiai[] = {

  // bushings/pulleys

  {    "3736.DAT", 44, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "4185.DAT", 30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "2983.DAT", 11, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "4265.DAT", 10, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "4265A.DAT", 10, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "4265B.DAT", 10, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "4265C.DAT", 10, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "13.DAT", 10, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "2736.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "6628.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},

  // axles

  {    "3704.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {   "32062.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "4519.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "6587.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3705.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  { "3705C01.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {   "32073.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {     "552.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3706.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3707.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  {    "3737.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},
  { "3737C01.DAT",  4, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "3708.DAT",  4, {{0, 0, 1},{0,1, 0},{-1,0,0}}},

  // axle joiner

  {   "6538A.DAT", 10, {{0, 0, 1},{0,1, 0},{ 0,0,1}}},

  // gears

  {   "32007.DAT", 32, {{1, 0, 0},{0,1, 0},{ 0,0,1}}, { 0, 32, 0}},
  {   "73071.DAT", 35, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "6573.DAT", 30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  { 9   "4019.DAT", 20, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "6542.DAT", 20, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "3648.DAT", 30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "60C01.DAT", 30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {   "3650A.DAT", 30, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "3649.DAT", 50, {{1, 0, 0},{0,1, 0},{ 0,0,1}}},
  {    "2855.DAT", 73, {{1, 0, 0},{0,0,-1},{ 0,1,0}}, { -8, 0, 0}},
  {    "6539.DAT", 20, {{1, 0, 0},{0,1, 0},{ 0,1,0}}},
  {    "LS00.DAT",  1, {{1, 0, 0},{0,0,-1},{ 0,1,0}}},
};

void
calc_crosses(
  LSL_band_constraint *k,
  LSL_band_constraint *m,
  int                 *layer,
  FILE *output)
{
  PRECISION xlk = k->end_line.loc.x - k->start_line.loc.x;
  PRECISION ylk = k->end_line.loc.y - k->start_line.loc.y;
  PRECISION xnm = m->end_line.loc.x - m->start_line.loc.x;
  PRECISION ynm = m->end_line.loc.y - m->start_line.loc.y;
  PRECISION xmk = m->start_line.loc.x - k->start_line.loc.x;
  PRECISION ymk = m->start_line.loc.y - k->start_line.loc.y;

  PRECISION det = xnm*ylk - ynm*xlk;
  if (fabs(det) < ACCY) {
    /* parallel lines */
  } else {
    PRECISION detinv = 1.0/det;
    PRECISION s = (xnm*ymk - ynm*xmk)*detinv;
    PRECISION t = (xlk*ymk - ylk*xmk)*detinv;
    if (s >= 0 && s <= 1.0 && t >= 0 && t <= 1.0) {
      PRECISION x = k->start_line.loc.x + xlk*s;
      PRECISION y = k->start_line.loc.y + ylk*s;
      if (k->n_crossings < 8) {
        k->crossings[k->n_crossings].loc.x = x;
        k->crossings[k->n_crossings].loc.y = y;
        k->crossings[k->n_crossings].loc.z = 0;
        k->n_crossings++;
      }
      if (m->n_crossings < 8) {
        m->crossings[m->n_crossings].loc.x = x;
        m->crossings[m->n_crossings].loc.y = y;
        m->crossings[m->n_crossings].loc.z = 0;
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
    first_x = k->end_angle.loc.x;
    first_y = k->end_angle.loc.y;
    last_x  = k->start_angle.loc.x;
    last_y  = k->start_angle.loc.y;
  } else {
    first_x = k->start_angle.loc.x;
    first_y = k->start_angle.loc.y;
    last_x  = k->end_angle.loc.x;
    last_y  = k->end_angle.loc.y;
  }

  dx = (first_x - k->part.loc.x);
  dy = (first_y - k->part.loc.y);
  dx /= k->radius;
  dy /= k->radius;

  if (dy > 0) {
    angle = acos(dx);
  } else {
    angle = -acos(dx);
  }

  k->s_angle = angle;

  n = type->scale * 2 * pi * k->radius + 0.5;

  if (type->fill == FIXED3) {
    PRECISION circ;
    if (angle < 0) {
      angle = - angle;
    }
    circ = angle*2*k->radius;
    k->n_steps = circ*type->scale+0.5;
  } else {

    // circumference
    for (i = 0; i < n; i++) {
      PRECISION f;
      f = i;
      f /= n;
      ta = angle - 2*pi*(1-f);
      dx = k->radius*cos(ta) + k->part.loc.x - last_x;
      dy = k->radius*sin(ta) + k->part.loc.y - last_y;
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

  xlk = l->part.loc.x - k->part.loc.x;
  ylk = l->part.loc.y - k->part.loc.y;

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
      c = -(rk + a*k->part.loc.x + b*k->part.loc.y);

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
        k->part.loc.x,
        k->part.loc.y,
        k->radius,
        &k->start_line.loc.x,
        &k->start_line.loc.y);
        k->start_line.loc.z = 0;
        k->start_angle = k->start_line;

      intersect_line_circle_2D(
        xo,
        yo,
        f,
        g,
        l->part.loc.x,
        l->part.loc.y,
        l->radius,
        &k->end_line.loc.x,
        &k->end_line.loc.y);
        k->end_line.loc.z = 0;
        l->end_angle = k->end_line;

      // this means we need our previous neighbor's end line and our
      // start line to know the arc
    }
  }
  return 0;
}

static void
rotate_point3(
  PRECISION t[3],
  PRECISION r[3],
  PRECISION m[3][3])
{
  t[0] = r[0]*m[0][0] + r[1]*m[0][1] + r[2]*m[0][2];
  t[1] = r[0]*m[1][0] + r[1]*m[1][1] + r[2]*m[1][2];
  t[2] = r[0]*m[2][0] + r[1]*m[2][1] + r[2]*m[2][2];
}

int draw_arc_line(
  band_attrib_t *type,
  LSL_band_constraint *f_constraint,
  int color,
  int draw_line,
  FILE *output,
  int   ghost,
  LSL_part_usage *absolute)
{
  PRECISION r;
  int       i,j,n;
  PRECISION L1,L2;
  PRECISION pi = 2*atan2(1,0);
  PRECISION angle, at;
  PRECISION f[3],s[3],rf[3],rs[3];
  PRECISION orient[3][3];
  PRECISION rot[3][3];
  PRECISION roffset[3];
  PRECISION dx,dy,dz;
  int steps;

  if (draw_line) {

    for (j = 0; j < f_constraint->n_crossings - 1; j++) {

      // determine the orientation of the part in the XY plane

      dx = f_constraint->crossings[j+1].loc.x - f_constraint->crossings[j].loc.x;
      dy = f_constraint->crossings[j+1].loc.y - f_constraint->crossings[j].loc.y;
      dz = f_constraint->crossings[j+1].loc.z - f_constraint->crossings[j].loc.z;

      L1 = sqrt(dx*dx + dy*dy);
      L2 = sqrt(dx*dx + dy*dy + dz*dz);

      if (L1 == 0) {

        f_constraint->start_line.orient[0][0] = 1;
        f_constraint->start_line.orient[1][0] = 0;
        f_constraint->start_line.orient[2][0] = 0;
        f_constraint->start_line.orient[0][1] = 0;
        f_constraint->start_line.orient[1][1] = 1;
        f_constraint->start_line.orient[2][1] = 0;
        f_constraint->start_line.orient[0][2] = 0;
        f_constraint->start_line.orient[1][2] = 0;
        f_constraint->start_line.orient[2][2] = 1;
      } else {

        f_constraint->start_line.orient[0][0] = dy/L1;
        f_constraint->start_line.orient[1][0] = -dx/L1;
        f_constraint->start_line.orient[2][0] = 0;
        f_constraint->start_line.orient[0][1] = dx/L2;
        f_constraint->start_line.orient[1][1] = dy/L2;
        f_constraint->start_line.orient[2][1] = dz/L2;
        f_constraint->start_line.orient[0][2] = -dx*dz/(L1*L2);
        f_constraint->start_line.orient[1][2] = -dy*dz/(L1*L2);
        f_constraint->start_line.orient[2][2] = L1/L2;
      }

      // rotate the orientation of tangent line and part_orient

      mat_mult(rot,f_constraint->start_line.orient,type->tangent.orient);

      // rotate the offset by

      rotate_point3(roffset,type->tangent.offset,rot);

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
        PRECISION loc[3];
        PRECISION floc[3];
        PRECISION trot[3][3];

        floc[0] = f_constraint->crossings[j].loc.x + dx * i / n - roffset[0];
        floc[1] = f_constraint->crossings[j].loc.y + dy * i / n - roffset[1];
        floc[2] = f_constraint->crossings[j].loc.z + dz * i / n - roffset[2];

        rotate_point3(loc, floc,absolute->orient);

        loc[0] += absolute->loc.x;
        loc[1] += absolute->loc.y;
        loc[2] += absolute->loc.z;

        if (type->fill == STRETCH) {
          PRECISION scale[3][3];
          PRECISION lrot[3][3];
          int i,j;

          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              lrot[i][j] = rot[i][j];
              scale[i][j] = 0;
            }
          }
          scale[0][0] = 1;
          scale[1][1] = L2;
          scale[2][2] = 1;
          mat_mult(rot,lrot,scale);
        }

        mat_mult(trot,absolute->orient,rot);

        output_line(
          output,
          ghost,
          color,
          loc[0],loc[1],loc[2],
          trot[0][0],
          trot[0][1],
          trot[0][2],
          trot[1][0],
          trot[1][1],
          trot[1][2],
          trot[2][0],
          trot[2][1],
          trot[2][2],
          type->tangent.part);
      }
    }
  }

  /* now for the arc */

  n = 2*pi*f_constraint->radius*type->scale + 0.5;
  angle = f_constraint->s_angle;

  f[0] = f_constraint->radius * cos(angle);
  f[1] = f_constraint->radius * sin(angle);
  f[2] = 0;

  if (type->fill == FIXED3) {
    steps = f_constraint->n_steps + 2;
  } else {
    steps = f_constraint->n_steps;
  }

  for (i = 1; i < steps; i++) {
    PRECISION loc[3];
    PRECISION floc[3];
    PRECISION trot[3][3];
    part_t   *part;

    if (type->fill == FIXED3) {
       if (i == 1) {
         part = &type->start_trans;
       } else if (i+1 == steps) {
         part = &type->end_trans;
       } else {
         part = &type->arc;
       }
    } else {
      part = &type->arc;
    }

    s[0] = f_constraint->radius * cos(angle + 2*pi*i/n);
    s[1] = f_constraint->radius * sin(angle + 2*pi*i/n);
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

    mat_mult(rot,orient,part->orient);

    // rotate the offset by

    rotate_point3(roffset,f_constraint->offset,rot);

    floc[0] = f[0] + f_constraint->part.loc.x - roffset[0];
    floc[1] = f[1] + f_constraint->part.loc.y - roffset[1];
    floc[2] = f[2] + f_constraint->part.loc.z - roffset[2];

    rotate_point3(loc, floc,absolute->orient);

    loc[0] += absolute->loc.x;
    loc[1] += absolute->loc.y;
    loc[2] += absolute->loc.z;

    mat_mult(trot,absolute->orient,rot);

    output_line(output,
      ghost,
      color,
      loc[0],loc[1],loc[2],
      trot[0][0],
      trot[0][1],
      trot[0][2],
      trot[1][0],
      trot[1][1],
      trot[1][2],
      trot[2][0],
      trot[2][1],
      trot[2][2],
      part->part);

    f[0] = s[0];
    f[1] = s[1];
    f[2] = s[2];
  }
  return 0;
}

static void
rotate_point(
  LSL_3D *loc,
  PRECISION m[3][3])
{
  PRECISION t[3];

  t[0] = loc->x*m[0][0] + loc->y*m[0][1] + loc->z*m[0][2];
  t[1] = loc->x*m[1][0] + loc->y*m[1][1] + loc->z*m[1][2];
  t[2] = loc->x*m[2][0] + loc->y*m[2][1] + loc->z*m[2][2];

  loc->x = t[0];
  loc->y = t[1];
  loc->z = t[2];
}

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
  LSL_part_usage absolute;
  PRECISION xangle,yangle,zangle;
  PRECISION xrot[3][3],yrot[3][3],zrot[3][3],trot[3][3];
  int layer = 0;
  band_attrib_t *band_type = NULL;
  LSL_band_constraint save_constraint;

  for (i = 0; i < sizeof(band_types)/sizeof(band_types[0]); i++) {
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
    constraints[i].radius     = 0;
    constraints[i].inside     = inside;
    constraints[i].cross      = cross;
    constraints[i].was_cross  = was_cross;
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
      int j;

      // search the constraints table

      for (j = 0; j < sizeof(radiai)/sizeof(radiai[0]); j++) {
        if (strcmp(constraints[i].part.type,radiai[j].type) == 0) {
          PRECISION torient[3][3];

          torient[0][0] = constraints[i].part.orient[0][0];
          torient[0][1] = constraints[i].part.orient[0][1];
          torient[0][2] = constraints[i].part.orient[0][2];
          torient[1][0] = constraints[i].part.orient[1][0];
          torient[1][1] = constraints[i].part.orient[1][1];
          torient[1][2] = constraints[i].part.orient[1][2];
          torient[2][0] = constraints[i].part.orient[2][0];
          torient[2][1] = constraints[i].part.orient[2][1];
          torient[2][2] = constraints[i].part.orient[2][2];

          mat_mult(constraints[i].part.orient,radiai[j].orient,torient);

          constraints[i].radius = radiai[j].radius;
          constraints[i].offset[0] = radiai[j].offset[0];
          constraints[i].offset[1] = radiai[j].offset[1];
          constraints[i].offset[2] = radiai[j].offset[2];
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

  /*****************************************************
   * adjust all the constraints so that the first
   * pully is the origin.
   *
   * rotate from the plane indicated by first pulley
   * into the x/y plane.
   *
   * hint: if rotated about the x and/or y plane, reverse the
   * rotation.
   *****************************************************/

  memcpy(&absolute,&constraints[first].part,sizeof(absolute));

  for (i = 0; i <= n_constraints; i++) {
    if (constraints[i].radius) {
      constraints[i].part.loc.x -= absolute.loc.x;
      constraints[i].part.loc.y -= absolute.loc.y;
      constraints[i].part.loc.z -= absolute.loc.z;
    }
  }

  xangle = asin(absolute.orient[1][2]);
  yangle = asin(-absolute.orient[0][2]);
  zangle = -asin(-absolute.orient[0][1]);

  xrot[0][0] = 1;
  xrot[0][1] = 0;
  xrot[0][2] = 0;
  xrot[1][0] = 0;
  xrot[1][1] = cos(xangle);
  xrot[1][2] = -sin(xangle);
  xrot[2][0] = 0;
  xrot[2][1] = sin(xangle);
  xrot[2][2] = cos(xangle);

  for (i = 0; i <= n_constraints; i++) {
    if (constraints[i].radius) {
      rotate_point(&constraints[i].part.loc,xrot);
    }
  }

  yrot[0][0] = cos(yangle);
  yrot[0][1] = 0;
  yrot[0][2] = sin(yangle);
  yrot[1][0] = 0;
  yrot[1][1] = 1;
  yrot[1][2] = 0;
  yrot[2][0] = -sin(yangle);
  yrot[2][1] = 0;
  yrot[2][2] = cos(yangle);

  for (i = 0; i <= n_constraints; i++) {
    if (constraints[i].radius) {
      rotate_point(&constraints[i].part.loc,yrot);
    }
  }

  zrot[0][0] =  cos(zangle);
  zrot[0][1] = -sin(zangle);
  zrot[0][2] = 0;
  zrot[1][0] =  sin(zangle);
  zrot[1][1] =  cos(zangle);
  zrot[1][2] = 0;
  zrot[2][0] = 0;
  zrot[2][1] = 0;
  zrot[2][2] = 1;

  for (i = 0; i <= n_constraints; i++) {
    if (constraints[i].radius) {
      rotate_point(&constraints[i].part.loc,zrot);
    }
  }


  for (i = 0; i < sizeof(radiai)/sizeof(radiai[0]); i++) {
    if (strcmp(constraints[0].part.type,radiai[i].type) == 0) {
      int j,k;
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          save_constraint.part.orient[j][k] = constraints[0].part.orient[j][k];
          constraints[0].part.orient[j][k] = radiai[i].orient[j][k];
        }
      }
    }
  }

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
  constraints[first].end_angle = constraints[last].end_angle;

  /* calculate intersections between band straight line segments, so
   * we can make the line segments go around each other.
   */

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      constraints[i].crossings[0] = constraints[i].start_line;
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
            constraints[i].crossings[j].loc.z =
              (constraints[i].layer-layer/2)*BAND_DIAM + layer_offset;
          }
        }
      }
    }
  }

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      constraints[i].crossings[constraints[i].n_crossings++] = constraints[i].end_line;
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

  for (i = 0; i < 3; i++) {
    int j;
    for (j = 0; j < 3; j++) {
      constraints[0].part.orient[i][j] = save_constraint.part.orient[i][j];
    }
  }

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
            &absolute);
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






