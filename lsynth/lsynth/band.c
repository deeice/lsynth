/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynth.h"
#include "band.h"
#include "tube.h"
#include <float.h>

#define SCALE 2
#define CHAIN_SCALE (1/16.0)
#define TREAD_SCALE (1/18.5)

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
  LSL_band_constraint *k,
  char *type,
  FILE *output)
{
  PRECISION first_x, first_y, last_x, last_y;
  PRECISION dx,dy;
  PRECISION angle,ta;
  int i;
  float n;
  PRECISION pi = 2*atan2(1,0);
  PRECISION delta;

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

  if (strcmp(type,RUBBER_BAND) == 0) {
    n = SCALE*(2*pi*k->radius);
    delta = 1;
  } else if (strcmp(type,CHAIN) == 0) {
    n = CHAIN_SCALE*2*pi*k->radius;
    delta = 8;
  } else if (strcmp(type,PLASTIC_TREAD) == 0) {
    n = CHAIN_SCALE * 2*pi*k->radius;
    delta = 8;
  } else if (strcmp(type,RUBBER_TREAD) == 0) {
    n = TREAD_SCALE * 2*pi*k->radius;
    delta = 4;
  }

  if (strcmp(type,RUBBER_TREAD) == 0) {
    PRECISION circ;
    if (angle < 0) {
      angle = pi/2 - angle;
    }
    circ = 2*angle*k->radius;
    k->n_steps = circ*TREAD_SCALE-1;
  } else {

    // circumference
    for (i = 0; i < n; i++) {
      PRECISION f;
      f = i;
      f /= n;
      ta = angle - 2*pi*(1-f);
      dx = k->radius*cos(ta) + k->part.loc.x - last_x;
      dy = k->radius*sin(ta) + k->part.loc.y - last_y;
      if (sqrt(dx*dx+dy*dy) < delta) {
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
  char *type,
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
  PRECISION part_orient[3][3];
  PRECISION rot[3][3];
  PRECISION offset[3];
  PRECISION roffset[3];
  PRECISION dx,dy,dz;
  char     *seg_type;

  part_orient[0][0] = 1;
  part_orient[0][1] = 0;
  part_orient[0][2] = 0;
  part_orient[1][0] = 0;
  part_orient[1][1] = 1;
  part_orient[1][2] = 0;
  part_orient[2][0] = 0;
  part_orient[2][1] = 0;
  part_orient[2][2] = 1;

  offset[0] = 0;
  offset[1] = 0;
  offset[2] = 0;

  if (strcmp(type,RUBBER_BAND) == 0) {
    seg_type = "LS31.dat";
  } else if (strcmp(type,CHAIN) == 0) {
    seg_type = "3711.DAT";
    part_orient[0][0] = 0;
    part_orient[0][1] = 1;
    part_orient[0][2] = 0;
    part_orient[1][0] = 0;
    part_orient[1][1] = 0;
    part_orient[1][2] = -1;
    part_orient[2][0] = -1;
    part_orient[2][1] = 0;
    part_orient[2][2] = 0;
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 16;
  } else if (strcmp(type,PLASTIC_TREAD) == 0) {
    seg_type = "3873.DAT";

    part_orient[0][0] = 0;
    part_orient[0][1] = 1;
    part_orient[0][2] = 0;
    part_orient[1][0] = 0;
    part_orient[1][1] = 0;
    part_orient[1][2] = -1;
    part_orient[2][0] = -1;
    part_orient[2][1] = 0;
    part_orient[2][2] = 0;
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 16;
  } else if (strcmp(type,RUBBER_TREAD) == 0) {
    seg_type = "680.DAT";
    part_orient[0][0] = 0;
    part_orient[0][1] = -1;
    part_orient[0][2] = 0;
    part_orient[1][0] = 1;
    part_orient[1][1] = 0;
    part_orient[1][2] = 0;
    part_orient[2][0] = 0;
    part_orient[2][1] = 0;
    part_orient[2][2] = 1;
    offset[0] = 0;
    offset[1] = 32;
    offset[2] = 0;
  }

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

      if (strcmp(type,RUBBER_BAND) == 0) {
        n = L2*SCALE;
      } else if (strcmp(type,CHAIN) == 0 || strcmp(type,PLASTIC_TREAD) == 0) {
        n = L2*CHAIN_SCALE;
      } else if (strcmp(type,RUBBER_TREAD) == 0) {
        n = L2*TREAD_SCALE;
      }

      // rotate the orientation of tangent line and part_orient

      mat_mult(rot,f_constraint->start_line.orient,part_orient);

      // rotate the offset by

      rotate_point3(roffset,offset,rot);

      for (i = 0; i < n; i++) {
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

        mat_mult(trot,absolute->orient,rot);

        fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          ghost ? "0 GHOST " : "",
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
          seg_type);
      }
    }
  }

  /* now for the arc */

  if (strcmp(type,RUBBER_BAND) == 0) {
    seg_type = "LS31.dat";
    n = 2*pi*f_constraint->radius*SCALE;
    angle = f_constraint->s_angle;
  } else if (strcmp(type,CHAIN) == 0) {
    seg_type = "3711.DAT";
    n = 2*pi*f_constraint->radius*CHAIN_SCALE;
    part_orient[0][0] = 0;
    part_orient[0][1] = -1;
    part_orient[0][2] = 0;
    part_orient[1][0] = 0;
    part_orient[1][1] = 0;
    part_orient[1][2] = 1;
    part_orient[2][0] = -1;
    part_orient[2][1] = 0;
    part_orient[2][2] = 0;
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;

    if (strcmp(f_constraint->part.type,"2855.DAT") == 0) {
      offset[0] = -8;
      offset[1] = 0;
      offset[2] = 0;
    }

    angle = f_constraint->s_angle;
  } else if (strcmp(type,PLASTIC_TREAD) == 0) {
    seg_type = "3873.DAT";
    n = 2*pi*f_constraint->radius*CHAIN_SCALE;
    part_orient[0][0] = 0;
    part_orient[0][1] = -1;
    part_orient[0][2] = 0;
    part_orient[1][0] = 0;
    part_orient[1][1] = 0;
    part_orient[1][2] = 1;
    part_orient[2][0] = -1;
    part_orient[2][1] = 0;
    part_orient[2][2] = 0;
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;

    if (strcmp(f_constraint->part.type,"2855.DAT") == 0) {
      offset[0] = -8;
      offset[1] = 0;
      offset[2] = 0;
    }

    angle = f_constraint->s_angle;
  } else if (strcmp(type,RUBBER_TREAD) == 0) {
    seg_type = "682.DAT";
    n = 2*pi*f_constraint->radius*CHAIN_SCALE;
    part_orient[0][0] = -0.342;
    part_orient[0][1] = 0.940;
    part_orient[0][2] = 0;
    part_orient[1][0] = -0.940;
    part_orient[1][1] = -0.342;
    part_orient[1][2] = 0;
    part_orient[2][0] = 0;
    part_orient[2][1] = 0;
    part_orient[2][2] = 1;

    offset[0] = 0;
    offset[1] = 32;
    offset[2] = 0;

    angle = f_constraint->s_angle;
  }

  f[0] = f_constraint->radius * cos(angle);
  f[1] = f_constraint->radius * sin(angle);
  f[2] = 0;

  for (i = 1; i < f_constraint->n_steps; i++) {
    PRECISION loc[3];
    PRECISION floc[3];
    PRECISION trot[3][3];

    s[0] = f_constraint->radius * cos(angle + 2*pi*i/n);
    s[1] = f_constraint->radius * sin(angle + 2*pi*i/n);
    s[2] = 0;

    //rotate(rf,f,f_constraint->part.orient);
    //rotate(rs,s,f_constraint->part.orient);

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

    mat_mult(rot,orient,part_orient);

    // rotate the offset by

    rotate_point3(roffset,offset,rot);

    floc[0] = f[0] + f_constraint->part.loc.x - roffset[0];
    floc[1] = f[1] + f_constraint->part.loc.y - roffset[1];
    floc[2] = f[2] + f_constraint->part.loc.z - roffset[2];

    rotate_point3(loc, floc,absolute->orient);

    loc[0] += absolute->loc.x;
    loc[1] += absolute->loc.y;
    loc[2] += absolute->loc.z;

    mat_mult(trot,absolute->orient,rot);

    fprintf(output,"%s1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
      ghost ? "0 GHOST " : "",
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
      seg_type);

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
  //PRECISION pi = 2*atan2(1,0);

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
      //cross ^= 1;
      //was_cross = 1;
      inside ^= 1;

    /* Pulleys/Pins */

    } else if (strcmp(constraints[i].part.type,"3736.DAT") == 0) {
      constraints[i].radius = 44;
    } else if (strcmp(constraints[i].part.type,"4185.DAT") == 0) {
      constraints[i].radius = 30;
    } else if (strcmp(constraints[i].part.type,"2983.DAT") == 0) {
      constraints[i].radius = 11;
    } else if (strcmp(constraints[i].part.type,"4265C.DAT") == 0 ||
               strcmp(constraints[i].part.type,"4265A.DAT") == 0 ||
               strcmp(constraints[i].part.type,"4265B.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3713.DAT") == 0) {
      constraints[i].radius = 10;
    } else if (strcmp(constraints[i].part.type,"2736.DAT") == 0 ||
               strcmp(constraints[i].part.type,"50.DAT") == 0) {
      constraints[i].radius = 4;

    /* Axles */
    } else if (strcmp(constraints[i].part.type,"3704.DAT") == 0 ||
               strcmp(constraints[i].part.type,"32062.DAT") == 0 ||
               strcmp(constraints[i].part.type,"4519.DAT") == 0 ||
               strcmp(constraints[i].part.type,"6587.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3705.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3705C01.DAT") == 0 ||
               strcmp(constraints[i].part.type,"32073.DAT") == 0 ||
               strcmp(constraints[i].part.type,"552.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3706.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3707.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3737.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3737C01.DAT") == 0 ||
               strcmp(constraints[i].part.type,"3708.DAT") == 0) {

      constraints[i].radius = 4;

    /* tread parts */

    } else if (strcmp(constraints[i].part.type,"32007.DAT") == 0) {
      constraints[i].radius = 30;

    /* gears */

    } else if (strcmp(constraints[i].part.type,"73071.DAT") == 0) {
      constraints[i].radius = 35;
    } else if (strcmp(constraints[i].part.type,"6573.DAT") == 0) {
      constraints[i].radius = 30;
    } else if (strcmp(constraints[i].part.type,"4019.DAT") == 0) {
      constraints[i].radius = 20;
    } else if (strcmp(constraints[i].part.type,"6542.DAT") == 0) {
      constraints[i].radius = 20;
    } else if (strcmp(constraints[i].part.type,"3648.DAT") == 0) {
      constraints[i].radius = 30;
    } else if (strcmp(constraints[i].part.type,"60C01.DAT") == 0) {
      constraints[i].radius = 30;
    } else if (strcmp(constraints[i].part.type,"3650A.DAT") == 0) {
      constraints[i].radius = 30;
    } else if (strcmp(constraints[i].part.type,"3649.DAT") == 0) {
      constraints[i].radius = 50;
    } else if (strcmp(constraints[i].part.type,"2855.DAT") == 0) {
      constraints[i].radius = 73;
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

#if 0
  {
    /* 0 1 2 no
       0 2 1 (closer)
       1 0 2 no
       1 2 0 no
       2 0 1 (closer)
       2 1 0
     */
    int i = 0, j = 2, k = 1;
    double sy = sqrt(absolute.orient[i][j]*absolute.orient[i][j] +
                     absolute.orient[i][k]*absolute.orient[i][k]);
    if (sy > 16*FLT_EPSILON) {
      xangle = atan2(absolute.orient[i][j], absolute.orient[i][k]);
      zangle = atan2(sy, absolute.orient[i][i]);
      yangle = atan2(absolute.orient[j][i], -absolute.orient[k][i]);
    } else {
      xangle = atan2(-absolute.orient[j][k], absolute.orient[j][j]);
      zangle = atan2(sy, absolute.orient[i][i]);
      yangle = 0;
    }
  }
#endif

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

  //mat_mult(absolute.orient,xrot,yrot);

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
          constraints[i].crossings[j].loc.z = (constraints[i].layer-layer/2)*BAND_DIAM + layer_offset;
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
      calc_angles(&constraints[i], type, output);
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
            type,
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
      if (j == n_constraints) {
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






