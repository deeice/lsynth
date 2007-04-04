/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynthcp.h"
#include "hose.h"
#include "curve.h"
#include "mathlib.h"

static void
mult_point(PRECISION r[3], PRECISION lhs[3], PRECISION rhs[3])
{
  PRECISION tt;

  r[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
  r[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
  r[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];

  tt = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

  r[0] /= tt;
  r[1] /= tt;
  r[2] /= tt;
}

static void
make_rotation_pp(
  PRECISION r[3][3],
  PRECISION up[3],
  PRECISION front[3])
{
  PRECISION side[3];

  mult_point(side,front,up);

  r[0][0] = up[0];
  r[1][0] = up[1];
  r[2][0] = up[2];
  r[0][1] = front[0];
  r[1][1] = front[1];
  r[2][1] = front[2];
  r[0][2] = side[0];
  r[1][2] = side[1];
  r[2][2] = side[2];
}

static void
rotate_point(
  PRECISION r[3],
  PRECISION m[3][3])
{
  PRECISION t[3];

  t[0] = r[0]*m[0][0] + r[1]*m[0][1] + r[2]*m[0][2];
  t[1] = r[0]*m[1][0] + r[1]*m[1][1] + r[2]*m[1][2];
  t[2] = r[0]*m[2][0] + r[1]*m[2][1] + r[2]*m[2][2];

  r[0] = t[0];
  r[1] = t[1];
  r[2] = t[2];
}

void
orient2(
  part_t       *start,
  part_t       *end,
  int           n_segments,
  part_t       *segments)
{
  PRECISION up[3];
  int i;

  up[0] = 1;
  up[1] = 0;
  up[2] = 0;

  rotate_point(up,start->orient);
  printf("%5.2f %5.2f %5.2f\n",up[0],up[1],up[2]);

  for (i = 0; i < n_segments-1; i++) {
    PRECISION r;
    PRECISION front[3];
    PRECISION t[3];

    front[0] = segments[i+1].offset[0] - segments[i].offset[0];
    front[1] = segments[i+1].offset[1] - segments[i].offset[1];
    front[2] = segments[i+1].offset[2] - segments[i].offset[2];

    r = sqrt(front[0]*front[0] + front[1]*front[1] + front[2]*front[2]);

    front[0] /= r;
    front[1] /= r;
    front[2] /= r;

    mult_point(t,front,up); // side
    mult_point(up,t,front); // side * front

    make_rotation_pp(segments[i].orient,up,front);
  }
}

PRECISION hose_length(
  int           n_segments,
  part_t       *segments)
{
  PRECISION length = 0;
  int i;

  for (i = 0; i < n_segments-1; i++) {
    PRECISION l[3];

    vectorsub3(l,segments[i+1].offset,segments[i].offset);

    length += vectorlen(l);
  }
  return length;
}

void
orient(
  part_t       *start,
  part_t       *end,
  int           n_segments,
  part_t       *segments)
{
  PRECISION start_up[3],end_up[3],up[3];
  PRECISION total_length = hose_length(n_segments,segments);
  PRECISION cur_length;
  int i;

  start_up[0] = 1;
  start_up[1] = 0;
  start_up[2] = 0;
  rotate_point(start_up,start->orient);

  end_up[0] = 1;
  end_up[1] = 0;
  end_up[2] = 0;
  rotate_point(end_up,end->orient);

  /* Up vector controls the twist
   *
   * Interpolate the up vector based on start up vector, and
   * end up vector, and how far we are down the hose's total
   * length
   */

  for (i = 0; i < n_segments-1; i++) {
    PRECISION r;
    PRECISION front[3];
    PRECISION t[3];

    cur_length = hose_length(i,segments);

    cur_length /= total_length;

    up[0] = start_up[0]*(1-cur_length) + end_up[0]*cur_length;
    up[1] = start_up[1]*(1-cur_length) + end_up[1]*cur_length;
    up[2] = start_up[2]*(1-cur_length) + end_up[2]*cur_length;

    r = sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
    up[0] /= r;
    up[1] /= r;
    up[2] /= r;

    front[0] = segments[i+1].offset[0] - segments[i].offset[0];
    front[1] = segments[i+1].offset[1] - segments[i].offset[1];
    front[2] = segments[i+1].offset[2] - segments[i].offset[2];

    r = sqrt(front[0]*front[0] + front[1]*front[1] + front[2]*front[2]);

    front[0] /= r;
    front[1] /= r;
    front[2] /= r;

    mult_point(t,front,up); // side
    mult_point(up,t,front); // side * front

    make_rotation_pp(segments[i].orient,up,front);
  }
}

int
synth_curve(
  part_t       *start,
  part_t       *end,
  part_t       *segments,
  int           n_segments,
  PRECISION     attrib,
  FILE         *output)
{
  PRECISION vector[3];
  PRECISION start_speed_v[3];
  PRECISION stop_speed_v[3];
  PRECISION time,time2,i_time,i_time_sum,i_time_sum2,n_time;
  PRECISION x,x2,y,y2,z,z2;
  PRECISION ptp,ratio;
  PRECISION ptp_sum;
  int i,j,n;

  vector[0] = 0;
  vector[1] = attrib;
  vector[2] = 0;

  for (j = 0; j < 3; j++)
  {
     x = 0.0;
     for (i = 0; i < 3; i++) {
        x += start->orient[j][i] * vector[i];
     }
     start_speed_v[j] = x;
  }

  vector[0] = 0;
  vector[1] = attrib;
  vector[2] = 0;

  for (j = 0; j < 3; j++)
  {
     x = 0.0;
     for (i = 0; i < 3; i++) {
        x -= end->orient[j][i] * vector[i];
     }
     stop_speed_v[j] = x;
  }

  // stop_speed_v[2] = -stop_speed_v[2];

  for (i = 0; i < n_segments; i++) {
    time  = (PRECISION) i/ (PRECISION) n_segments;

    segments[i].offset[0] =
      (1 - time) * (1 - time) * (1 - time) * start->offset[0] +
      (1 - time) * (1 - time) * 3 * time * (start->offset[0] - start_speed_v[0]) +
      (1 - time) * 3 * time * time * (end->offset[0] - stop_speed_v[0]) +
            time * time * time * end->offset[0];

    segments[i].offset[1] =
      (1 - time) * (1 - time) * (1 - time) * start->offset[1] +
      (1 - time) * (1 - time) * 3 * time * (start->offset[1] - start_speed_v[1]) +
      (1 - time) * 3 * time * time * (end->offset[1] - stop_speed_v[1]) +
            time * time * time * end->offset[1];
/*
=(1-$A8)^3*D$4 + (1-$A8)^2*3*$A8*(D$4-D$3) +  (1-$A8)*3*$A8^2*(D$5-D$6) + $A8^3*D$5
 */
    segments[i].offset[2] =
      (1 - time) * (1 - time) * (1 - time) * start->offset[2] +
      (1 - time) * (1 - time) * 3 * time * (start->offset[2] - start_speed_v[2]) +
      (1 - time) * 3 * time * time * (end->offset[2] - stop_speed_v[2]) +
            time * time * time * end->offset[2];
  }

  ptp_sum = 0;

  for (i = 0; i < n_segments - 1; i++) {
    x = segments[i + 1].offset[0] - segments[i].offset[0];
    y = segments[i + 1].offset[1] - segments[i].offset[1];
    z = segments[i + 1].offset[2] - segments[i].offset[2];

    ptp = sqrt(x*x + y*y + z*z);
    ptp_sum += ptp;
  }

  /* G8 */

  i_time_sum = 0;
  for (i = 0; i < n_segments - 1; i++) {
    time  = (PRECISION) (i+1)/ (PRECISION) n_segments;

    x = segments[i + 1].offset[0] - segments[i].offset[0];
    y = segments[i + 1].offset[1] - segments[i].offset[1];
    z = segments[i + 1].offset[2] - segments[i].offset[2];

    ptp = sqrt(x*x + y*y + z*z);

    ratio = ptp*n_segments/ptp_sum;
    if (ratio == 0) {
      ratio = 1e-20;
    }
    i_time = 1.0/(n_segments*ratio);
    i_time_sum += i_time;
  }

  /* H, I, J */
  n_time = 0;
  i_time_sum2 = 0;
  for (i = 0; i < n_segments - 1; i++) {
    PRECISION foo;

    x = segments[i + 1].offset[0] - segments[i].offset[0];
    y = segments[i + 1].offset[1] - segments[i].offset[1];
    z = segments[i + 1].offset[2] - segments[i].offset[2];

    ptp = sqrt(x * x + y * y + z * z);  /* E */
    ratio = ptp*n_segments/ptp_sum;     /* F */
    if (ratio == 0) {
      ratio = 1e-20;
    }
    i_time = 1.0/(n_segments*ratio);    /* G */
    i_time_sum2 += i_time;

    foo = 1.0/n_segments;
    foo /= ratio;
    foo /= i_time_sum;

    segments[i].offset[0] =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[0] +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[0] - start_speed_v[0]) +
      (1 - n_time) * 3 * n_time * n_time * (end->offset[0] - stop_speed_v[0]) +
       n_time * n_time * n_time * end->offset[0];

    segments[i].offset[1] =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[1] +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[1] - start_speed_v[1]) +
      (1 - n_time) * 3 * n_time * n_time * (end->offset[1] - stop_speed_v[1]) +
       n_time * n_time * n_time * end->offset[1];

    segments[i].offset[2] =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[2] +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[2] - start_speed_v[2]) +
      (1 - n_time) * 3 * n_time * n_time * (end->offset[2] - stop_speed_v[2]) +
       n_time * n_time * n_time * end->offset[2];

    n_time += foo;  /* H */
  }

  segments[i].offset[0] =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[0] +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[0] - start_speed_v[0]) +
    (1 - n_time) * 3 * n_time * n_time * (end->offset[0] - stop_speed_v[0]) +
     n_time * n_time * n_time * end->offset[0];

  segments[i].offset[1] =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[1] +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[1] - start_speed_v[1]) +
    (1 - n_time) * 3 * n_time * n_time * (end->offset[1] - stop_speed_v[1]) +
     n_time * n_time * n_time * end->offset[1];

  segments[i].offset[2] =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[2] +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[2] - start_speed_v[2]) +
    (1 - n_time) * 3 * n_time * n_time * (end->offset[2] - stop_speed_v[2]) +
     n_time * n_time * n_time * end->offset[2];

  // orient(n_segments, segments);

  return 0; /* it all worked */
}

