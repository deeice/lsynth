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

int
synth_curve(
  LSL_part_usage       *start,
  LSL_part_usage       *end,
  LSL_part_usage       *segments,
  int                   n_segments,
  LSL_tube_attributes  *attrib,
  FILE                 *output)
{
  PRECISION vector[3];
  PRECISION start_speed_v[3];
  PRECISION stop_speed_v[3];
  PRECISION time,time2,i_time,i_time_sum,i_time_sum2,n_time;
  PRECISION x,x2,y,y2,z,z2;
  PRECISION ptp,ratio;
  PRECISION ptp_sum;
  int i,j,n;
  PRECISION up[3];

  vector[0] = 0;
  vector[1] = attrib->exit_start;
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
  vector[1] = attrib->exit_stop;
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

    segments[i].loc.x =
      (1 - time) * (1 - time) * (1 - time) * start->loc.x +
      (1 - time) * (1 - time) * 3 * time * (start->loc.x - start_speed_v[0]) +
      (1 - time) * 3 * time * time * (end->loc.x - stop_speed_v[0]) +
            time * time * time * end->loc.x;

    segments[i].loc.y =
      (1 - time) * (1 - time) * (1 - time) * start->loc.y +
      (1 - time) * (1 - time) * 3 * time * (start->loc.y - start_speed_v[1]) +
      (1 - time) * 3 * time * time * (end->loc.y - stop_speed_v[1]) +
            time * time * time * end->loc.y;
/*
=(1-$A8)^3*D$4 + (1-$A8)^2*3*$A8*(D$4-D$3) +  (1-$A8)*3*$A8^2*(D$5-D$6) + $A8^3*D$5
 */
    segments[i].loc.z =
      (1 - time) * (1 - time) * (1 - time) * start->loc.z +
      (1 - time) * (1 - time) * 3 * time * (start->loc.z - start_speed_v[2]) +
      (1 - time) * 3 * time * time * (end->loc.z - stop_speed_v[2]) +
            time * time * time * end->loc.z;

#ifdef DEBUG
    if (output) {
      fprintf(output,"0 XYZ (%f %f %f %f)\n",
        (1 - time) * (1 - time) * (1 - time) * start->loc.z,
        (1 - time) * (1 - time) * 3 * time * (start->loc.z - start_speed_v[2]),
        (1 - time) * 3 * time * time * (end->loc.z - stop_speed_v[2]),
        time * time * time * end->loc.z);
    }
#endif
  }

  ptp_sum = 0;

  for (i = 0; i < n_segments - 1; i++) {
    x = segments[i + 1].loc.x - segments[i].loc.x;
    y = segments[i + 1].loc.y - segments[i].loc.y;
    z = segments[i + 1].loc.z - segments[i].loc.z;

    ptp = sqrt(x*x + y*y + z*z);
    ptp_sum += ptp;
  }

  /* G8 */

  i_time_sum = 0;
  for (i = 0; i < n_segments - 1; i++) {
    time  = (PRECISION) (i+1)/ (PRECISION) n_segments;

    x = segments[i + 1].loc.x - segments[i].loc.x;
    y = segments[i + 1].loc.y - segments[i].loc.y;
    z = segments[i + 1].loc.z - segments[i].loc.z;

    ptp = sqrt(x*x + y*y + z*z);

    ratio = ptp*n_segments/ptp_sum;
    i_time = 1.0/(n_segments*ratio);
    i_time_sum += i_time;
  }

  /* H, I, J */
  n_time = 0;
  i_time_sum2 = 0;
  for (i = 0; i < n_segments - 1; i++) {
    PRECISION foo;

    x = segments[i + 1].loc.x - segments[i].loc.x;
    y = segments[i + 1].loc.y - segments[i].loc.y;
    z = segments[i + 1].loc.z - segments[i].loc.z;

    ptp = sqrt(x * x + y * y + z * z);  /* E */
    ratio = ptp*n_segments/ptp_sum;     /* F */
    i_time = 1.0/(n_segments*ratio);    /* G */

    foo = 1.0/n_segments;
    foo /= ratio;
    foo /= i_time_sum;

#ifdef DEBUG
    if (output) {
      fprintf(output,"0 XYZ %f %f %f PTP %f ITIME %f ",segments[i].loc.x,segments[i].loc.y,segments[i].loc.z,ptp,i_time_sum2);
    }
#endif
    i_time_sum2 += i_time;

    segments[i].loc.x =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->loc.x +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->loc.x - start_speed_v[0]) +
      (1 - n_time) * 3 * n_time * n_time * (end->loc.x - stop_speed_v[0]) +
       n_time * n_time * n_time * end->loc.x;

    segments[i].loc.y =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->loc.y +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->loc.y - start_speed_v[1]) +
      (1 - n_time) * 3 * n_time * n_time * (end->loc.y - stop_speed_v[1]) +
       n_time * n_time * n_time * end->loc.y;

    segments[i].loc.z =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->loc.z +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->loc.z - start_speed_v[2]) +
      (1 - n_time) * 3 * n_time * n_time * (end->loc.z - stop_speed_v[2]) +
       n_time * n_time * n_time * end->loc.z;

    n_time += foo;  /* H */

#ifdef DEBUG
    if (output) {
      fprintf(output," XYZ %f %f %f \n",segments[i].loc.x,segments[i].loc.y,segments[i].loc.z);
    }
#endif
  }

  segments[i].loc.x =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->loc.x +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->loc.x - start_speed_v[0]) +
    (1 - n_time) * 3 * n_time * n_time * (end->loc.x - stop_speed_v[0]) +
     n_time * n_time * n_time * end->loc.x;

  segments[i].loc.y =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->loc.y +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->loc.y - start_speed_v[1]) +
    (1 - n_time) * 3 * n_time * n_time * (end->loc.y - stop_speed_v[1]) +
     n_time * n_time * n_time * end->loc.y;

  segments[i].loc.z =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->loc.z +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->loc.z - start_speed_v[2]) +
    (1 - n_time) * 3 * n_time * n_time * (end->loc.z - stop_speed_v[2]) +
     n_time * n_time * n_time * end->loc.z;

  up[0] = 1;
  up[1] = 0;
  up[2] = 0;

  rotate_point(up,start->orient);

  for (i = 0; i < n_segments-1; i++) {
    PRECISION r;
    PRECISION front[3];
    PRECISION t[3];

    front[0] = segments[i+1].loc.x - segments[i].loc.x;
    front[1] = segments[i+1].loc.y - segments[i].loc.y;
    front[2] = segments[i+1].loc.z - segments[i].loc.z;

    r = sqrt(front[0]*front[0] + front[1]*front[1] + front[2]*front[2]);

    front[0] /= r;
    front[1] /= r;
    front[2] /= r;

    mult_point(t,front,up);
    mult_point(up,t,front);

    make_rotation_pp(segments[i].orient,up,front);
  }
  return 0; /* it all worked */
}






