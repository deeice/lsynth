/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynthcp.h"
#include "hose.h"
#include "curve.h"
#include "mathlib.h"

/*
 * a 1x1 brick is 20 LDU wide and 24 LDU high
 *
 * hose length = 14 brick widths long = 280 LDU
 * number of ribs = 45
 * 6.2 LDU per rib
 *
 */

void orient(
  int     n_segments,
  part_t *segments)
{
  int i;

  for (i = 0; i < n_segments-1; i++) {

      PRECISION dx, dy, dz;
      PRECISION orient[3][3];
      PRECISION L1,L2;

      // determine the orientation of the part in the XY plane

      dx = segments[i+1].offset[0] - segments[i].offset[0];
      dy = segments[i+1].offset[1] - segments[i].offset[1];
      dz = segments[i+1].offset[2] - segments[i].offset[2];

      L1 = sqrt(dx*dx + dy*dy);
      L2 = sqrt(dx*dx + dy*dy + dz*dz);
      if (L1 == 0) {

        segments[i].orient[0][0] = 1;
        segments[i].orient[1][0] = 0;
        segments[i].orient[2][0] = 0;
        segments[i].orient[0][1] = 0;
        segments[i].orient[1][1] = 1;
        segments[i].orient[2][1] = 0;
        segments[i].orient[0][2] = 0;
        segments[i].orient[1][2] = 0;
        segments[i].orient[2][2] = 1;
      } else {
        segments[i].orient[0][0] =  dy/L1;  //  cos
        segments[i].orient[1][0] = -dx/L1;  //  sin
        segments[i].orient[2][0] =  0;
        segments[i].orient[0][1] =  dx/L2;  // -sin
        segments[i].orient[1][1] =  dy/L2;  //  cos
        segments[i].orient[2][1] =  dz/L2;
        segments[i].orient[0][2] = -dx*dz/(L1*L2);
        segments[i].orient[1][2] = -dy*dz/(L1*L2);
        segments[i].orient[2][2] =  L1/L2;
      }
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
  PRECISION up[3];

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

#ifdef DEBUG
    if (output) {
      fprintf(output,"0 XYZ (%f %f %f %f)\n",
        (1 - time) * (1 - time) * (1 - time) * start->offset[2],
        (1 - time) * (1 - time) * 3 * time * (start->offset[2] - start_speed_v[2]),
        (1 - time) * 3 * time * time * (end->offset[2] - stop_speed_v[2]),
        time * time * time * end->offset[2]);
    }
#endif
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
    i_time = 1.0/(n_segments*ratio);    /* G */

    foo = 1.0/n_segments;
    foo /= ratio;
    foo /= i_time_sum;

#ifdef DEBUG
    if (output) {
      fprintf(output,"0 XYZ %f %f %f PTP %f ITIME %f ",segments[i].offset[0],segments[i].offset[1],segments[i].offset[2],ptp,i_time_sum2);
    }
#endif
    i_time_sum2 += i_time;

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

#ifdef DEBUG
    if (output) {
      fprintf(output," XYZ %f %f %f \n",segments[i].offset[0],segments[i].offset[1],segments[i].offset[2]);
    }
#endif
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

  orient(n_segments, segments);

  return 0; /* it all worked */
}

