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

  tt = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
 
  if (tt > 0.0)
  {
    tt = sqrt(tt);

    r[0] /= tt;
    r[1] /= tt;
    r[2] /= tt;
  }
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





/***************************************************************/
// Quaternion q[4] is mostly a direction vector with a spin(roll) angle.
// ie. quaternion(x,y,z,spin)
// 
// Gotta compare quat FAQ conversions to ldglite conversions to
// see what to do with the weird ldraw coordinate system.  I think I 
// may have to reverse -y (and maybe z) to get right? handed system
// before using quaternions.  Make sure I copied all conversions from
// same place to avoid mixed systems.

/***************************************************************/

/***************************************************************/
int normalizequat(PRECISION q[4])
{
  PRECISION L;

  L = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);

  if (L == 0.0)
    return 0;  // Uh Oh, the magnitude is zero.

  q[0] /= L;
  q[1] /= L;
  q[2] /= L;
  q[3] /= L;

  return 1;
}

/***************************************************************/
int normalize(PRECISION v[3])
{
  PRECISION L;

  L = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  if (L == 0.0)
    return 0;  // Uh Oh, the magnitude is zero.

  v[0] /= L;
  v[1] /= L;
  v[2] /= L;

  return 1;
}

/***************************************************************/
/* returns det(m), m is 3x3 (3x4) matrix */
/***************************************************************/
PRECISION M3Det(PRECISION m[3][4]) /* Note argument type !            */
{
   return (m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) +
           m[1][0] * (m[2][1] * m[0][2] - m[0][1] * m[2][2]) +
           m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]));
}
/***************************************************************/
// Convert axis v[3] and angle a to quaternion q[4].
/***************************************************************/
void axisangle2quat(
  PRECISION v[3], 
  PRECISION a,
  PRECISION q[4])
{
  PRECISION sina,cosa;
  
  a = a / 2.0;
  sina = sin(a);
  cosa = cos(a);
  q[0] = v[0] * sina;
  q[1] = v[1] * sina;
  q[2] = v[2] * sina;
  q[3] = cosa;
}

/***************************************************************/
// Convert rotation matrix m[3][3] to quaternion q[4].
// Beware!  Quaternions cannot handle mirror matrices.
// This may be a problem because the ldraw coordinate system is mirrored?
/***************************************************************/
void matrix2quat(
  PRECISION m[3][3],
  PRECISION q[4])
{
  PRECISION T;
  PRECISION S;

  // NOTE:  First unmirror any mirrored matrix.
  // (Multiplying the rotation matrix with the sign of its determinant).
  PRECISION D = M3Det(m);
  if (D < 0.0)
  {
    // TODO: Unmirror
  }

#ifdef SPEEDY_WAY
#define max(a,b) ((a > b) ? (a) : (b))

  q[0] = sqrt( max( 0, 1 + m[0][0] - m[1][1] - m[2][2] ) ) / 2; 
  q[1] = sqrt( max( 0, 1 - m[0][0] + m[1][1] - m[2][2] ) ) / 2; 
  q[2] = sqrt( max( 0, 1 - m[0][0] - m[1][1] + m[2][2] ) ) / 2; 
  q[3] = sqrt( max( 0, 1 + m[0][0] + m[1][1] + m[2][2] ) ) / 2; 
  q[0] = _copysign( q[0], m[2][1] - m[1][2] );
  q[1] = _copysign( q[1], m[0][2] - m[2][0] );
  q[2] = _copysign( q[2], m[1][0] - m[0][1] );
  return;
#endif

  // Calculate the trace of the matrix T from the equation:
  T = 1 + m[0][0] + m[1][1] + m[2][2];
  if ( T > 0.00000001 )
  {
    S = sqrt(T) * 2;
    q[0] = ( m[2][1] - m[1][2] ) / S;
    q[1] = ( m[0][2] - m[2][0] ) / S;
    q[2] = ( m[1][0] - m[0][1] ) / S;
    q[3] = 0.25 * S;
  }
  // If the trace of the matrix is equal to zero then identify
  // which major diagonal element has the greatest value.
  // Depending on this, calculate the following:
  else if ( (m[0][0] > m[1][1]) && (m[0][0] > m[2][2]) )  {	// Column 0: 
    S  = sqrt( 1.0 + m[0][0] - m[1][1] - m[2][2] ) * 2;
    q[0] = 0.25 * S;
    q[1] = ( m[1][0] + m[0][1] ) / S;
    q[2] = ( m[0][2] + m[2][0] ) / S;
    q[3] = ( m[1][2] - m[2][1] ) / S; //    q[3] = ( m[2][1] - m[1][2] ) / S;
  } else if ( m[1][1] > m[2][2] ) {			// Column 1: 
    S  = sqrt( 1.0 + m[1][1] - m[0][0] - m[2][2] ) * 2;
    q[0] = ( m[1][0] + m[0][1] ) / S;
    q[1] = 0.25 * S;
    q[2] = ( m[2][1] + m[1][2] ) / S;
    q[3] = ( m[0][2] - m[2][0] ) / S;
  } else {						// Column 2:
    S  = sqrt( 1.0 + m[2][2] - m[0][0] - m[1][1] ) * 2;
    q[0] = ( m[0][2] + m[2][0] ) / S;
    q[1] = ( m[2][1] + m[1][2] ) / S;
    q[2] = 0.25 * S;
    q[3] = ( m[0][1] - m[1][0] ) / S; //    q[3] = ( m[1][0] - m[0][1] ) / S;
  }
}

/***************************************************************/
// Convert quaternion q[4] to rotation matrix m[3][3].
/***************************************************************/
void quat2matrix(
  PRECISION q[4],
  PRECISION m[3][3])
{
  PRECISION a,b,c,s;

  normalizequat( q );

  a = q[0];
  b = q[1];
  c = q[2];
  s = q[3];

  m[0][0]  = 1 - 2*b*b-2*c*c;
  m[0][1]  = 2*a*b - 2*s*c;
  m[0][2]  = 2*a*c + 2*s*b;
  m[1][0]  = 2*a*b+2*s*c;
  m[1][1]  = 1 - 2*a*a - 2*c*c;
  m[1][2]  = 2*b*c - 2*s*a;
  m[2][0]  = 2*a*c - 2*s*b;
  m[2][1]  = 2*b*c + 2*s*a;
  m[2][2] = 1 - 2*a*a - 2*b*b;
}

/***************************************************************/
// Convert quaternion q[4] to a rotation axis v[3] and angle a.
/***************************************************************/
 
PRECISION quat2axisangle(
  PRECISION q[4],
  PRECISION v[3],
  PRECISION *a)
{
  PRECISION cos_angle, sin_angle;

  normalizequat( q );

  cos_angle  = q[3];
  *a          = acos( cos_angle ) * 2;
  sin_angle  = sqrt( 1.0 - cos_angle * cos_angle );

  if ( fabs( sin_angle ) < 0.0005 )
    sin_angle = 1;

  v[0] = q[0] / sin_angle;
  v[1] = q[1] / sin_angle;
  v[2] = q[2] / sin_angle;

  return *a;
}

//**********************************************************************
PRECISION dotprod(PRECISION a[3], PRECISION b[3])
{
   return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

//**********************************************************************
// Calculate the the turn matrix M, turn axis T, and return the turn angle r.
// In case of emergency preload M and T with reasonable defaults.
// 
// NOTE: If I build each orient from previous instead of from start, 
// then I don't need to pass in or use MS.  M is preloaded with previous M.
// To do this I should pass previous v for a (instead of always using start_v).
//**********************************************************************
PRECISION get_turn_mat(PRECISION M[3][3], PRECISION MS[3][3], 
		       PRECISION a[3], PRECISION b[3], PRECISION T[3])
{
  PRECISION m[3][3];
  PRECISION q[4];
  PRECISION t[3];
  PRECISION r;

    //Dot product gives turn angle.  a.b=|a||b|cos(theta)
    //We normalized so |a|=|b|=1, which means theta = acos(a.b).
    r = dotprod(b, a);
    r = acos(r);

    //Cross product gives turn axis.
    mult_point(t, a, b);

    // NOTE: gotta check cross product.  
    // If magnitude is zero we have either a 180 or 0 degree turn.
    // 
    // Dot product can differentiate between 0 and 180 turn.
    // Dot product of acute is negative (dot product of 180 = -1).
    // Dot product of perpendicular is 0.
    // Dot product of obtuse is positive (dot product of 0 = 1).
    //
    // If we get a 0 degree turn, just reuse previous orient matrix.
    // If we get a 180 degree turn, reuse the previous turn axis.
    // (Any axis in the perpendicular plane can make the 180, 
    // but we want to maintain the up vector, so reuse the previous axis.)

    // I think I can reasonably expect to avoid the 180 degree turns
    // if we apply small incremental turns to the previous orient matrix
    // instead of always applying the whole turn from the start matrix.
    // But do we lose control of the up vector if we do that?

    if (normalize(t) == 0)
    {
      if (dotprod(b, a) < 0)
      {
	printf("***  180 degree turn!!!\n");
	// Just use the up vector (or previous t) which we preloaded into T.
	vectorcp(t, T);
      }
      else 
      {
	printf("***  0 degree turn!!!\n");
	// Just reuse the previous orient, which we preloaded into M.
	return r;
      }
    }

    // Convert to a quaternion, and convert that to a rotation matrix.
    axisangle2quat(t, r, q);
    normalizequat( q );
    quat2matrix(q, m);

    // Then combine the rotation from that matrix with the start orient.
    matrixmult3(M,m,MS);

    vectorcp(T, t);
    return r;
}

//**********************************************************************
// Convert start, end constraints and n points to oriented segments
// using quaternion.
//**********************************************************************
void
orientq(
  part_t       *start,
  part_t       *end,
  int           n_segments,
  part_t       *segments)
{
  PRECISION start_up[3],end_up[3],up[3];
  PRECISION start_v[3],end_v[3],v[3];
  PRECISION start_q[4],end_q[4],q[4];
  PRECISION m[3][3];
  PRECISION total_length = hose_length(n_segments,segments);
  PRECISION cur_length;
  PRECISION r;
  PRECISION front[3];
  PRECISION t[3];
  int i;

  // Get the start and end velocity vectors.
  // The basic HoseSeg.dat files all extend the hose along the -Y axis.
#if 1
  start_v[0] = end_v[0] = 0;
  start_v[1] = end_v[1] = -1; 
  start_v[2] = end_v[2] = 0;
#else
  start_v[0] = end_v[0] = 0;
  start_v[1] = end_v[1] = 1; // Try positive Y.  Maybe it'll fix offset.
  start_v[2] = end_v[2] = 0;
  // It does fix it.  Why?  Does some other matrix flip Y in hose.c?
  // Or is it just that I no longer need the 10 offset in lsynth.mpd?
  //
  // I substituted a constraint for one of the chain links and it points
  // backwards with this code.  So I want -Y and fix the lsynth.mpd offset.
  //
  // Run the original code and substitute a constraint to see if it flips Y.
  // Yes, they all point backwards.  
  // Why?  Is it just my LENGTH21 code in hose.c?  Or is it all hoses?
#endif
  rotate_point(start_v,start->orient);
  rotate_point(end_v,end->orient);
  normalize(start_v); // Normalize again,
  normalize(end_v);   // Just in case orient includes a stretch.

  // Get the up vectors.  The fin on LS00.dat is along -Z so rotate that.
  start_up[0] = end_up[0] = 0;
  start_up[1] = end_up[1] = 0;
  start_up[2] = end_up[2] = -1;
  rotate_point(start_up,start->orient);
  rotate_point(end_up,end->orient);
  normalize(start_up); // Normalize again,
  normalize(end_up);   // Just in case orient includes a stretch.

#define DEBUG_QUAT_MATH 1
#ifdef DEBUG_QUAT_MATH
  {
    PRECISION pi = 2*atan2(1,0);
    PRECISION degrees = 180.0 / pi;
      
    printf("V[S] = (%.2fx, %.2fy, %.2fz)\n", start_v[0],start_v[1],start_v[2]);
    printf("M[S] = (%.2f, %.2f, %.2f,   %.2f, %.2f, %.2f,   %.2f, %.2f, %.2f)\n",
	 start->orient[0][0],
	 start->orient[0][1],
	 start->orient[0][2],
	 start->orient[1][0],
	 start->orient[1][1],
	 start->orient[1][2],
	 start->orient[2][0],
	 start->orient[2][1],
	 start->orient[2][2]);
  }
#endif

  /* Up vector controls the twist
   *
   * Interpolate the up vector based on start up vector, and
   * end up vector, and how far we are down the hose's total
   * length
   */

  for (i = 0; i < n_segments-1; i++) {

    cur_length = hose_length(i,segments);

    cur_length /= total_length;

    up[0] = start_up[0]*(1-cur_length) + end_up[0]*cur_length;
    up[1] = start_up[1]*(1-cur_length) + end_up[1]*cur_length;
    up[2] = start_up[2]*(1-cur_length) + end_up[2]*cur_length;
    normalize(up);

    v[0] = segments[i+1].offset[0] - segments[i].offset[0];
    v[1] = segments[i+1].offset[1] - segments[i].offset[1];
    v[2] = segments[i+1].offset[2] - segments[i].offset[2];
    normalize(v);

    // Preload segments[i].orient with previous orient in case of 0 degree turn.
    // Preload turn axis t with previous t (or up) in case of 180 degree turn.
    if (i)
    {
      matrixcp(segments[i].orient, segments[i-1].orient);
      vectorcp(t, start_up); // Should reuse the previous t.  Did I save it?
    }
    else 
    {
      matrixcp(segments[i].orient, start->orient);
      vectorcp(t, start_up); // Default t to the up vector
    }

    r = get_turn_mat(segments[i].orient, start->orient, start_v, v, t);

#ifdef DEBUG_QUAT_MATH
  {
    PRECISION pi = 2*atan2(1,0);
    PRECISION degrees = 180.0 / pi;
      
    printf("V[%d] = (%.2fx, %.2fy, %.2fz) TURN = %.2fdeg (%.2f, %.2f, %.2f)\n",
	   i, v[0],v[1],v[2], r*degrees, t[0],t[1],t[2]);
#if 0
    printf("M[%d] = (%.2f, %.2f, %.2f,   %.2f, %.2f, %.2f,   %.2f, %.2f, %.2f)\n",
	   i,
	 segments[i].orient[0][0],
	 segments[i].orient[0][1],
	 segments[i].orient[0][2],
	 segments[i].orient[1][0],
	 segments[i].orient[1][1],
	 segments[i].orient[1][2],
	 segments[i].orient[2][0],
	 segments[i].orient[2][1],
	 segments[i].orient[2][2]);
#endif
  }
#endif

  }

#ifdef DEBUG_QUAT_MATH
  {
    PRECISION pi = 2*atan2(1,0);
    PRECISION degrees = 180.0 / pi;
      
    printf("V[E] = (%.2fx, %.2fy, %.2fz)\n", end_v[0],end_v[1],end_v[2]);
    printf("M[E] = (%.2f, %.2f, %.2f,   %.2f, %.2f, %.2f,   %.2f, %.2f, %.2f)\n",
	 end->orient[0][0],
	 end->orient[0][1],
	 end->orient[0][2],
	 end->orient[1][0],
	 end->orient[1][1],
	 end->orient[1][2],
	 end->orient[2][0],
	 end->orient[2][1],
	 end->orient[2][2]);
  }
#endif

  // Now we should be able to add the spin
  // Rotate start_up by segments[n].orient
  // Use dot product to get the extra spin angle between that and end_up.
  // Spread out the spin (as a pre-rotation around Y) over the length of hose.
  // Actually, we can calculate this extra spin ahead of time
  // and do the pre-rotation inside the for loop we already have.

#ifdef DEBUG_QUAT_MATH
  printf("U[S] = (%.2fx, %.2fy, %.2fz)\n", start_up[0],start_up[1],start_up[2]);
  i = n_segments-2;
  printf("M[%d] = (%.2f, %.2f, %.2f,   %.2f, %.2f, %.2f,   %.2f, %.2f, %.2f)\n",
	   i,
	 segments[i].orient[0][0],
	 segments[i].orient[0][1],
	 segments[i].orient[0][2],
	 segments[i].orient[1][0],
	 segments[i].orient[1][1],
	 segments[i].orient[1][2],
	 segments[i].orient[2][0],
	 segments[i].orient[2][1],
	 segments[i].orient[2][2]);
#endif

  // I know.  I can't use segments[n-2].orient because that includes
  // the rotation in start->orient.  If I want to use that I have to
  // pass it the unmodified -Z up vector.  Not the rotated start_up.
  start_up[0] = 0;
  start_up[1] = 0;
  start_up[2] = -1;
  //matrixcp(m, segments[n_segments-2].orient);
  matrixcp(m, start->orient); // Reuse start->orient for 0 degree twist.
  vectorcp(t, start_up); // Should I reuse the previous t.  Did I save it?
  r = get_turn_mat(m, start->orient, start_v, end_v, t);
  rotate_point(start_up,m);
  normalize(start_up); // Normalize again

  // The rotated u[S] should match the U[E], unless there is some extra twist.
  // It doesn't match, so something is wrong.

#ifdef DEBUG_QUAT_MATH
  printf("u[S] = (%.2fx, %.2fy, %.2fz)\n", start_up[0],start_up[1],start_up[2]);
  printf("U[E] = (%.2fx, %.2fy, %.2fz)\n", end_up[0],end_up[1],end_up[2]);
#endif

  //Dot product gives turn angle.  a.b=|a||b|cos(theta)
  //We normalized so |a|=|b|=1, which means theta = acos(a.b).
  r = dotprod(start_up, end_up);
  r = acos(r);

#ifdef DEBUG_QUAT_MATH
  {
    PRECISION pi = 2*atan2(1,0);
    PRECISION degrees = 180.0 / pi;
      
    //Cross product gives turn axis.
    mult_point(t, start_up, end_up);
    normalize(t);

    printf("u[E] = (%.2fx, %.2fy, %.2fz) TWIST = %.2fdeg (%.2f, %.2f, %.2f)\n",
	   start_up[0],start_up[1],start_up[2], r*degrees, t[0],t[1],t[2]);

    // I'm still getting a small TURN around the side axis (-X) 
    // that must be eliminated before we can isolate the extra TWIST.
    // Probably the difference between the last segment we calculated
    // using the bezier curve() fn and the orient of the end constraint.
    // Can we force the bezier fn to make the start and end segments
    // match the orient of the start and end constraints.  The original
    // code just replaces the orient of the first and last segments.
    // That's OK for very short segments, but could give a rather sharp
    // Turn for the 2nd and 2nd to last segments if they have any 
    // length to them.
  }
#endif

  // Ok now we have the extra up vector twist angle in r.
  // Interpolate it over the length of the hose.
  // I'm not sure this is quite right, but it goes from 0 to 100%.
  // Probably ought to draw a pretty picture to be absolutely sure.
  total_length = hose_length(n_segments-1,segments);
  if (r != 0.0)
   for (i = 1; i < n_segments; i++) {
    PRECISION m1[3][3];
    PRECISION angle;

    cur_length = hose_length(i,segments);
    cur_length /= total_length;
    angle = r * cur_length;
    angle = -angle;  // Rotate the opposite way.

#ifdef DEBUG_QUAT_MATH
    {
    PRECISION pi = 2*atan2(1,0);
    PRECISION degrees = 180.0 / pi;
      
    printf("TWIST[%d] at length %.2f of %.2f = %.2fdeg \n", 
	   i, cur_length, total_length, angle*degrees);
    }
#endif

    m1[0][0] = 1;
    m1[0][1] = 0;
    m1[0][2] = 0;
    m1[1][0] = 0;
    m1[1][1] = 1;
    m1[1][2] = 0;
    m1[2][0] = 0;
    m1[2][1] = 0;
    m1[2][2] = 1;

    m1[0][0] =   cos(angle);
    m1[0][2] =   sin(angle);
    m1[2][0] =  -sin(angle);
    m1[2][2] =   cos(angle);

    matrixmult(segments[i-1].orient, m1);
  }

}

