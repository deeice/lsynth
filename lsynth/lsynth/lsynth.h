/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef LSYNTH_H
#define LSYNTH_H

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include <string.h>

#ifndef PRECISION
#define PRECISION double
#endif

#define RIBBED_HOSE       "RIBBED_HOSE "
#define RUBBER_HOSE       "RUBBER_HOSE "
#define STRING            "STRING "
#define MINIFIG_CHAIN     "MINIFIG_CHAIN "
#define FLEXIBLE_HOSE     "FLEX_SYSTEM_HOSE_LD "
#define FLEX_CABLE        "FLEX_SYSTEM_CABLE "
#define RIGID_HOSE        "FLEX_SYSTEM_HOSE "
#define ELECTRIC_CABLE    "ELECTRIC_CABLE "
#define PNEUMATIC_HOSE    "PNEUMATIC_HOSE "
#define FLEXIBLE_AXLE     "FLEXIBLE_AXLE "
#define FIBER_OPTIC_CABLE "FIBER_OPTIC_CABLE "

#define ALL_HOSES RIBBED_HOSE,FLEX_CABLE,FLEXIBLE_HOSE,RIGID_HOSE,ELECTRIC_CABLE,PNEUMATIC_HOSE,FLEXIBLE_AXLE,FIBER_OPTIC_CABLE,RUBBER_HOSE,STRING,MINIFIG_CHAIN

#define RUBBER_BAND       "RUBBER_BAND "
#define RUBBER_BELT       "RUBBER_BELT "
#define CHAIN             "CHAIN "
#define PLASTIC_TREAD     "PLASTIC_TREAD "
#define RUBBER_TREAD      "RUBBER_TREAD "

#define ALL_BANDS RUBBER_BAND,RUBBER_BELT,CHAIN,PLASTIC_TREAD,RUBBER_TREAD

#define ACCY (1e-6)

typedef struct {
  PRECISION x,y,z;
} LSL_3D;

typedef struct {
  int color;
  char type[128];
  LSL_3D loc;
  PRECISION orient[3][3];
} LSL_part_usage;

typedef struct {
  LSL_part_usage part;
  PRECISION      radius;
  PRECISION      offset[3];
  int            inside;
  int            was_cross;
  int            cross;

  LSL_part_usage start_line;  /* from next to previous */
  LSL_part_usage end_line;
  LSL_part_usage crossings[8];
  int            n_crossings;
  int            layer;

  LSL_part_usage start_angle;
  LSL_part_usage end_angle;
  int            n_steps;
  PRECISION      s_angle;
} LSL_band_constraint;

extern PRECISION hose_res_angle;

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
  char            *type);


#endif