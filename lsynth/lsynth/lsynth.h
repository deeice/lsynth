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

#define TEST_HOSE         "TEST_HOSE "
#define RIBBED_HOSE       "RIBBED_HOSE "
#define RIBBED_TUBE       "RIBBED_TUBE "
#define FLEXIBLE_HOSE     "FLEX_SYSTEM_HOSE_LD "
#define FLEXIBLE_TUBE     "FLEX_SYSTEM_TUBE_LD "
#define FLEX_CABLE        "FLEX_SYSTEM_CABLE "
#define RIGID_HOSE        "FLEX_SYSTEM_HOSE "
#define RIGID_TUBE        "FLEX_SYSTEM_TUBE "
#define ELECTRIC_CABLE    "ELECTRIC_CABLE "
#define PNEUMATIC_TUBE    "PNEUMATIC_TUBE "
#define PNEUMATIC_HOSE    "PNEUMATIC_HOSE "
#define FLEXIBLE_AXLE     "FLEXIBLE_AXLE "
#define FIBER_OPTIC_CABLE "FIBER_OPTIC_CABLE "
#define FIBRE_OPTIC_CABLE "FIBRE_OPTIC_CABLE "

#define RUBBER_HOSE       "RUBBER_HOSE "

#define ALL_HOSES TEST_HOSE,RIBBED_HOSE,RIBBED_TUBE,FLEX_CABLE,FLEXIBLE_HOSE,FLEXIBLE_TUBE,RIGID_HOSE,RIGID_TUBE,ELECTRIC_CABLE,PNEUMATIC_HOSE,PNEUMATIC_TUBE,FLEXIBLE_AXLE,FIBER_OPTIC_CABLE,FIBRE_OPTIC_CABLE,RUBBER_HOSE

#define RUBBER_BAND       "RUBBER_BAND "
#define CHAIN             "CHAIN "
#define PLASTIC_TREAD     "PLASTIC_TREAD "
#define RUBBER_TREAD      "RUBBER_TREAD "

#define ALL_BANDS RUBBER_BAND,CHAIN,PLASTIC_TREAD,RUBBER_TREAD

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
  PRECISION exit_start, exit_stop;
} LSL_tube_attributes;

typedef struct {
  LSL_part_usage part;
  PRECISION      radius;
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

#endif
