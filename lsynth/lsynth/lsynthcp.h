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
#include "mathlib.h"

#define ACCY (1e-6)

typedef struct {
  char      type[128];
  PRECISION orient[3][3];
  PRECISION offset[3];
  int       attrib;
} part_t;

typedef struct {
  part_t    part;
  PRECISION radius;
  int       inside;
  int       was_cross;
  int       cross;
  PRECISION start_line[3];
  PRECISION end_line[3];
  PRECISION crossings[8][3];
  int       n_crossings;
  int       layer;

  PRECISION start_angle[3];
  PRECISION end_angle[3];
  int       n_steps;
  PRECISION s_angle;

  int       constraint_type_n;
} LSL_band_constraint;

extern PRECISION hose_res_angle;
extern PRECISION band_res;

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