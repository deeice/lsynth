/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef TUBE_H
#define TUBE_H

#ifdef _cplusplus
extern "C" {
#endif

int
synth_tube(
  char                 *type,
  int                   n_constraints,
  LSL_part_usage       *constraints,
  int                   ghost,
  int                   color,
  FILE *output);

void mat_mult( PRECISION res[3][3], PRECISION lft[3][3], PRECISION rht[3][3]);
void rotate( PRECISION res[3], PRECISION loc[3], PRECISION rot[3][3]);
#ifdef _cplusplus
};
#endif

#endif