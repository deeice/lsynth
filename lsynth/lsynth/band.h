/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef BAND_H
#define BAND_H

#ifdef _cplusplus
extern "C" {
#endif

int
synth_band(
  char                *type,
  int                  n_constraints,
  LSL_band_constraint *constraints,
  int                  color,
  FILE                *output,
  int                  ghost);

#ifdef _cplusplus
};
#endif

#endif