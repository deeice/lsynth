/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef CURVE_H
#define CURVE_H

#ifdef _cplusplus
extern "C" {
#endif

int
synth_curve(
  LSL_part_usage       *start,
  LSL_part_usage       *end,
  LSL_part_usage       *segments,
  int                   n_segments,
  LSL_tube_attributes  *attrib,
  FILE                 *output);

#ifdef _cplusplus
};
#endif

#endif