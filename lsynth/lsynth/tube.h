/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef TUBE_H
#define TUBE_H

#ifdef _cplusplus
extern "C" {
#endif

void list_hose_types(      void);
void list_hose_constraints(void);
void hose_ini(             void);
int ishosetype(      char *type);
int ishoseconstraint(char *type);

int
synth_hose(
  char   *type,
  int     n_constraints,
  part_t *constraints,
  int     ghost,
  int     color,
  FILE   *output);
#ifdef _cplusplus
};
#endif

#endif