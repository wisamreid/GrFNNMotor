#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _infot_reg(void);
extern void _intf6_reg(void);
extern void _intfsw_reg(void);
extern void _misc_reg(void);
extern void _nstim_reg(void);
extern void _staley_reg(void);
extern void _stats_reg(void);
extern void _vecst_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," infot.mod");
    fprintf(stderr," intf6.mod");
    fprintf(stderr," intfsw.mod");
    fprintf(stderr," misc.mod");
    fprintf(stderr," nstim.mod");
    fprintf(stderr," staley.mod");
    fprintf(stderr," stats.mod");
    fprintf(stderr," vecst.mod");
    fprintf(stderr, "\n");
  }
  _infot_reg();
  _intf6_reg();
  _intfsw_reg();
  _misc_reg();
  _nstim_reg();
  _staley_reg();
  _stats_reg();
  _vecst_reg();
}
