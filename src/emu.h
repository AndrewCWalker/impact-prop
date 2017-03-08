#ifndef EMU_H
#define EMU_H

void emuInit(struct RSM_struct *RSMdata);

void emu(struct RSM_struct *RSMdata, double *xstar, double *specw, double *ystar);

#endif
