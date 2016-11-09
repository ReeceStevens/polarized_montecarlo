#ifndef __scat_h__
#define __scat_h__

struct {
    double *s11;
    double *s12;
    double *s33;
    double *s43;
    double (*mu_s) (double); // Function defining mu_s at a given polarization state
} typedef scat_props_t;

#endif
