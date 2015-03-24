#include <stdio.h>
#include "common_funcs.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix prospect4_model(
        NumericVector params,
        NumericMatrix p4data){
    double N = params[1], Cab = params[2], Cw = params[3], Cm = params[4];
    int wl = p4data.nrow();
    NumericVector k(wl), theta(wl), Spec(wl);

    NumericMatrix::Column Cab_abs = p4data(_,0);
    NumericMatrix::Column Cw_abs = p4data(_,1);
    NumericMatrix::Column Cm_abs = p4data(_,2);
    NumericMatrix::Column tao1 = p4data(_,3);
    NumericMatrix::Column tao2 = p4data(_,4);
    NumericMatrix::Column rho1 = p4data(_,5);
    NumericMatrix::Column rho2 = p4data(_,6);
    NumericMatrix::Column x = p4data(_,7);
    NumericMatrix::Column y = p4data(_,8);

    k = (1.0/N) * (Cab * Cab_abs + Cw * Cw_abs + Cm * Cm_abs);

    for(int i=0; i<wl; i++) theta[i] = exp_int(k[i]);

    Spec = gpm(tao1, tao2, rho1, rho2, x, y, theta, N);
    return Spec;
}