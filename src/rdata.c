#include "include/rdata.h"

void freeReturnData(ReturnData rdata) {
    if(tsdata) free(tsdata);
    if(xdotdata) free(xdotdata);
    if(dxdotdpdata) free(dxdotdpdata);
    if(dydxdata) free(dydxdata);
    if(dydpdata) free(dydpdata);
    if(Jdata) free(Jdata);
    if(zdata) free(zdata);
    if(sigmazdata) free(sigmazdata);
    if(szdata) free(szdata);
    if(ssigmazdata) free(ssigmazdata);
    if(srzdata) free(srzdata);
    if(s2rzdata) free(s2rzdata);
    if(xdata) free(xdata);
    if(sxdata) free(sxdata);
    if(ydata) free(ydata);
    if(sigmaydata) free(sigmaydata);
    if(sydata) free(sydata);
    if(ssigmaydata) free(ssigmaydata);
    if(numstepsdata) free(numstepsdata);
    if(numrhsevalsdata) free(numrhsevalsdata);
    if(orderdata) free(orderdata);
    if(numstepsSdata) free(numstepsSdata);
    if(numrhsevalsSdata) free(numrhsevalsSdata);
    if(llhdata) free(llhdata);
    if(sllhdata) free(sllhdata);
    if(s2llhdata) free(s2llhdata);
    if(chi2data) free(chi2data);
    free(rdata);
}
