#include "include/rdata.h"
#include<include/rdata_accessors.h>

void freeReturnData(ReturnData *rdata) {
    if(tsdata) delete[] tsdata;
    if(xdotdata) delete[] xdotdata;
    if(dxdotdpdata) delete[] dxdotdpdata;
    if(dydxdata) delete[] dydxdata;
    if(dydpdata) delete[] dydpdata;
    if(Jdata) delete[] Jdata;
    if(zdata) delete[] zdata;
    if(sigmazdata) delete[] sigmazdata;
    if(szdata) delete[] szdata;
    if(ssigmazdata) delete[] ssigmazdata;
    if(srzdata) delete[] srzdata;
    if(s2rzdata) delete[] s2rzdata;
    if(xdata) delete[] xdata;
    if(xssdata) delete[] xssdata;
    if(sxdata) delete[] sxdata;
    if(ydata) delete[] ydata;
    if(sigmaydata) delete[] sigmaydata;
    if(sydata) delete[] sydata;
    if(ssigmaydata) delete[] ssigmaydata;
    if(numstepsdata) delete[] numstepsdata;
    if(numrhsevalsdata) delete[] numrhsevalsdata;
    if(orderdata) delete[] orderdata;
    if(numstepsSdata) delete[] numstepsSdata;
    if(numrhsevalsSdata) delete[] numrhsevalsSdata;
    if(ns_numstepsdata) delete[] ns_numstepsdata;
    if(ns_numlinstepsdata) delete[] ns_numlinstepsdata;
    if(ns_timedata) delete[] ns_timedata;
    if(newtdata) delete[] newtdata;
    if(llhdata) delete[] llhdata;
    if(sllhdata) delete[] sllhdata;
    if(s2llhdata) delete[] s2llhdata;
    if(chi2data) delete[] chi2data;
    delete rdata;
}
