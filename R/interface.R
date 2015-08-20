SetupKbac <- function(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen) {
    tmp <- .C("R_set_up_kbac_test", as.integer(nn), as.integer(qq), as.double(aa), as.double(mafUpper), as.double(as.vector(xdatIn)), as.double(as.vector(ydatIn)), as.double(as.vector(mafIn)), as.integer(xcol), as.integer(ylen));
}

DoKbac <- function(pval, sided, statistic) {
    tmp <- .C("R_do_kbac_test", pvalue = as.double(pval), alternative = as.integer(sided), kbac.stat=as.double(statistic));
    return(tmp);
}

RmKbac <- function() {
    tmp <- .C("R_clear_kbac_test");
}

KbacGetP <- function(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen, sided) {
    SetupKbac(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen);
    pvalue <- DoKbac(9.0, sided, 10000.0);
    RmKbac()
    return(pvalue);
}
