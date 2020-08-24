MWerror <- function(t1, dt1, t2, dt2, M=41.5, M.theory=19.25){
    #error transfer function to estimate molecular weight
    #t1, t2: diffusion time
    #dt1, dt2: standard error of diffusion time
    #M: refered protein weight
    #M.theory: theory weight of protein to estimate
    M.est <- (t1/t2)^3*M
    M.RE <- 3*sqrt((dt1/t1)^2+(dt2/t2)^2)
    M.err <- M.est*M.RE
    boundNum <- (M.est-40.25)/M.theory
    boundNum.err <- M.err/M.theory
    print(paste('Estimate weight:', round(M.est, 0), '+/-', round(M.err, 0), '(', round(M.RE*100, 0), '% ) kD'))
    print(paste('Bound number per substrate:', round(boundNum, 1), '+/-', round(boundNum.err, 1)))
}
