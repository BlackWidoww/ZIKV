###############################################################################
######     Zika Virus Model per Van Den Driessche and Watmough (2002)    ######
######      Created March 22, 2016 by Courtney Shelley                   ######
######   Last Modified April 14, 2016 by Courtney Shelley                ######
###############################################################################

rm(list=ls())
setwd("~/Dropbox/BarkerLab/ZIKV")
require(sensitivity); require(FME)


## Dummy Values to Run Model

betaHH = 0.5; betaV1H = 0.5; betaV2H = 0.5; betaHV1 = 0.5; betaHV2 = 0.5
SH = 1000; SV1 = 1000; SV2 = 1000
NH = 1000; NV1 = 1000; NV2 = 1000
KH = 1000; KV1 = 1000; KV2 = 1000;  k1 = NV1/KV1; k2 = NV2/KV2
eH = 0.02; eV1 = 0.02; eV2 = 0.02
dH = .01; dV1 = 0.01; dV2 = 0.01
bH = .01; bV1 = 0.01; bV2 = 0.01

gammaH = 0.01


Compute_R0 <- function(parms) {
  with(as.list(parms), {                 #  p := parameters
    
    FV <- matrix(c(0, betaHH*SH/NH^2, 
                   0,             bH), byrow = TRUE, nrow = 2)
    
    VV <- matrix(c(eH+dH,   betaHH/NH,
                     -eH, gammaH + dH), byrow = TRUE, nrow = 2)
    
    # FH = matrix of incoming horizontal infectious into each pool EH, IH, EV1, IV1, EV2, IV2.    
    FH <- matrix(c(0, betaHH*SH/NH^2, 0, betaV1H*SH/NH, 0, betaV2H*SH/NH, 
                   0,              0, 0,             0, 0,             0,
                   0, betaHV1*SV1/NH, 0,             0, 0,             0,
                   0,              0, 0,             0, 0,             0, 
                   0, betaHV2*SV2/NH, 0,             0, 0,             0, 
                   0,              0, 0,             0, 0,             0), 
                 byrow = TRUE, nrow = 6)
    
    
    # VH = matrix of all other movements into and out of pools such that dx = F - V.
    VH <- matrix(c(eH + dH,            0,            0,      0,            0,      0,
                       -eH, gammaH+dH-bH,            0,      0,            0,      0,
                         0,            0, eV1 + dV1*k1,      0,            0,      0,
                         0,            0,          eV1, dV1*k1,            0,      0,
                         0,            0,            0,      0, eV1 + dV2*k2,      0,
                         0,            0,            0,      0,         -eV2, dV2*k2), 
                 byrow = TRUE, nrow = 6)
    

    
    # Next generation matrices
    NGV <- FV %*% solve(VV)
    NGH <- FH %*% solve(VH)
    
    # Compute the eigenvalues of the next generation matrices
    EigenV <- eigen(NGV, only.values=TRUE)
    EigenH <- eigen(NGH, only.values=TRUE)
    
    # Compute the spectral radius of the next generation matrix
    # which is the maximum of the absolute values of its real-part eigenvalues
    R0V <- max(abs(EigenV$values))
    R0H <- max(abs(EigenH$values))
    
    as.vector(R0V + R0H)
  })
}




### Model Sensitivity Analysis ###

## 1. Generate Latin Hypercube Sample 

vars <- rbind(betaHH = c(0,1), betaV1H = c(0,1), betaV2H = c(0,1),
              betaHV1 = c(0,1), betaHV2 = c(0,1),
              eH = c(0, 1), eV1 = c(0,1), eV2 = c(0,1),
              dH = c(0,1), dV1 = c(0,1), dV2 = c(0,1),
              bH = c(0,1), bV1 = c(0,1), bV2 = c(0,1),
              gammaH = c(0,1))

factors <- c("betaHH", "betaV1H", "betaV2H", "betaHV1", "betaHV2", 
             "eH", "eV1", "eV2", "dH", "dV1", "dV2", "bH", "bV1", "bV2", "gammaH")


parRange <- data.frame(min = vars[,1], max = vars[,2])
rownames(parRange) <- factors

LHS <- data.frame(Latinhyper(parRange, 1000))
LHSplot <- pairs(LHS, main = "Latin Hypercube")

y <- vector()
for(i in 1:nrow(LHS)) {
  R0 <- Compute_R0(LHS[i,])
  y <- append(y, R0, after = length(y))
}

PCC <- pcc(X = LHS, y = y, rank = TRUE, nboot = 100, conf = 0.95)

## Plot

par(cex = 1.2)
prcc <- PCC$PRCC[,1]
barplot(prcc, ylim = c(-1,1), xaxt='n', ann=FALSE, yaxt='n',
        xlab = "Parameter Value", ylab = "PRCC")
abline(h=0, col = "red")
xes <- c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9, 9.1, 10.3, 11.5, 12.7, 13.9, 15.1, 16.3, 17.5)
segments(x0=xes, x1=xes, y0 = PCC$PRCC[,4], y1=PCC$PRCC[,5])        # vertical whiskers
segments(x0=xes-.1, x1=xes+.1, y0 = PCC$PRCC[,4], y1=PCC$PRCC[,4])  # bottom horizontal
segments(x0=xes-.1, x1=xes+.1, y0 = PCC$PRCC[,5], y1=PCC$PRCC[,5])  # top horizontal

axis(1, at = xes, labels = expression(beta[HH], beta[V1H], beta[V2H], beta[HV1], beta[HV2], 
                                      epsilon[H], epsilon[V1], epsilon[V2], 
                                      d[H], d[V1], d[V2], b[H], b[V1], b[V2], gamma[H]))


axis(2, at = seq(-0.8, 0.8, by = 0.2), label = c(-0.8,NA,NA,-0.2, NA,NA,0.4,NA,NA))
box()


