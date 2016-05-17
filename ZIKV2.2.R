##########################################################################
######                 Zika Virus Model 2.2                         ######
######         Working Copy to Play With Bite Rates                 ######
######                                                              ######
##########################################################################


### Host infection rates, beta_V*Host are modeled as:
###     beta_ij = (probability of vector i biting host j)*(transmission rate, t_i)*(bites per day) 
### In this way, a constant vector biting rate is spread between a preferred host and 
### a less desirable host. 

### Vector infection rates, beta_Host*V is modeled as:
###     beta_ji = (probability of vector i biting host j)*(transmission rate, t_H)*(bites per day)


rm(list=ls())
#Supporting package installation if necessary
#install.packages("sensitivity"); install.packages("FME")
require(sensitivity); require(FME)


## Dummy Values to Run Model

### Beta_ij          # i = vectors C and I, j = hosts H and D
# beta_ij = (prob. i bites j)*(transmission rate i)*(bites per day)
pCH = .8; pIH = .7        # probability of i biting human host, complement = p(i bites D)
tC = .01; tI = .001       # vector competence; transmission probability per vector i
rC = 20;  rI = 20         # total number of bites per day per vector

betaCH = pCH*tC*rC; betaCD = (1-pCH)*tC*rC
betaIH = pIH*tI*rI; betaID = (1-pIH)*tI*rI

### Beta_ji
# beta_ji = (prob. i bites j)*(transmission rate i)*(bites per day)
tH = .01                  # host competence; transmission probability per infected human host
betaHC = pCH*tH*rC; betaHI = pIH*tH*rI

betaHH = 0.5                # human-human transmission rate (sexual transmission)



SH = 1000; SD = 1000; SC = 1000; SI = 1000
NH = 1000; ND = 1000; NC = 1000; NI = 1000
KH = 1000; KD = 10000; KC = 1000; KI = 1000;  
kC = NC/KC; kI = NI/KI
eH = 0.02; eC = 0.02; eI = 0.02
dH = .01; dC = 0.01; dI = 0.01
bH = .01; bC = 0.01; bI = 0.01

gammaH = 0.01

parms <- c(betaHH, betaCH, betaIH, betaCD, betaID, betaHC, betaHI, 
           SH, SD, SC, SI, NH, ND, NC, NI, KH, KD, KC, KI, kC, kI,
           eH, eC, eI, dH, dC, dI, bH, bC, bI, gammaH)

Compute_R0 <- function(parms) {
  with(as.list(parms), {                 #  p := parameters
    
    FV <- matrix(c(0, betaHH*SH/NH^2, 
                   0,             bH), byrow = TRUE, nrow = 2)
    
    VV <- matrix(c(eH+dH,   betaHH/NH,
                     -eH, gammaH + dH), byrow = TRUE, nrow = 2)
    
    # FH = matrix of incoming horizontal infectious into each pool EH, IH, EC, IC, EI, II.    
    FH <- matrix(c(0, betaHH*SH/NH^2, 0, betaCH*SH/NH, 0, betaIH*SH/NH, 
                   0,              0, 0,             0, 0,             0,
                   0, betaHC*SC/NH, 0,             0, 0,             0,
                   0,              0, 0,             0, 0,             0, 
                   0, betaHI*SI/NH, 0,             0, 0,             0, 
                   0,              0, 0,             0, 0,             0), 
                 byrow = TRUE, nrow = 6)
    
    
    # VH = matrix of all other movements into and out of pools such that dx = F - V.
    VH <- matrix(c(eH + dH,            0,            0,      0,            0,      0,
                       -eH, gammaH+dH-bH,            0,      0,            0,      0,
                         0,            0, eC + dC*kC,      0,            0,      0,
                         0,            0,          eC, dC*kC,            0,      0,
                         0,            0,            0,      0, eC + dI*kI,      0,
                         0,            0,            0,      0,         -eI, dI*kI), 
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

Compute_R0(parms)


### Model Sensitivity Analysis ###

## 1. Generate Latin Hypercube Sample 

vars <- rbind(betaHH = c(0,1), betaCH = c(0,1), betaIH = c(0,1),
              betaHC = c(0,1), betaHI = c(0,1),
              eH = c(0, 1), eC = c(0,1), eI = c(0,1),
              dH = c(0,1), dC = c(0,1), dI = c(0,1),
              bH = c(0,1), bC = c(0,1), bI = c(0,1),
              gammaH = c(0,1))

factors <- c("betaHH", "betaCH", "betaIH", "betaHC", "betaHI", 
             "eH", "eC", "eI", "dH", "dC", "dI", "bH", "bC", "bI", "gammaH")


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

prcc <- PCC$PRCC[,1]
print(PCC$PRCC)
barplot(prcc)   #x axes are in the same order as print(PRCC$PRCC) above

## Publication quality plot, make sure x axes match pcc$PRCC print out
#par(cex = 1.2)
#prcc <- PCC$PRCC[,1]
#barplot(prcc, ylim = c(-1,1), xaxt='n', ann=FALSE, yaxt='n',
#        xlab = "Parameter Value", ylab = "PRCC")
#abline(h=0, col = "red")
#xes <- c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9, 9.1, 10.3, 11.5, 12.7, 13.9, 15.1, 16.3, 17.5)
#segments(x0=xes, x1=xes, y0 = PCC$PRCC[,4], y1=PCC$PRCC[,5])        # vertical whiskers
#segments(x0=xes-.1, x1=xes+.1, y0 = PCC$PRCC[,4], y1=PCC$PRCC[,4])  # bottom horizontal
#segments(x0=xes-.1, x1=xes+.1, y0 = PCC$PRCC[,5], y1=PCC$PRCC[,5])  # top horizontal

#axis(1, at = xes, labels = expression(beta[HH], beta[CH], beta[IH], beta[HC], beta[HI], 
#                                      epsilon[H], epsilon[C], epsilon[I], 
#                                      d[H], d[C], d[I], b[H], b[C], b[I], gamma[H]))


#axis(2, at = seq(-0.8, 0.8, by = 0.2), label = c(-0.8,NA,NA,-0.2, NA,NA,0.4,NA,NA))
#box()



