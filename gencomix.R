

################################################################################
################################################################################
################################################################################
# sfs_code:
./sfs_code 1 1 | ./sfs2vcf.pl --onlypoly

vcf --> input format

library(Rplinkseq)

setwd( "~/Dropbox/lab/hmm/" )
load.vcf( "data/test.vcf" )

# library(vcf2geno)
# readVCFToListByGene( "data/test.vcf" )
# readVCFToListByRange( "data/test.vcf", range=c(0,1000), annoType="", vcfColumn="", vcfInfo="" )


# setwd( "~/Dropbox/lab/hmm/" )
# library(VariantAnnotation)
# readVcf( "data/test.vcf", "hg00" )

scanVcf( "data/test.vcf" )
scanVcf("/home/ccampbell/cortex/vcf/out.recode.vcf" )


./sfs_code 3 1 -TS 0 0 1 -TE 2 -TJ 2 2 500 2 0 1 0.8 0.2 F 0.2 0.8 | ./sfs2vcf.pl --onlypoly > ~/Dropbox/lab/hmm/data/sfssim_02.vcf

"~/Dropbox/lab/hmm/data/sfssim_01.vcf"
################################################################################
################################################################################
################################################################################
# dataset generation:
S <- 50
N <- c(10,10)
########################################
# distances:
load("~/storage/jhap/annotations/Affy_6.0_SNP_na32_2012-0418.RData")
dvec <- annot$Position[annot$Chromosome==1&annot$Position<1600000][5+1:S]

# ref haplotypes:
refH1 <- lapply(1:N[1], function(x) sample( 0:1, S, replace=TRUE, prob=c(0.01,0.99) ) )
refH2 <- lapply(1:N[2], function(x) sample( 0:1, S, replace=TRUE, prob=c(0.99,0.01) ) )

# simple admixture:
tmp <- which(diff(dvec)>40000)
tmp <- cbind( c(1,tmp+1), c(tmp,S) )
am <- integer(S); h <- 0
for(i in 1:nrow(tmp))  {
    am[ tmp[i,1]:tmp[i,2] ] <- h
    if( h==0 ) { 
        h <- 1 
    } else h <- 0
}
am[50] <- 0
# admixture w/ gene conversion:
tmp0 <- which(diff(dvec)>100 & diff(dvec)<800)
tmp <- cbind( tmp0, tmp+1)
amgc <- am
# amgc[ 6 ] <- 1
# amgc[ 30:31 ] <- 0
# amgc[ 43:44 ] <- 0
amgc[ 36 ] <- 0
dvec[44 ] <- dvec[43]+500

sampHap <- list( amgc,am  )

nk <- length(refH1)+length(refH2) + length(sampHap)
sXh <- paste("h",1:nk,sep="")
i <- cumsum(c(length(refH1),length(refH2) , length(sampHap)))
names(refH1) <- sXh[ 1:i[1] ]
names(refH2) <- sXh[ (i[1]+1):i[2] ]
names(sampHap) <- sXh[ (i[2]+1):i[3] ]

plot(NULL,xlim=range(dvec), ylim=c(nk+1,1), xlab="Position (bp)", ylab="Haplotype", yaxt="n")
axis( 2, labels=sXh, at=1:nk+0.5, las=1)
axis( 4, labels=NA, at=1:nk+0.5 )
abline(v=dvec, lty=1, col="grey90" )
cnt <- 1
for(x in c(refH1,refH2,sampHap)) {
    xcol <- ifelse( x==0, 4,2 )
    abline(h=cnt+0.5,col="grey80")
    points( dvec, rep(cnt+0.5,S), pch=20, col=xcol )
    cnt <- cnt+1
}
abline(h=i[-3]+1,lwd=2)


writeHap <- function( hapList, pop="p1", datafile ) {
    for(x in seq_along(hapList)) {
        tmp <- c( pop,names(hapList)[x], hapList[[x]] )
        write.table( t(tmp),datafile, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)
    }
}
outfile <- paste("~/Dropbox/lab/hmm/data/sampleData_",paste(N,collapse="-"),"_cgchmm.txt",sep="")
write.table(t(as.matrix(c( sum(N)+2 , S ) )),outfile, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table( t(dvec),outfile, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)
writeHap( refH1, "p1", outfile)
writeHap( refH2, "p2", outfile)
writeHap( sampHap, "0", outfile)

dvec[17] <- dvec[16]+250
dvec[18] <- dvec[16]+500
obs <- rep(0,50)
obs[17] <- 1

################################################################################
################################################################################
################################################################################
# load in data:
#fname <- "~/Dropbox/lab/hmm/data/sampleData_cgchmm.txt"
#fname <- paste("~/Dropbox/lab/hmm/data/sampleData_",paste(N,collapse="-"),"_cgchmm.txt",sep="")
fname <- "~/Dropbox/lab/hmm/data/sampleData_10-10_cgchmm.txt"
#readHap <- function( fname="~/Dropbox/lab/hmm/data/sampleData_cgchmm.txt" ) {
tmp <- scan( fname, what=integer(), nlines=1, quiet=TRUE)
N <- tmp[1]; S <- tmp[2]
dvec <- scan( fname, what=double(), nlines=1, skip=1, quiet=TRUE)
hap <- list()
for(i in 1:N) {
    tmp <- scan( fname, what=character(), nlines=1, skip=1+i, quiet=TRUE)
    hap[[i]] <- as.integer(tmp[-c(1:2)])
    names(hap)[i] <- paste(tmp[1:2],collapse="-")
}
ncut <- c(4,4) # cut down hap in size
refH1 <- hap[ grep("^p1",names(hap))[1:ncut[1]] ]
names(refH1) <- sapply(strsplit(names(refH1),"-"),'[',2)

refH2 <- hap[ grep("^p2",names(hap))[1:ncut[2]] ]
#names(refH2) <- sapply(strsplit(names(refH2),"-"),'[',2)
names(refH2) <- paste("h",(ncut[1]+1):sum(ncut) ,sep="")
sampHap <- hap[ grep("^0",names(hap)) ]
#names(sampHap) <- sapply(strsplit(names(sampHap),"-"),'[',2)
names(sampHap) <- paste("h",(sum(ncut)+1):((sum(ncut))+length(sampHap)) ,sep="")

# refH1 <- refH1[1]
# refH2 <- refH2[1]; names(refH2)="h2"
# names(sampHap) <- c("h3","h4")
################################################################################
################################################################################

########################################
########################################
########################################
########################################
# gay/mcvean:
#dvec <- dvec * 1.3e-6
#dvec <- dvec / 1.3e-6
dvec <- dvec/10000


logsum <- function(x) {
    maxi <- which.max(x)
    maxl <- x[maxi]
    maxl + log1p(sum(exp(x[-maxi]-maxl))) 
}

S <- length(sampHap[[1]])# how many biallelic loci (SNPs, segregating sites). 2^S possible haplotypes
n <- length(sampHap) # how many haplotypes in the sample ###haploid individuals

### parameters:
# pop sizes:
n1 <- length( refH1 ) # European
n2 <- length( refH2 ) # African
#theta <- sum( 1/(1:(n-1)) )^-1
theta <- 1/1000 # (1 per kb)
rho <- 1/1000000 # crossover rate (1 per kb)
gam <- 1/1000 # gene conversion rate 
lam <- 1/500 # gc tract termination rate

# # random ordering:
# ord <- list()
# #for(i in 1:20) ord[[i]] <- sample(1:n,n)
# #ord[[20]] <- 1:n
# ord[[1]] <- 1:n
# piAp <- matrix( as.numeric(NA), ncol=length(ord), nrow=n, dimnames=list(names(sampHap),paste("perm",1:length(ord),sep="")) )
# for(o in seq_along(ord) ) {
#     cat("Permutation",o)
# #ord[[1]] <- 1:n

# piA <- numeric(n); names(piA) <- names(sampHap)[ ord[[o]] ]
# piA[1] <- log( (1/n)^2 )
p <- list()

k=1
#for(k in 2:length(sampHap) ) {
km1 <- n1+n2 #k-1
oindx <- k #ord[[o]][k] # "observed" haplotype within this perm
#hindx <- ord[[o]][1:km1] # haplotype order for this perm
obs <- sampHap[[k]] # copy of observed haplotype

transL <- list()
for(j in 1:(S-1) ) { # site along haplotype sequence
d <- diff( dvec[j:(j+1)] )

nk <- length(refH1)+length(refH2) # number of total haplotypes in both reference populations
sXh <- c(names(refH1),names(refH2) ) # haplotype (x) states
sXp <- c( rep("p1",length(refH1)), rep("p2",length(refH2)) ) # population states
sX <- paste(sXp,sXh,sep="-")

### Haplotype (X) crossover transition probabilities:
#transX <- matrix( numeric(1), nrow=km1,ncol=km1, dimnames=list(names(sampHap)[hindx],names(sampHap)[hindx]) )
transX <- matrix( numeric(1), nrow=nk,ncol=nk, dimnames=list(sXh,sXh) )
# if x != x' :
transX[,] <- 1/km1 * ( 1 - exp(-rho*d/km1) )
# if x == x' :
diag(transX) <- exp(-rho*d/km1) + 1/km1 * (1 - exp(-rho*d/km1))

### Gene conversion (G) transition probabilities:
nk <- length(refH1)+length(refH2) + 1
sGh <- c(rep("0",1), names(refH1),names(refH2) ) #,names(sampHap)[hindx]) # haplotypes
sGp <- c( "0", rep("p1",length(refH1)), rep("p2",length(refH2)) ) #, rep("0",km1) ) # populations
sG <- paste(sGp,sGh,sep="-")
transG <- matrix( numeric(1), nrow=nk,ncol=nk, dimnames=list(sGh,sGh) )
#transG <- matrix( numeric(1), nrow=k,ncol=k, dimnames=list(c("g0",names(sampHap)[hindx]),c("g0",names(sampHap)[hindx])) )
# 1: Pr(G[j+1]=0 | G[j] = 0):
a1 <- exp(-lam*d) * exp(-gam*d/km1) + (km1*lam*(1-exp(-d*(gam/km1+lam)))) / (km1*lam+gam)
# 2: Pr(G[j+1]=g | G[j] = 0):
a2 <- exp(-lam*d)/km1 * (1-exp(-gam*d/km1)) + ( lam*(exp(-d*(gam/km1+lam))-1) / (gam+km1*lam) + (1-exp(-d*lam))/km1 )
# 3: Pr(G[j+1]=0 | G[j] = g):
a3 <- ( km1*lam* (1-exp(-d*(lam+gam/km1))) ) / ( gam + lam*km1 )
# 5: Pr(G[j+1]=g' | G[j] = g):
a5 <- (1-exp(-d*lam))/km1 + (lam*(exp(-d*(lam+gam/km1))-1)) / (km1*lam+gam)
# 4: Pr(G[j+1]=g | G[j] = g):
a4 <- exp(-lam*d) + a5

transG["0","0"] <- a1
transG["0",-1] <- a2
transG[-1,"0"] <- a3
transG[-1,-1] <- a5
diag(transG)[-1] <- a4

### combined transition matrix:
states <- apply( expand.grid( rownames(transX), rownames(transG) ), 1, paste, collapse="-")
trans <- matrix( numeric(1), nrow=length(states),ncol=length(states), dimnames=list( states, states ))
sX <- sapply(strsplit(states,"-"),"[",1)
sG <- sapply(strsplit(states,"-"),"[",2)

for(i in 1:nrow(trans)) {
    for(q in 1:nrow(trans)) {
        trans[i,q] <- transX[sX[i],sX[q]] * transG[sG[i],sG[q]]
    }
}
    transL[[j]] <- trans
} # end j loop
# emissions:
L <- 1
emit <- c(
    match = (2*km1*L+theta) / (2*(km1*L+theta)), # mismatch
    mismatch = theta / ( 2*(km1*L+theta) ) # match
)

####################
# pre-compute a matrix of match/mismatches to use as indexes to emit:
emat <- t(sapply( c(refH1,refH2), function(x) x==obs ))
emat[emat==TRUE] <- 1
emat[emat==FALSE] <- 2

sXG <- sG
sXG[sG=="0"] <- sX[sG=="0"]
####################
# starting probablities for X and G:
sprob <- numeric(length(states))
g0i <- grep("0",sG)
g1i <- setdiff(1:length(states),g0i)
#eIndx <- ifelse( c(sapply(refH1,'[',1),sapply(refH2,'[',1))[sX] ==obs[1], 1,2)
eIndx <- emat[sXG,1]
sprob[g0i] <- log( 1/km1 * (lam*km1) / (lam*km1+gam) * emit[eIndx[g0i]] )
sprob[g1i] <- log( 1/km1 * gam / (km1*(lam*km1+gam)) * emit[eIndx[g1i]] )


#fwdc <- matrix( numeric(1), nrow=length(states), ncol=S, dimnames=list(states=sXG, obs=obs) )
#fwd:
fwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
fwd[,1] <- sprob
for(j in 2:S) {
    cnt <- 1
    for(state1 in states) { # hap loop
        lsum <- -Inf
        for(state0 in states) {
            tmp <- fwd[state0,j-1] + log( trans[state0,state1] )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
        }
        hname <- sXG[cnt]; cnt <- cnt+1
        fwd[state1,j] <- log( emit[ emat[hname,j] ] ) + lsum
        #fwdc[hname,j] <- fwdc[hname,j]+1
    }
}
#piA[k] <- logsum(fwd[,S])
Pxa <- logsum(fwd[,S])

#bwd:
bwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
bwd[,S] <- 0
for(j in (S-1):1) {
    for(state1 in states) {
        lsum <- -Inf; cnt <- 1
        for(state2 in states) {
            st2 <- sXG[cnt]; cnt <- cnt+1
            #st2 <- strsplit(state2,"-")[[1]][1]
            #eIndx <- ifelse( sapply(sampHap[st2],'[',j+1) == obs[j+1], 1,2)
            #eIndx <- ifelse( c(sapply(refH1,'[',j+1),sapply(refH2,'[',j+1))[st2] ==obs[j+1],1,2)
            eIndx <- emat[st2,j+1]
            tmp <- bwd[state2,j+1] +
            log( trans[state1,state2] * emit[eIndx] )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
        }
        bwd[state1,j] <- lsum
    }
} # seq loop j

Pxb <- bwd[1,1] + sprob[1]; cnt <- 2
if(length(states)>1) {
    for(state1 in states[-1]) {
        tmp <- bwd[state1,1] + fwd[cnt,1]
        if(tmp>-Inf) Pxb <- tmp + log(1+exp(Pxb-tmp))
        cnt <- cnt+1
    }
}
Pxa-Pxb
#piA[k]-Pxb
### posterior decoding:
pprob <- exp( fwd+bwd - Pxa )
path <- apply(pprob,2,max)
names(path) <- states[ apply(pprob,2,which.max) ]

#Lpac <- colSums(piAp)


#pdf("~/Dropbox/lab/hmm/results/paths_1-1-h9.pdf",width=11,height=6)
pdf( paste("~/Dropbox/lab/hmm/results/paths_bp_",n1,"-",n2,".pdf",sep=""),width=11,height=6)
split.screen( rbind(
                    c( 0, 1, 0.35, 1 ),
                    c( 0, 1, 0, 0.35 )
    ))
screen(2)
par( mar=c(5,5,0,1) )
plot(dvec, 1-colSums(pprob[ sG=="0",]), xlab="", ylim=c(0,1), las=1, type="n",pch=20, ylab="P(gene conv)")
abline(v=dvec, lty=1, col="grey90" )
points(dvec, 1-colSums(pprob[ sG=="0",]), type="o",pch=20 )
title(xlab="Position (bp)" )
screen(1)
np <- k # which haplotype to "observe"
nk <- n1+n2+1
par( mar=c(0,5,1,1) )
plot(NULL,xlim=range(dvec), ylim=c(nk+0.5,0), xlab="", ylab="Haplotype", yaxt="n",main="GenCo", xaxt="n")
axis( 2, labels=c("g0",sXh,"obs"), at=0:nk, las=1)
axis( 4, labels=NA, at=1:nk )
abline(v=dvec, lty=1, col="grey90" )
cnt <- 1
for(x in c(refH1,refH2,sampHap[k])) {
    xcol <- ifelse( x==0, 4,2 )
    abline(h=cnt,col="grey80")
    points( dvec, rep(cnt,S), pch=20, col=xcol )
    cnt <- cnt+1
}
abline(h=c(1,cumsum(c(length(refH1),length(refH2) ))+1)-0.5, lwd=2)
# draw path:
pm <- do.call("rbind",strsplit(names(path),"-"))
colnames(pm) <- c("X","G")
points( dvec, as.numeric(gsub("^h","",pm[,"X"])), type="b", cex=2)
points( dvec, as.numeric(gsub("^[a-z]","",pm[,"G"])), type="b", cex=1.5, pch=5)

close.screen(all=TRUE)

###########################################################
###########################################################
###########################################################
###########################################################
################################################################################
################################################################################
# modified hapmix / gay model:
logsum <- function(x) {
    maxi <- which.max(x)
    maxl <- x[maxi]
    maxl + log1p(sum(exp(x[-maxi]-maxl))) 
}

########################################
########################################

# pop sizes:
n1 <- length( refH1 ) # European
n2 <- length( refH2 ) # African
# free parameters:
T <- 7
u1 <- 0.8
# fixed parameters:
# rho1 <- 60000/n1 # per morgan, European
# rho2 <- 90000/n2 # per morgan, African
# theta1 <- 0.2/(0.2+n1)
# theta2 <- 0.2/(0.2+n2)
# theta3 <- 0.01
#theta <- sum( 1/(1:(n-1)) )^-1
#theta <- 1/1000 # (1 per kb)
#rho <- 1/1000 # crossover rate (1 per kb)
#gam <- 1/1000 # gene conversion rate 
#lam <- 1/500 # mean GC tract length
#
S <- length(sampHap[[1]])# how many biallelic loci (SNPs, segregating sites). 2^S possible haplotypes
n <- length(sampHap) # how many haplotypes in the sample ###haploid individuals


#Rprof()
k=1
obs <- sampHap[[k]] # copy of observed haplotype

transL <- list()
for(j in 1:(S-1) ) { # site along haplotype sequence
d <- diff( dvec[j:(j+1)] )

nk <- length(refH1)+length(refH2) # number of total haplotypes in both reference populations
sXh <- c(names(refH1),names(refH2) ) # haplotype (x) states
sXp <- c( rep("p1",length(refH1)), rep("p2",length(refH2)) ) # population states
sX <- paste(sXp,sXh,sep="-")

### Haplotype (X) crossover transition probabilities:
#transX <- matrix( numeric(1), nrow=nk,ncol=nk, dimnames=list(sX,sX) )
transX <- matrix( numeric(1), nrow=nk,ncol=nk, dimnames=list(sXh,sXh) )
# if x != x' and p != p':  (1-exp(-d*T)) * u1/n1
# if x != x' and p == p':  exp(-d*T) * (1-exp(-d*rho))/n1 + (1-exp(-d*T))*u1/n1
# if x == x' and p == p':  exp(-d*T) * exp(-d*rho) + exp(-d*T) * (1-exp(-d*rho))/n1 + (1-exp(-d*T))*u1/n1
for(f in 1:nk) { # transition from...
    for(t in 1:nk) { # ... to:
        if( sXp[t]=="p1" ) { u <- u1; nm <- n1 }
        if( sXp[t]=="p2" ) { u <- 1-u1; nm <- n2 }
        if( sXp[f] != sXp[t] ) { # ancestry switch
            #if( sXh[f] != sXh[t] ) tmp <- (1-exp(-d*T)) * u/nm # anc AND hap switch
            if( sXh[f] != sXh[t] ) tmp <- (1-exp(-d*rho*T)) * u/nm # anc AND hap switch
        } else { # no ancestry switch
            #if( sXh[f] != sXh[t] ) tmp <- exp(-d*T) * (1-exp(-d*rho))/nm + (1-exp(-d*T))*u/nm # hap switch, NO anc switch
            if( sXh[f] != sXh[t] ) tmp <- exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm # hap switch, NO anc switch
            if( sXh[f] == sXh[t] ) { # NO hap switch, NO anc switch
                #tmp <- exp(-d*T) * exp(-d*rho) + exp(-d*T) * (1-exp(-d*rho))/nm + (1-exp(-d*T))*u/nm
                tmp <- exp(-d*rho*T) * exp(-d*rho) + exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm
            }
        } # endif
        transX[f,t] <- tmp
    }
}
rowSums(transX)==1

tmp1 <- (1-exp(-d*rho*T)) * u/nm # anc AND hap switch
tmp2 <- exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm # hap switch, NO anc switch
tmp3 <- exp(-d*rho*T) * exp(-d*rho) + exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm

### Gene conversion (G) transition probabilities:
nk <- length(refH1)+length(refH2) + 1
sGh <- c(rep("0",1), names(refH1),names(refH2) ) #,names(sampHap)[hindx]) # haplotypes
sGp <- c( "0", rep("p1",length(refH1)), rep("p2",length(refH2)) ) #, rep("0",km1) ) # populations
sG <- paste(sGp,sGh,sep="-")
#transG <- matrix( numeric(1), nrow=nk,ncol=nk, dimnames=list(sG,sG) )
transG <- matrix( numeric(1), nrow=nk,ncol=nk, dimnames=list(sGh,sGh) )

# a1=1;a2=2;a3=3;a4=4;a5=5
for(f in 1:nk) { # transition from...
    for(t in 1:nk) { # ... to:
        u <- 1
        if( sGp[t]=="p1" ) {
            u <- u1
            nm <- n1
        }
        if( sGp[t]=="p2" ) {
            u <- 1-u1
            nm <- n2
        }
        # 3: Pr(G[j+1]=0 | G[j] = g):
        a3 <- lam*(n1+n2)* (1-exp(-d*(gam*T+lam*(n1+n2))/(n1+n2))) / (gam*T+lam*(n1+n2))
        # 5: Pr(G[j+1]=g' | G[j] = g):
        a5 <- (lam*u*(n1+n2) * (exp(d*(-gam*T/(n1+n2)-lam))-1)) / ((n1+n2)*nm*lam + gam*T*nm) + (u-u*exp(-d*lam))/nm
        # 1: Pr(G[j+1]=0 | G[j] = 0):
        a1 <- exp(-lam*d) * exp(-gam*T*d/(n1+n2)) + a3
        # 2: Pr(G[j+1]=g | G[j] = 0):
        a2 <- exp(-lam*d) * (1-exp(-gam*T*d/(n1+n2))) *u/nm + a5
        # 4: Pr(G[j+1]=g | G[j] = g):
        a4 <- exp(-lam*d) + a5
        if( sGh[f] == sGh[t] ) tmp <- a4 # also includes 0->0
        if( sGh[f]=="0" & sGh[t]=="0" ) tmp <- a1
        if( sGh[f]=="0" & sGh[t]!="0" ) tmp <- a2
        if( sGh[f]!="0" & sGh[t]=="0" ) tmp <- a3
        if( sGh[f]!="0" & sGh[t]!="0" & sGh[f]!=sGh[t] ) tmp <- a5
        transG[f,t] <- tmp
    }
}
# rowSums(transG)
# rowSums(transG)==1

### combined transition matrix:
#states <- apply( expand.grid( rownames(transX), rownames(transG) ), 1, paste, collapse="-")
states <- apply( expand.grid( sX, sG ), 1, paste, collapse="-")
trans <- matrix( numeric(1), nrow=length(states),ncol=length(states), dimnames=list( states, states ))
sP <- sapply(strsplit(states,"-"),"[",1)
sX <- sapply(strsplit(states,"-"),"[",2)
sG <- sapply(strsplit(states,"-"),"[",4)

for(i in 1:nrow(trans)) {
    for(q in 1:nrow(trans)) {
        trans[i,q] <- transX[ sX[i],sX[q]] * transG[ sG[i], sG[q] ]
    }
}
transL[[j]] <- trans

} # end j loop

L <- 1
emit <- c(
    match = (2*(n1+n2)*L+theta) / (2*((n1+n2)*L+theta)) , # match
    mismatch = theta / ( 2*((n1+n2)*L+theta) ) # mismatch
)

sXG <- sG
sXG[sG=="0"] <- sX[sG=="0"]
####################
# starting probablities for X and G:
sprob <- numeric(length(states))
g0i <- grep("0",sG)
g1i <- setdiff(1:length(states),g0i)
#eIndx <- ifelse( c(sapply(refH1,'[',1),sapply(refH2,'[',1))[sX] ==obs[1], 1,2)
eIndx <- emat[ sXG, 1 ]
# sprob[g0i] <- log( 1/(n1+n2) * (lam*(n1+n2)) / (lam*(n1+n2)+gam) * emit[eIndx[g0i]] )
# sprob[g1i] <- log( 1/(n1+n2) * gam / ((n1+n2)*(lam*(n1+n2)+gam)) * emit[eIndx[g1i]] )
# 
sprob[g0i] <- log( 1/(n1+n2) * 
                  (lam*(n1+n2)) / (lam*(n1+n2)+gam*T) 
                  * emit[eIndx[g0i]] )
sprob[g1i] <- log( 1/(n1+n2) * 
                  (gam*T) / ((n1+n2)*(lam*(n1+n2)+gam*T)) 
                  * emit[eIndx[g1i]] )


#fwd:
fwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
fwd[,1] <- sprob
for(j in 2:S) {
    cnt <- 1
    for(state1 in states) { # hap loop
        lsum <- -Inf
        for(state0 in states) {
            tmp <- fwd[state0,j-1] + log( trans[state0,state1] )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
        }
        hname <- sXG[cnt]; cnt <- cnt+1
        #fwd[state1,j] <- log( emit[ ifelse( c(sapply(refH1,'[',j),sapply(refH2,'[',j))[hname] ==obs[j],1,2) ] ) +lsum
        fwd[state1,j] <- log( emit[ emat[hname,j] ] ) + lsum
    }
}
Pxa <- logsum(fwd[,S])

#bwd:
bwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
bwd[,S] <- 0
for(j in (S-1):1) {
    for(state1 in states) {
        lsum <- -Inf
        cnt <- 1
        for(state2 in states) {
            st2 <- sXG[cnt]; cnt <- cnt+1
            #st2 <- strsplit(state2,"-")[[1]][2]
            #eIndx <- ifelse( c(sapply(refH1,'[',j+1),sapply(refH2,'[',j+1))[st2] ==obs[j+1],1,2)
            eIndx <- emat[st2,j+1]
            tmp <- bwd[state2,j+1] +
                log( trans[state1,state2] * emit[eIndx] )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
        }
        bwd[state1,j] <- lsum
    }
} # seq loop j

Pxb <- bwd[1,1] + sprob[1]; cnt <- 2
if(length(states)>1) {
    for(state1 in states[-1]) {
        tmp <- bwd[state1,1] + fwd[cnt,1]
        if(tmp>-Inf) Pxb <- tmp + log(1+exp(Pxb-tmp))
        cnt <- cnt+1
    }
}
Pxa-Pxb
### posterior decoding:
pprob <- exp( fwd+bwd - Pxa )
path <- apply(pprob,2,max)
names(path) <- states[ apply(pprob,2,which.max) ]

#Rprof(NULL); summaryRprof()

# rowMeans( piAp )
# apply( piAp, 2, prod )
# Lpac <- colSums(piAp)



split.screen( rbind(
                    c( 0, 1, 0.35, 1 ),
                    c( 0, 1, 0, 0.35 )
    ))
screen(2)
par( mar=c(5,5,0,1) )
plot(dvec, 1-colSums(pprob[ sG=="0",]), xlab="", ylim=c(0,1), las=1, type="n",pch=20, ylab="P(gene conv)")
abline(v=dvec, lty=1, col="grey90" )
points(dvec, 1-colSums(pprob[ sG=="0",]), type="o",pch=20 )
title(xlab="Position (bp)" )
screen(1)
np <- k # which haplotype to "observe"
nk <- n1+n2+1
par( mar=c(0,5,1,1) )
plot(NULL,xlim=range(dvec), ylim=c(nk+0.5,0), xlab="", ylab="Haplotype", yaxt="n",main="GenCo modified", xaxt="n")
axis( 2, labels=c("g0",sXh,"obs"), at=0:nk, las=1)
axis( 4, labels=NA, at=1:nk )
abline(v=dvec, lty=1, col="grey90" )
cnt <- 1
for(x in c(refH1,refH2,sampHap[k])) {
    xcol <- ifelse( x==0, 4,2 )
    abline(h=cnt,col="grey80")
    points( dvec, rep(cnt,S), pch=20, col=xcol )
    cnt <- cnt+1
}
abline(h=c(1,cumsum(c(length(refH1),length(refH2) ))+1)-0.5, lwd=2)
# draw path:
pm <- do.call("rbind",strsplit(names(path),"-"))
colnames(pm) <- c("PX","X","PG","G")
points( dvec, as.numeric(gsub("^h","",pm[,"X"])), type="b", cex=2)
points( dvec, as.numeric(gsub("^[a-z]","",pm[,"G"])), type="b", cex=1.5, pch=5)

close.screen(all=TRUE)


dev.off()


plot( diff(dvec), exp( -diff(dvec)*T ) )

T <- 7
rho <- 1/1000
lam <- 1/500
gam <- 1/1000
s <- c(1,10,100,1000,10000,100000,1000000)
for(i in 1:length(s)) {
    d <- dvec0/s[i]
    par(mfrow=c(2,2), mar=c(5,5,2,0.5) )
    plot(diff(d), 1-exp(-diff(d)*T), main=paste("d /",s[i]),pch=20 )
    plot(diff(d), 1-exp(-diff(d)*rho), main=paste("d /",s[i]),pch=20 )
    plot(diff(d), 1-exp(-diff(d)*lam), main=paste("d /",s[i]),pch=20 )
    plot(diff(d), 1-exp(-diff(d)*gam*T/(n1+n2)), main=paste("d /",s[i]),pch=20 )
    locator(1)
}


d <- dvec0
rho <- 1/1000
lam <- 1/500
gam <- 1/1000
T <- c(1,7,20,100,200,500,1000)
for(i in 1:length(T)) {
    plot(diff(d), exp(-diff(d)*T[i]) )
    locator(1)
}

rho <- c(1/10000, 1/1000, 1/100, 1/10,1,10,100)
for(i in 1:length(rho)) {
    plot(diff(d), exp(-diff(d)*rho[i]) )
    locator(1)
}

gam <- c(1/10000, 1/1000, 1/100, 1/10,1,10,100)
for(i in 1:length(gam)) {
    plot(diff(d), exp(-diff(d)*gam[i]) )
    locator(1)
}

lam <- c(1/10000, 1/1000, 1/500, 1/100, 1/10,1,10,100)
for(i in 1:length(lam)) {
    plot(diff(d), exp(-diff(d)*lam[i]) )
    locator(1)
}

lam*exp(-lam*d)

##############################

library(lattice)
levelplot( t(pprob) )

levelplot( t() )


plot(NULL, xlim=range(dvec), ylim=c(0,1) )
for(i in 1:nrow(pprob)) {
    lines( dvec, pprob[i,], col=i)
#    Sys.sleep(1)
}

plot(NULL, xlim=range(dvec), ylim=c(0,1), ylab="P", xlab="Position (bp)",las=1)
points( dvec, colSums(pprob[sP=="p1",]), type="o", pch=20, col=2)
points( dvec, colSums(pprob[sP=="p2",]), type="o", pch=20, col=4)
points( dvec, colSums(pprob[sG!="0",]), type="o", pch=20, col=3)
legend("left",c("P(pop1)","P(pop2)","P(gconv)"), pch=20,lty=1,col=c(2,4,3), bty="n")



plot(NULL, xlim=range(dvec), ylim=c(0,1), ylab="P", xlab="Position (bp)",las=1)
points( dvec, colSums(pprob[sP=="p1"&sG!="0",]), type="o", pch=20, col=2)
points( dvec, colSums(pprob[sP=="p2"&sG!="0",]), type="o", pch=20, col=2)
points( dvec, colSums(pprob[sG=="0",]), type="o", pch=20, col=2)


plot(NULL, xlim=range(dvec), ylim=c(0,1), ylab="P", xlab="Position (bp)",las=1)
points( dvec, pprob[1,], type="o", pch=20, col=2)
points( dvec, pprob[3,], type="o", pch=20, col=2)
points( dvec, pprob[5,], type="o", pch=20, col=2)

points( dvec, pprob[2,], type="o", pch=20, col=2)
points( dvec, pprob[4,], type="o", pch=20, col=2)
points( dvec, pprob[6,], type="o", pch=20, col=2)

points( dvec, colSums(pprob[sP=="p2",]), type="o", pch=20, col=4)
points( dvec, colSums(pprob[sG!="0",]), type="o", pch=20, col=3)
legend("left",c("P(pop1)","P(pop2)","P(gconv)"), pch=20,lty=1,col=c(2,4,3), bty="n")


########################################

gam <- 1/1000 # gene conversion rate 
T <- 7
d <- 1
n1 <- 1:10
n2 <- 3
exp(-gam*T*d/n1) * exp(-gam*T*d/n2) - exp((-gam*T*d*(n1+n2))/(n1*n2))
exp(-gam*T*d/n1) * exp(-gam*T*d/n2) == exp((-gam*T*d*(n1+n2))/(n1*n2))
-gam*T*d/n1 + -gam*T*d/n2 == -gam*T*d*(n1+n2)/(n1*n2)
exp(1/n1) * exp(1/n2) - exp((n1+n2)/(n1*n2))
# x=1
# y=2
# 1/x + 1/y; (x+y)/(x*y)
# exp(1/x) * exp(1/y); exp( (x+y)/(x*y) )
# exp( (n1+n2)/(n1*n2) )
# exp( 1/n1 ) * exp( 1/n2 )





########################################
########################################
#hapmix:

# pop sizes:
n1 <- 10 # European
n2 <- 10 # African
# free parameters:
T
u1
# fixed parameters:
#p1 <- 0.05 # miscopying, inferred through EM
#p2 <- 0.05 # miscopying, inferred through EM
rho1 <- 60000/n1 # per morgan
rho2 <- 90000/n2 # per morgan
theta1 <- 0.2/(0.2+n1)
theta2 <- 0.2/(0.2+n2)
theta3 <- 0.01



(1-exp(-d*T)) * u1 * 1/n1

exp(-d*T) * (1-exp(-d*rho1)) * 1/n1 + (1-exp(-d*T)) * u1 * 1/n1

exp(-d*T) * exp(-d*rho1) + exp(-d*T) * (1-exp(-d*rho1)) * 1/n1 + (1-exp(-d*T)) * u1 * 1/n1


# P( no switch in ancestry):
T <- 0:10
d <- 1
plot( T, exp(-d*T), type="o", xlab="T (generations since admixture)", ylab=expression(e^{-d*T}), main=expression(e^{-d*T}==P(no~ancestry~switch)), las=1)

# P( no recomb w/in pop):
rho1 <- 0:100 # 60000 / n1 * 0:10
k <- 10
plot( rho1, exp(-d*rho1), type="o", xlab=expression(rho) )
plot( rho1, exp(-d*rho1/k), type="o", xlab=expression(rho) )

del <- 0:10
T <- 1
plot( del, exp(-d*del*T), type="o", xlab=expression(delta) )

plot( del, exp(-d*del*T), type="o" )






################################################################################
################################################################################

refH1$h2[2] <- 1
refH1$h2[39] <- 1
obs[39] <- 0
xpath <- ifelse( obs,2,5)
# xpath[39] <- 2
xpath[35:47] <- 1
gpath <- rep(0,length(obs))
gpath[39] <- 4
sampHap[[k]] <- obs


np <- k # which haplotype to "observe"
nk <- n1+n2+1
par( mar=c(4,5,1,1) )
plot(NULL,xlim=range(dvec), ylim=c(nk+0.5,0), xlab="", ylab="Haplotype", yaxt="n",main="" ) #, xaxt="n")
title(xlab="Position (bp)",line=2 )
axis( 2, labels=c("g0",sXh,"obs"), at=0:nk, las=1)
axis( 4, labels=NA, at=1:nk )
cnt <- 1
for(x in c(refH1,refH2,sampHap[k])) {
    xcol <- ifelse( x==0, 4,2 )
    abline(h=cnt,col="grey80")
    points( dvec, rep(cnt,S), pch=20, col=xcol )
    cnt <- cnt+1
}
abline(h=c(1,cumsum(c(length(refH1),length(refH2) ))+1)-0.5, lwd=2)
# draw path:
pm <- do.call("rbind",strsplit(names(path),"-"))
colnames(pm) <- c("PX","X","PG","G")

points( dvec, xpath, type="b", cex=2)
points( dvec, gpath, type="b", cex=1.5, pch=5)


dvec[33:35]

obs <- c( rep(0,3), rep(1,4), rep(0,2), rep(1,5) )
obs[6] <- 0
refH1 <- list()
for(i in 1:3) refH1[[i]] <- rep(0,length(obs))
refH2 <- list()
for(i in 1:3) refH2[[i]] <- rep(1,length(obs))
# dvec <- sample(1:30,length(obs))
dvec <- sort(dvec)
xpath <- c(2,2,2, 5,5,5,5, 3,3, 6,6,4,4,4)
gpath <- rep(0,length(obs))
gpath[6] <- 3


pdf( "results/hmmSchematic.pdf", width=8,height=3)
par( mar=c(1,6,1,1) )
plot(NULL,xlim=range(dvec), ylim=c(nk+0.25,-0.25), xlab="", ylab="", yaxt="n",main="", xaxt="n", bty="n")
#title(ylab="Position (bp)",line=2.5 )
#axis( 4, labels=NA, at=1:nk )
cnt <- 1
for(x in c(refH1,refH2,list(obs))) {
    xcol <- ifelse( x==0, 4,2 )
    abline(h=cnt,col="grey80", lwd=6)
    points( dvec, rep(cnt,14), pch=20, col=xcol, cex=3 )
    cnt <- cnt+1
}
# axis( 2, labels=c("g0",sXh,"h*"), at=0:nk, las=1, tcl=0)
#axis( 2, labels=c("gene\nconversion\nstate","Reference\npopulation 1","Reference\npopulation 2","admixed\nhaplotype"), at=c(0,2,5,7), las=1, tcl=0)
text(x=par()$usr[1]-diff(par()$usr[1:2])*0.01, y=c(0,2,5,7), labels= c("gene\nconversion\nstate","Reference\npopulation 1","Reference\npopulation 2","admixed\nhaplotype"), xpd=TRUE, adj=1)
abline(h=c(1,cumsum(c(length(refH1),length(refH2) ))+1)-0.5, lwd=1, lty=3)
points( dvec, xpath, type="b", cex=3)
points( dvec, gpath, type="b", cex=2.5, pch=5)
dev.off()










