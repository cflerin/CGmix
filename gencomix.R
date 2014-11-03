
setwd("~/Dropbox/lab/hmm/")
################################################################################
################################################################################
################################################################################
if(0) {
# load in manually simulated data:
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
} # end skip
################################################################################
################################################################################

################################################################################
################################################################################
# load 1000G data:

# vcftools --gzvcf ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep admix/indivs-to-keep_full --chr 22 --from-bp 20000000 --to-bp 21000000 --phased --remove-indels --recode --out admix/CEU-YRI_full
# vcftools --gzvcf CEU-YRI_full.vcf.gz --ldhat --chr 22 --out CEU-YRI_full

start <- Sys.time()

sitefname <- "data/CEU-YRI_full.ldhat.sites"
locfname <- "data/CEU-YRI_full.ldhat.locs"

n1h <- 2 # 99
n2h <- 2 # 107

tmp <- scan( sitefname, what=character(), skip=1, quiet=TRUE)
hapnames <- gsub("^>","",tmp[ seq(1,length(tmp),by=2) ])
tmp <- tmp[ seq(2,length(tmp),by=2) ]
hap <- lapply( strsplit(tmp,""), as.integer )
names(hap)[1:(n1h*2)] <- paste("p1-",hapnames[1:(n1h*2)],sep="")
names(hap)[(n1h*2+1):length(hap)] <- paste("p2-",hapnames[(n1h*2+1):length(hap)],sep="")

# reduce the haplotype size:
hsize <- 10 # 11200
hap <- lapply(hap,function(x) x[1:hsize] )

# cut down hap in size
ncut <- c(n1h*2,n2h*2) 
# pop 1:
refH1 <- hap[ grep("^p1",names(hap))[1:ncut[1]] ]
#names(refH1) <- sapply(strsplit(names(refH1),"-"),'[',2)
names(refH1) <- paste("h",1:ncut[1] ,sep="")
# pop 2:
refH2 <- hap[ grep("^p2",names(hap))[1:ncut[2]] ]
names(refH2) <- paste("h",(ncut[1]+1):sum(ncut) ,sep="")

# chr locations:
dvec <- scan( locfname, what=numeric(), skip=1, quiet=TRUE)
dvec <- dvec * 1000
dvec <- dvec[1:hsize]

# create an admixed individual:
sampHap <- list()
S <- length(hap[[1]])
nbreaks <- 4
bp <- sort(sample(S-2,nbreaks)) +1

refH1mix <- sample(1:length(refH1),1)
refH2mix <- sample(1:length(refH2),1)

m <- cbind(
    refH1[[refH1mix]] ,
    refH2[[refH2mix]]
    )

bp2 <- diff(c(0,bp,S))
sw <- inverse.rle( list( lengths=bp2, values=rep(c(1,2),S)[1:length(bp2)] ) )  # switch vector
sampHap[[1]] <- sapply(1:nrow(m),function(x) m[x,sw[x]] )
mixcol <- ifelse(sw==1, 4,2)
 
# remove "parent" haplotypes:
refH1 <- refH1[-refH1mix]
refH2 <- refH2[-refH2mix]

################################################################################
################################################################################


################################################################################
logsum <- function(x) {
    maxi <- which.max(x)
    maxl <- x[maxi]
    maxl + log1p(sum(exp(x[-maxi]-maxl))) 
}

########################################
########################################
# parameters:

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
theta <- 1/1000 # (1 per kb)
rho <- 1/1000000 # crossover rate (1 per kb)
gam <- 1/10000 # gene conversion rate 
lam <- 1/500 # mean GC tract length
#
# S <- length(sampHap[[1]])# how many biallelic loci (SNPs, segregating sites). 2^S possible haplotypes
n <- length(sampHap) # how many haplotypes in the sample ###haploid individuals

k=1
obs <- sampHap[[k]] # copy of observed haplotype

param <- list(
    n1 = n1, n2 = n2,
    T = T,
    u1 = u1,
    rho = rho ,
    gam = gam, lam=lam ,
    theta = theta 
    )


################################################################################
L <- 1
emit <- c(
    match = (2*(n1+n2)*L+theta) / (2*((n1+n2)*L+theta)) , # match
    mismatch = theta / ( 2*((n1+n2)*L+theta) ) # mismatch
)

sXh <- c(names(refH1),names(refH2) ) # haplotype (x) states
sXp <- c( rep("p1",length(refH1)), rep("p2",length(refH2)) ) # population states
sX <- paste(sXp,sXh,sep="-")
#
sGh <- c("0",sXh)
sGp <- c("0",sXp)
#sG <- paste( c("0",sXp), c("0",sXh),sep="-")
sG <- paste( sGp, sGh, sep="-")
#states <- apply( expand.grid( sX, sG ), 1, paste, collapse="-")
st <- expand.grid( sX=sX, sG=sG, stringsAsFactors=FALSE )
st <- cbind( st, 
    do.call("rbind",strsplit(st$sX,"-")) ,
    do.call("rbind",strsplit(st$sG,"-")) ,
    apply(st,1,paste,collapse="-") ,
    stringsAsFactors=FALSE
    )
colnames(st) <- c("X","G","Xpop","Xhap","Gpop","Ghap","states")

ref <- do.call("rbind", c(refH1,refH2) )

states <- st$states

##############################
# how to determine which haplotype to use for the emissions match/mismatch?
# either use Xstate or G state? or
#sXG <- sG
#sXG[sG=="0"] <- sX[sG=="0"]
sXG <- st$Ghap
sXG[st$Ghap=="0"] <- st$Xhap[st$Ghap=="0"]
#fwd:
fwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
### starting prob:
for(i in 1:length(states) ) {
    if( ref[ st$Xhap[i] , 1 ] == obs[1] ) {
        e <- emit[1]
    } else {
        e <- emit[2]
    }
    if( st$Ghap[i] == "0" ) {
        fwd[i,1] <- log( 1/(n1+n2) * (lam*(n1+n2)) / (lam*(n1+n2)+gam*T) * e )
    }
    if( st$Ghap[i] != "0" ) {
        fwd[i,1] <- log( 1/(n1+n2) * (gam*T) / ((n1+n2)*(lam*(n1+n2)+gam*T)) * e )
    }
}
for(j in 2:S) {
    #d <- diff( dvec[j:(j+1)] )
    d <- diff( dvec[ (j-1) : (j)] )
    cnt <- 1
    newemit <- c()
    for(t in 1:length(states)) { # hap loop "to" this state
    #for(state1 in states) { # hap loop
        state1 <- states[t]
        lsum <- -Inf
        newtrans <- c()
        newtmp <- c()
        newprevfwd <- c()
        for(f in 1:length(states)) { # "from" this state
        #for(state0 in states) {
            state0 <- states[f]
            if( st$Xpop[t] == "p1" ) { # if "to" state is p1, use u1
                u <- u1
                nm <- n1
            } 
            if( st$Xpop[t] == "p2" ) { # if "to" state is p2, use u2,n2:
                u <- 1 - u1 
                nm <- n2
            }
            ##################
            ### haplotype transitions: 
            if( st$Xpop[f] != st$Xpop[t] ) { # anc AND hap switch:
                trX <- (1-exp(-d*rho*T)) * u/nm
            } else { # no ancestry switch
                if( st$Xhap[f] != st$Xhap[t] ) { # no anc switch, hap switch
                    trX <- exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm
                }
                if( st$Xhap[f] == st$Xhap[t] ) { # no anc switch, no hap switch
                    trX <- exp(-d*rho*T) * exp(-d*rho) + exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm
                }
            }
            ##################
            ### gene conversion transitions:
            if( st$Gpop[t] == "p1" ) {
                u <- u1
                nm <- n1
            }
            if( st$Gpop[t] == "p2" ) {
                u <- 1 - u1
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
            if( st$Ghap[f] == st$Ghap[t] )              trG <- a4 # erroneously includes 0 -> 0; overwritten on next line
            if( st$Ghap[f] == "0" & st$Ghap[t] == "0" ) trG <- a1
            if( st$Ghap[f] == "0" & st$Ghap[t] != "0" ) trG <- a2
            if( st$Ghap[f] != "0" & st$Ghap[t] == "0" ) trG <- a3
            if( st$Ghap[f] != "0" & st$Ghap[t] != "0" & st$Ghap[f] != st$Ghap[t] ) trG <- a5
            ###
            tmp <- fwd[ f , j-1 ] + log( trX * trG )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
            newtrans <- c(newtrans, trX * trG )
            newtmp <- c(newtmp,tmp)
            newprevfwd <- c( newprevfwd, fwd[ f , j-1 ] )
        } # end from/lsum loop
        ### emission prob:
        #if( ref[ st$Xhap[t] , j ] == obs[j] ) { # how should we determine emissions match/mismatch?? use gene conversion or transition state... they often differ!
        if( ref[ sXG[cnt] , j ] == obs[j] ) {
            e <- emit[1]
            newemit <- rbind(newemit,c(st$Xhap[t],"match") )
        } else {
            e <- emit[2]
            newemit <- rbind(newemit,c(st$Xhap[t],"mismatch") )
        }
        cnt <- cnt+1
        fwd[ t , j ] <- log( e ) + lsum
        #fwd[state1,j] <- log( emit[ emat[hname,j] ] ) + lsum
    } # end to loop
} # end site loop

Pxa <- logsum(fwd[,S])


all( fwd[,1] == fwd0[,1] )
cbind( oldtmp, newtmp, oldtmp==newtmp )
cbind( oldtrans, newtrans, oldtrans==newtrans )
cbind( oldprevfwd, newprevfwd, oldprevfwd==newprevfwd )
cbind( fwd0[,1], fwd[,1] )

########################################
########################################

cat("Calculating transition probabilities...\t"); start1 <- Sys.time()
# transition matrices for entire sequence:
transL <- list()
for(j in 1:(S-1) ) { # site along haplotype sequence
d <- diff( dvec[j:(j+1)] )

nk <- length(refH1)+length(refH2) # number of total haplotypes in both reference populations
sXh <- c(names(refH1),names(refH2) ) # haplotype (x) states
sXp <- c( rep("p1",length(refH1)), rep("p2",length(refH2)) ) # population states
sX <- paste(sXp,sXh,sep="-")

### Haplotype (X) crossover transition probabilities:
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

(1-exp(-d*rho*T)) * u/nm # anc AND hap switch

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
rm(trans)

cat("done in ",format(Sys.time()-start1),"\n")

##################################################
##################################################


L <- 1
emit <- c(
    match = (2*(n1+n2)*L+theta) / (2*((n1+n2)*L+theta)) , # match
    mismatch = theta / ( 2*((n1+n2)*L+theta) ) # mismatch
)

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
eIndx <- emat[ sXG, 1 ]
#eIndx <- emat[ st$Xhap, 1 ]
# sprob[g0i] <- log( 1/(n1+n2) * (lam*(n1+n2)) / (lam*(n1+n2)+gam) * emit[eIndx[g0i]] )
# sprob[g1i] <- log( 1/(n1+n2) * gam / ((n1+n2)*(lam*(n1+n2)+gam)) * emit[eIndx[g1i]] )
# 
sprob[g0i] <- log( 1/(n1+n2) * 
                  (lam*(n1+n2)) / (lam*(n1+n2)+gam*T) 
                  * emit[eIndx[g0i]] )
sprob[g1i] <- log( 1/(n1+n2) * 
                  (gam*T) / ((n1+n2)*(lam*(n1+n2)+gam*T)) 
                  * emit[eIndx[g1i]] )
### add noise:
#sprob <- sprob -rnorm( length(states), sd=0.0001 )


#sprob[8] <- sprob[8]+2
#sprob[9:length(sprob)] <- 0


cat("Calculating forward probabilities...\t"); start1 <- Sys.time()

#fwd:
fwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
fwd[,1] <- sprob
for(j in 2:S) {
    cnt <- 1
    oldemit <- c()
    for(state1 in states) { # hap loop
        lsum <- -Inf
        oldtrans <- c()
        oldtmp <- c()
        oldprevfwd <- c()
        for(state0 in states) {
            #tmp <- fwd[state0,j-1] + log( trans[state0,state1] )
            cat( "state1=",state1, "state0=", state0, transL[[j-1]][state0,state1], "\n" )
            tmp <- fwd[state0,j-1] + log( transL[[j-1]][state0,state1] )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
            oldtrans <- c(oldtrans, transL[[j-1]][state0,state1] )
            oldtmp <- c(oldtmp, tmp ) #transL[[j-1]][state0,state1] )
            oldprevfwd <- c( oldprevfwd, fwd[state0,j-1] )
        }
        hname <- sXG[cnt]; cnt <- cnt+1
        #fwd[state1,j] <- log( emit[ ifelse( c(sapply(refH1,'[',j),sapply(refH2,'[',j))[hname] ==obs[j],1,2) ] ) +lsum
        fwd[state1,j] <- log( emit[ emat[hname,j] ] ) + lsum
        oldemit <- rbind(oldemit, c( sXG[cnt-1], c("match","mismatch")[ emat[hname,j] ] ) )
    }
}
Pxa <- logsum(fwd[,S])

cat("done in ",format(Sys.time()-start1),"\n")
cat("Calculating backward probabilities...\t"); start1 <- Sys.time()

#bwd:
bwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
bwd[,S] <- 0
for(j in (S-1):1) {
    for(state1 in states) {
        lsum <- -Inf
        cnt <- 1
        for(state2 in states) {
            st2 <- sXG[cnt]; cnt <- cnt+1
            eIndx <- emat[st2,j+1]
            tmp <- bwd[state2,j+1] +
                log( transL[[j]][state1,state2] * emit[eIndx] )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
        }
        bwd[state1,j] <- lsum
    }
} # seq loop j

cat("done in ",format(Sys.time()-start1),"\n")

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
pprob <- ( fwd+bwd - Pxa )
for(i in 1:ncol(pprob)) {
    pprob[,i] <- pprob[,i]-max(pprob[,i])
}
pprob <- exp(pprob)
for(i in 1:ncol(pprob)) {
    pprob[,i] <- pprob[,i]/sum(pprob[,i])
}
path <- apply(pprob,2,max)
names(path) <- states[ apply(pprob,2,which.max) ]


cat("Calculating Viterbi probabilities...\t"); start1 <- Sys.time()

# viterbi:
vit <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
# initial state:
vit[,1] <- sprob
# recursion:
for(j in 2:S) {
    cnt <- 1
    for(state1 in states) {
        tmp <- vit[,j-1] + log(transL[[j-1]][states,state1])
        vmax <- max(tmp)
        #vit[state1,i] <- log(emit[state1,obs[i]]) + vmax
        hname <- sXG[cnt]; cnt <- cnt+1
        vit[state1,j] <- log( emit[ emat[hname,j] ] ) + vmax
    }
}
# termination:
vpath <- rep(NA,length(obs))
vpath[S] <- states[ which.max( vit[,length(obs)] ) ]
# traceback:
for(i in (length(obs)-1):1 ) {
    vpath[i] <- states[ which.max( vit[states,i] + log(transL[[i]][states,vpath[i+1]]) ) ]
}
cat("done in ",format(Sys.time()-start1),"\n")


path <- numeric(S)
names(path) <- vpath


emat[,47:50]
fwd[1:8,47:50]
bwd[1:8,47:50]
pprob[1:8,47:50]
transL[[48]][1:8,1:8]

print( fwd[1:8,11:13],digits=22 )
print( bwd[1:8,11:13],digits=22 )
print( pprob[1:8,11:13],digits=22 )


paste( paste(names(param),param,sep="=") ,collapse=";")

runtime <- Sys.time() - start
cat("Total runtime of ",format(runtime),"\n")

pdf( paste( "results/admix1000G_n",sum(ncut),"_S",S,".pdf",sep=""), width=15,height=6 )
##################################################
##################################################
# plot:

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
#######
screen(1)
np <- k # which haplotype to "observe"
nk <- n1+n2+1
par( mar=c(0,5,1,1) )
plot(NULL,xlim=range(dvec), ylim=c(nk+0.5,0), xlab="", ylab="Haplotype", yaxt="n",main="", xaxt="n")
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
xl <- dvec-diff(dvec)[c(1,1:(length(dvec)-1))]
xr <- dvec+diff(dvec)[c(1:(length(dvec)-1),length(dvec)-1)]
rect(xleft=xl,xright=xr,ybottom=cnt-0.75,ytop=cnt-0.5,border=NA,col=mixcol)
abline(h=c(1,cumsum(c(length(refH1),length(refH2) ))+1)-0.5, lwd=2)
# draw path:
pm <- do.call("rbind",strsplit(names(path),"-"))
colnames(pm) <- c("PX","X","PG","G")
points( dvec, as.numeric(gsub("^h","",pm[,"X"])), type="b", cex=2)
points( dvec, as.numeric(gsub("^[a-z]","",pm[,"G"])), type="b", cex=1.5, pch=5)
#
close.screen(all=TRUE)

dev.off()

pprobcol <- pprob[1:8,]
pprobcol <- bwd[1:8,]
u <- sort(unique(as.vector(bwd[1:8,])))
hcvec <- rainbow( length(u), alpha=0.6, start=0.5, end=1 )
for(i in 1:length(hcvec)) { pprobcol[pprobcol==u[i]] <- hcvec[i] }

#par(mar=c(5,5,1,8))
par( mar=c(5,4.5,1,8), oma=c(0,0,0,1) )
plot( NULL, xlim=c(0,ncol(pprob)), ylim=c(n1+n2+1.5,0), xlab="", ylab="Haplotype", yaxt="n" )
axis( 2, labels=c("g0",sXh,"obs"), at=0:nk, las=1)
axis( 4, labels=NA, at=1:nk )
cnt <- 1
for(x in c(refH1,refH2)) { #,sampHap[k])) {
    #xcol <- ifelse( x==0, 4,2 )
    xcol <- ifelse( x==0, "-","+" )
    abline(h=cnt,col="grey80")
    rect(xleft=c(1:ncol(pprob))-0.5,xright=c(1:ncol(pprob))+0.5, ybottom=cnt-0.5,ytop=cnt+0.5, col=pprobcol[cnt,], border="grey40" )
    #points( rep(cnt,S), pch=20, col=xcol )
    points( rep(cnt,S), pch=xcol )
    cnt <- cnt+1
}
abline(h=cnt,col="grey80")
points( rep(cnt,S), pch=ifelse(obs==0,"-","+") )
# draw path:
points( as.numeric(gsub("^h","",pm[,"X"])), type="b", cex=2,lwd=2)
# key:
yr <- par()$usr[3:4] #c(0,1)
left <- seq( yr[1], yr[2], length.out=length(hcvec) )
rect( xleft=diff(par()$usr[1:2])*1.01, xright=diff(par()$usr[1:2])*1.07, ybottom=left[-length(hcvec)], ytop=left[-1], xpd=TRUE, col=hcvec, border=NA )
# linear:
axl <- round( seq(min(u),max(u),length.out=4) ,3)
axl <- seq(min(u),max(u),length.out=4)
### axl <- pretty(u)
### axl[1] <- round(min(u),3)
### axl[length(axl)] <- round(max(u),3)
axR <- rangeScale( axl, yr[1], yr[2] )
axis(4,las=1,line=3.5,at=axR,labels=format(axl,digits=3) )
#mtext("G statistic",side=4,outer=TRUE)




axl <- pretty(log10(u))
axl <- c( ( trunc(min(log10(u))) : trunc(max(log10(u))) ), log10(max(u)) )
axl[1] <- log10(min(u))
#if( diff(axl[length(axl)+-1:0]) < 0.1 ) axl <- axl[-length(axl)]
axR <- rangeScale( axl, min(left), max(left) )
axis(4,las=1,line=2.5,at=axR,labels=10^axl )
#mtext("# of 50kb bins",side=4,outer=TRUE)





