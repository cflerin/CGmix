
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

n1h <- 3 #10 # 99
n2h <- 3 #10 # 107

tmp <- scan( sitefname, what=character(), skip=1, quiet=TRUE)
hapnames <- gsub("^>","",tmp[ seq(1,length(tmp),by=2) ])
tmp <- tmp[ seq(2,length(tmp),by=2) ]
hap <- lapply( strsplit(tmp,""), as.integer )
names(hap)[1:(n1h*2)] <- paste("p1-",hapnames[1:(n1h*2)],sep="")
names(hap)[(n1h*2+1):length(hap)] <- paste("p2-",hapnames[(n1h*2+1):length(hap)],sep="")

# reduce the haplotype size:
hsize <- 20 # 100 # 1000 # 20 # 11200
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
#bp <- sort(sample(S-2,nbreaks)) +1
#bp <- c(4,7,12) # 20
#bp <- c(69,117,193,947) # 1000
#bp <- c(13,68,70,106) # 200
bp <- c(9,12,14,21) # 200


refH1mix <- 11 #5 # sample(1:length(refH1),1)
refH2mix <- 15 #38 #5 # sample(1:length(refH2),1)

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


# determine haplotype (X) and gene conversion (G) transitions:
getXtrans <- function( to, from, st, param ) {
# from = index to $st for state we are transitioning from
# to   = index to $st for state we are transitioning to
# st   = matrix of states, labels and haplotype/gc labels
# params = parameters
    ### haplotype transitions: 
    if( st$Xpop[to] == "p1" ) { # if "to" state is p1, use u1
        u <- param$u1
        nm <- param$n1
    } 
    if( st$Xpop[to] == "p2" ) { # if "to" state is p2, use u2,n2:
        u <- 1 - param$u1 
        nm <- param$n2
    }
    if( st$Xpop[from] != st$Xpop[to] ) { # anc AND hap switch:
        trX <- with( param, (1-exp(-d*rho*T)) * u/nm )
    } else { # no ancestry switch
        if( st$Xhap[from] != st$Xhap[to] ) { # no anc switch, hap switch
            trX <- with( param, exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm )
        }
        if( st$Xhap[from] == st$Xhap[to] ) { # no anc switch, no hap switch
            trX <- with( param, exp(-d*rho*T) * exp(-d*rho) + exp(-d*rho*T) * (1-exp(-d*rho))/nm + (1-exp(-d*rho*T))*u/nm )
        }
    }
    return(trX)
}

getGtrans <- function( to, from, st, param ) {
    ### gene conversion transitions:
    u <- 0
    nm <- 0
    if( st$Gpop[to] == "p1" ) {
        u <- param$u1
        nm <- param$n1
    }
    if( st$Gpop[to] == "p2" ) {
        u <- 1 - param$u1
        nm <- param$n2
    }
    # 3: Pr(G[j+1]=0 | G[j] = g):
    a3 <- with( param, lam*(n1+n2)* (1-exp(-d*(gam*T+lam*(n1+n2))/(n1+n2))) / (gam*T+lam*(n1+n2)) )
    # 5: Pr(G[j+1]=g' | G[j] = g):
    a5 <- with( param, (lam*u*(n1+n2) * (exp(d*(-gam*T/(n1+n2)-lam))-1)) / ((n1+n2)*nm*lam + gam*T*nm) + (u-u*exp(-d*lam))/nm )
    # 1: Pr(G[j+1]=0 | G[j] = 0):
    a1 <- with( param, exp(-lam*d) * exp(-gam*T*d/(n1+n2)) + a3 )
    # 2: Pr(G[j+1]=g | G[j] = 0):
    a2 <- with( param, exp(-lam*d) * (1-exp(-gam*T*d/(n1+n2))) *u/nm + a5 )
    # 4: Pr(G[j+1]=g | G[j] = g):
    a4 <- with( param, exp(-lam*d) + a5 )
    if( st$Ghap[from] == st$Ghap[to] )              trG <- a4 # erroneously includes 0 -> 0; overwritten on next line
    if( st$Ghap[from] == "0" & st$Ghap[to] == "0" ) trG <- a1
    if( st$Ghap[from] == "0" & st$Ghap[to] != "0" ) trG <- a2
    if( st$Ghap[from] != "0" & st$Ghap[to] == "0" ) trG <- a3
    if( st$Ghap[from] != "0" & st$Ghap[to] != "0" & st$Ghap[from] != st$Ghap[to] ) trG <- a5
    return(trG)
}

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
sG <- paste( sGp, sGh, sep="-")
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

write.table( rbind(ref,obs), file=paste("src/admixSampleData_j",S,".sites",sep=""), sep="\t", col.names=FALSE,row.names=FALSE)
write.table( matrix(dvec,nrow=1), file=paste("src/admixSampleData_j",S,".locs",sep=""), sep="\t", col.names=FALSE,row.names=FALSE)
write.table( cbind( c(rownames(ref),"obs"), c(sXp,"p3") ), file=paste("src/admixSampleData_j",S,".hapnames",sep=""), sep="\t", col.names=FALSE,row.names=FALSE,quote=FALSE)
##############################
# how to determine which haplotype to use for the emissions match/mismatch?
# either use Xstate or G state? or
sXG <- st$Ghap
sXG[st$Ghap=="0"] <- st$Xhap[st$Ghap=="0"]

cat("Calculating forward probabilities...\t"); start1 <- Sys.time()
##############################
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
    d <- diff( dvec[ (j-1) : (j)] )
    cnt <- 1
    for(t in 1:length(states)) { # hap loop "to" this state
        state1 <- states[t]
        lsum <- -Inf
        for(f in 1:length(states)) { # "from" this state
            state0 <- states[f]
            ### haplotype transitions: 
            trX <- getXtrans( from=f, to=t, st=st, param=param )
            ### gene conversion transitions:
            trG <- getGtrans( from=f, to=t, st=st, param=param )
            ###
            tmp <- fwd[ f , j-1 ] + log( trX * trG )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
        } # end from/lsum loop
        ### emission prob:
        if( ref[ sXG[cnt] , j ] == obs[j] ) {
            e <- emit[1]
        } else {
            e <- emit[2]
        }
        cnt <- cnt+1
        fwd[ t , j ] <- log( e ) + lsum
    } # end to loop
} # end site loop
cat("done in ",format(Sys.time()-start1),"\n")

Pxa <- logsum(fwd[,S])

cat("Calculating backward probabilities...\t"); start1 <- Sys.time()
##############################
#bwd:
bwd <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
bwd[,S] <- 0
for(j in (S-1):1) {
    d <- diff( dvec[ (j) : (j+1)] )
    for(f in 1:length(states)) {
        state1 <- states[f]
        lsum <- -Inf
        cnt <- 1
        for(t in 1:length(states)) {
            state2 <- states[t]
            ### haplotype transitions: 
            trX <- getXtrans( from=f, to=t, st=st, param=param )
            ### gene conversion transitions:
            trG <- getGtrans( from=f, to=t, st=st, param=param )
            ###
            if( ref[ sXG[cnt] , j+1 ] == obs[j+1] ) {
                e <- emit[1]
            } else {
                e <- emit[2]
            }
            cnt <- cnt+1
            tmp <- bwd[t,j+1] + log( trX * trG * e )
            if(tmp>-Inf) lsum <- tmp + log(1+exp(lsum-tmp))
        }
        bwd[f,j] <- lsum
    }
} # seq loop j
cat("done in ",format(Sys.time()-start1),"\n")

Pxb <- bwd[1,1] + fwd[1,1]; cnt <- 2
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
# path <- apply(pprob,2,max)
# names(path) <- states[ apply(pprob,2,which.max) ]

cat("Calculating Viterbi probabilities...\t"); start1 <- Sys.time()
##############################
# viterbi:
vit <- matrix( as.numeric(NA), nrow=length(states), ncol=S, dimnames=list(states=states, obs=obs) )
# initial state:
vit[,1] <- fwd[,1]
# recursion:
for(j in 2:S) {
    d <- diff( dvec[ (j-1) : (j)] )
    cnt <- 1
    for(t in 1:length(states)) {
        state1 <- states[t]
        tmp <- numeric(length(states))
        for(f in 1:length(states)) { # "from" this state
            state0 <- states[f]
            ### haplotype transitions: 
            trX <- getXtrans( from=f, to=t, st=st, param=param )
            ### gene conversion transitions:
            trG <- getGtrans( from=f, to=t, st=st, param=param )
            ###
            tmp[f] <- vit[f,j-1] + log( trX * trG )
        } # end from loop
        vmax <- max(tmp)
        if( ref[ sXG[cnt] , j ] == obs[j] ) {
            e <- emit[1]
        } else {
            e <- emit[2]
        }
        cnt <- cnt+1
        vit[t,j] <- log( e ) + vmax
    } # end to loop
} # end site loop
# termination:
vpath <- rep(NA,length(obs))
vprob <- rep(NA,length(obs))
vpath[S] <- states[ which.max( vit[,length(obs)] ) ]
vprob[S] <- pprob[ which.max( vit[,length(obs)] ) , S ]
# traceback:
for(j in (S-1):1 ) {
    # find previous max state:
    t <- which( states == vpath[j+1] )
    state1 <- states[t]
    tmp <- numeric(length(states))
    for(f in 1:length(states)) { # "from" this state
        state0 <- states[f]
        ### haplotype transitions: 
        trX <- getXtrans( from=f, to=t, st=st, param=param )
        ### gene conversion transitions:
        trG <- getGtrans( from=f, to=t, st=st, param=param )
        ###
        tmp[f] <- vit[f,j] + log( trX * trG )
    } # end from loop
    vpath[j] <- states[ which.max(tmp) ]
    vprob[j] <- pprob[ which.max(tmp), j ]
    #vpath[j] <- states[ which.max( vit[states,j] + log(transL[[j]][states,vpath[j+1]]) ) ]
}
cat("done in ",format(Sys.time()-start1),"\n")

path <- numeric(S)
names(path) <- vpath

########################################
########################################
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
plot(dvec, 1-colSums(pprob[ st$Ghap=="0",]), xlab="", ylim=c(0,1), las=1, type="n",pch=20, ylab="P(gene conv)")
abline(v=dvec, lty=1, col="grey90" )
points(dvec, 1-colSums(pprob[ st$Ghap=="0",]), type="o",pch=20 )
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





