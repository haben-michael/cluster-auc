auc <- function(x,y,normalized=TRUE) {
    psi <- outer(x,y,function(a,b)(a<b)+1/2*(a==b))
    if(normalized)mean(psi) else sum(psi)
}


## population AUC variance estimator from Obuchowski '97
auc.obu <- function(x,y,alpha=.05) {
    
    I <- length(x)
    stopifnot(length(x)==length(y))
    m <- sapply(x,length)
    n <- sapply(y,length)
    I.10 <- sum(m>0)
    I.01 <- sum(n>0)
    M <- sum(m)
    N <- sum(n)

    theta.hat <- auc(unlist(x),unlist(y),normalized=FALSE) / (M*N)
    V.10 <- sapply(x,function(x.i)auc(x.i,unlist(y),normalized=FALSE))/N
    V.01 <- sapply(y,function(y.i)auc(unlist(x),y.i,normalized=FALSE))/M

    S.10 <- sum((V.10 - m*theta.hat)^2) * I.10 / ((I.10-1)*M)
    S.01 <- sum((V.01 - n*theta.hat)^2) * I.01 / ((I.01-1)*N)
    S.11 <- sum((V.10 - m*theta.hat)*(V.01 - n*theta.hat)) * I / (I-1)

    var.hat <- S.10/M + S.01/N + 2*S.11/(M*N)
    q <- qnorm(1-alpha/2)
    return(c(theta.hat=theta.hat,var.hat=var.hat,CI.lower=theta.hat-q*sqrt(var.hat),CI.upper=theta.hat+q*sqrt(var.hat)))

}


draw.ellipse <- function(mean.xy=c(0,0),cov.xy=diag(2),alpha=.05,resolution=1e2,...) {
    q <- qchisq(1-alpha,df=2)
    cov.xy.sqrt <- with(eigen(cov.xy), vectors%*%diag(sqrt(values))%*%t(vectors))
    theta <- seq(0,2*pi,len=resolution)
    circle <- rbind(cos(theta),sin(theta))
    ellipse <- mean.xy + cov.xy.sqrt%*%(sqrt(q)*circle)
    ## browser()
    lines(t(ellipse),...)
}

##check if point lies in ellipse given parametrically as center +
## sigma * circle(theta)
inside.ellipse <- function(pt,center,qform,radius) {
    pt <- as.matrix(pt,nrow=2)
    qform.inv.sqrt <- with(eigen(solve(qform)), vectors%*%diag(sqrt(values))%*%t(vectors))
    colSums((qform.inv.sqrt%*%(pt-center))^2) < radius^2
}

## sample according to normal model with exchangeable correlation
## structure, random effect, and location shift
rbinormal_old <- function(I,m,n,mean.x,mean.y,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1,var.y=1,mean.z=rep(0,I),var.z=0,plot=FALSE) {

    stopifnot(sum(m<=0)==0 && sum(n<=0)==0)
        
    data <- lapply(1:I, function(i) {
        ## print(i)
        ## m <- ms[i]; n <- ns[i]
        mean.xy <- c(rep(mean.x,m[i]),rep(mean.y,n[i]))
        vcov.xy <- rbind(cbind(matrix(cov.xx,m[i],m[i]),matrix(cov.xy,m[i],n[i])),
                         cbind(matrix(cov.xy,n[i],m[i]),matrix(cov.yy,n[i],n[i]))) +
            diag(c(rep(var.x-cov.xx,m[i]),rep(var.y-cov.yy,n[i])))
        xy <- rmvnorm(1,mean=mean.xy,sigma=vcov.xy) + rnorm(1,mean=mean.z[i],sd=sqrt(var.z))
        x <- xy[1:(m[i])]; y <- xy[(m[i]+1):(m[i]+n[i])]
        if(sum(is.na(unlist(y)))>0) browser()
        list(x=x,y=y)
    })
    x <- lapply(data,function(xy)xy$x)
    y <- lapply(data,function(xy)xy$y)

    if(plot) {
        plot(0,type='n',xlim=range(unlist(data)),ylim=c(0,I))
        invisible(sapply(1:I, function(i) points(c(x[[i]],y[[i]]),y=rep(i,m[i]+n[i]),col=rep(c(1,2),c(m[i],n[i])))))
        points(c(unlist(x),unlist(y)),y=rep(0,sum(m)+sum(n)),col=rep(c(1,2),c(sum(m),sum(n))))
    }

    return(list(x=x,y=y))
}


## Sample according to normal model with exchangeable correlation
## structure, random effect, and location shift update: improving
## efficiency, getting rid of Z--redundant to cov.xy. could probably
## scale variances and get rid of var.x and var.y.
rbinormal <- function(I,m,n,mean.x,mean.y,cov.xx=cov.xy,cov.xy=0,cov.yy=cov.xy,var.x=1,var.y=1,plot=FALSE) {

    stopifnot(sum(m<=0)==0 && sum(n<=0)==0)
    stopifnot(var.x>=cov.xx && var.y>=cov.yy && min(cov.xx,cov.yy)>=cov.xy)

    x <- matrix(rnorm(I*max(m),sd=sqrt(var.x-cov.xx)),I) + mean.x
    y <- matrix(rnorm(I*max(n),sd=sqrt(var.y-cov.yy)),I) + mean.y
    re <- rnorm(I,sd=sqrt(abs(cov.xy)))
    x <- x + sign(cov.xy)*re + rnorm(I,sd=sqrt(cov.xx-cov.xy))
    y <- y + sign(cov.xy)*re + rnorm(I,sd=sqrt(cov.yy-cov.xy))
    x <- lapply(1:I,function(i)x[i,1:m[i]])
    y <- lapply(1:I,function(i)y[i,1:n[i]])

    if(plot) {
        plot(0,type='n',xlim=range(c(unlist(x),unlist(y))),ylim=c(0,I))
        invisible(sapply(1:I, function(i) points(c(x[[i]],y[[i]]),y=rep(i,m[i]+n[i]),col=rep(c(1,2),c(m[i],n[i])))))
        points(c(unlist(x),unlist(y)),y=rep(0,sum(m)+sum(n)),col=rep(c(1,2),c(sum(m),sum(n))))
    }

    return(list(x=x,y=y))
}


auc.to.params.binormal <- function(theta.12,theta.11) {
    Delta <- sqrt(2)*qnorm(theta.12)
    rho <- 1 - 1/2*Delta^2 / qnorm(theta.11)^2
    return(c(Delta=Delta,rho=rho))
}







## binormal model with observations censored at +/- bound. assumes unit variances, zero mean control variance, etc.
rbinormal.censored <- function(I,m,n,bound=abs(qnorm(.2)),mean.x=0,mean.y,cov.xx=cov.xy,cov.xy=0,cov.yy=cov.xy,var.x=1,var.y=1,plot=FALSE) {
    stopifnot(sum(m<=0)==0 && sum(n<=0)==0)
    stopifnot(bound>0 & mean.x==0 & mean.y > mean.x & var.x==1 & var.y==1)    
    bnd <- bound
    Delta <- mean.y
    rho <- cov.xy
    xy <- lapply(1:I, function(i) sqrt(rho)*rnorm(1) + sqrt(1-rho)*rnorm(m[i]+n[i]))
    xy <- lapply(1:I, function(i) c(xy[[i]][1:m[i]],  xy[[i]][(m[i]+1):(m[i]+n[i])] + Delta))
    xy <- lapply(xy, function(row) {
        row[row < -bnd] <- -bnd
        row[row > bnd] <- bnd
        row
        })
    x <- lapply(1:I, function(i) xy[[i]][1:m[i]])
    y <- lapply(1:I, function(i) xy[[i]][(m[i]+1):(m[i]+n[i])])

    if(plot) {
        plot(0,type='n',xlim=range(c(unlist(x),unlist(y))),ylim=c(0,I))
        invisible(sapply(1:I, function(i) points(c(x[[i]],y[[i]]),y=rep(i,m[i]+n[i]),col=rep(c(1,2),c(m[i],n[i])))))
        points(c(unlist(x),unlist(y)),y=rep(0,sum(m)+sum(n)),col=rep(c(1,2),c(sum(m),sum(n))))
    }

    return(list(x=x,y=y))
}

# rho=0 case is for the population AUC
auc.binormal.censored <- function(bnd,delta,rho) {
    if(rho==0)return (
                  -integrate(function(x)dnorm(x)*pnorm(x-delta),-bnd,bnd)$val + 1/2*(pnorm(bnd)-pnorm(bnd-delta)-pnorm(-bnd-delta)+pnorm(bnd)*(pnorm(-bnd-delta)+pnorm(bnd-delta))+1)
              )
    f <- Vectorize(function(x,y)mvtnorm::dmvnorm(c(x,y),mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2)))
    F <- function(lower,upper)mvtnorm::pmvnorm(lower,upper,mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2))
    interior <- integrate(Vectorize(function(x)integrate(function(y)f(x,y),x,bnd)$val),-bnd,bnd)$val
    boundary  <- pnorm(-bnd) + 1-pnorm(bnd-delta) - F(c(-Inf,-Inf),c(-bnd,-bnd)) - F(c(bnd,bnd),c(Inf,Inf)) - F(c(-Inf,bnd),c(-bnd,Inf))
    diag <- F(c(-Inf,-Inf),c(-bnd,-bnd)) + F(c(bnd,bnd),c(Inf,Inf))
    as.numeric(interior + boundary + 1/2*diag)
}


auc.to.params.binormal.censored <- function(theta.12,theta.11,bound) {
    auc.pop <- theta.12; auc.indiv <- theta.11
    Delta.star <- uniroot(function(Delta)auc.binormal.censored(bound,Delta,rho=0) - theta.12,c(0,5),extendInt='yes')$root
    obj <- Vectorize(function(rho)abs(auc.binormal.censored(bound,Delta.star,rho)-theta.11))
    rho.star <- optimize(f=obj,interval=c(0,1))$min
return(c(Delta=Delta.star,rho=rho.star))
}


## get population and personalized AUC point estimates and variances
auc.cluster <- function(x,y,get.vcov=TRUE,alpha=.05) {
    stopifnot(length(x)==length(y))
    I <- length(x)
    psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
    m <- sapply(x,length); n <- sapply(y,length)
    M <- sum(m); N <- sum(n)
    phi <- diag(psi)/m/n
    theta.11.hat <- mean(diag(psi)/m/n)
    theta.12.hat <- (sum(psi)-sum(diag(psi)))/M/N * I/(I-1)
    out <- if(get.vcov ) {
        cov.hat <- mean(phi/mean(phi)*((colMeans(psi)+rowMeans(psi))/mean(psi) - m/mean(m)-n/mean(n))) * theta.11.hat*theta.12.hat
        var.11.hat <- var(phi)
        var.12.hat <- theta.12.hat^2*I*sum((1/theta.12.hat/M/N*(colSums(psi)+rowSums(psi)) - (m/M+n/N))^2) * (I/(I-1))^2
        q <- qnorm(1-alpha/2)
        list(theta.11.hat=theta.11.hat,theta.12.hat=theta.12.hat,vcov.hat=matrix(c(var.11.hat,cov.hat,cov.hat,var.12.hat),2),
             theta.11.CI=theta.11.hat+c(-1,1)*q*sqrt(var.11.hat)/sqrt(I),
             theta.12.CI=theta.12.hat+c(-1,1)*q*sqrt(var.12.hat)/sqrt(I),I=I,M=M,N=N)
    } else {
        c(theta.11.hat=theta.11.hat,theta.12.hat=theta.12.hat)
    }
    structure(out,class='auc.cluster')
}


plot.auc.cluster <- function(x,alpha=.05,resolution=1e2,add=FALSE,...) {
    auc0 <- x
    if(!add) {
        plot(0,type='n',xlim=c(0,1),ylim=c(0,1),asp=1,xlab=expression(paste(theta[11])),ylab=expression(paste(theta[12])),...)
        abline(h=1/2,v=1/2,lty=2)
        abline(a=0,b=1)
    }
    with(auc0, draw.ellipse(mean.xy=c(theta.11.hat,theta.12.hat), cov.xy=vcov.hat/I,alpha=alpha,resolution=resolution)  )
    with(auc0, points(theta.11.hat,theta.12.hat,pch=20,cex=.3))
}
