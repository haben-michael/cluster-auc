## 1. figures showing large indiv auc vs small pop auc and vice versa

## large indiv auc, small pop auc

require(mvtnorm)
source('auc.cluster.R')
I <- 10
max.cluster.size <- 20
m <- sample(1:max.cluster.size,I,replace=TRUE); n <- sample(1:max.cluster.size,I,replace=TRUE)
Delta <- .3
cc <- .99
data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=cc,cov.xy=cc,cov.yy=cc,var.x=1,var.y=1,plot=FALSE)
x <- data$x; y <- data$y
auc.cluster(x,y,get.vcov=FALSE)
save.image('211231a.RData')



## large pop auc, small indiv auc

require(mvtnorm)
source('auc.cluster.R')
I <- 10
max.cluster.size <- k <- 10
m <- rpois(I,1)+1
m <- m + rbinom(I,1,1/2)*(20-2*m)
n <- 20-m
D.corr <- .99 # controls how u-shaped m,n are
Sigma <- matrix(D.corr,k,k) + diag(k)*(1-D.corr)
D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
D.idx <- D.latent>0
m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
n <- k-m
Delta <- 0
mm <- (n*(n>m)-m*(n<=m))/5
data <- rbinormal(I=I,m=m,n=n,mean.x=mm,mean.y=mm+Delta,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1/2,var.y=1/2,plot=FALSE)
x <- data$x; y <- data$y
auc.cluster(x,y,get.vcov=FALSE)
save.image('211231b.RData')


## plotting routine after 1a and 1b have been run

rm(list=ls())
for(save.file in c('211231a','211231b')) {
    load(paste0(save.file,'.RData'))
    pch.control <- '|'
    pch.case <- '-'
    png(paste0(save.file,'.png'))
    plot(0,type='n',xlim=range(unlist(c(x,y))),ylim=c(0,I),xlab='predictor',ylab='clusters',xaxt='n',yaxt='n',mgp=c(1,1,0))   
    for(i in 1:I) points(x[[i]],y=rep(i,length(x[[i]])),pch=pch.control)
    for(i in 1:I) points(y[[i]],y=rep(i,length(y[[i]])),pch=pch.case)
    abline(h=1/2)
    points(unlist(c(x,y)),y=rep(0,sum(m)+sum(n)),pch=rep(c(pch.control,pch.case),c(sum(m),sum(n))))
    dev.off()
}







## 2. coverage simulation
## uncomment code labeled #1 for the binormal and #2 for the censored binormal models

require(mvtnorm)
require(parallel)
source('auc.cluster.R')

bound <- abs(qnorm(.2)) # 2
auc.to.params <- function(theta.12,theta.11)auc.to.params.binormal.censored(theta.12,theta.11,bound) #2
sampler <- rbinormal.censored #2
## auc.to.params <- auc.to.params.binormal # binormal #1
## sampler <- rbinormal #1

B <- 1e3
I <- 60
alpha <- .05
q <- qchisq(1-alpha,df=2)
ks <- 2:5
theta.12s <- c(.7,.8)
theta.11s <-sapply(theta.12s, function(theta.12)c(theta.12,theta.12+.1))
aucs <- do.call(rbind,lapply(1:length(theta.12s),function(i)cbind(theta.12=theta.12s[i],theta.11=theta.11s[,i])))
params <- expand.grid(list(D.corr=c(0,.1,.4,.8)))
params <- as.data.frame(cbind(aucs[rep(1:nrow(aucs),each=nrow(params)),], params[rep(1:nrow(params),nrow(aucs)),,drop=FALSE]))
params <- round(params,1)
rownames(params) <- NULL
by.param <- apply(params,1,function(r) {
    ## browser()
    print(r)
    k <- sample(ks,1)
    D.corr <- unname(r['D.corr']);  theta.12 <- unname(r['theta.12']);  theta.11 <- unname(r['theta.11'])
    Delta.rho <- auc.to.params(theta.12,theta.11)
    rho <- round(Delta.rho['rho'],9)
    Delta <- Delta.rho['Delta']
    ## browser()
    Sigma <- matrix(D.corr,k,k) + diag(k)*(1-D.corr)
    sim <- lapply(1:B, FUN=function(jj){
    ## sim <- mclapply(1:B, mc.cores=detectCores()-3,FUN=function(jj){
        D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
        D.idx <- D.latent>0
        m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
        n <- k-m
        data <- sampler(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xy=rho)
        hat <- auc.cluster(data$x,data$y)
        cover <-  
            inside.ellipse(pt=c(theta.11,theta.12),  center=c(hat$theta.11.hat,hat$theta.12.hat), qform=hat$vcov.hat/I, radius=sqrt(q))
        cover.11 <- cover.12 <- cover.12a <- NA
        list(theta.11.hat=hat$theta.11.hat,theta.12.hat=hat$theta.12.hat,cover=cover,cover.11=cover.11,cover.12=cover.12,cover.12a=cover.12a,vcov.hat=hat$vcov.hat)
    })
    theta.11.hats <- sapply(sim,function(ll)ll$theta.11.hat)
    theta.12.hats <- sapply(sim,function(ll)ll$theta.12.hat)
    coverage <- mean(sapply(sim,function(ll)ll$cover))
    coverage.11 <- mean(sapply(sim,function(ll)ll$cover.11))
    coverage.12 <- mean(sapply(sim,function(ll)ll$cover.12))
    coverage.12a <- mean(sapply(sim,function(ll)ll$cover.12a))
    vcov.hats <- simplify2array(lapply(sim,function(ll)ll$vcov.hat))
    list(coverage=coverage,coverage.11=coverage.11,coverage.12=coverage.12,coverage.12a=coverage.12a,bias.theta.11=mean(theta.11.hats)-theta.11, bias.theta.12=mean(theta.12.hats)-theta.12,vcov.hat=apply(vcov.hats,1:2,mean),vcov.mc=matrix(c(var(theta.11.hats),rep(cov(theta.11.hats,theta.12.hats),2),var(theta.12.hats))*I,2))
})


require(xtable)
out <- sapply(by.param, function(lst)with(lst,c(coverage=coverage,bias.theta.12=bias.theta.12,bias.theta.11=bias.theta.11,structure((vcov.hat-vcov.mc)[c(1,2,4)],names=c('vcov.11','vcov.12','vcov.22')))))
out <- cbind(params,t(out))
rownames(out) <- NULL
colnames(out) <- c('$\\theta_{12}$','$\\theta_{11}$','$\\rho_{MN}$','coverage','$\\theta_{12}$','$\\theta_{11}$','$\\Sigma_{11}$','$\\Sigma_{12}$','$\\Sigma_{22}$')
addtorow <- list()
addtorow$pos <- list(0,0,0,0)
addtorow$command <- 
    c('\\multicolumn{3}{|c||}{parameters} & coverage &\\multicolumn{5}{c|}{bias} \\\\\n',
      '\\hline\n',
      '\\hline\n',
                      '$\\theta_{12}$ & $\\theta_{11}$ & $\\rho_{MN}$ &  & $\\theta_{12}$ & $\\theta_{11}$  & $\\Sigma_{11}$ & $\\Sigma_{12}$ & $\\Sigma_{22}$ \\\\\n')
col.wid <- 'p{1.4cm}'
c('c',rep(c('|',rep(col.wid,3),'|'),3))
align.str <- paste(c('c',rep('p{1.4cm}|',9)),collapse='')
align.str <-  "|c|C{1.4cm}|C{1.4cm}|C{1.4cm}||C{1.4cm}||C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|"
## save.file <- 'figs/220105.tex' #1
save.file <- 'figs/220811.tex' #2
## sink(save.file)
print.xtable(xtable(out, align=align.str),include.rownames=FALSE,include.colnames=FALSE,sanitize.text.function=function(x)x, add.to.row=addtorow)
## sink() 
lines <- scan(save.file,what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
writeLines(lines[start.idx:end.idx],save.file)


## record simulation settings to be \input into the ms
cat(I,'%',file='input/sim_coverage_I.txt')
cat(format(B,big.mark=','),'%',file='input/sim_coverage_reps.txt')
cat(paste(ks,collapse=', '),'%',file='input/sim_coverage_ks.txt')
for(col in colnames(params)) 
    cat(paste(unique(params[,col]),collapse=', '),'%',file=paste0('input/sim_coverage_',col,'s.txt'))
    



## 3. data analysis

## 3a. read in and clean up nyc and boston data
data.dir <- '../data/'
## nyc data
filelist <- dir(data.dir)
filelist <- filelist[grep('^sqf[-0-9]+\\.csv',filelist)]
filelist <- paste0(data.dir,filelist)
snf <- lapply(filelist[], function(file) {
    with(read.csv(file),
         data.frame(duration=STOP_DURATION_MINUTES,race=SUSPECT_RACE_DESCRIPTION,cluster=ISSUING_OFFICER_COMMAND_CODE))
})
snf <- do.call(rbind,snf)
snf <- within(snf, {
    duration[duration > 300] <- NA
    race[race=='MALE'] <- NA
    race[race=='(null)'] <- NA
    race[race %in% c('AMER IND','AMERICAN INDIAN/ALASKAN N','AMERICAN INDIAN/ALASKAN NATIVE','MIDDLE EASTERN/SOUTHWEST','MIDDLE EASTERN/SOUTHWEST ASIAN')] <- 'other'#'OTHER'
    race[race %in% c('ASIAN/PAC.ISL','ASIAN / PACIFIC ISLANDER')] <- 'asian'
    race[race=='BLACK'] <- 'black.nonhisp'
    race[race=='BLACK HISPANIC'] <- 'black.hisp'
    race[race=='WHITE'] <- 'white.nonhisp'
    race[race=='WHITE HISPANIC'] <- 'white.hisp'
})
nyc <- snf
## boston data
years <- 2019:2021 # stop_duration coded differently for year<=2018
snfs <- lapply(years, function(year) {
    ## if(year==2017)browser()
    snf <- read.csv(paste0(data.dir,'fio_',year,'.csv'))
    snf.names <- read.csv(paste0(data.dir,'fio_names_',year,'.csv'))
    snf <- merge(snf,snf.names,by='fc_num')
    snf <- snf[order(snf$fc_num),]
    single.encounters <- with(rle(snf$fc_num), values[lengths==1])
    snf <- snf[snf$fc_num %in% single.encounters,]
    ## subset(snf,select=c(fc_num,race,ethnicity,stop_duration,contact_officer))
    with(snf, data.frame(duration=stop_duration,race.old=race,ethnicity=ethnicity,cluster=contact_officer))
})
snf <- do.call(rbind,snfs)
snf$duration[snf$duration=='NULL'] <- NA
snf$duration <- as.numeric(snf$duration)
snf$duration[snf$duration>300] <- NA
snf$race <- snf$race.old
snf <- within(snf, {
    race[race==''] <- NA
    race[race %in% c('Unknown','NULL')] <- NA
    race[race %in% c("American Indian or Alaskan Native","Native Hawaiian or Other Pacific Islander","Other")] <- 'other'
    race[race=='Asian'] <- 'asian'
    race[race=='Black' & ethnicity=='Hispanic Origin'] <- 'black.hisp'
    race[race=='Black' & ethnicity!='Hispanic Origin'] <- 'black.nonhisp'
    race[race=='White' & ethnicity=='Hispanic Origin'] <- 'white.hisp'
    race[race=='White' & ethnicity!='Hispanic Origin'] <- 'white.nonhisp'
})
boston <- snf



## 3b aggregate tables
summary.tables <- lapply(list(nyc=nyc,boston=boston), function(snf) {
    mean.durations <- aggregate(duration ~ race, mean, data=snf)
    sd.durations <- aggregate(duration ~ race, sd, data=snf)
    colnames(mean.durations)  <- c('group','mean.duration')
    colnames(sd.durations) <- c('group','sd.duration')
    counts <- table(snf$race)
    counts <- data.frame(group=names(counts),count=unname(as.integer(counts)))
    out <- merge(merge(mean.durations,sd.durations),counts)
    out$freq <- out$count / sum(out$count)
    out
})




## 3c AUC analysis
source('../auc.cluster.R')
statuses <- c(quote(race=='black.nonhisp' | race=='black.hisp'),
              quote(race=='black.nonhisp'),
              quote(race=='black.hisp'),
              quote(race=='white.nonhisp' | race=='white.hisp'),
              quote(race=='white.nonhisp'),
              quote(race=='white.hisp'),
              quote(race=='white.hisp' | race=='black.hisp'))
status.descriptions <- c('Black','Black non-Hispanic','Black Hispanic','White','White non-Hispanic','White Hispanic','Hispanic')
auc.estimates <- Map(function(snf,city.name) {
    sapply(statuses, function(status) {
        snf$status <- with(snf, eval(status))
        clusters <- split(snf,snf$cluster)
        x <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=status==FALSE,select=duration)))))
        y <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=status==TRUE,select=duration)))))
        omit.idx <- sapply(x,length)==0 | sapply(y,length)==0
        x <- x[!omit.idx]; y <- y[!omit.idx]
        auc0 <- auc.cluster(x=x,y=y) 
        contrast <- matrix(c(1,-1),nrow=2)
        alpha <- .05
        I <- length(x)
        z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
        if(paste(status,collapse='')==paste(statuses[[1]],collapse='')) {
            pdf(paste0('figs/220823_',city.name,'.pdf'))
            plot(auc0,resol=1e2)
            plot(auc0,alpha=.01,add=TRUE)
            dev.off()
        }
        with(auc0,structure(c(theta.12.hat,theta.12.CI,theta.11.hat,theta.11.CI,(1-pnorm(abs(z.stat)))*2,I,length(unlist(x)),length(unlist(y))),
                            names=c('theta.12.hat','theta.12.CI.lower','theta.12.CI.upper','theta.11.hat','theta.11.CI.lower','theta.11.CI.upper','pval','I','M','N')))
    })
}, list(nyc=nyc,boston=boston), c('nyc','boston'))
auc.estimates <- lapply(auc.estimates,t)
## auc.estimates <- lapply(auc.estimates,function(m)cbind(case=status.descriptions,as.data.frame(m)))
auc.estimates <- lapply(auc.estimates,function(m) {
    df <- cbind(case=status.descriptions,as.data.frame(m))
    rounded <- lapply(df[,grep('theta',colnames(df))],function(col)sprintf(col,fmt='%#.2f'))
    df$theta.12.str <- with(rounded, paste0(theta.12.hat,' [',theta.12.CI.lower,', ',theta.12.CI.upper,']'))
    df$theta.11.str <- with(rounded, paste0(theta.11.hat,' [',theta.11.CI.lower,', ',theta.11.CI.upper,']'))
    df$pval.str <- sprintf(df$pval,fmt='%#.2f')
    df
})


## 3d pretty print summaary data tables
require(xtable)
out <- lapply(summary.tables, function(x) {
    x <- x[order(x$group),]
    for(col in c('mean.duration','sd.duration','freq'))
        x[,col] <- sprintf(x[,col],fmt='%#.2f')
    x$mean.duration <- paste0(x$mean.duration,' (',x$sd.duration,')')
    subset(x,select=-sd.duration)
    })
out <- cbind(nyc=out[['nyc']],boston=subset(out[['boston']],select=-group))
colnames(out) <- c('group','mean duration (SD)','count','freq.','mean duration (SD)','count','freq.')
out$group <- c('Asian','Black Hispanic','Black non-Hispanic','other','White Hispanic','White non-Hispanic')
other.idx <- which(out$group=='other')
out <- out[c(1:(other.idx-1),(other.idx+1):nrow(out),other.idx),]
addtorow <- list()
addtorow$pos <- list(0,0,0,0)
addtorow$command <- c("& \\multicolumn{3}{c||}{NYC} & \\multicolumn{3}{c|}{Boston}\\\\\n",
                      '\\hline\n',
                      '\\hline\n',
                      paste0(paste0(colnames(out),collapse=' & '), '\\\\\n'))
align.str <-  "|c|l||C{3cm}|C{1cm}|C{1cm}||C{3cm}|C{1cm}|C{1cm}|"#|C{2cm}|C{2cm}|"
filename <- '220823.tex'
sink(filename)
print.xtable(xtable(out,align=align.str,digits=rep(0,ncol(out)+1)), add.to.row = addtorow,include.colnames=FALSE,include.rownames=FALSE,sanitize.text.function=function(x)x)
sink()
lines <- scan(filename,what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
writeLines(lines[start.idx:end.idx],filename)


## 3f record stats to be  \input to manuscript
for(city in c('nyc','boston')) {
    cat(format(nrow(get(city)),big.mark=','),'%',file=paste0('input/da_total_stops_',city,'.txt'))
    apply(auc.estimates[[city]],1,function(row) {
        cat(row['theta.12.str'],'%',sep='',file=tolower(paste0('input/da_',row['case'],'_',city,'_theta12.txt')))
        cat(row['theta.11.str'],'%',sep='',file=tolower(paste0('input/da_',row['case'],'_',city,'_theta11.txt')))
        cat(row['pval.str'],'%',sep='',file=tolower(paste0('input/da_',row['case'],'_',city,'_pval.txt')))
    })
}
with(auc.estimates$nyc, {
    fmt <- '%#.2f'
    cat(sprintf(fmt,theta.12.hat[case=='Black']),'%',file='input/da_black_nyc_theta12_mean.txt')
    cat(sprintf(fmt,theta.12.CI.lower[case=='Black']),'---',sprintf(fmt,theta.12.CI.upper[case=='Black']),'%',sep='',file='input/da_black_nyc_theta12_CI.txt')
    cat(sprintf(fmt,theta.11.hat[case=='Black']),'%',file='input/da_black_nyc_theta11_mean.txt')
    cat(sprintf(fmt,theta.11.CI.lower[case=='Black']),'---',sprintf(fmt,theta.11.CI.upper[case=='Black']),'%',sep='',file='input/da_black_nyc_theta11_CI.txt')
})



## 3e pretty print table of AUC estimate
require(xtable)
idx <- as.numeric(t(matrix(1:(2*nrow(auc.estimates[[1]])),ncol=2)))
out <- as.data.frame(do.call(rbind,auc.estimates)[idx,])
out$city <- c('NYC','Boston')
out$case.status <- rep(status.descriptions,each=2)
out <- subset(out,select=c(case.status,city,I,M,N,theta.12.str,theta.11.str,pval.str))
out <- within(out,case.status[duplicated(case.status)] <- '')
colnames(out) <- c('case group','data set','I','$\\Sigma M_i$','$\\Sigma N_i$','$\\theta_{12}$','$\\theta_{11}$','$H_{0}:\\theta_{12}=\\theta_{11}$')
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- c(    '\\hline\n')
align.str <-  "|c|l||c|c|c|c|c|c|c|"#C{3cm}|C{1cm}|C{1cm}||C{3cm}|C{1cm}|C{1cm}|"
## sink('220821.tex')
print.xtable(xtable(out,align=align.str,digits=rep(0,ncol(out)+1)),include.rownames=FALSE,add.to.row=addtorow,sanitize.text.function=function(x)x)
sink()
lines <- scan('220821.tex',what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
## writeLines(lines[start.idx:end.idx],'220821.tex')


