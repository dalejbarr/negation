nmc <- 9999L
library(dplyr)
library(tidyr)

con <- src_sqlite("raw_data/Experiment_Data.db")

tdata <- tbl(con, "TrialData") %>% collect(n = Inf) %>%
    select(RespID, SubjID, Cond)

pog <- tbl(con, "POG_exproffset") %>% collect(n = Inf) %>%
    filter(ms>=0, ms<=2500) %>%
    mutate(bin=floor((ms+25)/50)*50)

# the following lines divide the data into bins
# aggregated up to the Subject level
dat <- tdata %>% inner_join(pog, by="RespID") %>%
    group_by(SubjID, Cond, bin, ID) %>%
    summarize(Y=n()) %>%
        right_join(expand(tdata %>%
                              inner_join(pog, by="RespID"),
                          SubjID, Cond, bin, ID))
dat$Y <- with(dat, ifelse(is.na(Y), 0, Y))

dat_subj <- dat %>% group_by(SubjID, Cond, bin) %>%
    summarize(N=sum(Y)) %>% ungroup() %>% inner_join(dat) %>%
    mutate(p=ifelse(N==0, NA, Y/N))
getClust <- function(x) {
    ff <- x %>% mutate(cl=(p<.05)*sign(F))
    ff.runs <- rle(ff$cl)
    nruns <- length(ff.runs$lengths)
    clust.ix <- which(ff.runs$values!=0)
    if (length(clust.ix)) {
        res <- lapply(clust.ix, function(ix) {
            if (ix>1) {
                t0 <- sum(ff.runs$length[1:(ix-1)]) + 1
            } else {
                t0 <- 1
            }
            t1 <- t0 + ff.runs$lengths[ix] - 1
            csum <- sum(ff$F[t0:t1])
            data.frame(t0=t0, t1=t1, csum=csum)
        })
        res <- do.call("rbind", res)
    } else { # do something for zero case
        res <- data.frame(t0=NA, t1=NA, csum=0)
    }
    res
}

# generic t-test
ttest1 <- function(x) {
    lvls <- unique(x$Cond2)
    ff <- x %>% 
        spread(Cond2, eff)
    vec <- ff[[lvls[1]]] - ff[[lvls[2]]]
    vmean <- mean(vec)
    serr <- sd(vec)/sqrt(length(vec))
    tobs <- vmean / serr
    data_frame(tobs=tobs, pval=2*(1-pt(abs(tobs), 23)))    
}

flip1 <- function(x) {
    xx <- mutate(x, CondNew=Cond2)
    if (sample(c(TRUE, FALSE), 1)) {
        xx$CondNew <- rev(xx$CondNew)
    } else {}
    return(xx)
}

permuteOnce1 <- function(x, unitcond) {
    unitcond %>% group_by(UnitID) %>%
        do(flip1(.)) %>% ungroup() %>%
        inner_join(x, by=c("UnitID", "Cond2")) %>%
        select(-Cond2) %>% rename(Cond2=CondNew)
}

# generic: do one permutation run
do_once <- function(x, subjcond) {
    clust <- permuteOnce1(x, subjcond) %>%
        group_by(bin) %>% do(ttest1(.)) %>%
        rename(F=tobs, p=pval) %>%
        getClust()
    return(abs(max(clust$csum)))
}

# get pvalues for each cluster
get_clustp <- function(orig, nhd) {
    pval <- function(x) {
        data_frame(p=sum(abs(c(x$csum, nhd))>=abs(x$csum))/(length(nhd)+1))
    }
    orig$p <- orig %>% rowwise() %>% do(pval(.)) %>% `[[`("p")
    return(orig)
}

# compare unmentioned to critical
# do it by subject
dat_uvc <- dat_subj %>%
    group_by(SubjID, Cond, bin, ID) %>%
    summarize(Y=sum(Y)) %>% ungroup() %>%
    spread(ID, Y) %>%
    mutate(N=blnk+comp+negt+unme, comp_p=comp/N, unme_p=unme/N,
            eff=log((unme_p + .5) / (comp_p + .5))) %>%
    select(-blnk, -comp, -negt, -unme, -N)

dat_uvc2 <- dat_uvc %>%
    mutate(Cond2=ifelse(Cond=="NG", "NG", "GD")) %>%
    group_by(UnitID=SubjID, Cond2, bin) %>%
    summarize(eff=mean(eff)) %>% ungroup() %>%
    arrange(UnitID, bin, Cond2)

subjcond <- dat_uvc2 %>%
    select(UnitID, Cond2) %>%
    distinct()

clust_orig <- dat_uvc2 %>% group_by(bin) %>%
    do(ttest1(.)) %>% rename(F=tobs, p=pval) %>%
    getClust()

nhd <- replicate(nmc, do_once(dat_uvc2, subjcond))

clust_orig <- get_clustp(clust_orig, nhd)
print(clust_orig)

saveRDS(clust_orig, file="result/cluster_random_grounded_vs_ungrounded_by_subject.rds")
