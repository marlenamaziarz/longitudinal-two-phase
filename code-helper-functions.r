
main <- function(dc, d.full){
    set.seed(dc$seed)
    
    if(dc$ncc){ # NCC
        weights          <- get.weights.ncc(dc, d.full)      # same length as the number of obs in ncc sample
        weights.pert.mx  <- get.weights.pert.ncc(dc, d.full) # for var est 
        data.sample      <- d.full[d.full$in.ncc, ]
    }else{ # CCH
        weights          <- get.weights.cch(dc, d.full)      # same length as the number of obs in ncc sample
        weights.pert.mx  <- get.weights.pert.cch(dc, d.full) # for var est 
        data.sample      <- d.full[d.full$in.cch, ]
    }
    train.data.list <- list(data = data.frame(data.sample), 
                            w.pert.mx = weights.pert.mx,  # perturbed, matrix with P columns 
                            weights = weights)  # unperturbed
    
    test.data.list   <- train.data.list # training and testing is on the same set.
    test.data.list.s <- get.info.at.time.s(dc, data.list = test.data.list) # in this case LVCF, can do better by using BLUP etc. 
    
    risk.vector <- get.logistic.risk(dc, 
                                     train.data.list$data, 
                                     test.data.list.s$data.s, 
                                     w.fit = train.data.list$weights) # for est, null for cohort
    
    risk.mx <- matrix(NA, ncol = dc$P, nrow = dim(test.data.list.s$data.s[1]))
    for(i in 1:dc$P){
        risk.mx[,i] <- get.logistic.risk(dc, 
                                         train.data.list$data, 
                                         test.data.list.s$data.s, 
                                         w.fit = train.data.list$w.pert.mx[,i]) 
    }
    
    out.list      <- get.stats(dc = dc, 
                               risk = risk.vector, 
                               data.s = test.data.list.s$data.s, 
                               weights.s = test.data.list.s$weights.s)
    
    ix.col.na <- is.na(risk.mx[1,])
    sim.stats.var <- get.stats.var(dc,
                                   risk.mx = risk.mx[, !ix.col.na], 
                                   data.s = test.data.list.s$data.s, 
                                   weights.mx.s = test.data.list.s$w.pert.mx.s[, !ix.col.na])
    
    return(list(list(success = T, 
                     sim.stats = out.list$stats, 
                     sim.stats.se = sqrt(sim.stats.var),
                     n.case = out.list$n.case,
                     n.ctrl = out.list$n.ctrl,
                     n.cens = out.list$n.cens)))
}    






set.data.control <- function(si = NA, 
                             ti = NA,
                             P = 500,
                             ncc = NA, # has to be specified. T = NCC, F = CCH
                             stats = c( 'PE',
                                        'TPF(0.4)', 
                                        'FPF(0.3)',
                                        'AUC', 
                                        'PCF(0.2)', 
                                        'PNF(0.8)'),
                             results.digits = 3, 
                             seed = 1,
                             spline.s.df = 3,
                             meas.time.interval = 3){
    
    
    # control object integrity check
    # exit if error encountered
    if(is.na(si) | is.na(ti)){
        print('Error in the sim control: provide si and ti')
        return(NULL)
    }
    
    data.control <- list(si = si, 
                         ti = ti, 
                         P = P,
                         ncc = ncc,
                         pred.time = ti - si,
                         stats = stats, 
                         n.stats = length(stats),
                         results.digits = results.digits,
                         seed = seed, 
                         spline.s.df = spline.s.df,
                         meas.time.interval = meas.time.interval)
    return(data.control)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN SOURCES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run.sources <- function(path = NULL){
    
    source(paste(path, 'data-sim.r', sep = ''))
    source(paste(path, 'stats.r', sep = ''))
    source(paste(path, 'fitting.r', sep = ''))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FITTING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# given a full dataset, return data at time s (s defined in dc)
# only works if measurements are available at regular intervals, easy to modify to allow for meas. time variability
# this is implicitly conditional on survival up to time s,
# otherwise there would not be an observation at s.
get.info.at.time.s <- function(dc, data.list){
    
    # get the last observation for each subject available before time s, 
    data.list.up.to.s <- get.info.up.to.time.s(dc, data.list)
    ix.last           <- c(diff(data.list.up.to.s$data$sub.id) != 0, T)

    data.s            <- data.list.up.to.s$data.s[ix.last, ]
    data.s$meas.time  <- dc$si 
    
    w.pert.mx.s <- data.list.up.to.s$w.pert.mx.s[ix.last, ]
    
    weights.s <- data.list.up.to.s$weights.s[ix.last]
   
    rm(data.list.up.to.s, ix.last)
    return(list(data.s = data.s, w.pert.mx.s = w.pert.mx.s, weights.s = weights.s))
}






# given a full dataset, return data at time s (s defined in dc)
# only works if measurements are available at regular intervals, and thus also precisely at si.
# this is implicitly conditional on survival up to time s,
# otherwise there would not be an observation at s.
get.info.up.to.time.s <- function(dc, data.list){
    
    ix.event.before.s <- data.list$data$time <= dc$si
    ix.obs.before.s   <- data.list$data$meas.time <= dc$si 
    ix.up.to.s        <- !ix.event.before.s & ix.obs.before.s
    
    data.s      <- data.list$data[ix.up.to.s, ]
    w.pert.mx.s <- data.list$w.pert.mx[ix.up.to.s, ]
    
    weights.s <- data.list$weights[ix.up.to.s] 
    
    return(list(data.s = data.s, w.pert.mx.s = w.pert.mx.s, weights.s = weights.s))
}



# weights.s should be a vector of perturbation weights the same length as data.s
get.d.pred <- function(dc, risk, data.s, weights.s = NULL){
    

    ix.case <- (data.s$time <= dc$ti) & (data.s$status == 1)
    ix.ctrl <- (data.s$time > dc$ti)
    ix.cens <- (data.s$time <= dc$ti) & (data.s$status == 0)        
    risk    <- risk[!ix.cens]
    ix.case <- ix.case[!ix.cens]
    ix.ctrl <- ix.ctrl[!ix.cens] 


    # get censoring weights just for the non-censored individuals, sorted as in data.s        
    new.times.nc <- data.s$time[!ix.cens]
    w.cens <- get.censoring.weights(dc, data = data.s, weights = weights.s, 
                                    new.times = new.times.nc, 
                                    ix.ctrl = ix.ctrl, 
                                    time.scale = 'time')
    if(is.null(weights.s)){
        weights.s <- w.cens
    }else{
        weights.s <- weights.s[!ix.cens] 
        weights.s <- w.cens * weights.s
    }

    
    # sort everything according to decreasing risk.
    ix.sort     <- sort(risk, decreasing = T, index.return = T)$ix
    risk        <- risk[ix.sort]
    ix.case     <- ix.case[ix.sort]
    ix.ctrl     <- ix.ctrl[ix.sort]
    weights.s   <- weights.s[ix.sort]
    
    d.pred <- list(risk         = risk, 
                   risk.case    = risk[ix.case],
                   risk.ctrl    = risk[ix.ctrl],
                   ix.case      = ix.case,
                   ix.ctrl      = ix.ctrl,
                   weights      = weights.s,
                   weights.case = weights.s[ix.case],
                   weights.ctrl = weights.s[ix.ctrl])
    return(d.pred)
}


# data      = dataset on which survfit will be fit, data at time s.
# weights   = weights for each individual in data
# new.times = times at which the censoring distribution is to be estimated
# ix.ctrl   = index of controls, same length as new.times, ie. index with respect to new.times
# time.scale = legacy, we used to have the option of 't.star'.
get.censoring.weights <- function(dc, data, weights = NULL, new.times, 
                                  ix.ctrl, time.scale = 'time'){
    
    # this saves some time. the survival curve doesn't need to be estimated past ti.
    ix            <- data$time > (dc$ti + 0.1) 
    data$time[ix] <- (dc$ti + 0.1)
    
    cc  <- survfit(Surv(time, status == 0) ~ 1,
                   data = data, weights = weights, se.fit = F, type = 'kaplan-meier')
    
    
    new.times.case.ctrl          <- new.times
    new.times.case.ctrl[ix.ctrl] <- dc$ti

    new.times.case.ctrl.sorted.incr <- sort(new.times.case.ctrl, decreasing = F, method = 'shell')
    recover.original.order <- rank(new.times.case.ctrl)
    
    cens.weights <- summary(cc, times = new.times.case.ctrl.sorted.incr)$surv[recover.original.order]
    
    return(1/cens.weights)
}





# ix.ncc - index of NCC, length = n obs = N (so full set)
# data = full data, N (n obs) x n covariates 
#
# return NCC weights for those in the NCC subset
get.weights.ncc <- function(dc, d.full){
    
    d.base <- d.full[d.full$ix.first, ] # first visits
    
    m <- gam(in.ncc ~ s(time) + status, data = d.base, family = binomial)
    
    # est.weights for everyone (1 obs. per person)
    d.base$est.weights <- 1/fitted(m)
    
    # now, put them back in the logitudinal ncc sample.
    base.ncc.sub.id  <- d.base$sub.id[d.base$in.ncc]
    weights.ncc.base <- d.base$est.weights[d.base$in.ncc]
    long.ncc.sub.id  <- d.full$sub.id[d.full$in.ncc]
    
    # cleanup
    rm(d.base, m)
    
    start.row.next.person <- pmatch(base.ncc.sub.id, long.ncc.sub.id)
    # an inelegant hack, but makes the loop below cleaner.
    start.row.next.person <- c(start.row.next.person, (length(long.ncc.sub.id) + 1))
    
    weights.ncc.long <- vector(mode = 'numeric', length = length(long.ncc.sub.id))
    sub.id.unique    <- unique(sort(long.ncc.sub.id))
    
    for(i in 1:length(sub.id.unique)){
        weights.ncc.long[start.row.next.person[i]:(start.row.next.person[i+1]-1)] <- weights.ncc.base[i]
    }
    
    return(weights.ncc.long) # longitudinal sample for those in NCC
}








# ix.ncc - index of NCC, length = n obs = N (so full set)
# w.pert.mx - matrix of pert weights (N x P)
# data = full data, N (n obs) x n covariates 
#
# returns NCC pert weights for those in the NCC subset
get.weights.pert.ncc <- function(dc, d.full){
    
    ix.base <- d.full$ix.first
    d.base  <- d.full[ix.base, ] 
    # generate weights, P per person, 
    np <- sum(ix.base)*dc$P
    w.pert.mx.base <- matrix(rexp(np), ncol = dc$P)
    
    est.weights.pert.base.mx <- matrix(NA, nrow = dim(d.base)[1], ncol = dc$P)
    for(i in 1:dc$P){
        m <- gam(in.ncc ~ s(time) + status, 
                 data = d.base, family = binomial, weights = w.pert.mx.base[,i])
        
        # est.weights for everyone (1 obs. per person in the full set)
        est.weights.pert.base.mx[,i] <- w.pert.mx.base[,i]/fitted(m)
    }
    
    # now, put them back in the logitudinal ncc sample.
    base.ncc.sub.id  <- d.base$sub.id[d.base$in.ncc]
    est.weights.pert.base.ncc.mx <- est.weights.pert.base.mx[d.base$in.ncc, ]
    long.ncc.sub.id  <- d.full$sub.id[d.full$in.ncc]
    
    # cleanup
    rm(d.base, d.full, ix.base.ncc, ix.base, m, est.weights.pert.base.mx)
    
    start.row.next.person <- pmatch(base.ncc.sub.id, long.ncc.sub.id)
    # an inelegant hack, but makes the loop below cleaner.
    start.row.next.person <- c(start.row.next.person, (length(long.ncc.sub.id) + 1))
    
    w.pert.ncc.long.mx <- matrix(NA, nrow = length(long.ncc.sub.id), ncol = dc$P)
    sub.id.unique    <- unique(sort(long.ncc.sub.id))
    
    for(i in 1:length(sub.id.unique)){
        # check that this replacement does the right thing... checked, it works as expected.
        n.obs <- start.row.next.person[i+1]-start.row.next.person[i]
        w.pert.ncc.long.mx[start.row.next.person[i]:(start.row.next.person[i+1]-1), ] <- rep(est.weights.pert.base.ncc.mx[i, ], each = n.obs)
    }
    
    return(w.pert.ncc.long.mx)
    
}






# data should have 'in.cch'
# returns cch weights for those in the CCH
get.weights.cch <- function(dc, d.full){

    # unique subject ids
    sub.ids.1 <- d.full$sub.id[d.full$status == 1 & d.full$visit == 1]
    sub.ids.0 <- d.full$sub.id[d.full$status == 0 & d.full$visit == 1] 
    
    sub.ids.1.cch <- with(d.full, sub.id[status == 1 & visit == 1 & in.cch])
    sub.ids.0.cch <- with(d.full, sub.id[status == 0 & visit == 1 & in.cch]) 
    
    # these are inverse weights
    weight.1 <- length(sub.ids.1)/length(sub.ids.1.cch) # this should really be 1
    weight.0 <- length(sub.ids.0)/length(sub.ids.0.cch) # should be > 1
    
    weights.cch <- d.full$status * weight.1 + (1 - d.full$status) * weight.0
    weights.cch <- weights.cch[d.full$in.cch]
    return(weights.cch)
}





# d.full should have 'in.cch' and 'ix.first'
# output: matrix of perturbation weights that account for censoring (size n_cch (long) x P)
get.weights.pert.cch <- function(dc, d.full){
    
    # w.pert is a matrix of perturbation weights (long)
    sub.id.unique <- unique(d.full$sub.id)
    np <- length(sub.id.unique) * dc$P 
    w.pert.base <- matrix(rexp(np), ncol = dc$P) # pert weights are exp(1)
    
    # now, put weights in long format. long obs for a given individual have the same weigth for a given perturbation
    w.pert <- matrix(NA, nrow = dim(d.full)[1], ncol = dc$P)
    
    start.row.next.person <- with(d.full, c(pmatch(sub.id[ix.first], sub.id), length(sub.id) + 1)) # the last number makes the loop below easier
    
    for(i in 1:length(sub.id.unique)){
        n.obs <- start.row.next.person[i+1]-start.row.next.person[i]
        w.pert[start.row.next.person[i]:(start.row.next.person[i+1]-1), ] <- rep(w.pert.base[i, ], each = n.obs) # matrix n.obs x dc$P
    }
    
    # calc CCH weights
    # sum up weights in events and non-events, one weight per individual
    weights.1.num <- colSums(w.pert[d.full$status == 1 & d.full$visit == 1, ])
    weights.0.num <- colSums(w.pert[d.full$status == 0 & d.full$visit == 1, ])
    
    weights.1.den <- colSums(w.pert[d.full$in.cch & d.full$status == 1 & d.full$visit == 1, ])
    weights.0.den <- colSums(w.pert[d.full$in.cch & d.full$status == 0 & d.full$visit == 1, ])
    
    # inverse perturbed cch weights
    weights.1 <- weights.1.num/weights.1.den 
    weights.0 <- weights.0.num/weights.0.den
    
    w.pert.cch.naive <- w.pert[d.full$in.cch, ]         # individual pert weights in CCH
    status.cch       <- d.full$status[d.full$in.cch]    # event status in CCH
    
    # here we rescale the perturbation weights to account for CCH sampling.
    w.pert.cch.rescale <- matrix(NA, nrow = length(status.cch), ncol = dc$P)
    w.pert.cch.rescale[status.cch == 1, ] <- weights.1
    w.pert.cch.rescale[status.cch == 0, ] <- weights.0
    w.pert.cch <- w.pert.cch.naive * w.pert.cch.rescale
    return(w.pert.cch) 
}







# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# est risk
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
get.logistic.risk <- function(dc, training.data, data.s, w.fit = NULL){ 
    
    data.tmp <- data.frame(sub.id = training.data$sub.id, 
                           si     = training.data$meas.time, 
                           zi     = training.data$marker, 
                           xi     = training.data$time, 
                           di     = training.data$status) 
    
    glm.data <- preGLM.FUN(dc, data.tmp, marker.ix = 3) 
    
    glm.out  <- GLMt.FUN(dc             = dc, 
                         data           = glm.data,
                         pred.time.vec  = dc$si, 
                         tau            = dc$pred.time, 
                         y.vec          = data.s$marker,
                         w.fit          = w.fit,
                         training.data  = training.data)
    
    if(is.na(glm.out)){
        return(NA)
    }
    risk <- glm.out$predp.tau.s0.y0

    rm(data.tmp, glm.data, glm.out)
    return(risk)
}



#### the next three functions are used to estimate Pr(T<s+tau|T>=s,Z(s))
#### code is based on T. Cai's Kernel based IPW-GLM for longitudinal data
#### we use measurement as covariate in the model instead of kernel function 
preGLM.FUN <- function(dc, data, marker.ix, n.digits = 6){
    tau <- dc$pred.time
    ## unpack data 
    xi <- round(data$xi, digits = n.digits)   # surv time, ie. time
    si <- round(data$si, digits = n.digits)   # meas.time
    zi <- data[ , marker.ix, drop = FALSE]    # marker
    yi <- as.numeric(xi < si + tau)           # ie. I(time < meas.time + tau) = t* < tau
    
    # remove those who are censored and have been lost to followup before meas.time + tau
    ix.censored <- (data$di == 0) & (data$xi < si + tau) 
    
    working.dataset <- data.frame(sub.id = data$sub.id, 
                                  yi = yi, 
                                  si = si, 
                                  xi = xi, 
                                  di = data$di, 
                                  zi = zi)
    working.dataset.nc <- working.dataset[ix.censored != 1, ]
    
    rm(working.dataset, tau, xi, si, zi, yi, ix.censored)
    
    return(list(working.dataset.nc = working.dataset.nc, n.digits = n.digits))
}


GLMt.FUN <- function(dc, data, pred.time.vec, tau, y.vec, w.fit = NULL, training.data){ 
    
    # cenosored subjects removed.
    yi       <- data$working.dataset.nc$yi                   # status of event before s + t
    # xi     <- data$working.dataset.nc$xi                   # time
    # di     <- data$working.dataset.nc$di                   # status
    # si     <- data$working.dataset.nc$si                   # meas.time
    # zi     <- data$working.dataset.nc[,-(1:5),drop=FALSE]  # marker
    
    n <- length(yi)   
    bs.mx  <- ns(data$working.dataset.nc$si, df = dc$spline.s.df)
    colnames(bs.mx) = paste(rep('s', dc$spline.s.df), 1:dc$spline.s.df, sep = '')
    
    # matrix of covariates (spline function of measurement time, marker)
    bs.zi <- cbind(bs.mx, data$working.dataset.nc[,-(1:5),drop=FALSE])
    
    # design matrix with intercept, spline function of measurement time, and marker 
    ones <- rep.int(1, n)
    ui <- as.matrix(cbind(ones, bs.zi)) # rep(1, dim(zi)[1]), zi))  
    model.data <- as.data.frame(cbind(yi, ui))
    rm(yi)
    
    new.s.tau <- round(data$working.dataset.nc$si + tau, digits = data$n.digits)
    
    wgt.IPW <- IPW.FUN(est.time   = training.data$time[training.data$ix.first], 
                       est.status = training.data$status[training.data$ix.first], 
                       new.xi     = data$working.dataset.nc$xi,
                       new.s.tau  = new.s.tau, 
                       new.status = data$working.dataset.nc$di, 
                       weights    = w.fit[training.data$ix.first])
    
    
    if(!is.null(w.fit)){             # w.fit will be null for estimation in the cohort.
        ix.cens  <- (training.data$status == 0) & (training.data$time < training.data$meas.time + dc$pred.time)
        w.fit.nc <- w.fit[!ix.cens]
        weights <- wgt.IPW * w.fit.nc
    }else{
        weights <- wgt.IPW 
    }
    fmla <- as.formula(paste("yi ~ ", paste(colnames(bs.zi), collapse= "+")))
    fit  <- try(glm(fmla, weights = weights, data = model.data, family = "binomial"))
    if(class(fit) == 'try-error'){
        print('error in glm fit')
        return(NA)
    }
    
    beta            <- matrix(fit$coef, ncol = 1)
    predp.tau.s0.y0 <- NULL
    
    if(!is.null(pred.time.vec)){
        bs0 <- predict(bs.mx, pred.time.vec)
        ns0 <- length(pred.time.vec)
        ny0 <- length(y.vec)
        
        ui0 <- cbind(rep.int(1, ns0 * ny0),
                     matrix(rep(bs0, ny0), nrow = ny0*ns0, byrow = T),
                     rep(y.vec, each = ns0))
        
        predp.tau.s0.y0 <- g.logit(ui0 %*% beta) 
        rm(bs0, ns0, ny0, ui0)
    }
    return(list(beta = beta, predp.tau.s0.y0 = predp.tau.s0.y0))
}


# used as: IPW.FUN(xi.baseline, di.baseline, xi, round(si+tau, digits = data$n.digits), di)
IPW.FUN <- function(est.time, est.status, new.xi, new.s.tau, new.status, type = "kaplan", weights = NULL){
    numerator   <- new.status * (new.xi <= new.s.tau) + (new.xi >= new.s.tau) 
    min.time    <- pmin(new.xi, new.s.tau)
    ranked.min.time <- rank(min.time)
    denominator <- summary(survfit(Surv(est.time, est.status == 0) ~ 1, 
                                   se.fit = F, type = type, weights = weights), 
                           times = sort(min.time))$surv[ranked.min.time]
    result  <- numerator/denominator
    rm(numerator, min.time, ranked.min.time, denominator)
    return(result)
}


# in:
# vector.to.rep = vector to be replicated by row
# n.row = number of rows in the returned matrix made up of the vector vector.to.rep
# out:
# matrix n.rows x length(vector.to.rep)
VTM <- function(vector.to.rep, n.rows){
    matrix(vector.to.rep, ncol = length(vector.to.rep), nrow = n.rows, byrow = T)
}


g.logit <- function(x){
    return(1/(1 + exp(-x)))
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~


get.stats <- function(dc, risk, data.s, weights.s = NULL){
    
    # put all the info we need for getting stats into one data frame
    d.pred <- get.d.pred(dc, risk, data.s, weights.s = weights.s)
    
    n.case <- sum(d.pred$ix.case)
    n.ctrl <- sum(d.pred$ix.ctrl)
    n.cens <- dim(data.s)[1] - (n.case + n.ctrl)
    
    stats <- matrix(c(get.prediction.error(d.pred),
                      get.tpfp(d.pred, risk.thr = c(0.4, 0.3)),  
                      get.auc(d.pred),
                      get.pcf(d.pred, prop = 0.2),  
                      get.ppnf(d.pred, prop = 0.8)), 
                    ncol = 1)
    rownames(stats) <- dc$stats
    return(list(stats = stats, n.case = n.case, n.ctrl = n.ctrl, n.cens = n.cens))
}


# risk = vector of risks
# weights.mx.s = matrix of perturbation weights
get.stats.var <- function(dc, risk.mx, data.s, weights.mx.s){
    
    stats.mx <- matrix(NA, nrow = dc$n.stats, ncol = dc$P)
    
    for(i in 1:dc$P){
        
        # put all the info we need for getting stats into one data frame
        d.pred <- get.d.pred(dc, risk.mx[,i], data.s, weights.s = weights.mx.s[,i])
        
        stats.mx[,i] <- c(get.prediction.error(d.pred),
                          get.tpfp(d.pred, risk.thr = c(0.4, 0.3)), 
                          get.auc(d.pred),
                          get.pcf(d.pred, prop = 0.2), 
                          get.ppnf(d.pred, prop = 0.8)) 
    }
    stats.var <- apply(stats.mx, 1, var, na.rm = T)
    return(stats.var)
}


# ~~~~~~~~~~~~~~~~
# ~~~ get.tpfp ~~~
# ~~~~~~~~~~~~~~~~
# IN
# d.pred = estimated risks, ix of cases and controls
# risk.thr = risk threshold, could be a vector
# OUT
# matrix with two columns, TP and FP, estimated at each thtreshold in the d.pred dataset.
# if risk.thr is specified, returns a matrix of length(risk.thr) x 2
#       ie. TP and FP values for the thresholds in risk.thr
get.tpfp <- function(d.pred, risk.thr = NULL){
    
    if(is.null(risk.thr)){
        risk.thr.tp <- d.pred$risk
        risk.thr.fp <- d.pred$risk
    }else{ 
        if(length(risk.thr) == 2){
            risk.thr.tp <- risk.thr[1]
            risk.thr.fp <- risk.thr[2]
        }else{ # assumes single value if got to here
            risk.thr.tp <- risk.thr
            risk.thr.fp <- risk.thr
        }
    }
    
    # risk.thr.tp and risk.thr.fp have the same length
    tpfp <- matrix(NA, nrow = length(risk.thr.tp), ncol = 2)
    
    for(i in 1:length(risk.thr.tp)){ 
        tpfp[i, 1] <- sum(d.pred$weights.case[d.pred$risk.case > risk.thr.tp[i]])
        tpfp[i, 2] <- sum(d.pred$weights.ctrl[d.pred$risk.ctrl > risk.thr.fp[i]])
    }
    tpfp[, 1] <- tpfp[, 1] / sum(d.pred$weights.case)
    tpfp[, 2] <- tpfp[, 2] / sum(d.pred$weights.ctrl)  
    
    # Risk should be sorted, and so should tpfp. Just in case, this sorts tpfp by increasing FP
    if(!all(tpfp[,2]==sort(tpfp[,2], decreasing = F))){
         print('Error in get.tpfp() - TP FP not sorted in increasing order!')
         ix.sort <- sort(tpfp[,2], decreasing = F, index.return = T)$ix # changed april 7, 2014
         tpfp <- tpfp[ix.sort, ]
    }
    return(tpfp)  
}




# ~~~~~~~~~~~~~~~
# ~~~ get.auc ~~~
# ~~~~~~~~~~~~~~~
get.auc <- function(d.pred){
    tpfp <- get.tpfp(d.pred)
    
    # just a precaution, they should all be sorted already.
    ix.sort <- sort(tpfp[,1], decreasing = F, index.return = T)$ix
    tpfp <- tpfp[ix.sort,]
    
    x <- c(0, tpfp[, 2])
    y <- c(tpfp[, 1], 1)
    n <- length(x)
    dx <- x[-1] - x[-n]
    mid.y <- (y[-n] + y[-1])/2
    rm(x,y,n)
    return(sum(dx*mid.y))
}




# ~~~~~~~~~~~~~~~~
# ~~~ get.pcf ~~~
# ~~~~~~~~~~~~~~~~
# proportion of cases followed
# assumes that risk is sorted in a decreasing order.
get.pcf <- function(d.pred, prop){
    weighted.prop <- prop * (sum(d.pred$weights))
    
    population.q <- 0
    weights.q.cases <- 0
    count <- 1
    
    while(population.q <= weighted.prop){ 
        population.q <- population.q + d.pred$weights[count]
        if(d.pred$ix.case[count] == T){
            weights.q.cases <- weights.q.cases + d.pred$weights[count] 
        }
        count <- count + 1
    }
    pcf <- weights.q.cases/sum(d.pred$weights.case)
    
    return(pcf)
}





# ~~~~~~~~~~~~~~~~
# ~~~ get.ppnf ~~~
# ~~~~~~~~~~~~~~~~
# proportion of the population needed to be followed (to capture prop of the cases)
#  we count the number of cases at highest risk to get to prop of cases
#  then we look at how many of all subjects are above that cutoff
#
# assumes that the subjects are sorted by decreasing risk 
get.ppnf <- function(d.pred, prop){
    
    reached.prop.cases <- prop * sum(d.pred$weights.case) 
    
    sum.q.case <- 0
    sum.q.ctrl <- 0 
    
    # integer counters for cases, controls and the entire sample
    count.top.q.cases  <- 0
    count.top.q.ctrls  <- 0
    count.top.q.sample <- 0 
    
    while(sum.q.case <= reached.prop.cases){
        count.top.q.sample <- count.top.q.sample + 1
        if(d.pred$ix.case[count.top.q.sample] == T){
            count.top.q.cases <- count.top.q.cases + 1
            sum.q.case <- sum.q.case + d.pred$weights.case[count.top.q.cases]
        }else{
            count.top.q.ctrls <- count.top.q.ctrls + 1
            sum.q.ctrl <- sum.q.ctrl + d.pred$weights.ctrl[count.top.q.ctrls]
        }
    }
    sum.num <- sum.q.case + sum.q.ctrl
    sum.den <- sum(d.pred$weights.case) + sum(d.pred$weights.ctrl)
    return(sum.num/sum.den) 
}






get.prediction.error <- function(d.pred){
    mean.pe.case <- sum((1-d.pred$risk.case)^2 * d.pred$weights.case)
    mean.pe.ctrl <- sum((d.pred$risk.ctrl)^2 * d.pred$weights.ctrl)
    
    return((mean.pe.case + mean.pe.ctrl)/sum(d.pred$weights))
} 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions for getting nice output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mx is a matrix with columns of est and se, est, se etc.
get.nice.table <- function(mx, row.names = NULL, digits = 3){
    
    n.cols <- dim(mx)[2]
    mx <- formatC(mx, digits = digits, format = 'f')
    out <- NULL
    for(i in seq(1, n.cols, 2)){
        out <- cbind(out, paste(mx[,i], ' (', mx[,i+1], ')', sep = ''))
    }
    rownames(out) <- row.names
    return(out)
}






# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# summaries of various types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# get estimates of mean
# est.mx = n.stat x K
get.estimates <- function(est.mx){
    return(rowMeans(est.mx, na.rm = T))
}


# get standard deviations
# est.mx = n.stat x K
get.se.emp <- function(est.mx){
    return(apply(est.mx, 1, sd, na.rm = T))
}



# get standard errors
# var.mx = n.stat x K
get.se.pert <- function(var.mx){
    if(is.na(var.mx)){
        return(NA)
    }
    vars <- rowMeans(var.mx, na.rm = T)
    return(sqrt(vars))
}



