####################################################################################################
# Constraining values of regression coefficients
####################################################################################################
trans = function(c0){
  cc = c0
  
  #Make sure autoregression coefficients (phi and psi) remain stable
  iter = 1
  for(i in c("phi", rep("psi", length(names(c0)[grepl("psi", names(c0))])/2))){
    if(i == "phi"){
      v = paste0("phi", 1:2)
    }else if(i == "psi"){
      v = paste0(paste0("psi", iter), 1:2)
      iter = iter + 1
    }
    c1 = c0[v[1]]/(1 + abs(c0[v[1]]))
    c2 = c0[v[2]]/(1 + abs(c0[v[2]]))
    cc[v[1]] = c1 + c2
    cc[v[2]] = -1*c1*c2
  }
  
  #Keep probabilities between 0 and 1
  cc[grepl("p_", names(cc))] = exp(c0[grepl("p_", names(cc))])/(1 + exp(c0[grepl("p_", names(cc))]))
  
  cc[grepl("mu_d", names(cc))] = -exp(c0[grepl("mu_d", names(c0), ignore.case = T)])
  cc[grepl("mu_u", names(cc))] = exp(c0[grepl("mu_u", names(c0), ignore.case = T)])
  
  return(cc)
}

####################################################################################################
# Reverse the transformation for initial values
####################################################################################################
init = function(cc){
  c0 = cc
  
  iter = 1
  for(i in c("phi", rep("psi", length(names(cc)[grepl("psi", names(cc))])/2))){
    if(i == "phi"){
      v = paste0("phi", 1:2)
    }else if(i == "psi"){
      v = paste0(paste0("psi", iter), 1:2)
      iter = iter + 1
    }
    
    c2 = 1/2*(cc[v[1]] + sqrt(cc[v[1]]^2 + 4*cc[v[2]]))
    if(is.na(c2)){
      c2 = Inf
    }
    c1 = -cc[v[2]]/c2
    
    c0[v[1]] = ifelse(c1 < 0, c1/(1 + c1), c1/(1 - c1))
    c0[v[2]] = ifelse(c2 < 0, c2/(1 + c2), c2/(1 - c2))
  }
  
  c0[grepl("p_", names(cc))] = log(cc[grepl("p_", names(cc))]/(1 - cc[grepl("p_", names(cc))]))
  
  c0[grepl("mu_d", names(cc))] = log(-(cc[grepl("mu_d", names(c0), ignore.case = T)]))
  c0[grepl("mu_u", names(cc))] = log(cc[grepl("mu_u", names(c0), ignore.case = T)])
  
  c0[is.na(c0)] = 0
  
  return(c0)
}

#' State space model for markov switching mean
#' Creates a state space model in list form
#' yt = Ht %*% Bt + e_t
#' Bt = Mu_ts + Ft %*% B_{t=1} + u_t
#'
#' @param par Vector of named parameter values
#' @param yt Multivariate time series of data values
#' @param panelID Column name that identifies the cross section of the data
#' @param timeID Column name that identifies the date
#' @return List of space space matrices
#' @examples
#' SSmodel_ms(par, yt, panelID, timeID)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
SSmodel_ms = function(par, yt, panelID = NULL, timeID = NULL){
  #Get the parameters
  vars = dimnames(yt)[which(unlist(lapply(dimnames(yt), function(x){!is.null(x)})))][[1]]
  vars = vars[!vars %in% c(panelID, timeID)]
  phi = par[grepl("phi", names(par))]
  names(phi) = gsub("phi", "", names(phi))
  gamma = par[grepl("gamma", names(par))]
  names(gamma) = gsub("gamma", "", names(gamma))
  psi = par[grepl("psi", names(par))]
  names(psi) = gsub("psi", "", names(psi))
  sig = par[grepl("sigma", names(par))]
  names(sig) = gsub("sigma", "", names(sig))
  mu = par[grepl("mu", names(par))]
  names(mu) = gsub("mu_", "", names(mu))
  pr = par[grepl("p_", names(par))]
  names(pr) = gsub("p_", "", names(pr))
  
  #Build the transition equation matrix
  Fm = rbind(phi, c(1, 0))
  Fm = matrix(0, nrow = length(phi) + length(psi), ncol = length(phi) + length(psi))
  rownames(Fm) = c("ct0", "ct1", unlist(lapply(paste0("e", 1:length(vars)), function(x){paste0(x, 0:1)})))
  colnames(Fm) = c("ct1", "ct2", unlist(lapply(paste0("e", 1:length(vars)), function(x){paste0(x, 1:2)})))
  for(j in seq(1, nrow(Fm), 2)){
    Fm[j:(j + 1), j:(j + 1)] = t(matrix(c(c(phi, psi)[j:(j + 1)], 1, 0), nrow = 2, ncol = 2))
  }
  if(length(gamma) > length(vars)){
    Fm = cbind(Fm[, 1:2], matrix(0, ncol = max(table(substr(names(gamma)[nchar(names(gamma)) == 2], 1, 1))) - 1, nrow = nrow(Fm)), Fm[, 3:ncol(Fm)])
    Fm = rbind(Fm[1:2, ], matrix(0, ncol = ncol(Fm), nrow = max(table(substr(names(gamma)[nchar(names(gamma)) == 2], 1, 1))) - 1), Fm[3:nrow(Fm), ])
    colnames(Fm)[colnames(Fm) == ""] = paste0("ct", 3:(2 + length(colnames(Fm)[colnames(Fm) == ""])))
    rownames(Fm)[rownames(Fm) == ""] = paste0("ct", 2:(1 + length(rownames(Fm)[rownames(Fm) == ""])))
    diag(Fm[paste0("ct", 1:(max(table(substr(names(gamma)[nchar(names(gamma)) == 2], 1, 1))))), paste0("ct", 1:(max(table(substr(names(gamma)[nchar(names(gamma)) == 2], 1, 1)))))]) = 1
  }
  
  #Build the observation equation matrix
  Hm = matrix(gamma[1:length(vars)], ncol = 1)
  if(length(gamma) > length(vars)){
    mat = matrix(0, nrow = nrow(Hm), ncol = max(table(substr(names(gamma)[nchar(names(gamma)) == 2], 1, 1))))
    for(i in names(gamma[nchar(names(gamma)) == 2])){
      rc = as.numeric(strsplit(i, "")[[1]])
      mat[rc[1], rc[2]] = gamma[i]
    }
    Hm = cbind(Hm, mat)
  }else{
    Hm = cbind(Hm, matrix(0, ncol = 1, nrow = nrow(Hm)))
  }
  Hm = cbind(Hm, matrix(0, nrow = nrow(Hm), ncol = length(psi)))
  rownames(Hm) = vars
  colnames(Hm) = rownames(Fm)  
  diag(Hm[, grepl(paste0("e", paste0(1:length(vars), "0"), collapse = "|"), colnames(Hm))]) = 1
  
  #Build the transition equation covariance matrix
  #Set the dynamic common factor standard deviation to 1
  Qm = matrix(0, ncol = ncol(Fm), nrow = nrow(Fm))
  rownames(Qm) = rownames(Fm)
  colnames(Qm) = rownames(Qm)
  diag(Qm[rownames(Qm)[substr(rownames(Qm), nchar(rownames(Qm)), nchar(rownames(Qm))) == "0"], rownames(Qm)[substr(rownames(Qm), nchar(rownames(Qm)), nchar(rownames(Qm))) == "0"]]) = c(1, sig^2)
  
  #Build the observation eqution covariance matrix
  Rm = matrix(0, ncol = nrow(Hm), nrow = nrow(Hm))#diag(x = sig^2, nrow = length(sig), ncol = length(sig))
  colnames(Rm) = vars
  rownames(Rm) = vars
  
  #Build the Markov-switching mean matrix
  Mu_d = t(t(c(ifelse(length(mu[grepl("d", tolower(names(mu)))]) == 0, 0, mu[grepl("d", tolower(names(mu)))]), rep(0, ncol(Fm) - 1))))
  Mu_m = t(t(c(ifelse(length(mu[grepl("m", tolower(names(mu)))]) == 0, 0, mu[grepl("m", tolower(names(mu)))]), rep(0, ncol(Fm) - 1))))
  Mu_u = t(t(c(ifelse(length(mu[grepl("u", tolower(names(mu)))]) == 0, 0, mu[grepl("u", tolower(names(mu)))]), rep(0, ncol(Fm) - 1))))
  Mu = array(c(Mu_d = Mu_d, Mu_m = Mu_m, Mu_u = Mu_u), dim = c(nrow(Fm), 1, 3),
             dimnames = list(NULL, NULL, c("d", "m", "u")))
  
  #Observaton equation intercept matrix
  Am = matrix(0, nrow = 1, ncol = 1)
  
  #Steady state probabilities
  Tr_mat = matrix(0, nrow = 3, ncol = 3)
  rownames(Tr_mat) = colnames(Tr_mat) = unique(unlist(lapply(names(pr), function(x){strsplit(x, "")[[1]][2]})))
  for(j in names(pr)){
    Tr_mat[strsplit(j, "")[[1]][2], strsplit(j, "")[[1]][1]] = pr[j]
  }
  for(j in names(which(apply(Tr_mat, 2, sum) != 0))){
    row = rownames(Tr_mat)[!rownames(Tr_mat) %in% substr(names(pr)[substr(names(pr), 1, 1) == j], 2, 2)]
    row = row[!row %in% names(which(apply(Tr_mat, 2, sum) == 0))]
    Tr_mat[row, j] = 1 - sum(Tr_mat[names(which(apply(Tr_mat, 2, sum) != 0)), j])
  }
  
  #Initialize the filter for each state
  B0 = array(unlist(lapply(colnames(Tr_mat), function(x){ginv(diag(ncol(Fm)) - Fm) %*% Mu[,, x]})),
             dim = dim(Mu), dimnames = dimnames(Mu))
  VecP0 = ginv(diag(prod(dim(Fm))) - kronecker(Fm, Fm)) %*% matrix(as.vector(Qm), ncol = 1)
  P0 = array(unlist(lapply(dimnames(B0)[[3]], function(x){matrix(VecP0[, 1], ncol = ncol(Fm))})), 
             dim = c(nrow(Fm), ncol(Fm), dim(B0)[3]), dimnames = dimnames(B0))
  
  return(list(B0 = B0, P0 = P0, At = Am, Mu = Mu, Ht = Hm, Ft = Fm, Qt = Qm, Rt = Rm, Tr_mat = Tr_mat))
}

#' Markov switching model estimation by the Kim filter (Hamilton + Kalman filter)
#
#' Estimates a Markov switching mean state space model
#' The function checks for growth and unit variables in the data set provided.
#' It will log all growth variables and difference all unit root variables.
#' It also detects the frequency of the data if it is not provided.
#'
#' Model priors (initial values) 
#' Parameter naming convention
#'   DCF AR coefficients: phi1 and phi2 only
#'   Error MA coefficients: psi_i1 to psi_i2 only for each series i
#'   Error standard deviation: sigma_i only for each series i
#'   Observation coefficient on DCF:
#'     First gamma (gamma_i, ..., gamma_n) only 1 index number, not i0
#'     Any more gammas per equation: gamma_i1 to gamma_ik
#'   Markov switching growth rate: mu_d and mu_u
#'   Transition probabilities: p_dd, p_md (or p_mu), p_mm, p_md (or p_mu), p_uu, p_ud (or p_um)
#'
#'
#' @param y Multivariate time series of data values. May also be a data frame containing a date column
#' @param freq Seasonality of the data
#' @param panelID Column name that identifies the cross section of the data
#' @param timeID Column name that identifies the date
#' @param level Significance level of statistical tests (0.01, 0.05, 0.1)
#' @param formulas R formula describing the relationship between each data series, the unobserved dynamic common factor, and the erorr structure
#' @param optim_methods Vector of 1 to 3 optimization methods in order of preference ("NR", "BFGS", "CG", "BHHH", or "SANN")
#' @param maxit Maximum number of iterations for the optimization
#' @param maxtrials Maximum number of optimization trials to get convergence
#' @param theta Optional vector of named initial parameter guesses: DCF AR coefficients: phi1 and phi2 only; Error MA coefficients: psi_i1 to psi_i2 only for each series i; Error standard deviation: sigma_i only for each series i; Observation coefficient on DCF with first gamma (gamma_i, ..., gamma_n) only 1 index number, not i0 and any more gammas per equation: gamma_i1 to gamma_ik; Markov switching growth rate: mu_d and mu_u; Transition probabilities: p_dd, p_md (or p_mu), p_mm, p_md (or p_mu), p_uu, p_ud (or p_um)
#' @return List of estimation values including coefficients, convergence code, the panel and time ids, the variables in the data, the variable that were logged, and the variables that were differenced
#' @examples
#' ms_dcf_estim(y = DT[, c("date", "y")])
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
ms_dcf_estim = function(y, freq = NULL, panelID = NULL, timeID = NULL, level = 0.01,
                        formulas = c("y ~ c + e.l1 + e.l2"), theta = NULL,
                        optim_methods = c("BFGS", "CG", "NM"), maxit = 1000, maxtrials = 10, trace = F){
  
  dates = NULL
  if(is.null(freq)){
    #stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    y = data.table::as.data.table(y)
    `.SD` = data.table::.SD
    if(any(unlist(y[, lapply(.SD, function(x){class(x) == "Date"})]))){
      #Detect the frequency
      timeID = names(which(unlist(y[, lapply(.SD, function(x){class(x) == "Date"})])))
      if(length(timeID) > 1){
        stop("Too many date columns. Include only 1 date column or set the frequency manually.")
      }
      datediffs = unique(diff(y[, timeID, with = F][[1]]))
      freq = datediffs[which.max(tabulate(match(diff(y[, timeID, with = F][[1]]), datediffs)))]
      freq = c(365, 52, 12, 4, 1)[which.min(abs(freq - c(1, 7, 30, 90, 365)))]
      dates = y[, timeID, with = F][[1]]
      rm(datediffs)
    }else{
      stop("No date column detected. Include a date column or set the frequency.")
    }
  }else{
    if(!is.numeric(freq)){
      stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    }else if(!freq %in% c(1, 4, 12, 52, 365)){
      stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    }
  }
  if(level < 0.01 | level > 0.1){
    stop("level must be between 0.01 and 0.1.")
  }
  if(trace == F){
    trace = 0
  }else{
    trace = 2
  }

  
  if(is.null(panelID)){
    panelID = "panelid"
    y[, "panelid" := "panel"]
  }
  vars = colnames(y)[!colnames(y) %in% c(panelID, timeID)]
  
  objective = function(par, yt, panelID, timeID){
    par = trans(par)
    sp = SSmodel_ms(par, yt, panelID, timeID)
    toret = foreach::foreach(i = unique(yt[, c(panelID), with = F][[1]]), .packages = c("data.table")) %fun% {
      yti = t(yt[eval(parse(text = panelID)) == i, colnames(yt)[!colnames(yt) %in% c(panelID, timeID)], with = F])
      ans = kim_filter(sp$B0, sp$P0, sp$Mu, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
      return(ans$loglik)
    }
    return(mean(unlist(toret)))
  }
  
  #Check for growth variables and log them
  gr.test = colMeans(y[, lapply(.SD, function(x){
    d = x - shift(x, type = "lag", n = 1)
    return(t.test(d[!is.na(d)], mu = 0)$p.value)
  }), by = c(panelID), .SDcols = c(vars)][, c(vars), with = F])
  log.vars = names(gr.test)[which(gr.test <= level)]
  if(length(log.vars) > 0){
    y[, c(log.vars) := lapply(.SD, log), .SDcols = c(log.vars)]
  }
  
  #Unit root tests
  ur.vars = lapply(vars, function(x){
    ret = lapply(unique(y[, c(panelID), with = F][[1]]), function(z){
      tseries::adf.test(x = ts(y[eval(parse(text = panelID)) == z, c(x), with = F][[1]], freq = freq), alternative = "stationary")$p.value
    })
    names(ret) = unique(y[, c(panelID), with = F][[1]])
    return(ret)
  })
  names(ur.vars) = vars
  ur.vars = lapply(names(ur.vars), function(x){
    mean(unlist(ur.vars[[x]]))
  })
  names(ur.vars) = vars
  ur.vars = names(ur.vars)[ur.vars >= level]
  
  #First difference the relevant series
  yy_d = copy(y)
  if(length(ur.vars) > 0){
    yy_d[, c(ur.vars) := lapply(.SD, function(x){
      x - data.table::shift(x, type = "lag")
    }), by = c(panelID), .SDcols = c(ur.vars)]
  }
  yy_d = yy_d[complete.cases(yy_d), ]
  
  #Standardize the data
  yy_s = copy(yy_d)
  yy_s[, c(vars) := lapply(.SD, function(x){
    (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  }), by = c(panelID), .SDcols = vars]
  
  if(is.null(theta)){
    #Find the starting values
    y[, "C" := Matrix::rowMeans(y[, vars, with = F])]
    y[, "C" := lapply(.SD, smooth, twiceit = T), by = c(panelID), .SDcols = "C"]
    y[, "c" := C - data.table::shift(C, type = "lag", n = 1), by = c(panelID)]
    y[, "c" := (c - mean(c, na.rm = T))/sd(c, na.rm = T), by = c(panelID)]
    
    #Prior for the AR coefficients (phi)
    up = lapply(unique(y[, c(panelID), with = F][[1]]), function(x){
      c = y[eval(parse(text = panelID)) == x, ]$c
      return(forecast::Arima(c[which((c - mean(c, na.rm = T))/sd(c, na.rm = T) > 0.44)], order = c(2, 0, 0), include.mean = T))
    })
    mid = lapply(unique(y[, c(panelID), with = F][[1]]), function(x){
      c = y[eval(parse(text = panelID)) == x, ]$c
      return(forecast::Arima(c[which((c - mean(c, na.rm = T))/sd(c, na.rm = T) <= 0.44 & (c - mean(c, na.rm = T))/sd(c, na.rm = T) >= -0.44)], order = c(2, 0, 0), include.mean = F))
    })
    down = lapply(unique(y[, c(panelID), with = F][[1]]), function(x){
      c = y[eval(parse(text = panelID)) == x, ]$c
      return(forecast::Arima(c[which((c - mean(c, na.rm = T))/sd(c, na.rm = T) < -0.44)], order = c(2, 0, 0), include.mean = T))
    })
    names(up) = names(mid) = names(down) = unique(y[, c(panelID), with = F][[1]])
    
    theta = colMeans(do.call("rbind", lapply(names(up), function(x){
      c = y[eval(parse(text = panelID)) == x, ]$c
      coeff = data.table::data.table(cbind(up = up[[x]]$coef, mid = c(mid[[x]]$coef, 0), down = down[[x]]$coef), keep.rownames = T)
      coeff[, "m" := up*length(c[c > quantile(c, probs = 0.67, na.rm = T)])/length(c) + 
              mid*length(c[c >= quantile(c, probs = 0.33, na.rm = T) & c <= quantile(c, probs = 0.67, na.rm = T)])/length(c) + 
              down*length(c[c < quantile(c, probs = 0.33, na.rm = T)])/length(c)]
      ret = coeff$m
      names(ret) = coeff$rn
      return(ret)
    })))[c("ar1", "ar2")]
    theta = c(theta, colMeans(do.call("rbind", lapply(names(up), function(x){
      c = y[eval(parse(text = panelID)) == x, ]$c
      coeff = data.table::data.table(cbind(up = up[[x]]$coef, mid = c(mid[[x]]$coef, 0), down = down[[x]]$coef), keep.rownames = T)[rn == "intercept", ]
      ret = unlist(c(coeff[, c("up", "mid", "down"), with = F]))
      return(ret)
    })))[c("up", "down")])
    names(theta) = gsub("ar", "phi", names(theta))
    names(theta)[names(theta) == "up"] = "mu_u"
    names(theta)[names(theta) == "down"] = "mu_d"
    
    #Prior for the Markov-switching means and probabilities
    y[, "cd" := C - data.table::shift(C, type = "lag", n = 1), by = c(panelID)]
    y[, "S_d" := ifelse((cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) < -0.44, 1, 0)]
    y[, "S_m" := ifelse((cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) >= -0.44 & 
                                (cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) <= 0.44, 1, 0)]
    y[, "S_u" := ifelse((cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) > 0.44, 1, 0)]
    y[, "S_dd" := ifelse(S_d == 1 & data.table::shift(S_d, type = "lead", n = 1) == 1, 1, 0)]
    y[, "S_mm" := ifelse(S_m == 1 & data.table::shift(S_m, type = "lead", n = 1) == 1, 1, 0)]
    y[, "S_uu" := ifelse(S_u == 1 & data.table::shift(S_u, type = "lead", n = 1) == 1, 1, 0)]
    
    theta = c(theta, 
              p_dd = 0.90*sum(y$S_dd, na.rm = T)/sum(y$S_d, na.rm = T),
              p_dm = (1 - 0.9*sum(y$S_dd, na.rm = T)/sum(y$S_d, na.rm = T))/2,
              p_md = (1 - 0.9*sum(y$S_mm, na.rm = T)/sum(y$S_m, na.rm = T))/2,
              p_mm = 0.9*sum(y$S_mm, na.rm = T)/sum(y$S_m, na.rm = T),
              p_ud = (1 - 0.9*sum(y$S_uu, na.rm = T)/sum(y$S_u, na.rm = T))/2,
              p_uu = 0.9*sum(y$S_uu, na.rm = T)/sum(y$S_u, na.rm = T))
    
    #Define the equations by series
    if(length(formulas) < length(vars)){
      formulas = rep(formulas, length(vars))[1:length(vars)]
    }
    if(is.null(names(formulas))){
      names(formulas) = vars
    }
    
    #Prior for the MA coefficients (psi) as wella s the gamma and sigma values 
    theta = c(theta, unlist(lapply(vars, function(z){
      #Average the series
      ytemp = Matrix::rowMeans(yy_s[, c(z), with = F])
      ts = data.table::data.table(y = yy_s[, c(z), with = F], c = y[!is.na(c), ]$c)
      colnames(ts) = c("y", "c")
      #Use up to 1 quarter lags of the DCF for each series
      miter = ifelse(freq == 12, 3, ifelse(freq == 4, 1, 0))
      if(miter > 0){
        for(m in 1:miter){
          ts[, paste0("c.l", m) := shift(c, type = "lag", n = m)]
        }
      }
      ts = ts[complete.cases(ts), ]
      #Step 1: Regress the series on the DCF and any lags
      lm.vars = trimws(strsplit(formulas[z], "\\~|\\+")[[1]])
      lm.vars = lm.vars[lm.vars %in% colnames(ts)]
      fit = lm(as.formula(paste("y ~", paste(lm.vars[2:length(lm.vars)], collapse = " + "))), data = ts)
      
      #Get the errors
      ts = cbind(ts, e = fit$residuals)
      ts[, "e.l1" := shift(e, type = "lag", n = 1)]
      ts[, "e.l2" := shift(e, type = "lag", n = 2)]
      ts = ts[complete.cases(ts), ]
      #Step 2: Regress the series on the DCF and the errors from step 1
      fit = lm(as.formula(paste(formulas[z], "- 1")), data = ts[, colnames(ts)[colnames(ts) != "e"], with = F])
      
      #Get the estimated coefficients
      idx = which(vars == z)
      coeff = fit$coefficients
      names(coeff)[names(coeff) == "c"] = paste0("gamma", idx)
      names(coeff)[grepl("c\\.l", names(coeff))] = paste0(paste0("gamma", idx), 1:(length(names(coeff)[grepl("c\\.l", names(coeff))])))
      names(coeff)[grepl("e\\.l", names(coeff))] = paste0(paste0("psi", idx), 1:length(names(coeff)[grepl("e\\.l", names(coeff))]))
      
      #Get the sigma parameter
      ts$e.fit = as.matrix(ts[, colnames(ts)[grepl("e\\.", colnames(ts))], with = F]) %*% as.matrix(coeff[grepl("psi", names(coeff))])
      ts[, "eps" := e - e.fit]
      coeff = c(coeff, sd(ts$eps, na.rm = T))
      names(coeff)[length(coeff)] = paste0("sigma", idx)
      return(coeff)
    })))
    rm(up, mid, down)
  }
  y = y[, c(panelID, timeID, vars), with = F]
  
  #Estimate the model
  cl = parallel::makeCluster(min(c(length(unique(yy_s[, c(panelID), with = F][[1]])), parallel::detectCores())))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  out = tryCatch(maxLik::maxLik(logLik = objective, start = init(theta), method = optim_methods[1],
                                finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), #init = init,
                                yt = yy_s, panelID = panelID, timeID = timeID),
                 error = function(err){
                   tryCatch(maxLik::maxLik(logLik = objective, start = init(theta), method = optim_methods[min(c(2, length(optim_methods)))],
                                           finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), #init = init,
                                           yt = yy_s, panelID = panelID, timeID = timeID),
                            error = function(err){
                              tryCatch(maxLik::maxLik(logLik = objective, start = init(theta), method = optim_methods[min(c(3, length(optim_methods)))],
                                                      finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), #init = init,
                                                      yt = yy_s, panelID = panelID, timeID = timeID),
                                       error = function(err){NULL})
                            })
                 })
  #Attempt to get convergence if it failed the first time
  if(!is.null(out)){
    trials = 1
    while(out$code != 0 & trials < maxtrials){
      out2 = tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[1],
                                     finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), #init = init,
                                     yt = yy_s, panelID = panelID, timeID = timeID),
                      error = function(err){
                        tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(2, length(optim_methods)))],
                                                finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), #init = init,
                                                yt = yy_s, panelID = panelID, timeID = timeID),
                                 error = function(err){
                                   tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(3, length(optim_methods)))],
                                                           finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), #init = init,
                                                           yt = yy_s, panelID = panelID, timeID = timeID),
                                            error = function(err){NULL})
                                 })
                      })
      #End the loop if no parameters changed or estimation failed
      if(!is.null(out2) & !all(coef(out) == coef(out2))){
        out = out2
        trials = trials + 1
      }else{
        break
      }
    }
  }
  snow::stopCluster(cl)
  return(list(coef = trans(coef(out)), convergence = out$code, loglik = out$maximum, panelID = panelID, timeID = timeID,
              vars = vars, log.vars = log.vars, ur.vars = ur.vars))
}

#' Kim filter (Hamilton + Kalman filter) an estimated model from ms_dcf_estim
#' @param y Multivariate time series of data values. May also contain a date column.
#' @param model Structural time series model estimated using ms_dcf_estim.
#' @param plot Logial, whether to plot the output or not.
#' @return List of data tables containing the filtered and smoothed series.
#' @examples
#' ms_dcf_filter(y = DT[, c("date", "y")], model = model)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
ms_dcf_filter = function(y, model, plot = F){
  
  #Get the data
  if(any(grepl(model$panelID, colnames(y)))){
    y = y[, c(model$panelID, model$timeID, model$vars), with = F]
  }else{
    y[, "panelid" := "panel"]
  }
  
  #Log the relevant variables
  if(length(model$log.vars) > 0){
    y[, c(model$log.vars) := lapply(.SD, log), .SDcols = c(model$log.vars)]
  }
  
  #Difference the relevant variables
  yy_d = copy(y)
  if(length(model$ur.vars) > 0){
    yy_d[, c(model$ur.vars) := lapply(.SD, function(x){
      x - data.table::shift(x, type = "lag")
    }), by = c(model$panelID), .SDcols = c(model$ur.vars)]
  }
  yy_d = yy_d[complete.cases(yy_d), ]
  
  #Standardize the data
  yy_s = copy(yy_d)
  yy_s[, c(model$vars) := lapply(.SD, function(x){
    (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  }), by = c(model$panelID), .SDcols = c(model$vars)]
  
  #Filter the data using the estimated model
  sp = SSmodel_ms(model$coef, yy_s, model$panelID, model$timeID)
  cl = parallel::makeCluster(min(c(length(unique(yy_s[, c(model$panelID), with = F][[1]])), parallel::detectCores())))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  uc = foreach::foreach(i = unique(yy_s[, c(model$panelID), with = F][[1]]), .packages = c("data.table", "MASS"), .export = c("SSmodel_ms", "kim_filter", "kim_smoother", "ss_prob", "v_prob")) %fun% {
    yti = t(yy_s[eval(parse(text = model$panelID)) == i, colnames(yy_s)[!colnames(yy_s) %in% c(model$panelID, model$timeID)], with = F])
    ans = kim_filter(sp$B0, sp$P0, sp$Mu, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
    
    #Get the steady state Kalman gain approximation
    K = ans$K[,, dim(ans$K)[3]]
    W = ginv(diag(nrow(K)) - (diag(nrow(K)) - K %*% sp$Ht) %*% sp$Ft) %*% K
    means = colMeans(as.matrix(yy_d[eval(parse(text = model$panelID)) == i, c(model$vars), with = F]))
    sds = matrixStats::colSds(as.matrix(yy_d[eval(parse(text = model$panelID)) == i, c(model$vars), with = F]))
    
    #W(1) is the first row of W
    d = W[1, ] %*% means
    
    #Get the intercept terms
    D = means - sp$Ht %*% matrix(c(rep(d, length(colnames(sp$Ft)[grepl("ct", colnames(sp$Ft))])), rep(0, length(colnames(sp$Ft)[grepl("e", colnames(sp$Ft))]))), ncol = 1)
    
    #Initialize first element of the dynamic common factor
    Y1 = t(y[eval(parse(text = model$panelID)) == i, ][eval(parse(text = model$timeID)) == min(eval(parse(text = model$timeID)), na.rm = T), c(model$vars), with = F])
    initC = function(par){
      return(sum((Y1 - D - sp$Ht %*% matrix(c(par, rep(0, length(colnames(sp$Ft)[grepl("e", colnames(sp$Ft))])))))^2))
    }
    C10 = optim(par = rep(mean(Y1)/mean(model$coef[grepl("gamma", names(model$coef))]), length(colnames(sp$Ft)[grepl("ct", colnames(sp$Ft))])), 
                fn = initC, method = "BFGS", control = list(trace = F))$par[1]
    dcf_loc = which(rownames(sp$Ft) == "ct0")
    
    toret = list()
    for(k in c("filter", "smooth")){
      Ctt = rep(C10, ncol(yti) + 1)
      if(k == "smooth"){
        ans = kim_smoother(B_tlss = ans$B_tlss, B_tts = ans$B_tts, B_tt = ans$B_tt, P_tlss = ans$P_tlss, P_tts = ans$P_tts, 
                          Pr_tls = ans$Pr_tls, Pr_tts = ans$Pr_tts, Ft = sp$Ft, Mu = sp$Mu, Tr_mat = sp$Tr_mat)
      }
        
      #Build the rest of the DCF series
      #Rescale the first differenced DCF series to match that of the actual series
      ctt = c(0, (t(ans$B_tt[, dcf_loc]))*mean(sds)/sd(ans$B_tt[, dcf_loc]))
      mutt = c(0, t(ans$Mu_t[, dcf_loc])*mean(sds)/sd(ans$B_tt[, dcf_loc]))
      for(j in 2:length(Ctt)){
        #First element of dCtt is C_22 - C_11 = dCtt_2
        Ctt[j] = ctt[j] + Ctt[j - 1] + c(d)
      }
      Ctt = data.table::data.table(panelID = i, date = y[eval(parse(text = model$panelID)) == i, c(model$timeID), with = F][[1]], Ctt = Ctt, ctt = ctt, mutt = mutt)
      colnames(Ctt) = c(model$panelID, model$timeID, "Ctt", "ctt", "mutt")
      
      #Build the probability series
      prob = data.table::data.table(panelID = i, date = yy_s[eval(parse(text = model$panelID)) == i, c(model$timeID), with = F][[1]], data.table::data.table(ans$Pr_tts))
      colnames(prob) = c(model$panelID, model$timeID, paste0("Pr_", dimnames(sp$Mu)[[3]]))
      toret[[k]] = merge(Ctt, prob, by = c(model$timeID, model$panelID), all = T)
    }
    return(toret)
  }
  snow::stopCluster(cl)
  names(uc) = unique(yy_s[, c(model$panelID), with = F][[1]])
  
  if(plot == T){
    for(j in c("filter", "smooth")){
      for(i in unique(y[, c(model$panelID), with = F][[1]])){  
        toplot1 = melt(y[eval(parse(text = model$panelID)) == i, ], id.vars = c(model$panelID, model$timeID))
        toplot1 = merge(toplot1, toplot1[eval(parse(text = model$timeID)) == min(eval(parse(text = model$timeID))), c("variable", "value"), with = F], by = c("variable"), all.x = T, all.y = F, suffixes = c("", "_start"))
        toplot1[, "value2" := value - (value_start - mean(toplot1[eval(parse(text = model$timeID)) == min(eval(parse(text = model$timeID))), ]$value))*9/10, by = "variable"]
        g1 = ggplot2::ggplot(toplot1) + 
          ggplot2::ggtitle(paste(i, "Data Series"), subtitle = "Levels") + 
          ggplot2::scale_y_continuous(name = "Log Levels (Intercept Adjusted)") + 
          ggplot2::scale_x_date(name = "") + 
          ggplot2::geom_line(ggplot2::aes(x = eval(parse(text = model$timeID)), y = value2, group = variable, color = variable)) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + ggplot2::guides(color = ggplot2::guide_legend(title = NULL))
        
        toplot2 = melt(yy_s[eval(parse(text = model$panelID)) == i, ], id.vars = c(model$panelID, model$timeID))
        g2 = ggplot2::ggplot(toplot2) + 
          ggplot2::ggtitle(paste(i, "Data Series"), subtitle = "Differenced & Standardized") + 
          ggplot2::scale_y_continuous(name = "Differences") + 
          ggplot2::scale_x_date(name = "") + 
          ggplot2::geom_hline(yintercept = 0, color = "black") + 
          ggplot2::geom_line(ggplot2::aes(x = eval(parse(text = model$timeID)), y = value, group = variable, color = variable)) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + ggplot2::guides(color = ggplot2::guide_legend(title = NULL))
        
        toplot3 = melt(uc[[i]][[j]], id.vars = c(model$panelID, model$timeID))
        d_range1 = range(toplot3[variable == "Ctt", ]$value, na.rm = T)
        p_range1 = range(toplot3[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], ]$value, na.rm = T)
        toplot3[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], "value" := (value - p_range1[1])/diff(p_range1) * diff(d_range1) + d_range1[1], by = "variable"]
        g3 = ggplot2::ggplot() +  
          ggplot2::ggtitle(paste(i, "Dynamic Common Factor"), subtitle = "Differenced") + 
          ggplot2::scale_y_continuous(name = "Dynamic Common Factor", limits = range(toplot3[variable == "Ctt", ]$value, na.rm = T), 
                                      sec.axis = ggplot2::sec_axis(name = "Probability of State (%)", ~((. - d_range1[1])/diff(d_range1) * diff(p_range1) + p_range1[1]) * 100)) + 
          ggplot2::scale_x_date(name = "") +
          ggplot2::geom_ribbon(data = toplot3[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_d|Pr_u", colnames(uc[[i]][[j]]), ignore.case = T)], ], 
                               ggplot2::aes(x = eval(parse(text = model$timeID)), ymin = d_range1[1], ymax = value, group = variable, fill = variable), alpha = 0.5) + 
          ggplot2::geom_hline(yintercept = 0, color = "grey") + 
          ggplot2::geom_line(data = toplot3[variable == "Ctt", ], 
                             ggplot2::aes(x = eval(parse(text = model$timeID)), y = value, group = variable, color = variable), color = "black") + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + 
          ggplot2::guides(color = ggplot2::guide_legend(title = NULL), fill = ggplot2::guide_legend(title = NULL))
        
        toplot4 = melt(uc[[i]][[j]], id.vars = c(model$panelID, model$timeID))
        d_range2 = range(toplot4[variable == "ctt", ]$value, na.rm = T)
        p_range2 = range(toplot4[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], ]$value, na.rm = T)
        toplot4[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], "value" := (value - p_range2[1])/diff(p_range2) * diff(d_range2) + d_range2[1], by = "variable"]
        g4 = ggplot2::ggplot() +  
          ggplot2::ggtitle(paste(i, "Dynamic Common Factor"), subtitle = "Differenced") + 
          ggplot2::scale_y_continuous(name = "Dynamic Common Factor", limits = range(toplot4[variable == "ctt", ]$value, na.rm = T), 
                             sec.axis = ggplot2::sec_axis(name = "Probability of State (%)", ~((. - d_range2[1])/diff(d_range2) * diff(p_range2) + p_range2[1]) * 100)) + 
          ggplot2::scale_x_date(name = "") +
          ggplot2::geom_ribbon(data = toplot4[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_d|Pr_u", colnames(uc[[i]][[j]]), ignore.case = T)], ], 
                               ggplot2::aes(x = eval(parse(text = model$timeID)), ymin = d_range2[1], ymax = value, group = variable, fill = variable), alpha = 0.5) + 
          ggplot2::geom_hline(yintercept = 0, color = "grey") + 
          ggplot2::geom_line(data = toplot4[variable == "ctt", ], 
                             ggplot2::aes(x = eval(parse(text = model$timeID)), y = value, group = variable, color = variable), color = "black") + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + 
          ggplot2::guides(color = ggplot2::guide_legend(title = NULL), fill = ggplot2::guide_legend(title = NULL))
        
        gridExtra::grid.arrange(g1, g2, g3, g4, layout_matrix = matrix(c(1, 3, 2, 4), nrow = 2),
                                top = grid::textGrob(j, gp = grid::gpar(fontsize = 20, font = 3)))
        Sys.sleep(0.1)
      }
    }
  }
  
  names(uc) = "US"
  uc[["CAN"]] = copy(uc[["US"]])
  uc[["CAN"]]$filter[, "panelid" := "CAN"]
  uc[["CAN"]]$smooth[, "panelid" := "CAN"]
  
  ret = list(filter = do.call("rbind", lapply(names(uc), function(x){
    uc[[x]]$filter
    })), 
    smooth = do.call("rbind", lapply(names(uc), function(x){
      uc[[x]]$smooth
    }))
  )
  
  return(ret)
}

########
#Call these to build the package
#devtools::document()
#devtools::build_vignettes()
#devtools::install()
#library(projectmap)
#git config remote.origin.url git@github.com:opendoor-labs/MarkovSwitchingDCF.git
