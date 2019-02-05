#' Constraining values of regression coefficients
#' @export
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
  
  #Keep probabilities between 0 and 1 and sum of probabilities between 0 and 1
  states = unique(substr(gsub("p_", "", names(cc)[grepl("p_", names(cc))]), 1, 1))
  ss = paste0("p_", paste0(states, states))
  cc[names(cc) %in% ss] = 1/(1 + exp(-c0[names(c0) %in% ss]))
  if(length(states) > 2){
    for(s in states){
      cc[names(cc)[grepl("p_", names(cc)) & !names(cc) %in% ss & substr(names(cc), 3, 3) == s]] =  (1 - cc[names(cc) %in% ss & substr(names(cc), 3, 3) == s])/(1 + exp(-c0[names(c0)[grepl("p_", names(c0)) & !names(c0) %in% ss & substr(names(c0), 3, 3) == s]]))
    }
  }
  
  #Force mu_d to be negative and mu_u to be positive
  cc[grepl("mu_d", names(cc))] = -exp(c0[grepl("mu_d", names(c0), ignore.case = T)])
  cc[grepl("mu_u", names(cc))] = exp(c0[grepl("mu_u", names(c0), ignore.case = T)])
  
  #Force sd multiples and sigma parameters to be positive
  cc[grepl("sd_", names(cc))] = exp(c0[grepl("sd_", names(c0))])
  cc[grepl("sigma", names(cc))] = exp(c0[grepl("sigma", names(c0))])
  
  return(cc)
}

#' Reverse the transformation for initial values
#' @export
init_trans = function(cc){
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
  
  #Keep probabilities between 0 and 1 and sum of probabilities between 0 and 1
  states = unique(substr(gsub("p_", "", names(cc)[grepl("p_", names(cc))]), 1, 1))
  ss = paste0("p_", paste0(states, states))
  c0[names(c0) %in% ss] = log(cc[names(cc) %in% ss]/(1 - cc[names(cc) %in% ss]))
  if(length(states) > 2){
    for(s in states){
      c0[names(c0)[grepl("p_", names(c0)) & !names(c0) %in% ss & substr(names(c0), 3, 3) == s]] = log(cc[names(c0)[grepl("p_", names(cc)) & !names(cc) %in% ss & substr(names(cc), 3, 3) == s]]/(1 - sum(cc[grepl("p_", names(cc)) & substr(names(cc), 3, 3) == s])))
    }
  }
  
  #Force mu_d to be negative and mu_u to be positive
  c0[grepl("mu_d", names(c0))] = log(-(cc[grepl("mu_d", names(c0), ignore.case = T)]))
  c0[grepl("mu_u", names(c0))] = log(cc[grepl("mu_u", names(c0), ignore.case = T)])
  
  #Force sd multiples and sigma parameters to be positive
  c0[grepl("sd_", names(c0))] = log(cc[grepl("sd_", names(cc))])
  c0[grepl("sigma", names(c0))] = log(cc[grepl("sigma", names(cc))])
  
  c0[is.na(c0)] = 0
  
  return(c0)
}

#' State space model for markov switching mean
#' Creates a state space model in list form
#' yt = At + Ht %*% Bt + e_t
#' Bt = D_t + Ft %*% B_{t=1} + u_t
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
SSmodel_ms = function(par, yt, n_states, ms_var, panelID = NULL, timeID = NULL, init = NULL){
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
  sd = par[grepl("sd", names(par))]
  names(sd) = gsub("sd_", "", names(sd))
  pr = par[grepl("p_", names(par))]
  names(pr) = gsub("p_", "", names(pr))
  if(length(sd) == 0 & n_states > 1 & is.finite(n_states)){
    sd = rep(1, n_states)
    names(sd) = substr(names(pr[substr(names(pr), 1, 1) == substr(names(pr), 2, 2)]), 1, 1)
  }
  
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
  if(is.infinite(n_states)){
    Fm = cbind(matrix(0, nrow = nrow(Fm), ncol = 1), Fm)
    Fm = rbind(matrix(0, nrow = 1, ncol = ncol(Fm)), Fm)
    colnames(Fm)[1] = "mt1"
    rownames(Fm)[1] = "mt0"
    Fm[c("mt0", "ct0"), "mt1"] = 1
  }
  if(n_states == 1 | is.infinite(n_states)){
    Fm = array(c(Fm), dim = c(nrow(Fm), ncol(Fm), 1),
               dimnames = list(rownames(Fm), colnames(Fm), c("m")))
  }else if(n_states == 2){
    Fm = array(c(Fm, Fm), dim = c(nrow(Fm), ncol(Fm), n_states),
               dimnames = list(rownames(Fm), colnames(Fm), c("d", "u")))
  }else if(n_states == 3){
    Fm = array(c(Fm, Fm, Fm), dim = c(nrow(Fm), ncol(Fm), n_states),
               dimnames = list(rownames(Fm), colnames(Fm), c("d", "m", "u")))
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
  if(is.infinite(n_states)){
    Hm = cbind(matrix(0, nrow = nrow(Hm), ncol = 1), Hm)
  }
  colnames(Hm) = rownames(Fm)  
  if(is.matrix(Hm[, grepl(paste0("e", paste0(1:length(vars), "0"), collapse = "|"), colnames(Hm))])){
    diag(Hm[, grepl(paste0("e", paste0(1:length(vars), "0"), collapse = "|"), colnames(Hm))]) = 1
  }else{
    Hm[, grepl(paste0("e", paste0(1:length(vars), "0"), collapse = "|"), colnames(Hm))] = 1
  }
  if(n_states == 1 | is.infinite(n_states)){
    Hm = array(c(Hm), dim = c(nrow(Hm), ncol(Hm), 1),
               dimnames = list(rownames(Hm), colnames(Hm), c("m")))
  }else if(n_states == 2){
    Hm = array(c(Hm, Hm), dim = c(nrow(Hm), ncol(Hm), n_states),
               dimnames = list(rownames(Hm), colnames(Hm), c("d", "u")))
  }else if(n_states == 3){
    Hm = array(c(Hm, Hm, Hm), dim = c(nrow(Hm), ncol(Hm), n_states),
               dimnames = list(rownames(Hm), colnames(Hm), c("d", "m", "u")))
  }
  
  #Build the transition equation covariance matrix
  #Set the dynamic common factor standard deviation to 1
  Qm = matrix(0, ncol = ncol(Fm), nrow = nrow(Fm))
  rownames(Qm) = rownames(Fm)
  colnames(Qm) = rownames(Qm)
  diag(Qm[c("ct0", paste0("e", 1:nrow(Hm), 0)), c("ct0", paste0("e", 1:nrow(Hm), 0))]) = c(1, sig[names(sig) %in% as.character(1:nrow(Hm))]^2)
  if(is.infinite(n_states)){
    Qm["mt0", "mt0"] = sig["M"]
  }
  if(n_states == 1 | is.infinite(n_states)){
    Qm = array(c(Qm), dim = c(nrow(Qm), ncol(Qm), 1),
               dimnames = list(rownames(Qm), colnames(Qm), c("m")))
  }else if(n_states == 2){
    Qm = array(c(sd["d"]*Qm, sd["u"]*Qm), dim = c(nrow(Qm), ncol(Qm), n_states),
               dimnames = list(rownames(Qm), colnames(Qm), c("d", "u")))
  }else if(n_states == 3){
    Qm = array(c(sd["d"]*Qm, Qm, sd["u"]*Qm), dim = c(nrow(Qm), ncol(Qm), n_states),
               dimnames = list(rownames(Qm), colnames(Qm), c("d", "m", "u")))
  }
  
  #Build the observation eqution covariance matrix
  Rm = matrix(0, ncol = nrow(Hm), nrow = nrow(Hm))#diag(x = sig^2, nrow = length(sig), ncol = length(sig))
  colnames(Rm) = vars
  rownames(Rm) = vars
  if(n_states == 1 | is.infinite(n_states)){
    Rm = array(c(Rm), dim = c(nrow(Rm), ncol(Rm), 1),
               dimnames = list(rownames(Rm), colnames(Rm), c("m")))
  }else if(n_states == 2){
    Rm = array(c(Rm, Rm), dim = c(nrow(Rm), ncol(Rm), n_states),
               dimnames = list(rownames(Rm), colnames(Rm), c("d", "u")))
  }else if(n_states == 3){
    Rm = array(c(Rm, Rm, Rm), dim = c(nrow(Rm), ncol(Rm), n_states),
               dimnames = list(rownames(Rm), colnames(Rm), c("d", "m", "u")))
  }
  
  #Tranistion equation intercept matrix
  #The Markov-switching mean matrix
  Dm = matrix(c(1, rep(0, nrow(Fm) - 1)), nrow = nrow(Fm), ncol = 1)
  if(n_states == 1 | is.infinite(n_states)){
    Dm = array(c(Dm/Inf), dim = c(nrow(Fm), 1, 1),
               dimnames = list(rownames(Fm), NULL, c("m")))
  }else if(n_states == 2){
    Dm = array(c(mu[grepl("d", names(mu))]*Dm, mu[grepl("u", names(mu))]*Dm), dim = c(nrow(Fm), 1, n_states),
               dimnames = list(rownames(Fm), NULL, c("d", "u")))
  }else if(n_states == 3){
    Dm = array(c(mu[grepl("d", names(mu))]*Dm, Dm/Inf, mu[grepl("u", names(mu))]*Dm), dim = c(nrow(Fm), 1, n_states),
               dimnames = list(rownames(Fm), NULL, c("d", "m", "u")))
  }
  
  #Observaton equation intercept matrix
  Am = matrix(0, nrow = nrow(Hm), ncol = 1)
  if(n_states == 1 | is.infinite(n_states)){
    Am = array(c(Am), dim = c(nrow(Am), ncol(Am), 1),
               dimnames = list(rownames(Hm), NULL, c("m")))
  }else if(n_states == 2){
    Am = array(c(Am, Am), dim = c(nrow(Am), ncol(Am), n_states),
               dimnames = list(rownames(Hm), NULL, c("d", "u")))
  }else if(n_states == 3){
    Am = array(c(Am, Am, Am), dim = c(nrow(Am), ncol(Am), n_states),
               dimnames = list(rownames(Hm), NULL, c("d", "m", "u")))
  }
  
  #Steady state probabilities
  Tr_mat = matrix(NA, nrow = ifelse(is.finite(n_states), n_states, 1), ncol = ifelse(is.finite(n_states), n_states, 1))
  rownames(Tr_mat) = colnames(Tr_mat) = dimnames(Am)[[3]]
  for(j in names(pr)){
    Tr_mat[strsplit(j, "")[[1]][2], strsplit(j, "")[[1]][1]] = pr[j]
  }
  for(j in 1:ncol(Tr_mat)){
    Tr_mat[which(is.na(Tr_mat[, j])), j] = 1 - sum(Tr_mat[, j], na.rm = T)
  }
  
  #Initialize the filter for each state
  if(is.null(init)){
    B0 = array(unlist(lapply(colnames(Tr_mat), function(x){MASS::ginv(diag(ncol(Fm[,, x])) - Fm[,, x]) %*% as.matrix(Dm[,, x])})),
             dim = dim(Dm), dimnames = dimnames(Dm))
  }else{
    B0 = init[["B0"]]
  }
  if(is.null(init)){
    VecP0 = array(unlist(lapply(colnames(Tr_mat), function(x){MASS::ginv(diag(prod(dim(Fm[,, x]))) - kronecker(Fm[,, x], Fm[,, x])) %*% matrix(as.vector(Qm[,, x]), ncol = 1)})),
                         dim = dim(Dm), dimnames = list(NULL, NULL, colnames(Tr_mat)))    
    P0 = array(unlist(lapply(colnames(Tr_mat), function(x){matrix(VecP0[, 1, x], ncol = ncol(Fm))})), 
               dim = c(nrow(Fm), ncol(Fm), dim(B0)[3]), dimnames = dimnames(Qm))
    for(i in 1:dim(P0)[3]){
      if(any(diag(P0[,,i]) < 0)){
        diag(P0[,,i]) = 1
      }
    }
  }else{
    P0 = init[["P0"]]
  }
    
  return(list(B0 = B0, P0 = P0, At = Am, Dt = Dm, Ht = Hm, Ft = Fm, Qt = Qm, Rt = Rm, Tr_mat = Tr_mat))
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
#' @param prior "estimate", "uninformative" or vector of named prior parameter guesses: DCF AR coefficients: phi1 and phi2 only; Error MA coefficients: psi_i1 to psi_i2 only for each series i; Error standard deviation: sigma_i only for each series i; Observation coefficient on DCF with first gamma (gamma_i, ..., gamma_n) only 1 index number, not i0 and any more gammas per equation: gamma_i1 to gamma_ik; Markov switching growth rate: mu_d and mu_u; Transition probabilities: p_dd, p_md (or p_mu), p_mm, p_md (or p_mu), p_uu, p_ud (or p_um)
#' @param log.vars Character vector of variables to be logged
#' @param ur.vars Character vector of unit root variables to be differenced
#' @param n_states Number of states to include in the Markov switching model
#' @param detect.lag.length Logical, detect lag length of the dynamic common factor to include in each observation equation using the cross correlation function up to a max of 3
#' @return List of estimation values including coefficients, convergence code, the panel and time ids, the variables in the data, the variable that were logged, and the variables that were differenced
#' @examples
#' ms_dcf_estim(y = DT[, c("date", "y")])
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
ms_dcf_estim = function(y, freq = NULL, panelID = NULL, timeID = NULL, level = 0.01, detect.lag.length = F, 
                        log.vars = NULL, ur.vars = NULL, n_states = 2, ms_var = F,
                        formulas = c("y ~ c + e.l1 + e.l2"), prior = "estimate",
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
  if(!n_states %in% c(1, 2, 3, Inf)){
    stop("n_states must be 1, 2, 3, or Inf.")
  }
  if(is.infinite(n_states) & ms_var == T){
    warning("Infinite state variance not allowed. Model will contain a constant variance structure.")
    ms_var = T
  }
  
  if(is.null(panelID)){
    panelID = "panelid"
    y[, "panelid" := "panel"]
  }
  if(!is.logical(ms_var)){
    stop("ms_var must be T or F.")
  }
  vars = colnames(y)[!colnames(y) %in% c(panelID, timeID)]
  
  #Check for growth variables and log them
  if(is.null(log.vars)){
    gr.test = colMeans(y[, lapply(.SD, function(x){
      d = x - shift(x, type = "lag", n = 1)
      return(t.test(d[!is.na(d)], mu = 0)$p.value)
    }), by = c(panelID), .SDcols = c(vars)][, c(vars), with = F])
    log.vars = names(gr.test)[which(gr.test <= level)]
    if(length(log.vars) > 0){
      y[, c(log.vars) := lapply(.SD, log), .SDcols = c(log.vars)]
    }
  }
  
  #Unit root tests
  if(is.null(ur.vars)){
    ur.vars = lapply(vars, function(x){
      ret = lapply(unique(y[, c(panelID), with = F][[1]]), function(z){
        tseries::adf.test(x = ts(y[eval(parse(text = panelID)) == z &!is.na(eval(parse(text = paste0("`", x, "`")))), c(x), with = F][[1]], freq = freq), alternative = "stationary")$p.value
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
  }
  
  #First difference the relevant series
  yy_d = copy(y)
  if(length(ur.vars) > 0){
    yy_d[, c(ur.vars) := lapply(.SD, function(x){
      x - data.table::shift(x, type = "lag")
    }), by = c(panelID), .SDcols = c(ur.vars)]
  }
  yy_d = yy_d[2:.N, ]
  
  #Standardize the data
  yy_s = copy(yy_d)
  yy_s[, c(vars) := lapply(.SD, function(x){
    (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  }), by = c(panelID), .SDcols = vars]
  
  theta = NULL
  if(all(prior %in% c("estimate", "uninformative"))){
    #Find the starting values
    y[, "C" := Matrix::rowMeans(y[, c(vars), with = F])]
    y[, "C" := imputeTS::na.kalman(C), by = c(panelID)]
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
    names(mid) = unique(y[, c(panelID), with = F][[1]])
    down = lapply(unique(y[, c(panelID), with = F][[1]]), function(x){
      c = y[eval(parse(text = panelID)) == x, ]$c
      return(forecast::Arima(c[which((c - mean(c, na.rm = T))/sd(c, na.rm = T) < -0.44)], order = c(2, 0, 0), include.mean = T))
    })
    names(up) = names(down) = unique(y[, c(panelID), with = F][[1]])
    
    theta = colMeans(do.call("rbind", lapply(names(up), function(x){
      c = y[eval(parse(text = panelID)) == x, ]$c
      coeff = data.table::data.table(cbind(up = up[[x]]$coef, mid = c(mid[[x]]$coef, 0), down = down[[x]]$coef), keep.rownames = T)
      coeff[, "m" := up*length(c[which((c - mean(c, na.rm = T))/sd(c, na.rm = T) > 0.44)])/length(c) + 
              mid*length(c[which((c - mean(c, na.rm = T))/sd(c, na.rm = T) <= 0.44 & (c - mean(c, na.rm = T))/sd(c, na.rm = T) >= -0.44)])/length(c) + 
              down*length(c[which((c - mean(c, na.rm = T))/sd(c, na.rm = T) < -0.44)])/length(c)]
      if(n_states > 1 & is.finite(n_states)){
        return(c(ar = coeff[grepl("ar", rn), ]$m, mu_u = coeff[rn == "intercept", ]$up, mu_d = coeff[rn == "intercept", ]$down))
      }else{
        return(c(ar = coeff[grepl("ar", rn), ]$m))
      }
    })))
    names(theta) = gsub("ar", "phi", names(theta))
    
    #Prior for the Markov-switching means and probabilities
    # y[, "cd" := C - data.table::shift(C, type = "lag", n = 1), by = c(panelID)]
    # y[, "S_d" := ifelse((cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) < ifelse(n_states == 3, -0.44, 0), 1, 0)]
    # y[, "S_u" := ifelse((cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) > ifelse(n_states == 3, 0.44, 0), 1, 0)]
    # y[, "S_dd" := ifelse(S_d == 1 & data.table::shift(S_d, type = "lead", n = 1) == 1, 1, 0)]
    # y[, "S_uu" := ifelse(S_u == 1 & data.table::shift(S_u, type = "lead", n = 1) == 1, 1, 0)]
    # if(n_states == 3){
    #   y[, "S_m" := ifelse((cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) >= -0.44 &
    #                         (cd - mean(cd, na.rm = T))/sd(cd, na.rm = T) <= 0.44, 1, 0)]
    #   y[, "S_mm" := ifelse(S_m == 1 & data.table::shift(S_m, type = "lead", n = 1) == 1, 1, 0)]
    # }
    
    if(n_states > 1 & is.finite(n_states)){
      theta = c(theta, 
              p_dd = 0.95, #sum(y$S_dd, na.rm = T)/sum(y$S_d, na.rm = T),
              p_uu = 0.95)#sum(y$S_uu, na.rm = T)/sum(y$S_u, na.rm = T))
    }
    if(n_states == 3){
      theta = c(theta, 
                p_dm = (1 - 0.95)/2, #(1 - sum(y$S_dd, na.rm = T)/sum(y$S_d, na.rm = T))/2,
                p_md = (1 - 0.95)/2, #(1 - sum(y$S_mm, na.rm = T)/sum(y$S_m, na.rm = T))/2,
                p_mm = 0.95, #sum(y$S_mm, na.rm = T)/sum(y$S_m, na.rm = T),
                p_ud = (1 - 0.95)/2) #(1 - sum(y$S_uu, na.rm = T)/sum(y$S_u, na.rm = T))/2)
    }
    
    #Define the equations by series
    if(length(formulas) < length(vars)){
      formulas = rep(formulas, length(vars))[1:length(vars)]
      names(formulas) = vars
      
      if(detect.lag.length == T){
        c.lags = lapply(unique(yy_s[, c(panelID), with = F][[1]]), function(x){
          c.lag = unlist(lapply(vars, function(v){
            max.lag = max(c(signif(floor(freq/4),1), 1))
            ccf = ccf(x = y[!is.na(c), ]$c, y = yy_s[, c(v), with = F][[1]], na.action = na.pass, lag.max = max.lag, plot = F)
            ccf = data.table(lag = ccf$lag, value = ccf$acf, low = qnorm((1 - 0.95)/2)/sqrt(nrow(y[complete.cases(y), ])/2), up = -qnorm((1 - 0.95)/2)/sqrt(nrow(y[complete.cases(y), ])/2))
            ccf = ccf[lag < 0 & (value > up | value < low), ]
            return(max(abs(ccf$lag)))
          }))
          names(c.lag) = vars
          return(c.lag)
        })
        names(c.lags) = unique(yy_s[, c(panelID), with = F][[1]])
        c.lags = floor(colMeans(do.call("rbind", c.lags)))
        
        for(v in vars){
          if(c.lags[v] > 0){
            seq = signif(round(seq(0, max(c(signif(floor(freq/4),1), 1)), max(c(signif(floor(freq/4),1), 1))/3)), 1)
            formulas[v] = gsub("c +", paste("c +", paste(paste0("c.l", seq[seq > 0]), collapse = " + "), ""), formulas[v])
          }
        }
      }
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
     
      #Create lags
      for(m in 1:3){
        ts[, paste0("c.l", m) := shift(c, type = "lag", n = m)]
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
    if(ms_var == T & n_states > 1 & is.finite(n_states)){
      theta = c(theta, sd_ = c(1.5, 0.5))
      names(theta)[grepl("sd_", names(theta))] = paste0("sd_", c("d", "u"))
    }else if(is.infinite(n_states)){
      theta["sigmaM"] = 1
    }
    suppressWarnings(rm(up, mid, down))
  }
  if(all(prior == "uninformative")){
    theta[grepl("phi|psi", names(theta))] = 0
    theta[names(theta) %in% paste0("gamma", 1:length(vars))] = 1
    theta[grepl("gamma", names(theta)) & !names(theta) %in% paste0("gamma", 1:length(vars))] = 0
    theta[grepl("sig", names(theta))] = 1
    theta[grepl("mu_d", names(theta))] = -1.5
    theta[grepl("mu_u", names(theta))] = 1.5
    if(is.infinite(n_states)){
      theta["sigmaM"] = 1
    }
  }
  if(is.null(theta) & is.numeric(prior)){
    theta = prior
  }
  y = y[, c(panelID, timeID, vars), with = F]
  
  cl = parallel::makeCluster(min(c(length(unique(yy_s[, c(panelID), with = F][[1]])), parallel::detectCores())))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  
  #Get initial values for the filter
  sp = SSmodel_ms(theta, yy_s, n_states, ms_var, panelID, timeID)
  init = foreach::foreach(i = unique(y[, c(panelID), with = F][[1]]), .packages = c("data.table"), .export = c("SSmodel_ms", "kim_filter", "kim_smoother")) %fun% {
    yti = t(yy_s[eval(parse(text = panelID)) == i, colnames(yy_s)[!colnames(yy_s) %in% c(panelID, timeID)], with = F])
    ret = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
    smooth = kim_smoother(ret$B_tlss, ret$B_tts, ret$B_tt, ret$P_tlss, ret$P_tts, ret$Pr_tls, ret$Pr_tts, 
                       sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, Tr_mat = sp$Tr_mat)
    
    ret = list(B0 = array(unlist(lapply(1:dim(ret$B_tts)[3], function(x){matrix(ret$B_tts[,1,x], ncol = 1)})),
                          dim = dim(sp$B0)),
               P0 = array(unlist(lapply(1:length(ret$P_tts), function(x){ret$P_tts[[x]][,,1]})),
                          dim = dim(sp$P0)))
   return(ret)
  }
  names(init) = unique(yy_s[, c(panelID), with = F][[1]])
  #init = NULL
  
  objective = function(par, n_states, ms_var, panelID, timeID, init = NULL, na_locs = NULL){
    par = trans(par)
    yy_s = get("yy_s")
    toret = foreach::foreach(i = unique(yy_s[, c(panelID), with = F][[1]]), .packages = c("data.table"), .export = c("SSmodel_ms", "kim_filter", "kim_smoother")) %fun% {
      sp = SSmodel_ms(par, yy_s, n_states, ms_var, panelID, timeID, init = init[[i]])
      yti = t(yy_s[eval(parse(text = panelID)) == i, colnames(yy_s)[!colnames(yy_s) %in% c(panelID, timeID)], with = F])
      ans = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
      if(!is.null(na_locs)){
        # ans = kim_smoother(ans$B_tlss, ans$B_tts, ans$B_tt, ans$P_tlss, ans$P_tts, ans$Pr_tls, ans$Pr_tts,
        #                    sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, Tr_mat = sp$Tr_mat)
        fc = do.call("cbind", lapply(1:nrow(ans$B_tt), function(x){ans$H_tt[,,x] %*% ans$B_tt[x, ]}))
        rownames(fc) = rownames(yti)
        for(x in rownames(fc)){
          yti[x, na_locs[[i]][[x]]] = fc[x, na_locs[[i]][[x]]]
        }
        ans = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
        ret = list(loglik = ans$loglik, yy_s = data.table(t(yti)))
        ret$yy_s[, "temp" := i]
        colnames(ret$yy_s)[colnames(ret$yy_s) == "temp"] = panelID
      }else{
        ret = list(loglik = ans$loglik)
      }
      return(ret)
    }
    if(!is.null(na_locs)){
      assign("yy_s", do.call("rbind", lapply(toret, function(x){x$yy_s})), .GlobalEnv)
    }
    return(mean(unlist(lapply(toret, function(x){x$loglik}))))
  }
  
  #Estimate the model
  if(any(is.na(yy_s[, c(vars), with = F]))){
    na_locs = lapply(unique(yy_s[, c(panelID), with = F][[1]]), function(x){
      ret = lapply(vars, function(v){
        yy_s[eval(parse(text = panelID)) == x, c(v), with = F][is.na(eval(parse(text = paste0("`", v, "`")))), which = T]
      })
      names(ret) = vars
      return(ret)
    })
    names(na_locs) = unique(yy_s[, c(panelID), with = F][[1]])
  }else{
    na_locs = NULL
  }
  
  theta = init_trans(theta)
  out = tryCatch(maxLik::maxLik(logLik = objective, start = theta, method = optim_methods[1],
                                finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit),
                                na_locs = na_locs, n_states = n_states, ms_var = ms_var, panelID = panelID, timeID = timeID, init = init),
                 error = function(err){
                   tryCatch(maxLik::maxLik(logLik = objective, start = theta, method = optim_methods[min(c(2, length(optim_methods)))],
                                           finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit),
                                           na_locs = na_locs, n_states = n_states, ms_var = ms_var, panelID = panelID, timeID = timeID, init = init),
                            error = function(err){
                              tryCatch(maxLik::maxLik(logLik = objective, start = theta, method = optim_methods[min(c(3, length(optim_methods)))],
                                                      finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit),
                                                      na_locs = na_locs, n_states = n_states, ms_var = ms_var, panelID = panelID, timeID = timeID, init = init),
                                       error = function(err){NULL})
                            })
                 })
  
  #Attempt to get convergence if it failed the first time
  if(!is.null(out)){
    trials = 1
    while(out$code != 0 & trials < maxtrials){
      out2 = tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[1],
                                     finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), 
                                     na_locs = na_locs, n_states = n_states, ms_var = ms_var, panelID = panelID, timeID = timeID, init = init),
                      error = function(err){
                        tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(2, length(optim_methods)))],
                                                finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit),
                                                na_locs = na_locs, n_states = n_states, ms_var = ms_var, panelID = panelID, timeID = timeID, init = init),
                                 error = function(err){
                                   tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(3, length(optim_methods)))],
                                                           finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit),
                                                           na_locs = na_locs, n_states = n_states, ms_var = ms_var, panelID = panelID, timeID = timeID, init = init),
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
  return(list(coef = trans(coef(out)), coef_init = trans(theta), convergence = out$code, loglik = out$maximum, panelID = panelID, timeID = timeID,
              vars = vars, log.vars = log.vars, ur.vars = ur.vars, n_states = n_states, ms_var = ms_var, 
              formulas = formulas, prior = ifelse(is.numeric(prior), "given", prior)))
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
  yy_d = yy_d[2:.N, ]
  
  #Standardize the data
  yy_s = copy(yy_d)
  yy_s[, c(model$vars) := lapply(.SD, function(x){
    (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  }), by = c(model$panelID), .SDcols = c(model$vars)]
  
  if(any(is.na(yy_s[, c(model$vars), with = F]))){
    na_locs = lapply(unique(yy_s[, c(model$panelID), with = F][[1]]), function(x){
      ret = lapply(model$vars, function(v){
        yy_s[eval(parse(text = model$panelID)) == x, c(v), with = F][is.na(eval(parse(text = paste0("`", v, "`")))), which = T]
      })
      names(ret) = model$vars
      return(ret)
    })
    names(na_locs) = unique(yy_s[, c(model$panelID), with = F][[1]])
  }else{
    na_locs = NULL
  }
  
  #Filter the data using the estimated model
  sp = SSmodel_ms(model$coef, yy_s, model$n_states, model$ms_var, model$panelID, model$timeID)
  cl = parallel::makeCluster(min(c(length(unique(yy_s[, c(model$panelID), with = F][[1]])), parallel::detectCores())))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  uc = foreach::foreach(i = unique(yy_s[, c(model$panelID), with = F][[1]]), .packages = c("data.table", "MASS"), .export = c("SSmodel_ms", "kim_filter", "kim_smoother")) %fun% {
    yti = t(yy_s[eval(parse(text = model$panelID)) == i, colnames(yy_s)[!colnames(yy_s) %in% c(model$panelID, model$timeID)], with = F])
    init = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
    init = kim_smoother(init$B_tlss, init$B_tts, init$B_tt, init$P_tlss, init$P_tts, init$Pr_tls, init$Pr_tts, 
                        sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat)
    init = list(B0 = array(unlist(lapply(1:dim(init$B_tts)[3], function(x){matrix(init$B_tts[,1,x], ncol = 1)})),
                           dim = dim(sp$B0)),
                P0 = array(unlist(lapply(1:length(init$P_tts), function(x){init$P_tts[[x]][,,1]})),
                           dim = dim(sp$P0)))
    sp2 = SSmodel_ms(model$coef, yy_s[eval(parse(text = model$panelID)) == i, ], model$n_states, model$ms_var, model$panelID, model$timeID, init = init)
    ans = kim_filter(sp2$B0, sp2$P0, sp2$At, sp2$Dt, sp2$Ft, sp2$Ht, sp2$Qt, sp2$Rt, sp2$Tr_mat, yti)
    if(!is.null(na_locs)){
      fc = do.call("cbind", lapply(1:nrow(ans$B_tt), function(x){ans$H_tt[,,x] %*% ans$B_tt[x, ]}))
      rownames(fc) = rownames(yti)
      for(x in rownames(fc)){
        yti[x, na_locs[[i]][[x]]] = fc[x, na_locs[[i]][[x]]]
      }
      ans = kim_filter(sp2$B0, sp2$P0, sp2$At, sp2$Dt, sp2$Ft, sp2$Ht, sp2$Qt, sp2$Rt, sp2$Tr_mat, yti)
    }
    
    #Get the steady state Kalman gain approximation
    Ht = Reduce("+", lapply(1:dim(ans$H_tt)[3], function(x){ans$H_tt[,, x]}))/dim(ans$H_tt)[3]
    Ft = Reduce("+", lapply(1:dim(ans$F_tt)[3], function(x){ans$F_tt[,, x]}))/dim(ans$F_tt)[3]
    rownames(Ht) = rownames(yti)
    colnames(Ht) = colnames(Ft) = colnames(sp2$Ft[,, 1])
    rownames(Ft) = rownames(sp2$Ft[,, 1])
    K = matrix(ans$K[,, dim(ans$K)[3]], nrow = dim(ans$K)[1], ncol = dim(ans$K)[2])
    W = MASS::ginv(diag(nrow(K)) - (diag(nrow(K)) - K %*% Ht) %*% Ft) %*% K
    means = unlist(yy_d[eval(parse(text = model$panelID)) == i, lapply(.SD, mean, na.rm = T), .SDcols = c(model$vars)])
    sds = unlist(yy_d[eval(parse(text = model$panelID)) == i, lapply(.SD, sd, na.rm = T), .SDcols = c(model$vars)])
    
    #W(1) is the first row of W
    d = W[1, ] %*% means
    
    #Get the intercept terms
    D = means - Ht[, grepl("ct", colnames(Ht))] %*% matrix(rep(d, ncol(Ht[, grepl("ct", colnames(Ht))])))
    
    #Initialize first element of the dynamic common factor
    Y1 = as.matrix(yti[, 1])
    nonna_idx = which(!is.na(Y1))
    initC = function(par){
      return(sum((as.matrix(Y1[nonna_idx, ]) - as.matrix(D[nonna_idx, ]) - as.matrix(Ht[nonna_idx, grepl("ct", colnames(Ht))]) %*% matrix(par))^2))
    }
    C10 = optim(par = rep(mean(Y1)/mean(model$coef[grepl("gamma", names(model$coef))]), length(colnames(Ft)[grepl("ct", colnames(Ft))])), 
                fn = initC, method = "BFGS", control = list(trace = F))$par[1]
    dcf_loc = which(rownames(Ft) == "ct0")
    
    toret = list()
    for(k in c("filter", "smooth")){
      Ctt = rep(C10, ncol(yti) + 1)
      if(k == "smooth"){
        ans = kim_smoother(ans$B_tlss, ans$B_tts, ans$B_tt, ans$P_tlss, ans$P_tts, ans$Pr_tls, ans$Pr_tts, 
                           sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat)
      }
      
      #Build the rest of the DCF series
      #Rescale the first differenced DCF series to match that of the actual series
      ctt = c(0, (t(ans$B_tt[, dcf_loc]))*mean(sds)/sd(ans$B_tt[, dcf_loc]))
      if(is.finite(model$n_states)){
        mutt = c(0, t(ans$D_tt[dcf_loc,,])*mean(sds)/sd(ans$B_tt[, dcf_loc]))
      }else{
        m_loc = which(rownames(Ft) == "mt0")
        mutt = c(0, t(ans$B_tt[, m_loc])*mean(sds)/sd(ans$B_tt[, m_loc]))
      }
      for(j in 2:length(Ctt)){
        #First element of dCtt is C_22 - C_11 = dCtt_2
        Ctt[j] = ctt[j] + Ctt[j - 1] + c(d)
      }
      Ctt = data.table::data.table(panelID = i, date = y[eval(parse(text = model$panelID)) == i, c(model$timeID), with = F][[1]], DCF = Ctt, d.DCF = ctt, Mu = mutt)
      colnames(Ctt) = c(model$panelID, model$timeID, "DCF", "d.DCF", "Mu")
      
      #Build the probability series
      prob = data.table::data.table(panelID = i, date = yy_s[eval(parse(text = model$panelID)) == i, c(model$timeID), with = F][[1]], data.table::data.table(ans$Pr_tts))
      colnames(prob) = c(model$panelID, model$timeID, paste0("Pr_", dimnames(sp2$Dt)[[3]]))
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
        toplot1 = merge(toplot1, toplot1[, .(value_start = mean(value, na.rm = T)), by = c(model$panelID, "variable")], by = c("variable"), all.x = T, all.y = F)
        toplot1[, "value2" := value - (value_start - mean(toplot1[eval(parse(text = model$timeID)) == min(eval(parse(text = model$timeID))), ]$value, na.rm = T))*9/10, by = "variable"]
        g1 = ggplot2::ggplot(toplot1) + 
          ggplot2::ggtitle(paste(ifelse(i == "panel" & model$panelID == "panelid", "", i), "Data Series"), subtitle = "Levels") + 
          ggplot2::scale_y_continuous(name = "Levels (Intercept Adjusted)") + 
          ggplot2::scale_x_date(name = "") + 
          ggplot2::geom_line(ggplot2::aes(x = eval(parse(text = model$timeID)), y = value2, group = variable, color = variable)) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + ggplot2::guides(color = ggplot2::guide_legend(title = NULL))
        
        toplot2 = melt(yy_s[eval(parse(text = model$panelID)) == i, ], id.vars = c(model$panelID, model$timeID))
        g2 = ggplot2::ggplot(toplot2) + 
          ggplot2::ggtitle(paste(ifelse(i == "panel" & model$panelID == "panelid", "", i), "Data Series"), subtitle = "Differenced & Standardized") + 
          ggplot2::scale_y_continuous(name = "Differences") + 
          ggplot2::scale_x_date(name = "") + 
          ggplot2::geom_hline(yintercept = 0, color = "black") + 
          ggplot2::geom_line(ggplot2::aes(x = eval(parse(text = model$timeID)), y = value, group = variable, color = variable)) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + ggplot2::guides(color = ggplot2::guide_legend(title = NULL))
        
        toplot3 = melt(uc[[i]][[j]], id.vars = c(model$panelID, model$timeID))
        d_range1 = range(toplot3[variable == "DCF", ]$value, na.rm = T)
        p_range1 = range(toplot3[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], ]$value, na.rm = T)
        toplot3[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], "value" := (value - p_range1[1])/diff(p_range1) * diff(d_range1) + d_range1[1], by = "variable"]
        g3 = ggplot2::ggplot() +  
          ggplot2::ggtitle(paste(ifelse(i == "panel" & model$panelID == "panelid", "", i), "Dynamic Common Factor"), subtitle = "Levels") + 
          ggplot2::scale_x_date(name = "") +
          ggplot2::geom_hline(yintercept = 0, color = "grey") + 
          ggplot2::geom_line(data = toplot3[variable == "DCF", ], 
                             ggplot2::aes(x = eval(parse(text = model$timeID)), y = value, group = variable, color = variable)) + 
          ggplot2::scale_color_manual(values = c("black")) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + 
          ggplot2::guides(color = ggplot2::guide_legend(title = NULL), fill = ggplot2::guide_legend(title = NULL))
        if(model$n_states == 1 | is.infinite(model$n_states)){
          g3 = g3 + ggplot2::scale_y_continuous(name = "Dynamic Common Factor", limits = range(toplot3[variable == "DCF", ]$value, na.rm = T))
        }else{
          g3 = g3 + 
            ggplot2::scale_y_continuous(name = "Dynamic Common Factor", limits = range(toplot3[variable == "DCF", ]$value, na.rm = T), 
                                        sec.axis = ggplot2::sec_axis(name = "Probability of State (%)", ~((. - d_range1[1])/diff(d_range1) * diff(p_range1) + p_range1[1]) * 100)) + 
            ggplot2::geom_ribbon(data = toplot3[variable %in% colnames(uc[[i]][[j]])[grepl(ifelse(model$n_states == 2, "Pr_d", "Pr_d|Pr_u"), colnames(uc[[i]][[j]]), ignore.case = T)], ], 
                                 ggplot2::aes(x = eval(parse(text = model$timeID)), ymin = d_range1[1], ymax = value, group = variable, fill = variable), alpha = 0.5) + 
            ggplot2::scale_fill_manual(values = c("red", "green"))
        }
        
        toplot4 = melt(uc[[i]][[j]], id.vars = c(model$panelID, model$timeID))
        d_range2 = range(toplot4[variable %in% c("d.DCF", "Mu"), ]$value, na.rm = T)
        p_range2 = range(toplot4[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], ]$value, na.rm = T)
        toplot4[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], "value" := (value - p_range2[1])/diff(p_range2) * diff(d_range2) + d_range2[1], by = "variable"]
        g4 = ggplot2::ggplot() +  
          ggplot2::ggtitle(paste(ifelse(i == "panel" & model$panelID == "panelid", "", i), "Dynamic Common Factor"), subtitle = "Differenced") + 
          ggplot2::scale_x_date(name = "") +
          ggplot2::geom_hline(yintercept = 0, color = "grey") + 
          ggplot2::geom_line(data = toplot4[variable %in% c("d.DCF", "Mu"), ], 
                             ggplot2::aes(x = eval(parse(text = model$timeID)), y = value, group = variable, color = variable)) + 
          ggplot2::scale_color_manual(values = c("black", "blue")) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + 
          ggplot2::guides(color = ggplot2::guide_legend(title = NULL), fill = ggplot2::guide_legend(title = NULL))
        if(model$n_states == 1 | is.infinite(model$n_states)){
          g4 = g4 + ggplot2::scale_y_continuous(name = "Dynamic Common Factor", limits = range(toplot4[variable %in% c("d.DCF", "Mu"), ]$value, na.rm = T)) 
        }else{
          g4 = g4 + ggplot2::scale_y_continuous(name = "Dynamic Common Factor", limits = range(toplot4[variable %in% c("d.DCF", "Mu"), ]$value, na.rm = T), 
                                                sec.axis = ggplot2::sec_axis(name = "Probability of State (%)", ~((. - d_range2[1])/diff(d_range2) * diff(p_range2) + p_range2[1]) * 100)) + 
            ggplot2::geom_ribbon(data = toplot4[variable %in% colnames(uc[[i]][[j]])[grepl(ifelse(model$n_states == 2, "Pr_d", "Pr_d|Pr_u"), colnames(uc[[i]][[j]]), ignore.case = T)], ], 
                                 ggplot2::aes(x = eval(parse(text = model$timeID)), ymin = d_range2[1], ymax = value, group = variable, fill = variable), alpha = 0.5) + 
            ggplot2::scale_fill_manual(values = c("red", "green"))
        }
        
        gridExtra::grid.arrange(g1, g2, g3, g4, layout_matrix = matrix(c(1, 3, 2, 4), nrow = 2),
                                top = grid::textGrob(j, gp = grid::gpar(fontsize = 20, font = 3)))
        Sys.sleep(0.1)
      }
    }
  }
  
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
