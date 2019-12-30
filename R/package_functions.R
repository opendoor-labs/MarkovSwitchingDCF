#' State space model for markov switching mean
#' Creates a state space model in list form
#' yt = At + Ht * Bt + e_t
#' Bt = D_t + Ft * B_tl + u_t
#'
#' @param par Vector of named parameter values
#' @param yt Multivariate time series of data values
#' @param panelID Column name that identifies the cross section of the data
#' @param timeID Column name that identifies the date
#' @param ms_var, Logical, T for Markow switching variance, default is F
#' @param n_states Number of states to include in the Markov switching model
#' @return List of space space matrices
#' @examples
#' SSmodel_ms(par = theta, yt = yt, n_states = 2, panelID = "panel", timeID = "date")
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
SSmodel_ms = function(par, yt, n_states, ms_var = F, panelID = NULL, timeID = NULL){
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
  if(n_states == 1 | is.infinite(n_states)){
    states = "m"
  }else{
    states = sort(unique(substr(names(pr), 1, 1)))
  }
  
  #Steady state probabilities
  Tr_mat = matrix(NA, nrow = ifelse(is.finite(n_states), n_states, 1), ncol = ifelse(is.finite(n_states), n_states, 1))
  rownames(Tr_mat) = colnames(Tr_mat) = states
  for (j in names(pr)){
    Tr_mat[strsplit(j, "")[[1]][2], strsplit(j, "")[[1]][1]] = pr[j]
  }
  for(j in 1:ncol(Tr_mat)){
    Tr_mat[which(is.na(Tr_mat[, j])), j] = 1 - sum(Tr_mat[, j], na.rm = T)
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
  Fm = array(Fm, dim = c(nrow(Fm), ncol(Fm), ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(rownames(Fm), colnames(Fm), states))
  
  #Build the observation equation matrix
  Hm = matrix(gamma[as.character(1:length(vars))], ncol = 1, nrow = length(vars))
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
  if(is.matrix(Hm[, grepl(paste0("e", paste0(1:length(vars), "0"), collapse = "|"), colnames(Hm))])) {
    diag(Hm[, grepl(paste0("e", paste0(1:length(vars), "0"), collapse = "|"), colnames(Hm))]) = 1
  }else{
    Hm[, grepl(paste0("e", paste0(1:length(vars), "0"), collapse = "|"), colnames(Hm))] = 1
  }
  Hm = array(Hm, dim = c(nrow(Hm), ncol(Hm), ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(rownames(Hm), colnames(Hm), states))
  
  #Build the transition equation covariance matrix
  #Set the dynamic common factor standard deviation to 1
  Qm = matrix(0, ncol = ncol(Fm), nrow = nrow(Fm))
  rownames(Qm) = rownames(Fm)
  colnames(Qm) = rownames(Qm)
  diag(Qm[c("ct0", paste0("e", 1:nrow(Hm), 0)), c("ct0", paste0("e", 1:nrow(Hm), 0))]) = c(1, sig[names(sig) %in% as.character(1:nrow(Hm))]^2)
  if (is.infinite(n_states)) {
    Qm["mt0", "mt0"] = sig["M"]
  }
  Qm = array(Qm, dim = c(nrow(Qm), ncol(Qm), ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(rownames(Qm), colnames(Qm), states))
  if(ms_var == T){
    for(i in names(sd)){
      Qm[, , i] = sd[i] * Qm[, , i]
    }
  }
  
  #Build the observation eqution covariance matrix
  Rm = matrix(0, ncol = nrow(Hm), nrow = nrow(Hm))
  colnames(Rm) = vars
  rownames(Rm) = vars
  Rm = array(Rm, dim = c(nrow(Rm), ncol(Rm), ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(rownames(Rm), colnames(Rm), states))
  
  #Tranistion equation intercept matrix
  #The Markov-switching mean matrix
  Dm = matrix(0, nrow = nrow(Fm), ncol = 1)
  Dm = array(Dm, dim = c(nrow(Dm), 1, ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(rownames(Fm), NULL, states))
  if(length(mu) > 0){
    for(i in names(mu)){
      Dm["ct0", , i] = mu[i]
    }
  }
  
  #Observaton equation intercept matrix
  Am = matrix(0, nrow = nrow(Hm), ncol = 1)
  Am = array(Am, dim = c(nrow(Am), ncol(Am), ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(vars, NULL, states))
  
  #Initialize the filter for each state
  B0 = matrix(0, nrow(Fm), 1)
  B0 = array(B0, dim = c(nrow(B0), ncol(B0), ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(rownames(Fm), NULL, states))
  P0 = diag(nrow(Fm))
  P0 = array(P0, dim = c(nrow(P0), ncol(P0), ifelse(is.infinite(n_states), 1, n_states)), dimnames = list(rownames(B0), colnames(B0), states))
  
  return(list(B0 = B0, P0 = P0, At = Am, Dt = Dm, Ht = Hm, Ft = Fm, Qt = Qm, Rt = Rm, Tr_mat = Tr_mat))
}

#' Transform data for estimation
#'
#' @param y Multivariate time series of data values. May also be a data frame containing a date column
#' @param model Estimated model from ms_dcf_estim
#' @param freq Seasonality of the data
#' @param panelID Column name that identifies the cross section of the data
#' @param timeID Column name that identifies the date
#' @param level Significance level of statistical tests (0.01, 0.05, 0.1)
#' @param detect.formula Logical, detect lag length of the dynamic common factor to include in each observation equation using the cross correlation function up to a max of 3
#' @param detect.growth Logical, detect which variables to log
#' @param detect.diff Logical, detect which variables to difference
#' @param log.vars Character vector of variables to be logged
#' @param diff.vars Character vector of unit root variables to be differenced
#' @param diff.lag Integer, number of lags to use for differencing
#' @return list containing differenced data (yy_d) and standardized data (yy_s)
#' @examples
#' data_trans(y = y, panelID = "panel", timeID = "date")
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
data_trans = function(y, model = NULL, log.vars = NULL, diff.vars = NULL, freq = NULL,
                      detect.diff = F, detect.growth = F, level = 0.01, diff.lag = 1,
                      panelID = NULL, timeID = NULL){
  y = copy(y)
  if(!is.null(model)){
    vars = model$vars
    log.vars = model$log.vars
    diff.vars = model$diff.vars
    freq = model$freq
    detect.diff = model$detect.diff
    detect.growth = model$detect.growth
    level = model$level
    panelID = model$panelID
    timeID = model$timeID
  }else{
    vars = colnames(y)[!colnames(y) %in% c(panelID, timeID)]
  }

  #Check for growth variables and log them
  if(detect.growth == T){
    gr.test = colMeans(y[, lapply(.SD, function(x){
      d = x - shift(x, type = "lag", n = diff.lag)
      return(t.test(d[!is.na(d)], mu = 0)$p.value)
    }), by = c(panelID), .SDcols = c(vars)][, c(vars), with = F])
    log.vars = names(gr.test)[which(gr.test <= level)]
  }
  if(length(log.vars) > 0){
    if(any(log.vars %in% colnames(y))){
      y[, c(log.vars[log.vars %in% colnames(y)]) := lapply(.SD, log),
        .SDcols = c(log.vars[log.vars %in% colnames(y)])]
    }
  }

  #Check for nonstationary variables to difference
  if(detect.diff == T){
    diff.vars = lapply(vars, function(x){
      ret = lapply(unique(y[, c(panelID), with = F][[1]]), function(z){
        tseries::adf.test(x = ts(y[eval(parse(text = panelID)) == z &!is.na(eval(parse(text = paste0("`", x, "`")))), c(x), with = F][[1]], freq = freq), alternative = "stationary")$p.value
      })
      names(ret) = unique(y[, c(panelID), with = F][[1]])
      return(ret)
    })
    names(diff.vars) = vars
    diff.vars = lapply(names(diff.vars), function(x){
      mean(unlist(diff.vars[[x]]))
    })
    names(diff.vars) = vars
    diff.vars = names(diff.vars)[diff.vars >= level]
  }

  #Difference the relevant variables
  yy_d = copy(y)
  if(length(diff.vars) > 0){
    if(any(diff.vars %in% colnames(yy_d))){
      yy_d[, c(diff.vars[diff.vars %in% colnames(yy_d)]) := lapply(.SD, function(x){
        x - data.table::shift(x, type = "lag", n = diff.lag)
      }), by = c(panelID), .SDcols = c(diff.vars[diff.vars %in% colnames(yy_d)])]
    }
  }
  yy_d = yy_d[(diff.lag + 1):.N, ]

  #Standardize the data
  yy_s = copy(yy_d)
  yy_s[, c(vars[vars %in% colnames(yy_s)]) := lapply(.SD, function(x){
    (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  }), by = c(panelID), .SDcols = c(vars[vars %in% colnames(yy_s)])]
  
  yy_d = yy_d[, c(panelID, timeID, vars), with = F]
  yy_s = yy_s[, c(panelID, timeID, vars), with = F]
  setcolorder(yy_d, c(panelID, timeID, vars))
  setcolorder(yy_s, c(panelID, timeID, vars))
  return(list(yy_d = yy_d, yy_s = yy_s))
}

#' Set the priors for estimation
#' 
#' @param yy_s Multivariate time series of standardized data values from data_trans
#' @param panelID Column name that identifies the cross section of the data
#' @param timeID Column name that identifies the date
#' @param formulas R formula describing the relationship between each data series, the unobserved dynamic common factor, and the erorr structure
#' @param prior "estimate", "uninformative" or vector of named prior parameter guesses: DCF AR coefficients: phi1 and phi2 only; Error MA coefficients: psi_i1 to psi_i2 only for each series i; Error standard deviation: sigma_i only for each series i; Observation coefficient on DCF with first gamma (gamma_i, ..., gamma_n) only 1 index number, not i0 and any more gammas per equation: gamma_i1 to gamma_ik; Markov switching growth rate: mu_d and mu_u; Transition probabilities: p_dd, p_md (or p_mu), p_mm, p_md (or p_mu), p_uu, p_ud (or p_um)
#' @param n_states Number of states to include in the Markov switching model
#' @param ms_var, Logical, T for Markow switching variance, default is F
#' @param level Significance level of statistical tests (0.01, 0.05, 0.1)
#' @param detect.formula Logical, detect lag length of the dynamic common factor to include in each observation equation using the cross correlation function up to a max of 3
#' @return vector of initial coefficient values
#' @examples
#' set_priors(yy_s = yy_s, prior = prior, panelID = "panel", timeID = "date")
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
set_priors = function(yy_s, prior, panelID, timeID, n_states = 2, ms_var = F, detect.formula = F,
                      level = 0.01, formulas = c("y ~ c + e.l1 + e.l2")){
  yy_s = copy(yy_s)
  vars = colnames(yy_s)[!colnames(yy_s) %in% c(panelID, timeID)]
  
  #Set the priors
  theta = NULL
  if(!all(is.numeric(prior))){
    yy_s[complete.cases(yy_s), "c" := psych::fa(yy_s[, c(vars), with = F], nfactors = 1, rotate = "none", fm = "pa")$scores]
    yy_s[, "c" := imputeTS::na_kalman(c), by = c(panelID)]
    
    #Prior for the AR coefficients (phi)
    up = lapply(unique(yy_s[, c(panelID), with = F][[1]]), function(x){
      c = yy_s[eval(parse(text = panelID)) == x, ]$c
      return(forecast::Arima(c[c > 0.44], order = c(2, 0, 0), include.mean = T))
    })
    mid = lapply(unique(yy_s[, c(panelID), with = F][[1]]), function(x){
      c = yy_s[eval(parse(text = panelID)) == x, ]$c
      return(forecast::Arima(c[c <= 0.44 & c >= -0.44], order = c(2, 0, 0), include.mean = F))
    })
    names(mid) = unique(yy_s[, c(panelID), with = F][[1]])
    down = lapply(unique(yy_s[, c(panelID), with = F][[1]]), function(x){
      c = yy_s[eval(parse(text = panelID)) == x, ]$c
      return(forecast::Arima(c[c < -0.44], order = c(2, 0, 0), include.mean = T))
    })
    names(up) = names(down) = unique(yy_s[, c(panelID), with = F][[1]])
    
    theta = colMeans(do.call("rbind", lapply(names(up), function(x){
      c = yy_s[eval(parse(text = panelID)) == x, ]$c
      coeff = data.table::data.table(cbind(up = up[[x]]$coef, mid = c(mid[[x]]$coef, 0), down = down[[x]]$coef), keep.rownames = T)
      coeff[, "m" := up*length(c[c > 0.44])/length(c) + 
              mid*length(c[c <= 0.44 & c >= -0.44])/length(c) + 
              down*length(c[c < -0.44])/length(c)]
      if(n_states > 1 & is.finite(n_states)){
        return(c(ar = coeff[grepl("ar", rn), ]$m, mu_u = coeff[rn == "intercept", ]$up, mu_d = coeff[rn == "intercept", ]$down))
      }else{
        return(c(ar = coeff[grepl("ar", rn), ]$m))
      }
    })))
    names(theta) = gsub("ar", "phi", names(theta))
    
    #Prior for the Markov-switching means and probabilities
    yy_s[, "S_d" := ifelse(c < ifelse(n_states == 3, -0.44, 0), 1, 0)]
    yy_s[, "S_u" := ifelse(c > ifelse(n_states == 3, 0.44, 0), 1, 0)]
    yy_s[, "S_dd" := ifelse(S_d == 1 & data.table::shift(S_d, type = "lead", n = 1) == 1, 1, 0)]
    yy_s[, "S_uu" := ifelse(S_u == 1 & data.table::shift(S_u, type = "lead", n = 1) == 1, 1, 0)]
    if(n_states == 3){
      yy_s[, "S_m" := ifelse(c >= -0.44 & c <= 0.44, 1, 0)]
      yy_s[, "S_mm" := ifelse(S_m == 1 & data.table::shift(S_m, type = "lead", n = 1) == 1, 1, 0)]
    }
    
    if(n_states > 1 & is.finite(n_states)){
      theta = c(theta, 
                p_dd = max(0.5, min(c(0.95, sum(yy_s$S_dd, na.rm = T)/sum(yy_s$S_d, na.rm = T)))),
                p_uu = max(0.5, min(c(0.95, sum(yy_s$S_uu, na.rm = T)/sum(yy_s$S_u, na.rm = T)))))
    }
    if(n_states == 3){
      theta = c(theta, 
                p_dm = (1 - sum(yy_s$S_dd, na.rm = T)/sum(yy_s$S_d, na.rm = T))/2,
                p_md = (1 - sum(yy_s$S_mm, na.rm = T)/sum(yy_s$S_m, na.rm = T))/2,
                p_mm = sum(yy_s$S_mm, na.rm = T)/sum(yy_s$S_m, na.rm = T),
                p_ud = (1 - sum(yy_s$S_uu, na.rm = T)/sum(yy_s$S_u, na.rm = T))/2)
    }
    
    #Define the equations by series
    if(length(formulas) < length(vars)){
      formulas = rep(formulas, length(vars))[1:length(vars)]
    }
    names(formulas) = vars
    
    if(detect.formula == T){
      c.lags = lapply(unique(yy_s[, c(panelID), with = F][[1]]), function(x){
        c.lag = unlist(lapply(vars, function(v){
          
          #max lags is the based on the number of paramters to be estimated per equation
          max.lag = max(c(floor(nrow(yy_s)/30 - (length(which(gregexpr("e\\.", formulas[v])[[1]] > 0)) + 1)), 1))
          ccf = TSA::prewhiten(x = yy_s[!is.na(c), ]$c, y =  yy_s[, c(v), with = F][[1]], plot = F)$ccf
          
          df = nrow(yy_s[complete.cases(yy_s), ]) - 2
          critical.t = qt(level/2, df, lower.tail = F)
          critical.r = sqrt((critical.t^2)/((critical.t^2) + df))
          
          ccf = data.table(lag = ccf$lag, value = ccf$acf, low = -abs(critical.r), up = abs(critical.r))
          ccf = ccf[lag < 0 & (value > up | value < low), ]
          if(nrow(ccf) == 0){
            return(0)
          }else{
            return(min(c(max.lag, max(abs(ccf$lag)))))
          }
        }))
        names(c.lag) = vars
        return(c.lag)
      })
      names(c.lags) = unique(yy_s[, c(panelID), with = F][[1]])
      c.lags = floor(colMeans(do.call("rbind", c.lags)))
      c.lags = sapply(names(c.lags), function(x){max(c.lags)})
      
      for(v in vars){
        if(c.lags[v] > 0){
          seq = seq(0, c.lags[v], 1)
          formulas[v] = gsub("c +", paste("c +", paste(paste0("c.l", seq[seq > 0]), collapse = " + "), ""), formulas[v])
        }
      }
    }
    
    theta = c(theta, unlist(lapply(vars, function(z){
      #Average the series
      ytemp = Matrix::rowMeans(yy_s[, c(z), with = F])
      ts = data.table::data.table(y = yy_s[!is.na(eval(parse(text = z))) & !is.na(c), c(z, "c"), with = F])
      colnames(ts) = c("y", "c")
      
      #Create lags
      for(m in 1:length(which(gregexpr("c\\.", formulas[z])[[1]] > 0))){
        ts[, paste0("c.l", m) := shift(c, type = "lag", n = m)]
      }
      ts = ts[complete.cases(ts), ]
      
      #Step 1: Regress the series on the DCF and any lags
      lm.vars = trimws(strsplit(formulas[z], "\\~|\\+")[[1]])
      lm.vars = lm.vars[lm.vars %in% colnames(ts)]
      fit = lm(as.formula(paste("y ~", paste(lm.vars[2:length(lm.vars)], collapse = " + "), "-1")), data = ts)
      
      #Get the errors
      ts = cbind(ts, e = fit$residuals)
      ts[, "e.l1" := shift(e, type = "lag", n = 1)]
      ts[, "e.l2" := shift(e, type = "lag", n = 2)]
      ts[, "fit" := fit$fitted.values]
      ts = ts[complete.cases(ts), ]
      
      #Step 2: Regress the series on the DCF and the errors from step 1
      fit2 = lm("y ~ offset(fit) + e.l1 + e.l2 - 1", data = ts[, colnames(ts)[colnames(ts) != "e"], with = F])
      
      #Get the estimated coefficients
      idx = which(vars == z)
      coeff = c(fit$coefficients, fit2$coefficients)
      names(coeff)[names(coeff) == "c"] = paste0("gamma", idx)
      names(coeff)[grepl("c\\.l", names(coeff))] = paste0(paste0("gamma", idx), 1:(length(names(coeff)[grepl("c\\.l", names(coeff))])))
      names(coeff)[grepl("e\\.l", names(coeff))] = paste0(paste0("psi", idx), 1:length(names(coeff)[grepl("e\\.l", names(coeff))]))
      
      #Get the sigma parameter
      coeff = c(coeff, summary(fit)$sigma)
      names(coeff)[length(coeff)] = paste0("sigma", idx)
      return(coeff)
    })))
    if(!is.null(prior)){
      if(prior == "uninformative"){
        theta[names(theta) %in% paste0("gamma", 1:length(vars))] = 1
        theta[!names(theta) %in% paste0("gamma", 1:length(vars)) & grepl("gamma", names(theta))] = 0
        theta[grepl("psi", names(theta))] = 0
        theta[grepl("sigma", names(theta))] = 1
      }
    }
    if(ms_var == T & n_states > 1 & is.finite(n_states)){
      theta = c(theta, sd_ = c(1.5, 0.5))
      names(theta)[grepl("sd_", names(theta))] = paste0("sd_", c("d", "u"))
    }else if(is.infinite(n_states)){
      theta["sigmaM"] = 1
    }
    suppressWarnings(rm(up, mid, down))
  }else{
    theta = prior
  }
  return(theta)
}

#' Set the constraints for estimation
#' 
#' @param yy_s Multivariate time series of standardized data values from data_trans
#' @param panelID Column name that identifies the cross section of the data
#' @param timeID Column name that identifies the date
#' @param theta Named ector of initial coefficient values 
#' @param n_states Number of states to include in the Markov switching model
#' @return vector of initial coefficient values
#' @examples
#' set_constraintsy(yy_s = yy_s, theta = theta, n_states = n_states, panelID = panelID, timeID = timeID)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
set_constraints = function(yy_s, theta, n_states, panelID, timeID){
  yy_s = copy(yy_s)
  
  #Set the constraints
  nseries = ncol(yy_s[, colnames(yy_s)[!colnames(yy_s) %in% c(panelID, timeID)], with = F])
  nrow = 2 + 2*nseries + 
    length(theta[grepl("sigma", names(theta))]) +
    ifelse(is.finite(n_states) & n_states > 1, 2*length(theta[grepl("p_", names(theta))]), 0) + 
    ifelse(is.finite(n_states) & n_states > 1, length(theta[grepl("mu_", names(theta))]), 0) +
    ifelse(is.finite(n_states) & n_states == 2, 4, ifelse(is.finite(n_states) & n_states == 3, 6, 0))
  ineqA = matrix(0, nrow = nrow, ncol = length(theta)) 
  ineqB = matrix(0, nrow = nrow(ineqA), ncol = 1)
  colnames(ineqA) = names(theta)
  #-1 < sum(phi) < 1
  ineqA[1:2, c("phi1", "phi2")] = c(1, -1)
  ineqB[1:2, ] = 1
  #-1 < sum(psi_i) < 1 & sigma > 0
  rn = 3
  for(i in 1:nseries){
    ineqA[c(rn, rn + 1), colnames(ineqA)[grepl(paste0("psi", i), colnames(ineqA))]] = c(1, -1)
    ineqA[rn + 2, colnames(ineqA)[grepl(paste0("sigma", i), colnames(ineqA))]] = 1
    ineqB[c(rn, rn + 1), ] = 1
    rn = rn + 3
  }
  if(is.finite(n_states) & n_states > 1){
    #0 < p < 1
    for(i in 1:length(theta[grepl("p_", names(theta))])){
      ineqA[c(rn, rn + 1), names(theta)[grepl("p_", names(theta))][i]] = c(1, -1)
      ineqB[rn + 1, ] = 1
      rn = rn + 2
    }
    #0 < sum(p) < 1
    ineqA[c(rn, rn + 1), colnames(ineqA)[grepl("p_d", colnames(ineqA))]] = c(1, -1)
    ineqB[c(rn + 1), ] = 1
    rn = rn + 2
    ineqA[c(rn, rn + 1), colnames(ineqA)[grepl("p_u", colnames(ineqA))]] = c(1, -1)
    ineqB[c(rn + 1), ] = 1
    rn = rn + 2
    if(n_states == 3){
      ineqA[c(rn, rn + 1), colnames(ineqA)[grepl("p_m", colnames(ineqA))]] = c(1, -1)
      ineqB[c(rn + 1), ] = 1
      rn = rn + 2
    }
    
    #mu_d < 0 & 0 < mu_u
    ineqA[rn, "mu_d"] = -1
    ineqA[rn + 1, "mu_u"] = 1
    rn = rn + 2
  }
  if("sigmaM" %in% names(theta)){
    ineqA[rn, which(colnames(ineqA) == "sigmaM")] = 1
    rn = rn + 1
  }
  if("sigmaV" %in% names(theta)){
    ineqA[rn , which(colnames(ineqA) == "sigmaV")] = 1
  }
  
  #Make sure initial guesses are within the constraints
  if(any(ineqA %*% as.matrix(theta) + ineqB < 0)){
    wh = which(any(ineqA %*% as.matrix(theta) + ineqB < 0))
    for(j in wh){
      wh_vars = names(which(ineqA[j, ] != 0))
      for(k in wh_vars){
        if(grepl("phi|psi|gmma", k)){
          theta[k] = 0
        }else if(grepl("sigma", k)){
          theta[k] = 1
        }else if(k == "mu_d"){
          theta[k] = -1.5
        }else if(k == "mu_u"){
          theta[k] = 1.5
        }
      }
    }
  }
  constraints = list(ineqA = ineqA, ineqB = ineqB)
  rm(ineqA, ineqB)
  return(constraints)
}

#' Detect data frequency
#' 
#' @param y Multivariate time series of data values. May also be a data frame containing a date column
#' @param timeID Column name that identifies the date
#' @return Integer
#' @examples
#' detect_frequency(y = y, timeID = timeID)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
detect_frequency = function(y, timeID){
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
    return(freq)
  }else{
    stop("No date column detected. Include a date column or set the frequency.")
  }
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
#' @param y Multivariate time series of data values. May also be a data frame containing a date column
#' @param freq Seasonality of the data
#' @param panelID Column name that identifies the cross section of the data
#' @param timeID Column name that identifies the date
#' @param level Significance level of statistical tests (0.01, 0.05, 0.1)
#' @param formulas R formula describing the relationship between each data series, the unobserved dynamic common factor, and the erorr structure
#' @param optim_methods Vector of 1 to 3 optimization methods in order of preference ("NR", "BFGS", "CG", "BHHH", or "SANN")
#' @param maxit Maximum number of iterations for the optimization
#' @param prior "estimate", "uninformative" or vector of named prior parameter guesses: DCF AR coefficients: phi1 and phi2 only; Error MA coefficients: psi_i1 to psi_i2 only for each series i; Error standard deviation: sigma_i only for each series i; Observation coefficient on DCF with first gamma (gamma_i, ..., gamma_n) only 1 index number, not i0 and any more gammas per equation: gamma_i1 to gamma_ik; Markov switching growth rate: mu_d and mu_u; Transition probabilities: p_dd, p_md (or p_mu), p_mm, p_md (or p_mu), p_uu, p_ud (or p_um)
#' @param log.vars Character vector of variables to be logged
#' @param diff.vars Character vector of unit root variables to be differenced
#' @param diff.lag Integer, number of lags to use for differencing
#' @param n_states Number of states to include in the Markov switching model
#' @param ms_var, Logical, T for Markow switching variance, default is F
#' @param detect.formula Logical, detect lag length of the dynamic common factor to include in each observation equation using the cross correlation function up to a max of 3
#' @param detect.growth Logical, detect which variables to log
#' @param detect.diff Logical, detect which variables to difference
#' @param weighted Logical, use weighted maximum likelihood. Weights are the rescaled inverse of the determinant of the forecast error covariance matrix for each observation.
#' @return List of estimated values including coefficients, convergence code, the panel and time ids, the variables in the data, the variable that were logged, and the variables that were differenced
#' @examples
#' ms_dcf_estim(y = DT[, c("date", "y")])
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
ms_dcf_estim = function(y, freq = NULL, panelID = NULL, timeID = NULL, level = 0.01, 
                        log.vars = NULL, diff.vars = NULL,  detect.formula = F,
                        detect.growth = F, detect.diff = F, diff.lag = 1,
                        n_states = 2, ms_var = F, prior = NULL,
                        formulas = c("y ~ c + e.l1 + e.l2"), weighted = F,
                        optim_methods = c("BFGS", "CG", "NM"), maxit = 1000, trace = F){  

  y = copy(y)
  y = data.table::as.data.table(y)
  if(is.null(freq)){
    freq = detect_frequency(y = y, timeID = timeID)
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
  if(!is.logical(weighted)){
    stop("weighted must be T, F.")
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
  if(!is.numeric(diff.lag)){
    stop("diff.lag must be numeric.")
  }
  if(!is.logical(detect.growth)){
    stop("detect.growth must be T, F.")
  }
  if(!is.logical(detect.diff)){
    stop("detect.diff must be T, F.")
  }
  if(!is.null(prior)){
    if(!prior %in% c("estimate", "uninformative") | !all(is.numeric(prior))){
      stop("prior must be NULL or a named numeric verctor.")
    }
  }
  vars = colnames(y)[!colnames(y) %in% c(panelID, timeID)]
  
  #Transform the data
  data = data_trans(y = y, log.vars = log.vars, diff.vars = diff.vars, freq = freq, diff.lag = diff.lag,
                    detect.diff = detect.diff, detect.growth = detect.growth, level = level,
                    panelID = panelID, timeID = timeID)
  yy_d = data$yy_d
  yy_s = data$yy_s
  rm(data)
  
  #Set the priors
  theta = set_priors(yy_s = yy_s, prior = prior, panelID = panelID, timeID = timeID,
                     n_states = n_states, ms_var = ms_var, 
                     formulas = formulas, detect.formula = detect.formula)
  
  #Set the constraints
  constraints = set_constraints(yy_s = yy_s, n_states = n_states, theta = theta, panelID = panelID, timeID = timeID)
  
  #Find location of NA values
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
  
  #Define the objective function
  objective = function(par, n_states, ms_var, panelID, timeID, na_locs = NULL, weighted){
    yy_s = get("yy_s")
    toret = foreach::foreach(i = unique(yy_s[, c(panelID), with = F][[1]]), .packages = c("data.table"), .export = c("SSmodel_ms", "kim_filter", "kim_smoother")) %fun% {
      sp = SSmodel_ms(par, yy_s, n_states, ms_var, panelID, timeID)
      yti = t(yy_s[eval(parse(text = panelID)) == i, colnames(yy_s)[!colnames(yy_s) %in% c(panelID, timeID)], with = F])
      ans = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti, weighted)
      if(!is.null(na_locs)){
        fc = do.call("cbind", lapply(1:nrow(ans$B_tt), function(x){ans$H_tt[,,x] %*% ans$B_tt[x, ]}))
        rownames(fc) = rownames(yti)
        for(x in rownames(fc)){
          yti[x, na_locs[[i]][[x]]] = fc[x, na_locs[[i]][[x]]]
        }
        ans = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti, weighted)
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
  
  #Initiate the clusters for parallel computing
  cl = parallel::makeCluster(min(c(length(unique(yy_s[, c(panelID), with = F][[1]])), parallel::detectCores())))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  
  #Get initial values for the filter
  sp = SSmodel_ms(theta, yy_s, n_states, ms_var, panelID, timeID)
  
  #Estimate the model
  for(m in optim_methods){
    out = tryCatch(maxLik::maxLik(logLik = objective, start = theta, method = m,
                                  finalHessian = F, hess = NULL, control = list(printLevel = trace, iterlim = maxit), constraints = constraints,
                                  na_locs = na_locs, n_states = n_states, ms_var = ms_var, panelID = panelID, timeID = timeID, weighted = weighted),
                   error = function(err){NULL})
    if(!is.null(out)){
      break
    }
  }
  
  snow::stopCluster(cl)
  return(list(coef = out$estimate, prior = theta, convergence = out$code, loglik = out$maximum, panelID = panelID, timeID = timeID,
              vars = vars, log.vars = log.vars, diff.vars = diff.vars, n_states = n_states, ms_var = ms_var,  level = level, freq = freq,
              diff.lag = diff.lag, formulas = formulas, detect.diff = detect.diff, detect.growth = detect.growth, detect.formula = detect.formula))
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
  y = copy(y)
  data = data_trans(y = y, model = model)
  
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
  
  #Get the unobserved components
  uc = foreach::foreach(i = unique(yy_s[, c(model$panelID), with = F][[1]]), .packages = c("data.table", "MASS"), .export = c("SSmodel_ms", "kim_filter", "kim_smoother")) %fun% {
    yti = t(yy_s[eval(parse(text = model$panelID)) == i, colnames(yy_s)[!colnames(yy_s) %in% c(model$panelID, model$timeID)], with = F])
    ans = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
    if(!is.null(na_locs)){
      fc = do.call("cbind", lapply(1:nrow(ans$B_tt), function(x){ans$H_tt[,,x] %*% ans$B_tt[x, ]}))
      rownames(fc) = rownames(yti)
      for(x in rownames(fc)){
        yti[x, na_locs[[i]][[x]]] = fc[x, na_locs[[i]][[x]]]
      }
      ans = kim_filter(sp$B0, sp$P0, sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat, yti)
    }
    
    #Get the steady state Kalman gain approximation
    Ht = matrix(Reduce("+", lapply(1:dim(ans$H_tt)[3], function(x){ans$H_tt[,, x]}))/dim(ans$H_tt)[3],
                nrow = dim(sp$Ht)[1], ncol = dim(sp$Ht)[2])
    Ft = matrix(Reduce("+", lapply(1:dim(ans$F_tt)[3], function(x){ans$F_tt[,, x]}))/dim(ans$F_tt)[3],
                nrow = dim(sp$Ft)[1], ncol = dim(sp$Ft)[2])
    rownames(Ht) = rownames(yti)
    colnames(Ht) = rownames(Ft) = rownames(sp$Ft[,, 1])
    colnames(Ft) = colnames(sp$Ft[,, 1])
    
    #Get the steady state Kalman gain approximation
    Ht = matrix(Reduce("+", lapply(1:dim(ans$H_tt)[3], function(x){ans$H_tt[,, x]}))/dim(ans$H_tt)[3],
                nrow = dim(sp$Ht)[1], ncol = dim(sp$Ht)[2]) 
    Ft = matrix(Reduce("+", lapply(1:dim(ans$F_tt)[3], function(x){ans$F_tt[,, x]}))/dim(ans$F_tt)[3],
                nrow = dim(sp$Ft)[1], ncol = dim(sp$Ft)[2])
    rownames(Ht) = rownames(yti)
    colnames(Ht) = rownames(Ft) = rownames(sp2$Ft[,, 1])
    colnames(Ft) = colnames(sp2$Ft[,, 1])
    
    dcf_loc = which(rownames(Ft) == "ct0")
    means = unlist(yy_s[eval(parse(text = model$panelID)) == i, lapply(.SD, mean, na.rm = T), .SDcols = c(model$vars)])
    sds = unlist(yy_s[eval(parse(text = model$panelID)) == i, lapply(.SD, sd, na.rm = T), .SDcols = c(model$vars)])
    
    #Check fo convergence of Kt to steady state K
    K = matrix(ans$K[,, dim(ans$K)[3]], nrow(ans$K), ncol(ans$K))
    W = MASS::ginv(diag(nrow(K)) - (diag(nrow(K)) - K %*% Ht) %*% Ft) %*% K
    #W(1) is the first row of W
    d = W[1, ] %*% as.matrix(means)
    
    #Get the intercept terms
    temp = matrix(Ht[, grepl("ct", colnames(Ht))], nrow = nrow(Ht))
    D = means - temp %*% matrix(rep(d, ncol(temp)))
    
    #Initialize first element of the dynamic common factor
    Y1 = as.matrix(yti[, 1])
    nonna_idx = which(!is.na(Y1))
    initC = function(par){
      return(sum((as.matrix(Y1[nonna_idx, ]) - as.matrix(D[nonna_idx, ]) - matrix(Ht[nonna_idx, grepl("ct", colnames(Ht))], nrow = nrow(temp)) %*% matrix(par))^2))
    }
    C10 = optim(par = rep(mean(Y1)/mean(model$coef[grepl("gamma", names(model$coef))]), length(colnames(Ft)[grepl("ct", colnames(Ft))])), 
                fn = initC, method = "BFGS", control = list(trace = F))$par[1]
    
    toret = list()
    for(k in c("filter", "smooth")){
      Ctt = rep(C10, ncol(yti) + 1)
      if(k == "smooth"){
        ans = kim_smoother(ans$B_tlss, ans$B_tts, ans$B_tt, ans$P_tlss, ans$P_tts, ans$Pr_tls, ans$Pr_tts, 
                           sp$At, sp$Dt, sp$Ft, sp$Ht, sp$Qt, sp$Rt, sp$Tr_mat)
        if(any(is.na(ans$B_tt))){
          break
        }
      }
      
      #Build the rest of the DCF series
      #Rescale the first differenced DCF series to match that of the actual series
      ctt = c(0, (t(ans$B_tt[, dcf_loc])))
      if(is.finite(model$n_states)){
        mutt = c(0, t(ans$D_tt[dcf_loc,,]))
      }else{
        m_loc = which(rownames(Ft) == "mt0")
        mutt = c(0, t(ans$B_tt[, m_loc]))
      }
      vartt = c(NA, ans$Q_tt[dcf_loc, dcf_loc, ])
      
      for(j in 2:length(Ctt)){
        #First element of dCtt is C_22 - C_11 = dCtt_2
        Ctt[j] = ctt[j] + Ctt[j - 1] + c(d)
      }
      ctt = ctt + c(d)
      Ctt = data.table::data.table(panelID = i, date = y[(model$diff.lag):.N, ][eval(parse(text = model$panelID)) == i, c(model$timeID), with = F][[1]], DCF = Ctt, d.DCF = ctt, Mu = mutt, Var = vartt)
      colnames(Ctt) = c(model$panelID, model$timeID, "DCF", "d.DCF", "Mu", "Var")
      
      #Build the probability series
      prob = data.table::data.table(panelID = i, date = yy_s[eval(parse(text = model$panelID)) == i, c(model$timeID), with = F][[1]], data.table::data.table(ans$Pr_tts))
      colnames(prob) = c(model$panelID, model$timeID, paste0("Pr_", dimnames(sp$Dt)[[3]]))
      toret[[k]] = merge(Ctt, prob, by = c(model$timeID, model$panelID), all = T)
    }
    return(toret)
  }
  snow::stopCluster(cl)
  names(uc) = unique(yy_s[, c(model$panelID), with = F][[1]])
  
  #Plot the outputs
  if(plot == T){
    for(j in c("filter", "smooth")){
      for(i in unique(y[, c(model$panelID), with = F][[1]])){  
        toplot1 = melt(y[eval(parse(text = model$panelID)) == i, ], id.vars = c(model$panelID, model$timeID))
        toplot1[, "value2" := (value - min(value, na.rm = T))/(diff(range(value, na.rm = T))), by = c("variable")]
        g1 = ggplot2::ggplot(toplot1) + 
          ggplot2::ggtitle(paste(ifelse(i == "panel" & model$panelID == "panelid", "", i), "Data Series"), subtitle = "Levels") + 
          ggplot2::scale_y_continuous(name = "Levels (Rescaled)") + 
          ggplot2::scale_x_date(name = "") + 
          ggplot2::geom_line(ggplot2::aes(x = eval(parse(text = model$timeID)), y = value2, group = variable, color = variable)) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + ggplot2::guides(color = ggplot2::guide_legend(title = NULL))
        
        toplot2 = melt(yy_s[eval(parse(text = model$panelID)) == i, ], id.vars = c(model$panelID, model$timeID))
        g2 = ggplot2::ggplot(toplot2) + 
          ggplot2::ggtitle(paste(ifelse(i == "panel" & model$panelID == "panelid", "", i), "Data Series"), subtitle = "Differenced & Standardized") + 
          ggplot2::scale_y_continuous(name = "Value") + 
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
          g3 = g3 + ggplot2::scale_y_continuous(name = "DCF Value", limits = range(toplot3[variable == "DCF", ]$value, na.rm = T))
        }else{
          g3 = g3 + 
            ggplot2::scale_y_continuous(name = "DCF Value", limits = range(toplot3[variable == "DCF", ]$value, na.rm = T), 
                                        sec.axis = ggplot2::sec_axis(name = "State Probability (%)", ~((. - d_range1[1])/diff(d_range1) * diff(p_range1) + p_range1[1]) * 100)) + 
            ggplot2::geom_ribbon(data = toplot3[variable %in% colnames(uc[[i]][[j]])[grepl(ifelse(model$n_states == 2, "Pr_d", "Pr_d|Pr_u"), colnames(uc[[i]][[j]]), ignore.case = T)], ], 
                                 ggplot2::aes(x = eval(parse(text = model$timeID)), ymin = d_range1[1], ymax = value, group = variable, fill = variable), alpha = 0.5) + 
            ggplot2::scale_fill_manual(values = c("red", "green"))
        }
        
        toplot4 = melt(uc[[i]][[j]], id.vars = c(model$panelID, model$timeID))
        d_range2 = range(toplot4[variable %in% c("d.DCF"), ]$value, na.rm = T)
        p_range2 = range(toplot4[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], ]$value, na.rm = T)
        toplot4[variable %in% colnames(uc[[i]][[j]])[grepl("Pr_", colnames(uc[[i]][[j]]))], "value" := (value - p_range2[1])/diff(p_range2) * diff(d_range2) + d_range2[1], by = "variable"]
        g4 = ggplot2::ggplot() +  
          ggplot2::ggtitle(paste(ifelse(i == "panel" & model$panelID == "panelid", "", i), "Dynamic Common Factor"), subtitle = "Differenced & Standardized") + 
          ggplot2::scale_x_date(name = "") +
          ggplot2::geom_hline(yintercept = 0, color = "grey") + 
          ggplot2::geom_line(data = toplot4[variable %in% c("d.DCF"), ], 
                             ggplot2::aes(x = eval(parse(text = model$timeID)), y = value, group = variable, color = variable)) + 
          ggplot2::scale_color_manual(values = c("black")) + 
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") + 
          ggplot2::guides(color = ggplot2::guide_legend(title = NULL), fill = ggplot2::guide_legend(title = NULL))
        if(model$n_states == 1 | is.infinite(model$n_states)){
          g4 = g4 + ggplot2::scale_y_continuous(name = "DCF Value", limits = range(toplot4[variable %in% c("d.DCF"), ]$value, na.rm = T)) 
        }else{
          g4 = g4 + ggplot2::scale_y_continuous(name = "DCF Value", limits = range(toplot4[variable %in% c("d.DCF"), ]$value, na.rm = T), 
                                                sec.axis = ggplot2::sec_axis(name = "State Probability (%)", ~((. - d_range2[1])/diff(d_range2) * diff(p_range2) + p_range2[1]) * 100)) + 
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
  
  #Get the filtered and smoothed values to return
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
#library(MarkovSwitchingDCF)
#git config remote.origin.url git@github.com:opendoor-labs/MarkovSwitchingDCF.git

