###############################################################################
################## Functions estimating stability measures #######################
###############################################################################



# variable inclusion frequency
vif <- function(resampling_coefs, global_coefs) {
  vif <- apply(resampling_coefs, 2, function(x)
    mean(x != 0))
  names(vif) <- names(global_coefs)
  
  return(vif)
}


# relative conditional bias
rcb <- function(resampling_coefs, global_coefs) {
  nom <- apply(resampling_coefs, 2, mean)
  denom <-
    global_coefs * apply(resampling_coefs, 2, function(x)
      mean(x != 0))
  
  rcb <- (nom / denom - 1) * 100
  
  names(rcb) <- names(global_coefs)
  
  return(rcb)
}

# root mean squared difference ratio
rmsdr <- function(resampling_coefs,
                  global_coefs,
                  global_se) {
  global_coef_large <-
    matrix(
      global_coefs,
      nrow = nrow(resampling_coefs),
      ncol = length(global_coefs),
      byrow = T
    )
  
  nom <-
    sqrt(apply((resampling_coefs - global_coef_large) ^ 2, 2, mean))
  denom <- global_se
  
  rmsdr <- nom / denom
  
  names(rmsdr) <- names(global_coefs)
  
  return(rmsdr)
}

# summary matrix
matrix_summary <-
  function(resampling_coefs,
           global_coefs,
           global_se,
           selected_coefs,
           selected_se = NULL) {
    mat <- round(
      cbind(
        "Estimate, global" = global_coefs,
        "Standard error, global" = global_se,
        "Variable inclusion frequency (%)" = vif(resampling_coefs, global_coefs) *
          100,
        "Estimate, selected" = selected_coefs,
        "Standard error, selected" = if (is.null(selected_se))
          rep(NA, length.out = selected_coefs)
        else
          selected_se,
        "RMSD ratio" = rmsdr(resampling_coefs, global_coefs, global_se),
        "Relative conditional bias (%)" = rcb(resampling_coefs, global_coefs)
      )
      ,
      3
    )
    mat[order(vif(resampling_coefs, global_coefs), decreasing = T), ]
  }