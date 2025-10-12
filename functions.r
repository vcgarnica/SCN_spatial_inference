## Updated Functions for Model Evaluation and Contrasts
##
##' This module provides tools to evaluate and compare ASReml models, including
##' AIC/BIC calculation (`icREML`) and custom contrast estimation for seed treatment
##' and cultivar effects (`st_contrasts`, `cult_contrasts_ss`, `cult_contrasts_ms`).
##'
##' The original `icREML` function (by Ari Verbyla) is used for
##' calculating AIC and BIC for a list of fitted models.
##'
##' Additional functions were added to estimate biologically meaningful contrasts
##' from fitted ASReml models:
##' \itemize{
##'   \item \code{st_contrasts}: Compares seed treatment effects (e.g., Fluopyram vs Base)
##'   \item \code{cult_contrasts_ss}: Computes contrasts between Peking and PI88788 cultivar groups for single-site models
##'   \item \code{cult_contrasts_ms}: Same as above, but for multi-site models
##'   \item \code{cult_contrasts_ms}: Same as above, but for multi-site models
##'   \item \code{cult_contrasts_ms}: Same as above, but for multi-site models   
##' }
##' Each contrast function optionally transforms the log-scale contrast to percent efficacy.
##'
##' @param fm A \code{list} of asreml fitted model objects for `icREML`
##' @param model A single fitted asreml model for contrast functions
##' @param scale Scaling factor for log-determinant stability (default = 1)
##' @param transf Logical; if TRUE, contrasts are back-transformed to percent efficacy
##' @param year Required for `cult_contrasts_ss` to determine cultivar groups
##' @return For `icREML`: a data.frame summarizing AIC/BIC and log-likelihood metrics.
##'         For contrast functions: a data.frame of estimated contrasts and confidence intervals.
##'
##' @author Original AIC/BIC function: Ari Verbyla (2019). DOI: https://doi.org/10.1111/anzs.12254
##'         Additional functions and modifications: Vinicius Garnica
##

library(TPSbits)  
library(tidymodels)

icREML = function(fm, scale=1) {
    if(!is.list(fm)) stop(" Models need to be in a list\n")
    if(is.null(names(fm))) namesfm = paste("model", 1:length(fm))
    else namesfm = names(fm)
    require(asreml)
    asreml.options(Cfixed = TRUE, gammaPar=FALSE)
    fm = lapply(fm, function(el) {
        if(is.null(el$Cfixed)) {
            out = update(el, maxit=50, trace = FALSE) }
        else out = el
        out})
    logl = lapply(fm, function(el) el$loglik)
    summ = lapply(fm, function(el) summary(el, coef=TRUE)$coef.fixed)
    which.X0 = lapply(summ, function(el) !is.na(el[, "z.ratio"]))
    p.0 = lapply(which.X0, function(el) sum(el))
    Cfixed = lapply(fm, function(el) el$Cfixed)
    logdet = lapply(1:length(fm), function(el, Cfixed, which.X0, scale) {
        log(prod(svd(as.matrix(scale*Cfixed[[el]][which.X0[[el]], which.X0[[el]]]))$d))
    }, Cfixed, which.X0, scale)
    vparam = lapply(fm, function(el) summary(el)$varcomp)
    q.0 = lapply(vparam, function(el) sum(!(el$bound == "F" | el$bound == "B")))
    b.0 = lapply(vparam, function(el) sum(el$bound == "F" | el$bound == "B"))
    logl = lapply(1:length(fm), function(el, logl, logdet, p.0) {
        logl[[el]] - logdet[[el]]/2}, logl, logdet,p.0)
    aic = unlist(lapply(1:length(fm), function(el, logl, p.0, q.0) {
        -2*logl[[el]] + 2*(p.0[[el]] + q.0[[el]])}, logl, p.0, q.0))
    bic = unlist(lapply(1:length(fm), function(el, logl, p.0, q.0, fm) {
        -2*logl[[el]] + log(fm[[el]]$nedf+p.0[[el]])*(p.0[[el]] + q.0[[el]])},
        logl, p.0, q.0, fm))
    results = data.frame(model=namesfm, loglik = unlist(logl), p=unlist(p.0),
                          q=unlist(q.0), b = unlist(b.0), AIC = aic, BIC = bic, logdet=unlist(logdet))%>%
      mutate(across(where(is.numeric), ~ round(.x, 1)))
    row.names(results) = 1:dim(results)[1]
    invisible(results)
}


round_df = function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ round(., 1)))
}

round_transf_df = function(df) {
  df %>%
    mutate(
      across(c(predicted.value, standard.error, lower.Confidence.limit, upper.Confidence.limit), 
             ~ round(., 2)),
      across(c(Bt, Bt_lower_CI, Bt_upper_CI), 
             ~ round(., 3)),  # More precision for ratios
      across(c(Eff, Eff_lower_CI, Eff_upper_CI), 
             ~ round(., 1)),   # Percentages get 1 decimal
      across(c(Reference), 
             ~ round(., 3))    # 3 decimals for reference
    )
}

st_contrasts = function(model, transf = FALSE) {
  
  wald_info = wald(model, denDF = 'numeric',trace=FALSE)$Wald
  
  # 1. Seed Treatment Main Effect Contrast ----------------------------
  pred = predict(model, classify = "seed_trt",trace=FALSE, vcov = TRUE)
  pred_df = pred$pvals
  vcov_mat = as.matrix(pred$vcov)
  
  seedtrt_labels = levels(pred_df$seed_trt)
  L_seedtrt = setNames(rep(0, length(seedtrt_labels)), seedtrt_labels)
  L_seedtrt[seedtrt_labels[1]] = -1
  L_seedtrt[seedtrt_labels[2]] = 1
  L_seedtrt = matrix(L_seedtrt, nrow = 1)
  
  seedtrt_contrast = predictPlus.asreml(
    classify = "seed_trt",
    linear.transformation = L_seedtrt,
    asreml.obj = model,
    wald.tab = wald_info,
    tables = "none",
    trace = FALSE
  )$predictions
  
  ref = as.numeric((L_seedtrt * (L_seedtrt < 0) / sum(L_seedtrt[L_seedtrt < 0])) %*% pred_df$predicted.value)
  
  if (transf) {
    
      bt_diff = exp(seedtrt_contrast$predicted.value)
      bt_lower_ci = exp(seedtrt_contrast$lower.Confidence.limit)
      bt_upper_ci = exp(seedtrt_contrast$upper.Confidence.limit)
      
      eff = (1 - bt_diff) * 100
      eff_lower = (1 - bt_upper_ci) * 100  # upper CI on log scale → lower efficacy
      eff_upper = (1 - bt_lower_ci) * 100  # lower CI on log scale → upper efficacy
    
    out = data.frame(
      Contrast = "Fluopyram vs Base",
      seedtrt_contrast,
      Bt = bt_diff,
      Bt_lower_CI = bt_lower_ci,
      Bt_upper_CI = bt_upper_ci,
      Eff = eff,
      Eff_lower_CI = eff_lower,
      Eff_upper_CI = eff_upper,
      Reference = exp(ref)
    )  %>%
      round_transf_df()
  } else {
    out = data.frame(
      Contrast = "Fluopyram vs Base",
      seedtrt_contrast,
      Reference = ref
    ) %>%
      round_df()
  }
  
  # Calculate variance of the cultivar contrast
  var_seed = as.numeric(L_seedtrt %*% vcov_mat %*% t(L_seedtrt))
  var_seed = data.frame(Contrast = "Fluopyram vs Base",var=var_seed)
  
  return(list(contrast = out %>% dplyr::select(-Combination,-est.status),var_contrast=var_seed))
}

cult_contrasts_ss = function(model, year, transf = FALSE) {
  
  contrast_groups_2022 = list(
    Peking = c("P27A26PR", "P34A59PR"),
    PI88788 = c("P28A83PR", "P30A46PR", "P34A65PR")
  )
  
  contrast_groups_2023 = list(
    Peking = c("P27A26PR", "P31A48PR", "P34A59PR"),  
    PI88788 = c("P28A83PR", "P30A46PR", "P34A65PR") 
  )
  
  # Select contrast group by year
  if (year == 2022) {
    contrast_groups = contrast_groups_2022
  } else if (year == 2023) {
    contrast_groups = contrast_groups_2023
  } else {
    stop("Unsupported year. Only 2022 and 2023 are valid.")
  }
  
  wald_info = wald(model, denDF = 'numeric',trace=FALSE)$Wald
  
  # 1. Cultivar Main Effect Contrast ----------------------------
  pred_cult = predict(model, classify = "cult",trace=FALSE, vcov = TRUE)
  pred_cult_df = pred_cult$pvals
  vcov_cult_mat = as.matrix(pred_cult$vcov)
  
  cult_labels = levels(pred_cult_df$cult)
  L_cult = setNames(rep(0, length(cult_labels)), cult_labels)
  L_cult[contrast_groups$Peking] = 1/length(contrast_groups$Peking)
  L_cult[contrast_groups$PI88788] = -1/length(contrast_groups$PI88788)
  L_cult = matrix(L_cult, nrow = 1)
  
  cult_contrast = predictPlus.asreml(
    classify = "cult",
    linear.transformation = L_cult,
    asreml.obj = model,
    wald.tab = wald_info,
    tables = "none",
    trace = FALSE
  )$predictions
  
  # Calculate variance of the cultivar contrast
  var_cult = as.numeric(L_cult %*% vcov_cult_mat %*% t(L_cult))
  var_cult = data.frame(Contrast = "Peking vs PI88788",var=var_cult)
  
  ref_cult = as.numeric((L_cult * (L_cult < 0) / sum(L_cult[L_cult < 0])) %*% pred_cult_df$predicted.value)
  
  # 2. Seed Treatment × Cultivar Interaction Contrast -----------
  pred_in = predict(model, classify = "cult:seed_trt",trace=FALSE, vcov = TRUE)
  pred_in_df = pred_in$pvals
  vcov_in_mat = as.matrix(pred_in$vcov)
  
  interaction_labels = with(pred_in_df, paste(cult, seed_trt, sep = ":"))
  
  L_int = setNames(rep(0, length(interaction_labels)), interaction_labels)
  
  for (cultivar in contrast_groups$Peking) {
    L_int[paste0(cultivar, ":Base")] = 1 / length(contrast_groups$Peking)
  }
  for (cultivar in contrast_groups$PI88788) {
    L_int[paste0(cultivar, ":Fluopyram")] = -1 / length(contrast_groups$PI88788)
  }
  
  L_int = matrix(L_int[match(interaction_labels, names(L_int))], nrow = 1)
  
  
  int_contrast = predictPlus.asreml(
    classify = "cult:seed_trt",
    linear.transformation = L_int,
    asreml.obj = model,
    wald.tab = wald_info,
    tables = "none",
    trace = FALSE
  )$predictions
  
  # Calculate variance of the interaction contrast
  var_int = as.numeric(L_int %*% vcov_in_mat %*% t(L_int))
  var_int = data.frame(Contrast = "Peking + Base vs PI88788 + Fluopyram",var=var_int)
  
  ref_in = as.numeric((L_int * (L_int < 0) / sum(L_int[L_int < 0])) %*% pred_in_df$predicted.value)

  if (transf) {
    bt_diff_cult = exp(cult_contrast$predicted.value)
    bt_lower_ci_cult = exp(cult_contrast$lower.Confidence.limit)
    bt_upper_ci_cult = exp(cult_contrast$upper.Confidence.limit)
    
    eff_cult = (1 - bt_diff_cult) * 100
    eff_lower_cult = (1 - bt_upper_ci_cult) * 100
    eff_upper_cult = (1 - bt_lower_ci_cult) * 100
    
    cult_df = data.frame(
      Contrast = "Peking vs PI88788",
      cult_contrast,
      Bt = bt_diff_cult,
      Bt_lower_CI = bt_lower_ci_cult,
      Bt_upper_CI = bt_upper_ci_cult,
      Eff = eff_cult,
      Eff_lower_CI = eff_lower_cult,
      Eff_upper_CI = eff_upper_cult,
      Reference = exp(ref_cult)) %>%
      round_transf_df()
    
    bt_diff_int = exp(int_contrast$predicted.value)
    bt_lower_ci_int = exp(int_contrast$lower.Confidence.limit)
    bt_upper_ci_int = exp(int_contrast$upper.Confidence.limit)
    
    eff_int = (1 - bt_diff_int) * 100
    eff_lower_int = (1 - bt_upper_ci_int) * 100
    eff_upper_int = (1 - bt_lower_ci_int) * 100
    
    int_df = data.frame(
      Contrast = "Peking + Base vs PI88788 + Fluopyram",
      int_contrast,
      Bt = bt_diff_int,
      Bt_lower_CI = bt_lower_ci_int,
      Bt_upper_CI = bt_upper_ci_int,
      Eff = eff_int,
      Eff_lower_CI = eff_lower_int,
      Eff_upper_CI = eff_upper_int,
      Reference = exp(ref_in)) %>%
      round_transf_df()
    
  } else {
    cult_df = data.frame(
      Contrast = "Peking vs PI88788",
      cult_contrast,
      Reference = ref_cult
    )%>%
      round_df()
    int_df = data.frame(
      Contrast = "Peking + Base vs PI88788 + Fluopyram",
      int_contrast,
      Reference = ref_in
    )%>%
      round_df()
  }
  
  return(list(contrast = bind_rows(cult_df, int_df) %>% dplyr::select(-Combination,-est.status),var_contrast= bind_rows(var_cult, var_int)))
}


cult_contrasts_ms = function(model, transf=FALSE) {
  
  contrast_groups = list(
    Peking = c("P27A26PR", "P34A59PR", "P31A48PR"),  
    PI88788 = c("P28A83PR", "P30A46PR", "P34A65PR") 
  )
  
  wald_info = wald(model, denDF = 'numeric',trace=FALSE)$Wald
  
  # 1. Cultivar Main Effect Contrast ----------------------------
  pred_cult = predict(model, classify = "cult",trace=FALSE, vcov = TRUE)
  pred_cult_df = pred_cult$pvals
  vcov_cult_mat = as.matrix(pred_cult$vcov)
  
  cult_labels = levels(pred_cult_df$cult)
  L_cult = setNames(rep(0, length(cult_labels)), cult_labels)
  L_cult[contrast_groups$Peking] = 1/length(contrast_groups$Peking)
  L_cult[contrast_groups$PI88788] = -1/length(contrast_groups$PI88788)
  L_cult = matrix(L_cult, nrow = 1)
  
  cult_contrast = predictPlus.asreml(
    classify = "cult",
    linear.transformation = L_cult,
    asreml.obj = model,
    wald.tab = wald_info,
    tables = "none",
    trace = FALSE
  )$predictions
  
  # Calculate variance of the cultivar contrast
  var_cult = as.numeric(L_cult %*% vcov_cult_mat %*% t(L_cult))
  var_cult = data.frame(Contrast = "Peking vs PI88788",var=var_cult)
  
  ref_cult = as.numeric((L_cult * (L_cult < 0) / sum(L_cult[L_cult < 0])) %*% pred_cult_df$predicted.value)
  
  # 2. Seed Treatment × Cultivar Interaction Contrast -----------
  pred_in = predict(model, classify = "cult:seed_trt",trace=FALSE, vcov = TRUE)
  pred_in_df = pred_in$pvals
  vcov_in_mat = as.matrix(pred_in$vcov)
  
  interaction_labels = with(pred_in_df, paste(cult, seed_trt, sep = ":"))
  
  L_int = setNames(rep(0, length(interaction_labels)), interaction_labels)
  
  for (cultivar in contrast_groups$Peking) {
    L_int[paste0(cultivar, ":Base")] = 1 / length(contrast_groups$Peking)
  }
  for (cultivar in contrast_groups$PI88788) {
    L_int[paste0(cultivar, ":Fluopyram")] = -1 / length(contrast_groups$PI88788)
  }
  
  L_int = matrix(L_int[match(interaction_labels, names(L_int))], nrow = 1)
  
  
  int_contrast = predictPlus.asreml(
    classify = "cult:seed_trt",
    linear.transformation = L_int,
    asreml.obj = model,
    wald.tab = wald_info,
    tables = "none",
    trace = FALSE
  )$predictions
  
  # Calculate variance of the interaction contrast
  var_int = as.numeric(L_int %*% vcov_in_mat %*% t(L_int))
  var_int = data.frame(Contrast = "Peking + Base vs PI88788 + Fluopyram",var=var_int)
  
  
  ref_in = as.numeric((L_int * (L_int < 0) / sum(L_int[L_int < 0])) %*% pred_in_df$predicted.value)
  
  
  if (transf) {
    bt_diff_cult = exp(cult_contrast$predicted.value)
    bt_lower_ci_cult = exp(cult_contrast$lower.Confidence.limit)
    bt_upper_ci_cult = exp(cult_contrast$upper.Confidence.limit)
    
    eff_cult = (1 - bt_diff_cult) * 100
    eff_lower_cult = (1 - bt_upper_ci_cult) * 100
    eff_upper_cult = (1 - bt_lower_ci_cult) * 100
    
    cult_df = data.frame(
      Contrast = "Peking vs PI88788",
      cult_contrast,
      Bt = bt_diff_cult,
      Bt_lower_CI = bt_lower_ci_cult,
      Bt_upper_CI = bt_upper_ci_cult,
      Eff = eff_cult,
      Eff_lower_CI = eff_lower_cult,
      Eff_upper_CI = eff_upper_cult,
      Reference = exp(ref_cult)) %>%
      round_transf_df()
    
    bt_diff_int = exp(int_contrast$predicted.value)
    bt_lower_ci_int = exp(int_contrast$lower.Confidence.limit)
    bt_upper_ci_int = exp(int_contrast$upper.Confidence.limit)
    
    eff_int = (1 - bt_diff_int) * 100
    eff_lower_int = (1 - bt_upper_ci_int) * 100
    eff_upper_int = (1 - bt_lower_ci_int) * 100
    
    int_df = data.frame(
      Contrast = "Peking + Base vs PI88788 + Fluopyram",
      int_contrast,
      Bt = bt_diff_int,
      Bt_lower_CI = bt_lower_ci_int,
      Bt_upper_CI = bt_upper_ci_int,
      Eff = eff_int,
      Eff_lower_CI = eff_lower_int,
      Eff_upper_CI = eff_upper_int,
      Reference = exp(ref_in)) %>%
      round_transf_df()
    
  } else {
    cult_df = data.frame(
      Contrast = "Peking vs PI88788",
      cult_contrast,
      Reference = ref_cult
    )%>%
      round_df()
    int_df = data.frame(
      Contrast = "Peking + Base vs PI88788 + Fluopyram",
      int_contrast,
      Reference = ref_in
    )%>%
      round_df()
  }
  
  return(list(contrast = bind_rows(cult_df, int_df) %>% dplyr::select(-Combination,-est.status),var_contrast= bind_rows(var_cult, var_int)))
}


spatial_matrix = function(data, x.coord.name = "col", y.coord.name = "row", at.name = "site") {
  data = data.frame(data)
  data[[at.name]] = as.factor(data[[at.name]])
  at.levels = levels(data[[at.name]])
  
  x.coord = data[[x.coord.name]]
  y.coord = data[[y.coord.name]]
  at = data[[at.name]]
  
  dat0 = data.frame(x.coord, y.coord, at)
  colnames(dat0) = c(x.coord.name, y.coord.name, at.name)
  dat0[[at.name]] = as.factor(dat0[[at.name]])
  
  data1 = split(dat0, dat0[[at.name]])
  
  multires = lapply(data1, function(dxy) {
    nseg = c(length(unique(dxy[["col"]])), length(unique(dxy[["row"]])) / 2)
    TPXZg = tpsmmb("col", "row", dxy, nsegments = nseg, asreml = "grp")
    
    all = TPXZg$data[, TPXZg$grp$All]
    
    fC = TPXZg$data[, TPXZg$grp$TP.R.1_fcol]
    fR = TPXZg$data[, TPXZg$grp$TP.C.1_frow]
    fC.R = TPXZg$data[, TPXZg$grp$TP.R.2_fcol]
    C.fR = TPXZg$data[, TPXZg$grp$TP.C.2_frow]
    fC.fR = TPXZg$data[, TPXZg$grp$TP_fcol_frow]
    rest = TPXZg$data[, min(c(which(colnames(TPXZg$data) =="TP.col"), 
                               which(colnames(TPXZg$data) == "TP.row"))):(min(c(TPXZg$grp$TP.R.1_fcol,TPXZg$grp$TP.C.1_frow)) - 1)]
    list(fC, fR, fC.R, C.fR, fC.fR, rest,all)
  })
  
  nrows = sapply(data1, nrow)
  end = cumsum(nrows)
  start = c(1, head(end, -1) + 1)
  
  nColList = lapply(1:7, function(k) {
    unique(unlist(lapply(multires, function(x) colnames(x[[k]]))))
  })
  
  Zl = lapply(1:7, function(k) {
    Z = matrix(0, nrow = nrow(dat0), ncol = length(nColList[[k]]))
    colnames(Z) = nColList[[k]]
    for (j in seq_along(multires)) {
      prov = multires[[j]][[k]]
      Z[start[j]:end[j], colnames(prov)] = as.matrix(prov)
    }
    attr(Z, "variables") = c(x.coord.name, y.coord.name)
    Z
  })
  names(Zl) = c("fC", "fR", "fC.R", "C.fR", "fC.fR","rest","all")
  
  dataToreturn = list()
  dataToreturn$data=
    cbind(
      data,
      data.frame(Zl$rest),
      data.frame(Zl$fR),
      data.frame(Zl$fC.R),
      data.frame(Zl$fC),
      data.frame(Zl$C.fR),
      data.frame(Zl$fC.fR)
    )
  
  col_names = colnames(dataToreturn$data)
  TP.C.1_frow  = grep("^TP\\.C\\.1_frow_", col_names)
  TP.C.2_frow  = grep("^TP\\.C\\.2_frow_", col_names)
  TP.R.1_fcol  = grep("^TP\\.R\\.1_fcol_", col_names)
  TP.R.2_fcol  = grep("^TP\\.R\\.2_fcol_", col_names)
  TP_fcol_frow = grep("^TP_fcol_frow_",  col_names)
  
  grp = list(
    TP.C.1_frow  = TP.C.1_frow,
    TP.C.2_frow  = TP.C.2_frow,
    TP.R.1_fcol  = TP.R.1_fcol,
    TP.R.2_fcol  = TP.R.2_fcol,
    TP_fcol_frow = TP_fcol_frow,
    All = c(TP.C.1_frow, TP.C.2_frow, TP.R.1_fcol, TP.R.2_fcol, TP_fcol_frow)
  )
  dataToreturn$grp = grp
  
  return(dataToreturn)
}



pretty_kable = function(df, align_left_cols = 1, font_size = 12) {
  kable(
    df,
    align = c(rep("l", align_left_cols), rep("r", ncol(df) - align_left_cols)),
    format.args = list(big.mark = ",")
  ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE,
      position = "center",
      font_size = font_size
    )
}

extract_varcomps = function(model) {
  vc = summary(model)$varcomp %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Component") %>%
    dplyr::rename(
      Variance = component,
      SE = std.error
    ) %>%
    dplyr::mutate(
      Formatted = dplyr::case_when(
        Variance > 1000 ~ sprintf("%.1f (%.1f)", Variance, SE),
        Variance < 1 ~ sprintf("%.2f (%.2f)", Variance, SE),
        TRUE ~ sprintf("%.2f (%.2f)", Variance, SE)
      )
    )
  return(vc)
}


extract_logSCN_effect = function(model) {
  coefs = summary(model, coef = TRUE)$coef.fixed
  
  intercept = coefs["(Intercept)", "solution"]
  
  resis_names = c("resis_Peking:log_centered_initial", "resis_PI88788:log_centered_initial")
  
  results = lapply(resis_names, function(coef_name) {
    if (coef_name %in% rownames(coefs)) {
      est = round(coefs[coef_name, "solution"],1)
      se  = round(coefs[coef_name, "std error"],1)
      pct = round((est / intercept) * 100,1)
      tibble(
        Resistance = sub("resis_", "", sub(":log_centered", "", coef_name)),
        Estimate = est,
        StdError = se,
        PercentEffect = pct
      )
    } else {
      tibble(
        Resistance = sub("resis_", "", sub(":log_centered", "", coef_name)),
        Estimate = NA_real_,
        StdError = NA_real_,
        PercentEffect = NA_real_
      )
    }
  }) %>% bind_rows()
}


