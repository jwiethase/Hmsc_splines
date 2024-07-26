prepareMGCVsplines <- function(df, colname, n_splines = 2, smooth = 'tp'){
      #' See ?smooth.terms for a list of available spline functions
      #' fx=FALSE: Penalized regression splines: Number of basis functions (and thus the smoothness) 
      #' is not fixed and can be chosen as part of model fitting
      #' This uses n_knots-2 representative subsets of the data to create spline basis functions.
      #' The last spline in s(x) is always perfectly correlated with x. 
      require(mgcv)
      expr_str <- paste0("smoothCon(s(", colname, ", bs = '", smooth, "', k = ", n_splines+1, ", fx = FALSE), data = df, absorb.cons = TRUE)[[1]]$X")
      base_functions <- eval(parse(text = expr_str))
      
      num_columns <- if (n_splines == 1) { n_splines + 1 } else { n_splines }
      
      for(i in 1:num_columns){
            df[, paste0(colname, "_spline", i)] <- base_functions[, i]
      }
      return(df)
}

constructGradient_JHW <- function (hM, focalVariables, ngrid = 25, coordinates = list()) {
      if (!require("mgcv", quietly = TRUE)) {
            stop("The mgcv package is required but is not installed.")
      }
      
      if (!require("dplyr", quietly = TRUE)) {
            stop("The dplyr package is required but is not installed.")
      }
      
      vars = all.vars(hM$XFormula)
      
      # Identify indices of focal and non-focal variables
      focal_indices = which(vars %in% focalVariables)
      non_focals = setdiff(1:length(vars), focal_indices)
      
      # Create a sequence for the primary focal variable
      if(length(focalVariables) > 1){
            base_variable = sub("_spline[0-9]+$", "", focalVariables[1]) # Removing spline suffix to get base variable
      }
      v.focal = hM$XData[, base_variable, drop = FALSE]
      grid_values = seq(min(v.focal, na.rm = TRUE), max(v.focal, na.rm = TRUE), length.out = ngrid)
      
      # Prepare a new data frame for the gradient
      XDataNew = data.frame(base_variable = grid_values)
      names(XDataNew) <- base_variable
      
      if(length(focalVariables) > 1){
            n_splines = length(focalVariables)
            # Apply splines transformations to the grid
            expr_str <- paste0("smoothCon(s(", base_variable, ", bs = 'tp', k = ", n_splines+1,", fx = FALSE), data = XDataNew, absorb.cons = TRUE)[[1]]$X")
            spline_bases <- eval(parse(text = expr_str))
            
            num_columns <- if (n_splines == 1) { n_splines + 1 } else { n_splines }
            
            for(i in 1:num_columns){
                  XDataNew[, paste0(base_variable, "_spline", i)] <- spline_bases[, i]
            }
      }
      
      # Handling non-focal variables by averaging
      # for (i in non_focals) {
      #       XDataNew[, vars[i]] = mean(hM$XData[, vars[i], drop = FALSE][[1]], na.rm = TRUE)
      # }
      # 
      
      # Set the values of the non-focal variable to most likely value, given the value of focal variable, 
      # based on a linear relationship
      for (i in non_focals) {
            lm_df <- data.frame(non_focal = hM$XData[, vars[i]], focal = v.focal[[1]])
            mylm = lm(non_focal ~ focal, data = lm_df)
            yy = predict(mylm, newdata = data.frame(focal = grid_values))
            XDataNew[, vars[i]] = yy
      }
      
      dfPiNew = matrix("new_unit", ngrid, hM$nr)
      colnames(dfPiNew) = hM$rLNames
      dfPiNew = as.data.frame(dfPiNew, stringsAsFactors = TRUE)
      rLNew = vector("list", hM$nr)
      names(rLNew) = hM$rLNames
      for (r in seq_len(hM$nr)) {
            rLname <- hM$rLNames[r]
            coord <- coordinates[[rLname]]
            if (!is.null(coord) && !(is.numeric(coord) || coord %in% 
                                     c("c", "i"))) 
                  stop("'coordinates' must be 'c', 'i' or numeric in random level ", 
                       sQuote(rLname))
            rL1 = hM$rL[[r]]
            xydata = rL1$s
            if (!is.null(xydata)) {
                  if (is(xydata, "Spatial")) {
                        if (!is.null(coord) && coord == "i") 
                              centre <- rep(Inf, NCOL(coordinates(xydata)))
                        else if (is.numeric(coord)) 
                              centre <- coord
                        else centre <- as.data.frame(t(colMeans(coordinates(xydata))))
                        rownames(centre) <- "new_unit"
                        coordinates(centre) <- colnames(centre)
                        proj4string(centre) <- proj4string(xydata)
                        xydata <- rbind(xydata, centre)
                  }
                  else {
                        if (!is.null(coord) && coord == "i") 
                              centre <- rep(Inf, NCOL(xydata))
                        else if (is.numeric(coord)) 
                              centre <- coord
                        else centre <- colMeans(xydata)
                        xydata = rbind(xydata, new_unit = centre)
                  }
            }
            rL1$s = xydata
            distMat = rL1$distMat
            if (!is.null(distMat)) {
                  units1 = c(rownames(distMat), "new_unit")
                  if (is.numeric(coord)) {
                        stop("numeric coordinates are not enabled for 'distMat'")
                  }
                  else if (!is.null(coord) && coord == "i") {
                        newdist <- rep(Inf, NCOL(distMat))
                  }
                  else {
                        newdist <- distMat^2/2
                        newdist <- sweep(newdist, 2L, colMeans(newdist), 
                                         check.margin = FALSE)
                        newdist <- sweep(newdist, 1L, rowMeans(newdist), 
                                         check.margin = FALSE)
                        newdist <- sqrt(diag(-newdist))
                  }
                  distMat1 = cbind(distMat, newdist)
                  distMat1 = rbind(distMat1, c(newdist, 0))
                  rownames(distMat1) = units1
                  colnames(distMat1) = units1
                  rL1$distMat = distMat1
            }
            rL1$pi = c(rL1$pi, "new_unit")
            rL1$N = rL1$N + 1
            rLNew[[r]] = rL1
      }
      
      Gradient = list(XDataNew = XDataNew, studyDesignNew = dfPiNew, rLNew = rLNew)
      return(Gradient)
}

calculate_quantiles <- function(data, probs) {
      apply(data, 1, quantile, probs = probs, na.rm = TRUE)
}

get_measure_data <- function(measure, hM, pred, index) {
      switch(measure,
             "E" = list(data = lapply(pred, calculate_pielou), ylabel = "Pielou's Evenness"),
             "S" = list(data = lapply(pred, rowSums), 
                        ylabel = ifelse(all(hM$distr[,1] == 2), "Species richness", 
                                        ifelse(all(hM$distr[, 1] == 3), "Total count", "Summed response"))),
             "Y" = list(data = abind::abind(pred, along = 3)[, index, ], ylabel = hM$spNames[[index]]),
             "T" = list(data = abind::abind(lapply(pred, function(a) (exp(a) %*% hM$Tr)/matrix(rep(rowSums(exp(a)), hM$nt), ncol = hM$nt)), along = 3)[, index, ], ylabel = hM$trNames[[index]])
      )
}

calculate_pielou <- function(data) {
      if (!require("vegan", quietly = TRUE)) {
            stop("The vegan package is required but is not installed.")
      }
      data[is.na(data)] <- 0
      H <- diversity(data)
      pielou <- H/log(specnumber(data))
      return(pielou)
}

plotGradient_JHW <- function (hM, Gradient, pred, measure, xlabel = NULL, ylabel = NULL, index = 1,
                              q = c(0.025, 0.5, 0.975), showData = TRUE, 
                              showPosteriorSupport = TRUE, back_transform_Y = NULL, unscale_df = NULL, ...) 
{
      if (!require("tidyverse", quietly = TRUE)) {
            stop("The tidyverse package is required but is not installed.")
      }
      if (!require("abind", quietly = TRUE)) {
            stop("The abind package is required but is not installed.")
      }
      
      if (is.null(xlabel)) {
            xlabel <- colnames(Gradient$XDataNew)[1]
      }
      
      new_x_data <- Gradient$XDataNew[1]
      
      if(!is.null(unscale_df)){
            var <- names(new_x_data)
            new_x_data[[var]] <- unscale(new_x_data[[var]], 
                                                     scaling_parameters[[paste0(var, "_mean")]], 
                                                     scaling_parameters[[paste0(var, "_sd")]])
      }

      if (!is.null(back_transform_Y)) {
            transform_func <- switch(back_transform_Y, 
                                     sqrt = function(x) x^2, 
                                     log = function(x) exp(x), 
                                     identity)
            hM$Y <- transform_func(hM$Y)
      }
      
      # Calculate the direction and strength of the effect of the gradient 
      measure_data <- get_measure_data(measure, hM, pred, index)
      data <- abind::abind(measure_data$data, along = 2)
      Pr <- mean(data[NROW(new_x_data), ] > data[1, ])

      qpred <- calculate_quantiles(data, q)
      
      if (is.null(ylabel)) ylabel <- measure_data$ylabel
      
      plot_data <- data.frame(new_x_data = new_x_data[[1]], median = qpred[2, ], lo = qpred[1, ], hi = qpred[3, ])
      
      pl <- ggplot(plot_data, aes(x = new_x_data, y = median)) + 
            geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.5, fill = 'blue') + 
            geom_line(lwd = 1) + 
            xlab(xlabel) + 
            ylab(ylabel) + 
            theme_bw()
      
      if (showData) {
            pX <- hM$XData[, names(new_x_data)]
            
            if(!is.null(unscale_df)){
                  pX <- unscale(pX, scaling_parameters[[paste0(var, "_mean")]], 
                                scaling_parameters[[paste0(var, "_sd")]])
            }
            
            pY <- switch(measure,
                         "E" = calculate_pielou(hM$Y),  # The observed evenness values
                         "S" = rowSums(hM$Y, na.rm = TRUE),
                         "Y" = hM$Y[, index],
                         "T" = (hM$Y %*% hM$Tr) / matrix(rep(rowSums(hM$Y), hM$nt), ncol = hM$nt))
            
            dataToPlot <- data.frame(pX = pX, pY = as.numeric(pY))
            pl <- pl + 
                  geom_point(data = dataToPlot, aes(x = pX, y = pY), alpha = 0.5)
      }
      
      if (showPosteriorSupport && !is.factor(new_x_data)) {
            direction <- ifelse(Pr < 0.5, "decrease", "increase")
            Pr_value <- ifelse(Pr < 0.5, 1 - Pr, Pr)
            pl <- pl + 
                  ggtitle(sprintf("Probability of %s %s with %s = %.2f", 
                                  ylabel, direction, xlabel, Pr_value)) +
                  theme(plot.title = element_text(size = 10))
      }
      
      return(pl)
}

bind0 = function(...) {
      if (!require("abind", quietly = TRUE)) {
            stop("The abind package is required but is not installed.")
      }
      abind(..., along = 0)
}

