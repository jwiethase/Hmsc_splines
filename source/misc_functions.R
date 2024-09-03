prepareMGCVsplines <- function(df, colname, n_splines = 2){
      #' See ?smooth.terms for a list of available spline functions
      #' fx=FALSE: Penalized regression splines: Number of basis functions (and thus the smoothness) 
      #' is not fixed and can be chosen as part of model fitting
      #' This uses n_knots-2 representative subsets of the data to create spline basis functions.
      #' The last spline in s(x) is always perfectly correlated with x. 
      require(mgcv)
      expr_str <- paste0("smoothCon(s(", colname, ", bs = 'tp', k = ", n_splines+1, ", fx = FALSE), data = df, absorb.cons = TRUE)[[1]]$X")
      base_functions <- eval(parse(text = expr_str))
      
      num_columns <- if (n_splines == 1) { n_splines + 1 } else { n_splines }
      
      for(i in 1:num_columns){
            df[, paste0(colname, "_spline", i)] <- base_functions[, i]
      }
      return(df)
}

constructGradient_JHW <- function (hM, focalVariables, ngrid = NULL, coordinates = list(),
                                   type = "marginal") {
      require(mgcv)  
      require(dplyr) 
      require(sp)
      
      vars = all.vars(hM$XFormula)
      
      # Identify indices of focal and non-focal variables
      focal_indices = which(vars %in% focalVariables)
      non_focals = setdiff(1:length(vars), focal_indices)
      
      # Grab the base variable name
      if(length(focalVariables) > 1){
            base_variable = sub("_spline[0-9]+$", "", focalVariables[1]) # Removing spline suffix to get base variable
      } else {
            base_variable = focalVariables
      }
      
      # Create a sequence for the primary focal variable
      v.focal = hM$XData[, base_variable, drop = FALSE]
      
      if(is.null(ngrid)){
            # If ngrid is not specified, use the original values for the prediction
            grid_values = v.focal
      } else {
            grid_values = seq(min(v.focal, na.rm = TRUE), max(v.focal, na.rm = TRUE), length.out = ngrid)
      }
      
      # Prepare a new data frame for the gradient
      XDataNew = data.frame(base_variable = grid_values)
      names(XDataNew) <- base_variable
      
      # If splines are being used, we need to create prediction data using the original spline basis functions
      if(length(focalVariables) > 1){
            # Get the number of splines
            n_splines = length(which(vars %in% focalVariables))
            
            # Reconstruct the original spline basis functions, and check if they match to the original ones used for model fitting
            smooth_construct <- eval(parse(text = paste0("smoothCon(s(", base_variable, ", bs = 'tp', k = ", n_splines+1,", fx = FALSE), data = hM$XData, absorb.cons = TRUE)[[1]]")))
            if(!identical(smooth_construct$X[, 1], hM$XData[, paste0(base_variable, "_spline", 1)])){
                  stop("Spline basis function for predictions is not identical to spline basis functions used for model fitting.")
            }
            
            # Create extended spline basis functions across the new prediction gradient
            spline_bases <- eval(parse(text = paste0("PredictMat(smooth_construct, data.frame(", base_variable, " = grid_values))")))
            
            num_columns <- if (n_splines == 1) { n_splines + 1 } else { n_splines }
            
            for(i in 1:num_columns){
                  XDataNew[, paste0(base_variable, "_spline", i)] <- spline_bases[, i]
            }
      }
      
      if(type == "total"){
            # Sets the values of the non-focal variable to most likely value, given the value of focal variable, 
            # based on a linear relationship. This is essentially the change in Y with the focal variable 
            # under an environment shift, it does not isolate the effect
            for (i in non_focals) {
                  lm_df <- data.frame(non_focal = hM$XData[, vars[i]], focal = v.focal[[1]])
                  mylm = lm(non_focal ~ focal, data = lm_df)
                  yy = predict(mylm, newdata = data.frame(focal = grid_values))
                  XDataNew[, vars[i]] = yy
            }
      }
      
      if(type == "marginal"){
            # Handling non-focal variables by averaging. This gets the marginal effect, i.e. the isolated
            # predicted Y change with the gradient of the focal variable, all else being constant at typical
            # conditions
            for (i in non_focals) {
                  XDataNew[, vars[i]] = mean(hM$XData[, vars[i], drop = FALSE][[1]], na.rm = TRUE)
            }
      }
      
      if(is.null(ngrid)){ngrid_mat = length(hM[["ranLevels"]][["site_ID"]][["pi"]])} else {ngrid_mat = ngrid} 
      dfPiNew = matrix("new_unit", ngrid_mat, hM$nr)
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
                        ylabel = ifelse(all(hM$distr[,1] == 2), "Richness", 
                                        ifelse(all(hM$distr[, 1] == 3), "Total count", "Summed response"))),
             "Y" = list(data = abind::abind(pred, along = 3)[, index, ], ylabel = hM$spNames[[index]]),
             "T" = list(data = abind::abind(lapply(pred, function(a) (a %*% hM$Tr)/matrix(rep(rowSums(a), hM$nt), ncol = hM$nt)), along = 3)[, index, ], ylabel = hM$trNames[[index]])
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

plotGradient_JHW <- function(hM, Gradient, pred, measure, xlabel = NULL, ylabel = NULL, index = 1,
                             q = c(0.025, 0.5, 0.975), showData = TRUE, supportLevel = 0.9, 
                             showPosteriorSupport = TRUE, back_transform_Y = NULL, back_transform_T = FALSE, ...) 
{
      require(tidyverse)
      require(abind)
      
      if (is.null(xlabel)) {
            xlabel <- colnames(Gradient$XDataNew)[1]
      }
      
      new_x_data <- Gradient$XDataNew[1]
      
      if (!is.null(back_transform_Y)) {
            transform_func <- switch(back_transform_Y, 
                                     sqrt = function(x) x^2, 
                                     log = function(x) exp(x), 
                                     identity)
            hM$Y <- transform_func(hM$Y)
      }
      
      if(back_transform_T){
            hM$Tr[, index] <- exp(hM$Tr[, index])
      }
      
      # Calculate the direction and strength of the effect of the gradient 
      measure_data <- get_measure_data(measure, hM, pred, index) # For each site, get all posterior estimates
      data <- abind::abind(measure_data$data, along = 2)  # Rearrange so each site is a row, and each column has the posterior samples
      
      # Get the proportion of posterior samples where the predicted value at the end of the gradient is greater 
      # than the predicted value at the beginning of the gradient.
      # This is the likelihood of an effect in a specific direction.
      Pr <- mean(
            data[NROW(new_x_data), ] >        # The last row
                  data[1, ], na.rm = T)  # The first row
      qpred <- calculate_quantiles(data, q)
      
      if (is.null(ylabel)) ylabel <- measure_data$ylabel
      
      plot_data <- data.frame(new_x_data = new_x_data[[1]], 
                              median = qpred[2, ], 
                              lo = qpred[1, ], 
                              hi = qpred[3, ])
      
      pl <- ggplot(plot_data, aes(x = new_x_data, y = median)) + 
            geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.5, fill = 'blue') + 
            geom_line(lwd = 1) + 
            xlab(xlabel) + 
            ylab(ylabel) + 
            theme_bw()
      
      if (showData) {
            pX <- hM$XData[, names(new_x_data)]
            
            Y_data <- hM$Y
            Y_data[is.na(Y_data)] <- 0
            
            pY <- switch(measure,
                         "E" = calculate_pielou(Y_data),  # The observed evenness values
                         "S" = rowSums(Y_data, na.rm = TRUE),
                         "Y" = Y_data[, index],
                         "T" = (Y_data %*% hM$Tr) / matrix(rep(rowSums(Y_data), hM$nt), ncol = hM$nt))
            
            if(measure == 'T'){
                  dataToPlot <- data.frame(pX = pX, pY = as.numeric(pY[, index]))
            } else {
                  dataToPlot <- data.frame(pX = pX, pY = as.numeric(pY))
            }
            
            pl <- pl + 
                  geom_point(data = dataToPlot, aes(x = pX, y = pY), alpha = 0.5)
      }
      
      if (showPosteriorSupport && !is.factor(new_x_data)) {
            direction <- ifelse(Pr < 0.5, "decrease", "increase")
            Pr_value <- round(ifelse(Pr < 0.5, 1 - Pr, Pr), digits = 2)
            pl <- pl + 
                  labs(title = paste0("P(", direction, ") = ", Pr_value)) +
                  theme(plot.title = element_text(size = 10,
                                                  color = ifelse(Pr_value >= supportLevel, "darkgreen", "black"),
                                                  face = ifelse(Pr_value >= supportLevel, "bold", "plain")))
      }
      return(pl)
}

bind0 = function(...) {
      if (!require("abind", quietly = TRUE)) {
            stop("The abind package is required but is not installed.")
      }
      abind(..., along = 0)
}

