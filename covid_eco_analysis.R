#..........................................................................................
###  ECOLOGICAL STUDY OF FACTORS ASSOCIATED WITH COVID-19 TRANSMISSIBILITY AND SEVERITY ###
#..........................................................................................

#..........................................................................................
## ------ R CODE TO PREPARE DATA, VISUALISE PATTERNS AND FIT STATISTICAL MODELS -------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Oct 2020)
                                          # (unless otherwise noted)

                                          # francesco.checchi@lshtm.ac.uk 


#..........................................................................................
### Preparatory steps
#..........................................................................................

  #...................................      
  ## Install or load required R packages
    
    # List of required packages
    x1 <- c("conflicted", "countrycode", "data.table", "dismo", "GGally", "ggpubr", "gtools", "lubridate", "MASS", 
      "mice", "patchwork", "randomForest", "randomForestExplainer", "RColorBrewer", "readxl", "scales", "tidyverse", "hrbrthemes")
    
    # Install any packages not yet installed
    x2 <- x1 %in% row.names(installed.packages())
    if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
    
    # Load all packages    
    lapply(x1, library, character.only = TRUE)


  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
    # Set font
    windowsFonts(Arial=windowsFont("Arial"))

    # Set working directory to where this file is stored
    current_path = rstudioapi::getActiveDocumentContext()$path 
    setwd(dirname(current_path ))
    print( getwd() )
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    

#.........................................................................................
### Specify parameters
#.........................................................................................    

  #...................................   
  ## Completeness threshold for accepting variables into regression
  threshold_compl <- 60 # at least this percentage must be non-missing
     
  #...................................   
  ## Completeness threshold for accepting variables into regression
  k_folds <- NA # Number of folds for cross-validation; NA = leave-one-out CV

  
      
#.........................................................................................
### Bespoke functions
#.........................................................................................    

  #...................................  
  ## Function to perform K-fold cross-validation and, if desired, plot predictions vs. observations for each fold
    # f_fit = model fit object; f_k_folds = see Parameters; f_plot = TRUE or FALSE (do you want a plot?)
    f_cv <- function(f_fit, f_k_folds, f_plot) {
      
      # access dataset used for fit
      f_obs <- fit$data

      # select observations for which all the formula variables are non-missing
      f_obs <- f_obs[complete.cases(f_obs[, all.vars(formula(f_fit)) ] ), ] 

      # determine dependent variable
      f_dep <- all.vars(formula(f_fit))[1]
    
      # determine number of folds if f_k_folds = NA (i.e. LOOCV case)
      if (is.na(f_k_folds) == TRUE) { x1 <- nrow(f_obs) }
      if (is.na(f_k_folds) == FALSE) { x1 <- f_k_folds }

      # shuffle dataset
      f_obs <- f_obs[sample(nrow(f_obs), nrow(f_obs), replace=FALSE), ]
    
      # split data into K folds
        # remove a few rows so as to come up with a n row divisible by K
        f_obs <- f_obs[1:(floor(nrow(f_obs)/x1) * x1), ]
        # split
        folds <- split(f_obs, (0:(nrow(f_obs)-1) %/% (nrow(f_obs)/x1)))
     
      # fit model on all the unfolded sets and track square residuals of model fit on each fold, as well as predictions for each fold  
        # vector to hold squared residuals and other statistics
        errors <- c()
        aic <- c()
        # vector to hold observations and predictions for each fold
        observations <- c()
        predictions <- c()
        
      for (i in 1:length(folds) ) {	
        # progress statement
        print(paste("now on fold ", i, " of ", length(folds), sep = "") )
        # fit on all data but the fold
        data_now <- do.call("rbind", folds[-i])
        cv.fit <- update(f_fit, formula=formula(f_fit),  data = data_now)
        # calculate squared residual of model when predicting fold data, add to other errors for single test observations
        x1 <- predict(cv.fit, newdata = folds[[i]], type = "response", allow.new.levels = TRUE)
        # update output
        observations <- c(observations, folds[[i]][, f_dep])
        predictions <- c(predictions, x1)
        errors <- c(errors , (folds[[i]][, f_dep] -  x1 )^2 )
        aic <- c(aic, AIC(cv.fit))

      }
      
      # return RMSE across all folds
      print("mean RMSE across all folds:")
      print(mean(errors, na.rm = TRUE)^0.5)
      
      # return AIC across all folds
      print("mean AIC across all folds:")
      print(mean(aic, na.rm = TRUE))

      # if plot is desired...
      if (f_plot == TRUE) {
        # prepare data
        x1 <- as.data.frame(cbind(observations, predictions))
        colnames(x1) <- c("observations", "predictions")
          # if on log scale, back-transform to linear
          if (grepl("ln", f_dep) ) {x1[, c("observations", "predictions")] <- exp(x1[, c("observations", "predictions")]) }

        # plot
        plot <- ggplot(x1) +
          geom_point(aes(x = observations, y = predictions), size=2, colour = brewer_pal(palette = "Dark2")(2)[1]) + 
          theme_bw() +
          scale_x_continuous("observed") +
          scale_y_continuous("predicted") +  
          geom_abline(intercept = 0, slope = 1, colour = brewer_pal(palette = "Dark2")(2)[2] ) +
          theme(axis.title = element_text(colour="grey20")) +
         ggtitle(paste("accuracy of predictions on cross-validation; model to predict ", all.vars(formula(f_fit))[1] , sep="") ) +
         theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") )
        
        print("plot shows accuracy of predictions on cross-validation")
        print(plot)
        
        # return plot data
        invisible(x1)
        
      }
          
    }
    
       
     
  #...................................
  ## Function to examine model diagnostics of OLS fit
    f_diag_ols <- function(f_fit) {
      
      # is the model mixed (i.e. does it include a random effect)?  
      f_lmm <- grepl("|", deparse1(formula(f_fit)) )
      
      # prepare data to plot
      if (f_lmm == TRUE) { x1 <- data.frame(resid(f_fit), fitted(f_fit), sd(resid(f_fit)) ) }
      if (f_lmm == FALSE) { x1 <- data.frame(residuals(f_fit), fitted(f_fit), sd(residuals(f_fit)) ) }
      colnames(x1) <- c("residuals", "fitted_values", "sd_residuals")
      
      # normal distribution of residuals
      plot <- ggplot(x1) +
        geom_histogram(aes(x = residuals), colour = brewer_pal(palette = "Greens")(9)[9], 
          fill = brewer_pal(palette = "Greens")(9)[4] ) +
        theme_bw() +
        scale_x_continuous("residuals" ) +
        ggtitle(paste("distribution of residuals for ", all.vars(formula(f_fit))[1] , sep="")) +
        theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") ) +
        stat_function(fun = dnorm, args = list(mean = mean(x1$residuals), sd = sd(x1$residuals) ) )
      print(plot)
      
      # homoskedasticity
      plot <- ggplot(x1) +
        geom_point(aes(x = fitted_values, y = residuals), size = 2, colour = brewer_pal(palette = "Dark2")(2)[1] ) +
        theme_bw() +
        scale_x_continuous("fitted values" ) +
        scale_y_continuous("residuals" ) +
        geom_abline(intercept = 0, slope = 0, linetype = "dotted", colour = brewer_pal(palette = "Dark2")(2)[2] ) +
        ggtitle(paste("residuals vs. fitted values for ", all.vars(formula(f_fit))[1] , sep="")) +
        theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") )
      print(plot)
      
      # quantile-quantile plot
      plot <- ggplot(x1, aes(sample = residuals / sd_residuals)) +
        theme_bw() +        
        stat_qq(size = 2, colour = brewer_pal(palette = "Dark2")(2)[1]) + 
        stat_qq_line(colour = brewer_pal(palette = "Dark2")(2)[2]) +
        ggtitle(paste("quantile-quantile plot of residuals for ", all.vars(formula(f_fit))[1] , sep="") ) +
        theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") ) +
        scale_x_continuous("theoretical quantiles" ) +
        scale_y_continuous("standardised residual quantiles")
      print(plot)      

    }  
   

  #...................................   
  ## Function to plot histograms of variables
  f_hist <- function(f_var, f_data, f_lims) {
      
    plot <- ggplot(f_data)
        
      # if the variable has >= 20 unique values...
        if (length(unique(na.omit(f_data[, f_var]))) >= 20) {
          plot <- plot + geom_histogram(aes(x = as.numeric(f_data[, f_var]) ), 
            color="seagreen", fill="seagreen3", alpha = 0.5 ) +
            theme_bw() + xlab(f_var) + scale_x_continuous(expand = c(0, 0), limits = f_lims )
        }
   
      # otherwise...
        if (length(unique(na.omit(f_data[, f_var]))) < 20) {
          plot <- plot + geom_histogram(aes(x = as.numeric(f_data[, f_var]) ), stat="count", 
            color="seagreen", fill="seagreen3", alpha = 0.5) +
            theme_bw() + xlab(f_var) + scale_x_continuous(expand = c(0, 0), limits = f_lims )
        }
          
      print(plot)
    }


  #...................................
  ## Function to fit a simple linear model and display clean results
  f_lm <- function(f_dep, f_preds, f_data) {
    # write the model formula
    form <- as.formula( paste(f_dep, " ~ ", paste(f_preds, collapse= " + "), sep="")  )
    
    # fit linear model
    fit <- lm(form, data = f_data)
    
    # return fit
    return(fit)
  }

    
  # #...................................
  # ## Function to fit a mixed linear model and display clean results
  #   f_lmm <- function(f_dep, f_preds, f_reff, f_data) {
  #     # write the model formula
  #       form <- as.formula( paste(f_dep, " ~ ", paste(f_preds, collapse= " + "), " + ", f_reff, sep="")  )
  #     # fit linear mixed model
  #       fit <- lmer(form, data = f_data, REML = FALSE)
  #     # return fit
  #       return(fit)
  #   }


  # #...................................      
  # ## Function to prepare data for linear regression
  # f_lm_prep <- function(f_outcome, f_model, f_vars, f_df, f_threshold_compl, 
  #   f_collinear_vars, f_impute, f_region) {
  #   # Select data
  #     # identify variables
  #     x1 <- f_vars[! is.na(f_vars[, f_model]), "variable"]
  #       # select only variables whose completeness is at least the threshold value
  #       x1 <- f_vars[f_vars$variable %in% x1 & f_vars$percent_w_data >= f_threshold_compl, "variable"]
  #       # exclude variables that are causing a lot of collinearity, and don't help to identify clear biological mechanisms for associations
  #       x1 <- x1[! x1 %in% f_collinear_vars]
  # 
  #   # If analysis is only being done for a given WHO region, select only data from this region
  #   if (! is.na(f_region) ) {f_df <- subset(f_df, region == f_region)}
  #           
  #   # Select complete cases or perform imputation 
  #   if (f_impute == FALSE) {
  #     # select only complete cases
  #     x2 <- f_df[complete.cases(f_df[, c(f_outcome, x1)]), c(f_outcome, x1)]
  #   }
  #   
  #   if (f_impute == TRUE) {
  #     # select cases for which the outcome is complete
  #     x3 <- f_df[complete.cases(f_df[, f_outcome]), c(f_outcome, x1)]
  #     # perform imputation by proximity method for missing independent variable values
  #     x2 <- rfImpute(as.formula(paste(f_outcome, "~", ".", sep = "")), x3, iter = 50)
  #   }
  # 
  # # Output prepared dataframe            
  # return(x2)
  # }

  
  #...................................
  ## Function to write both full and reduced final linear models
  f_lm_write <- function(f_filename) {
      
    # Write full model
    write.table("Full linear regression model:", f_filename, sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(tidy(fit, conf.int = TRUE), f_filename, sep = ",", row.names = FALSE, append = TRUE)
    write.table("--------------------", f_filename, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table(glance(fit), f_filename, sep = ",", row.names = FALSE, append = TRUE)
      
    write.table("+++++++++++++++++++++++++++++++++++++", f_filename, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
      
    # Append reduced model to the same file
    write.table("Reduced linear regression model with influential interactions:", f_filename, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table(tidy(fit_red, conf.int = TRUE), f_filename, sep = ",", row.names = FALSE, append = TRUE)
    write.table("--------------------", f_filename, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table(glance(fit_red), f_filename, sep = ",", row.names = FALSE, append = TRUE)
    
  }  

  
  #...................................      
  ## Function to grow random forest and produce relevant output
  f_rf <- function(f_outcome, f_model, f_df, f_effect_modifiers, f_filename) {
    
    # Grow random forest
    print("now growing forest...")
    fit <- randomForest(as.formula(paste(f_outcome, "~", ".", sep = "")), data = f_df, na.action = na.fail, localImp = TRUE)
    
    # Calculate R-squared
    rsq <- 1 - sum((f_df[, f_outcome] - predict(fit) )^2) / sum((f_df[, f_outcome] - mean(f_df[, f_outcome]))^2)

    # Plot of prediction accuracy
    x1 <- as.data.frame(cbind(f_df[, f_outcome], predict(fit)))
    colnames(x1) <- c("observations", "predictions")
    plot <- ggplot(x1) +
      geom_point(aes(x = observations, y = predictions), size=2, colour = brewer_pal(palette = "Dark2")(2)[1]) + 
      theme_bw() +
      scale_x_continuous("observed") +
      scale_y_continuous("predicted") +  
      geom_abline(intercept = 0, slope = 1, colour = brewer_pal(palette = "Dark2")(2)[2] ) +
      theme(axis.title = element_text(colour="grey20")) +
      ggtitle(paste("accuracy of random forest to predict ", all.vars(formula(fit))[1] , sep="") ) +
      theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") )
    print(plot)
        
    # Generate variable importance metrics
    print("now evaluating variable importance...")
    out <- measure_importance(fit)
      # evaluate importance of candidate interactions
      print("now evaluating interactions...")
      out_em <- min_depth_interactions(fit, vars = f_effect_modifiers[f_effect_modifiers$model == f_model, "effect_modifier"] )
      # select only interactions of interest
      x4 <- c()
      for (i in 1:nrow(out_em) ) {
        for (j in 1:nrow(f_effect_modifiers[f_effect_modifiers$model == f_model, ]) ) {
          if( length(intersect(
              c(out_em[i, "variable"], as.character(out_em[i, "root_variable"])) ,
              c(subset(f_effect_modifiers, model == f_model)[j, c("exposure", "effect_modifier")])
            )  
          ) == 2) {x4 <- c(x4, i )}
        }
      }
      out_em <- out_em[x4, ]

    # Write output to file
    write.table(out, f_filename, sep = ",", row.names = FALSE)
    write.table("--------------------", f_filename, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table(out_em, f_filename, sep = ",", row.names = FALSE, append = TRUE)
    write.table(c("--------------------", "r-squared:", rsq), f_filename, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    print("done")
    
    # Return random forest
    return(fit)
  }
    


  #...................................      
  ## Function to prepare data for random forest
  f_rf_prep <- function(f_outcome, f_model, f_vars, f_df, f_threshold_compl, 
    f_collinear_vars, f_impute, f_region) {
    # Select data
      # identify variables
      x1 <- f_vars[! is.na(f_vars[, f_model]), "variable"]
        # select only continuous variables
        x1 <- x1[-grep("_cat", x1)]
        # select only variables whose completeness is at least the threshold value
        x1 <- f_vars[f_vars$variable %in% x1 & f_vars$percent_w_data >= f_threshold_compl, "variable"]
        # exclude variables that are causing a lot of collinearity, and don't help to identify clear biological mechanisms for associations
        x1 <- x1[! x1 %in% f_collinear_vars]
  
    # If analysis is only being done for a given WHO region, select only data from this region
    if (! is.na(f_region) ) {f_df <- subset(f_df, region == f_region)}
            
    # Select complete cases or perform imputation 
    if (f_impute == FALSE) {
      # select only complete cases
      x2 <- f_df[complete.cases(f_df[, c(f_outcome, x1)]), c(f_outcome, x1)]
    }
    
    if (f_impute == TRUE) {
      # select cases for which the outcome is complete
      x3 <- f_df[complete.cases(f_df[, f_outcome]), c(f_outcome, x1)]
      # perform imputation by proximity method for missing independent variable values
      x2 <- rfImpute(as.formula(paste(f_outcome, "~", ".", sep = "")), x3, iter = 50)
    }
  
  # Output prepared dataframe            
  return(x2)
  }

                
  #...................................      
  ## Function to do univariate analysis (plots and statistical associations) for a given outcome
  f_univar <- function(f_outcome, f_model, f_vars, f_df) {
    # Continuous variables
      # identify continuous versions of independent variables
      x1 <- f_vars[f_vars[, f_model] %in% c("exposure", "confounder"), "variable"]
      x1 <- x1[- grep("_cat", x1)]
      
      # name variables (for labelling)
      names(x1) <- f_vars[f_vars$variable %in% x1, "variable_long"]
      
      # visualise correlations
        # define colours for WHO regions
        cols <- brewer.pal(length(unique(df$region)), "Set1")
        names(cols) <- sort(unique(df$region))

        # create scatter plots
        for (i in 1:length(x1) ) {
          x2 <- f_df[, c(f_outcome, x1[i], "region")]
          colnames(x2) <- c("outcome", "x_var", "region")
          if(grepl("ln", f_outcome) ) {x2[, "outcome"] <- exp(x2[, "outcome"]) }
          assign(paste("plot_", i, sep = ""), ggplot(x2, aes(y = outcome, x = x_var, colour = region) ) +
            geom_point(size = 2, alpha = 0.5) +
            theme_bw() +
            labs(title = names(x1[i])) +
            theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
              plot.title = element_text(hjust = 0.5, size = 10)) +
            scale_colour_manual(values = cols) +
            scale_x_continuous(limits = c(NA, quantile(x2$x_var, 0.95, na.rm = TRUE)))
          )
        }
      
        # identify plot names
        x2 <- paste("plot_", c(1:length(x1)), sep = "" )
      
        # arrange plots in a grid and save
        plot <- ggarrange(plotlist = mget(x2), ncol = 3, nrow = ceiling(length(x1) / 3),
          common.legend = TRUE, legend = "right", align = "v")
        plot
        ggsave(paste(f_outcome, "_univar_cont.png", sep=""), height = 33, width = 25, units = "cm", dpi = "print")
          
      # compute univariate associations in an OLS model
        # set up output
        out <- data.frame(matrix(NA, nrow = length(x1), ncol = 7))  
        colnames(out) <- c("variable", "variable_long", "n_observations", "coefficient", "95_CI", "p_value", "r_squared")  
        
        # fit all univariate models and extract statistics
        for (i in 1:length(x1) ) {
            fit <- lm(as.formula(paste(f_outcome, " ~ ", x1[i], sep = "")), data = f_df)
            out[i, "variable"] <- x1[i]
            out[i, "variable_long"] <- names(x1[i])
            out[i, "n_observations"] <- nrow(f_df[! is.na(f_df[, f_outcome]) & ! is.na(f_df[, x1[i]]), ])
            
            x2 <- coef(summary(fit))

            out[i, "coefficient"] <- round(x2[2, 1], digits = 3)
            out[i, "95_CI"] <- paste(round(x2[2, 1] - 1.96*x2[2, 2], digits = 3), 
              " to ", round(x2[2, 1] + 1.96*x2[2, 2], digits = 3), sep = "")
            out[i, "p_value"] <- round(x2[2,4], digits = 3)
            out[i, "r_squared"] <- round(summary(fit)$r.squared, digits = 2)
            
        }  
        
        # save output
        write.csv(out, paste("out_", f_outcome, "_univar_cont.csv", sep = ""), row.names = FALSE)
      
    # Univariate models with categorical variables
      # identify categorical versions of independent variables
      x1 <- f_vars[f_vars[, f_model] %in% c("exposure", "confounder"), "variable"]
      x1 <- x1[grep("_cat", x1)]
      
      # name variables (for labelling)
      names(x1) <- f_vars[f_vars$variable %in% x1, "variable_long"]
      
      # visualise correlations
        # create box plots
        for (i in 1:length(x1) ) {
          x2 <- f_df[, c(f_outcome, x1[i], "region")]
          colnames(x2) <- c("outcome", "x_var", "region")
          if(grepl("ln", f_outcome) ) {x2[, "outcome"] <- exp(x2[, "outcome"]) }
          assign(paste("plot_", i, sep = ""), ggplot(x2, aes(y = outcome, x = x_var) ) +
            geom_boxplot(lwd = 0.5, fill = palette_cb[sample.int(length(palette_cb), 1)], alpha = 0.5 ) +
            theme_bw() +
            labs(title = names(x1[i])) +
            theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
              axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(hjust = 0.5, size = 10))
          )
        }
      
        # identify plot names
        x2 <- paste("plot_", c(1:length(x1)), sep = "" )
      
        # arrange plots in a grid and save
        plot <- ggarrange(plotlist = mget(x2), ncol = 3, nrow = ceiling(length(x1) / 3),
          common.legend = TRUE, legend = "right", align = "hv")
        plot
        ggsave(paste(f_outcome, "_univar_cat.png", sep=""), height = 33, width = 25, units = "cm", dpi = "print")
          
      # compute univariate associations in an OLS model
        # set up output
        out <- data.frame(matrix(NA, nrow = sum(sapply(f_df[, x1], nlevels) ), ncol = 8))
        colnames(out) <- c("variable", "variable_long", "n_observations", "category", 
          "coefficient", "95_CI", "p_value", "r_squared")  
        
        # fit all univariate models and extract statistics
        x2 <- 1
        for (i in 1:length(x1) ) {
            fit <- lm(as.formula(paste(f_outcome, " ~ ", x1[i], sep = "")), data = f_df)

            x3 <- c( x2:(x2 + nlevels(f_df[, x1[i] ]) - 1) )
            out[x3, "variable"] <- x1[i]
            out[x3, "variable_long"] <- names(x1[i])
            out[x3, "n_observations"] <- table(f_df[! is.na(f_df[, f_outcome]) & ! is.na(f_df[, x1[i]]), x1[i] ])
            out[x3, "category"] <- levels(f_df[, x1[i] ])
            
            x4 <- as.data.frame(coef(summary(fit)))
            colnames(x4) <- c("coefficient", "se", "t", "p_value")
            x4 <- x4[-1, ]
            x4[, "category"] <- gsub(x1[i], "", row.names(x4) )
            
            x4[, "95_CI"] <- paste(round(x4[, "coefficient"] - 1.96 * x4[, "se"], digits = 3), " to ",
              round(x4[, "coefficient"] + 1.96 * x4[, "se"], digits = 3), sep = "")
            x4[, c("coefficient", "p_value")] <- round(x4[, c("coefficient", "p_value")], digits = 3)
            
            x5 <- which(out[, "category"] %in% x4$category & out[, "variable"] == x1[i])
            out[x5, c("coefficient", "95_CI", "p_value")] <- x4[, c("coefficient", "95_CI", "p_value")]
            x5 <- which(out[, "category"] == levels(f_df[, x1[i]])[1] & out[, "variable"] == x1[i])
            out[x5, c("coefficient", "95_CI", "p_value")] <- "ref."
            out[x3, "r_squared"] <- rep(round(summary(fit)$r.squared, digits = 2), nlevels(f_df[, x1[i] ]) )
            
            x2 <- x2 + nlevels(f_df[, x1[i] ])
        }  
        
        # save output
        write.csv(out[, c("variable", "variable_long", "category", "n_observations", 
          "coefficient", "95_CI", "p_value", "r_squared")], 
          paste("out_", f_outcome, "_univar_cat.csv", sep = ""), row.names = FALSE)
  }    
    

      
                  
#.........................................................................................
### Reading in required files
#.........................................................................................
    
  #...................................      
  ## Variable dictionary
  dict <- read_excel("covid_eco_analysis.xlsx", sheet = "dictionary")
    # remove tibble nonsense
    dict <- as.data.frame(dict)

  #...................................      
  ## Read in all the datasets needed
    # which datasets
    dfs <- as.data.frame(table(dict[, c("worksheet", "use")]) )
    dfs <- subset(dfs, use == "Y" & Freq > 0)[, "worksheet"]
    
    # for each dataset...
    for (i in dfs) {
      # read in
      assign(i, read_excel("covid_eco_analysis.xlsx", sheet = i) )
        # remove tibble nonsense
        assign(i, as.data.frame(get(i)) )
      # only keep needed columns
      x1 <- dict[dict$worksheet == i & dict$use == "Y", "variable"]
      x2 <- get(i)[, x1]
      assign(i, x2)
    }
    

#.........................................................................................                            
### Preparing data for analysis
#.........................................................................................    
    
  #...................................      
  ## Country age distributions

    # Retain only columns and observations needed, and format these as needed
      # set aside world population for age standardisation
      age_dist_world <- age_dist[age_dist$type == "World", ]
      age_dist_world <- age_dist_world[, colnames(age_dist_world) %in% paste("age", c(0:100), sep = "_") ]
         
      # all countries
      age_dist <- age_dist[age_dist$type == "Country/Area", ]
      age_dist <- age_dist[, ! colnames(age_dist) %in% c("type")]
        # format as needed
        colnames(age_dist) <- c("country", paste("age", c(0:100), sep = "_"))
        age_dist[, "country_iso"] <- countryname(age_dist[, "country"], "iso3c")
        age_dist <- subset(age_dist, ! is.na(country_iso) )
        age_dist <- age_dist[, c("country", "country_iso", paste("age", c(0:100), sep = "_"))]
        # calculate total population
        age_dist[, "pop_total"] <- rowSums(age_dist[, paste("age", c(0:100), sep = "_")])
        # calculate mean age
        age_dist[, "age_mean"] <- apply(age_dist[, paste("age", c(0:100), sep = "_")], 1, 
          function(x) {weighted.mean(c(0:100), x)} )
        
  #...................................      
  ## Country COVID data
    # Add country codes
    covid_data[, "country_iso"] <- countryname(covid_data[, "country"], "iso3c")
        
    # Separate out case/death age data and keep for later analysis
     # cases
      cases <- subset(covid_data, is.na(age_min) == FALSE & is.na(n_cases) == FALSE)
      cases <- cases[, c("region", "country", "country_iso", "n_cases", "age_min", "age_max")]
        
      # deaths  
      deaths <- subset(covid_data, is.na(age_min) == FALSE & is.na(n_deaths) == FALSE)
      deaths <- deaths[, c("region", "country", "country_iso", "n_deaths", "age_min", "age_max")] 
    
    # Separate out all-age data
      covid_data <- subset(covid_data, is.na(age_min))
      # make sure all values are integers
      covid_data[, "n_cases"] <- as.integer(covid_data[, "n_cases"])
      covid_data[, "n_deaths"] <- as.integer(covid_data[, "n_deaths"])

  #...................................      
  ## Average Rt values per country
    rt_imperial <- aggregate(rt_imperial[, "rt_median"], by=list(country = rt_imperial$country), FUN = mean)
    colnames(rt_imperial) <- c("country", "rt_mean")

  #...................................      
  ## Other dataset-specific changes
    # names
    colnames(pf_malaria) <- c("country", "prev_pf")  
    colnames(pv_malaria) <- c("country", "prev_pv")  
    colnames(helminths) <- c("country", "prev_helminths")  
    colnames(schisto) <- c("country", "prev_schisto")  
    colnames(filaria) <- c("country", "prev_filaria")
    
    # malaria: NAs mean zero
    pf_malaria[is.na(pf_malaria$prev_pf), "prev_pf"] <- 0
    pv_malaria[is.na(pv_malaria$prev_pv), "prev_pv"] <- 0

    # express proportions or percents as 0-100
    pf_malaria[, "prev_pf"] <- pf_malaria[, "prev_pf"] * 100
    pv_malaria[, "prev_pv"] <- pv_malaria[, "prev_pv"] * 100
    helminths[, "prev_helminths"] <- helminths[, "prev_helminths"] * 100
    schisto[, "prev_schisto"] <- schisto[, "prev_schisto"] * 100
    filaria[, "prev_filaria"] <- filaria[, "prev_filaria"] * 100
    comorbidities[, "percent_increased_risk_age_std"] <- comorbidities[, "percent_increased_risk_age_std"] * 100
    
  #...................................      
  ## Add country codes to all datasets
    for (i in dfs) {
      x1 <- get(i)
      x1[, "country_iso"] <- countryname(x1[, "country"], "iso3c")
      assign(i, x1)
    }
            
      
  #...................................      
  ## Merge datasets into one
    # Start with age distribution dataset, as most complete in terms of geographical entities
    df <- age_dist[, c("country", "country_iso", "pop_total", "age_mean")]
    
    # Merge in other datasets
    for (i in dfs[dfs != "age_dist"]) {
      x1 <- get(i)
      df <- merge(df, x1[, ! colnames(x1) %in% c("country")], by = "country_iso", all.x = TRUE )
    }
    
    # Attribute missing WHO regions
    df[df$country %in% c("Aruba", "Curaçao", "Guadeloupe", "French Guiana", "Martinique",  "Puerto Rico"), "region"] <- "PAHO"
    df[df$country %in% c("Guam", "New Caledonia", "French Polynesia", "China, Macao SAR", "China, Taiwan Province of China") , "region"] <- "WPRO"
    df[df$country %in% c("Réunion", "Mayotte" , "Western Sahara"), "region"] <- "AFRO"
                        

#..........................................................................................
### Managing and age-standardising case/death age data, by country
#..........................................................................................
    
  #...................................      
  ## Prepare data
    
    # Force maximum age_max to be 100
    cases <- cases[order(cases[, "country"], cases[, "age_max"]), ]
    x1 <- cumsum(table(cases$country))
    cases[x1, "age_max"] <- 100
    
    deaths <- deaths[order(deaths[, "country"], deaths[, "age_max"]), ]
    x1 <- cumsum(table(deaths$country))
    deaths[x1, "age_max"] <- 100
    
    # Make sure that age_max is always 1 less than next category's age_min
    for (i in 1:nrow(cases) ) {
      if (cases[i, "age_max"] < 100) { cases[i, "age_max"] <- cases[i+1, "age_min"] - 1 }     
    }  

    for (i in 1:nrow(deaths) ) {
      if (deaths[i, "age_max"] < 100) { deaths[i, "age_max"] <- deaths[i+1, "age_min"] - 1 }     
    }  
            
    # Calculate population within case/deaths age categories
    cases[, "pop"] <- NA
    for (i in 1:nrow(cases)) {
      x1 <- paste("age", c(cases[i, "age_min"] : cases[i, "age_max"]), sep = "_" )
      cases[i, "pop"] <- sum(age_dist[age_dist$country_iso == cases[i, "country_iso"], x1])
    }

    deaths[, "pop"] <- NA
    for (i in 1:nrow(deaths)) {
      x1 <- paste("age", c(deaths[i, "age_min"] : deaths[i, "age_max"]), sep = "_" )
      deaths[i, "pop"] <- sum(age_dist[age_dist$country_iso == deaths[i, "country_iso"], x1])
    }
    
    # Make sure cases and deaths are all numeric - if datum is "< x" assume mid-point between 0 and x
    x1 <- grep("<", cases$n_cases)
    if (length(x1) > 0) { cases[x1, "n_cases"] <- ceiling(as.numeric(gsub("[^0-9.-]", "", cases[x1, "n_cases"]) ) / 2) }
    cases[, "n_cases"] <- as.integer(cases[, "n_cases"])

    x1 <- grep("<", deaths$n_deaths)
    if (length(x1) > 0) { deaths[x1, "n_deaths"] <- ceiling(as.numeric(gsub("[^0-9.-]", "", deaths[x1, "n_deaths"]) ) / 2) }
    deaths[, "n_deaths"] <- as.integer(deaths[, "n_deaths"])
    
    # Compute per capita incidence/deaths within each age category (per 1000)
    cases[, "cases_rate"] <- cases[, "n_cases"] / cases[, "pop"]    
    deaths[, "deaths_rate"] <- deaths[, "n_deaths"] / deaths[, "pop"]    

              
  #...................................      
  ## Calculate age-standardised median age of cases and deaths (assuming world's population distribution)
    
    # Compute the world's population within each age category
      # cases
      cases[, "pop_world"] <- NA
      for (i in 1:nrow(cases)) {
        x1 <- paste("age", c(cases[i, "age_min"] : cases[i, "age_max"]), sep = "_" )
        cases[i, "pop_world"] <- sum(unlist(age_dist_world[x1]) )
      }
      
      # deaths  
      deaths[, "pop_world"] <- NA
      for (i in 1:nrow(deaths)) {
        x1 <- paste("age", c(deaths[i, "age_min"] : deaths[i, "age_max"]), sep = "_" )
        deaths[i, "pop_world"] <- sum(unlist(age_dist_world[x1]) )
      }
      
    # Compute age-standardised cases and deaths
      cases[, "n_cases_world"] <- cases[, "cases_rate"] * cases[, "pop_world"]
      deaths[, "n_deaths_world"] <- deaths[, "deaths_rate"] * deaths[, "pop_world"]
      
    # Compute median age of cases and deaths from linear interpolation of cumulative age-standardised cases / deaths
      # cases
      cases <- cases[order(cases[, "country"], cases[, "age_max"]), ]
      out_cases <- as.data.frame(matrix(NA, ncol = 5, nrow = length(unique(cases$country)) ) )
      colnames(out_cases) <- c("country", "country_iso", "median_age_cases", "quartile_age_cases_lower", "quartile_age_cases_upper")
      out_cases[, c("country", "country_iso")] <- unique(cases[, c("country", "country_iso")] )
      for (i in out_cases$country ) {
        # select country data
        x1 <- subset(cases, country == i)[, c("age_max", "n_cases_world")]
        # cumulative proportions of cases
        x1[, "cum"] <- cumsum(x1[, "n_cases_world"]) / sum(x1[, "n_cases_world"])
          # add 0 starting values
          x1 <- rbind(c(-1, 0, 0), x1)
        # linear interpolation across all ages from 0 to 100
        x2 <- approx(x1[, c("age_max", "cum")], xout = c(-1:100) )
        x2 <- as.data.frame(x2)[-1, ]
        # median and IQR of age
        x3 <- c(which.min(abs(0.50 - x2$y) ) , which.min(abs(0.25 - x2$y) ) , which.min(abs(0.75 - x2$y) ) )
        out_cases[out_cases$country == i, c("median_age_cases", "quartile_age_cases_lower", "quartile_age_cases_upper")] <- x2$x[x3]
      }
      
      out_cases[, "cases_table_format"] <- paste(out_cases[, "median_age_cases"], " (", out_cases[, "quartile_age_cases_lower"], " to ",
        out_cases[, "quartile_age_cases_upper"], ")", sep ="")
 
      # deaths
      deaths <- deaths[order(deaths[, "country"], deaths[, "age_max"]), ]
      out_deaths <- as.data.frame(matrix(NA, ncol = 5, nrow = length(unique(deaths$country)) ) )
      colnames(out_deaths) <- c("country", "country_iso", "median_age_deaths", "quartile_age_deaths_lower", "quartile_age_deaths_upper")
      out_deaths[, c("country", "country_iso")] <- unique(deaths[, c("country", "country_iso")] )
      for (i in out_deaths$country ) {
        # select country data
        x1 <- subset(deaths, country == i)[, c("age_max", "n_deaths_world")]
        # cumulative proportions of deaths
        x1[, "cum"] <- cumsum(x1[, "n_deaths_world"]) / sum(x1[, "n_deaths_world"])
          # add 0 starting values
          x1 <- rbind(c(-1, 0, 0), x1)
        # linear interpolation across all ages from 0 to 100
        x2 <- approx(x1[, c("age_max", "cum")], xout = c(-1:100) )
        x2 <- as.data.frame(x2)[-1, ]
        # median and IQR of age
        x3 <- c(which.min(abs(0.50 - x2$y) ) , which.min(abs(0.25 - x2$y) ) , which.min(abs(0.75 - x2$y) ) )
        out_deaths[out_deaths$country == i, c("median_age_deaths", "quartile_age_deaths_lower", "quartile_age_deaths_upper")] <- x2$x[x3]
      }
      
      out_deaths[, "deaths_table_format"] <- paste(out_deaths[, "median_age_deaths"], " (", out_deaths[, "quartile_age_deaths_lower"], " to ",
        out_deaths[, "quartile_age_deaths_upper"], ")", sep ="")
    
    # Write output and merge into main database    
    out <- merge(out_cases, out_deaths, by= c("country_iso", "country"), all= TRUE)
    write.csv(out, "out_median_age_std.csv", row.names = FALSE)
    
    df <- merge(df, out[, c("country_iso", "median_age_cases", "median_age_deaths")], by = "country_iso", all.x = TRUE)
    
          
  #...................................      
  ## Calculate crude and age-standardised case-fatality ratio (CFR)
    
    # Crude
    df[, "cfr_crude"] <- df[, "n_deaths"] / df[, "n_cases"]
      # add 0.0001 to CFR of countries with CFR = 0
      df[, "cfr_crude"] <- ifelse(df[, "cfr_crude"] == 0, 0.0001, df[, "cfr_crude"]) 
      # eliminate CFRs based on < 50 cases
      df[, "cfr_crude"] <- ifelse(df[, "n_cases"] < 50, NA, df[, "cfr_crude"]) 
    
    # Age-standardised - version 1 (age-specific cases if country X had world's population structure x age-specific CFRs in country X)
      # sort cases and deaths, and merge to select countries with matching case/death age categories
      cases <- cases[order(cases[, "country"], cases[, "age_max"]), ]
      deaths <- deaths[order(deaths[, "country"], deaths[, "age_max"]), ]
      cfr <- merge(cases, deaths, by = c("country_iso", "country", "age_min", "age_max"))
      cfr <- cfr[, c("country", "country_iso", "age_max", "n_deaths", "n_cases", "n_cases_world")]
      
      # calculate age-specific CFRs
      cfr[, "cfr_age"] <- cfr[, "n_deaths"] / cfr[, "n_cases"]
      
      # calculate age-standardised CFR
      cfr[, "n_deaths_world"] <- cfr[, "n_cases_world"] * cfr[, "cfr_age"]
      cfr <- aggregate(cfr[, c("n_deaths_world", "n_cases_world")], by = list(country_iso = cfr$country_iso), FUN = sum)
      cfr[, "cfr_std1"] <- cfr[, "n_deaths_world"] / cfr[, "n_cases_world"]
      
      # merge into the main database
      cfr <- merge(cfr, unique(df[, c("country_iso", "country")]), by = "country_iso", all.x = TRUE)
      df <- merge(df, cfr[, c("country_iso", "cfr_std1")], by = "country_iso", all.x = TRUE)
 
    # Age-standardised - version 2 (basis: age-specific cases in South Korea x observed age-specific CFRs from country X)
      # sort cases and deaths, and merge to select countries with matching case/death age categories
      cases <- cases[order(cases[, "country"], cases[, "age_max"]), ]
      deaths <- deaths[order(deaths[, "country"], deaths[, "age_max"]), ]
      cfr <- merge(cases, deaths, by = c("country_iso", "country", "age_min", "age_max"))
      cfr <- cfr[, c("country", "country_iso", "age_max", "n_deaths", "n_cases", "n_cases_world")]
      
      # calculate age-specific CFRs
      cfr[, "cfr_age"] <- cfr[, "n_deaths"] / cfr[, "n_cases"]

      # estimate cumulative cases by year cohort in South Korea, using linear interpolation
        # select S. Korea data
        x1 <- subset(cases, country_iso == "KOR")[, c("age_max", "n_cases")]
        # cumulative numbers of cases
        x1[, "cum"] <- cumsum(x1[, "n_cases"])
          # add 0 starting values
          x1 <- rbind(c(-1, 0, 0), x1)
        # linear interpolation across all ages from 0 to 100
        x2 <- approx(x1[, c("age_max", "cum")], xout = c(-1:100) )
        x2 <- as.data.frame(x2)[-1, ]
        colnames(x2) <- c("age_max", "n_cases_kor_cum")

      # for each country's age category, find S Korea's caseload
      cfr <- merge(cfr, x2, by = "age_max")
      cfr <- cfr[order(cfr[, "country"], cfr[, "age_max"]), ]
      cfr[, "n_cases_kor"] <- NA
      for (i in cfr$country) {
        x1 <- subset(cfr, country == i)[, "n_cases_kor_cum"]
        x2 <- c(x1[1], diff(x1))
        cfr[cfr$country == i, "n_cases_kor"] <- x2
      }
      
      # calculate CFR for each country based on S. Korea caseload
      cfr[, "n_deaths_kor"] <- cfr[, "cfr_age"] * cfr[, "n_cases_kor"]
      cfr <- aggregate(cfr[, c("n_deaths_kor", "n_cases_kor")], by = list(country_iso = cfr$country_iso), FUN = sum)
      cfr[, "cfr_std2"] <- cfr[, "n_deaths_kor"] / cfr[, "n_cases_kor"]
      
      # also calculate CFR ratio (country vs S Korea)
      cfr[, "cfr_std2_ratio"] <- cfr[, "cfr_std2"] / cfr[cfr$country_iso == "KOR", "cfr_std2"]
      
      # merge into the main database
      cfr <- merge(cfr, unique(df[, c("country_iso", "country")]), by = "country_iso", all.x = TRUE)
      df <- merge(df, cfr[, c("country_iso", "cfr_std2", "cfr_std2_ratio")], by = "country_iso", all.x = TRUE)
    
    # Save all CFR estimates
    x1 <- df[, c("country", "country_iso", "cfr_crude", "cfr_std1", "cfr_std2", "cfr_std2_ratio")]
    x1[, c("cfr_crude", "cfr_std1", "cfr_std2")] <- round(x1[, c("cfr_crude", "cfr_std1", "cfr_std2")] * 100, digits = 1)
    x1[, "cfr_std2_ratio"] <- round(x1[, "cfr_std2_ratio"], digits = 2)
    write.csv(x1, "out_cfr.csv", row.names = FALSE)
 
         
#..........................................................................................
### Preparing data for statistical analysis
#..........................................................................................

  #...................................      
  ## Create additional variables and modify others
    # Testing rate per 1000
    df[, "test_rate"] <- df[, "n_tests"] / df[, "pop_total"]
      
    # Region as factor
    df[, "region"] <- as.factor(df[, "region"])
    
  #...................................      
  ## Apply transformations to the outcomes, if appropriate
    # Mean R_t
    f_hist("rt_mean", df, c(NA,NA))
      # no transformation
    
    # Median age of cases (age-standardised)
    f_hist("median_age_cases", df, c(NA,NA))
      # no transformation
    
    # Median age of deaths (age-standardised)
    f_hist("median_age_deaths", df, c(NA,NA))
      # no transformation
    
    # Crude CFR
    f_hist("cfr_crude", df, c(NA,NA))
      # natural log
      df[, "cfr_crude_ln"] <- log(df[, "cfr_crude"])
      f_hist("cfr_crude_ln", df, c(NA,NA))
        # set outliers to NA
        df[! is.na(df$cfr_crude_ln) & df$cfr_crude_ln < (-8), "cfr_crude_ln"] <- NA
      
    # Age-standardised CFR 1
    f_hist("cfr_std1", df, c(NA,NA))
      # natural log
      df[, "cfr_std1_ln"] <- log(df[, "cfr_std1"])
      f_hist("cfr_std1_ln", df, c(NA,NA))

    # Age-standardised CFR 2
    f_hist("cfr_std2", df, c(NA,NA))
      # no transformation
      

  #...................................      
  ## Categorise independent variables, as an alternative to their continuous form, and set reference categories
    
    # Mean age of population
    f_hist("age_mean", df, c(NA, NA))
    df[, "age_mean_cat"] <- cut(df[, "age_mean"], c(0, 25, 35, 100), 
      labels = c("<25 yo", "25 to 34.9yo", ">= 35yo"), include.lowest = TRUE)
    table(df$age_mean_cat)
    df[, "age_mean_cat"] <- relevel(df[, "age_mean_cat"], "<25 yo" )
        
    # Prevalence of filaria
    f_hist("prev_filaria", df, c(NA, NA))
    f_hist("prev_filaria", df, c(NA, 1))
    df[, "prev_filaria_cat"] <- cut(df[, "prev_filaria"], c(0, 0.00001, 2.5, 100), 
      labels = c("0%", "0.1 to 2.4%", ">= 2.5%"), include.lowest = TRUE)
    table(df$prev_filaria_cat)
    df[, "prev_filaria_cat"] <- relevel(df[, "prev_filaria_cat"], "0%" )
      
    # Prevalence of helminths
    f_hist("prev_helminths", df, c(NA, NA))
    df[, "prev_helminths_cat"] <- cut(df[, "prev_helminths"], c(0, 0.00001, 10.0, 20.0, 100), 
      labels = c("0%", "0.1 to 9.9%", "10.0% to 19.9%", ">= 20.0%"), include.lowest = TRUE)
    table(df$prev_helminths_cat)
    df[, "prev_helminths_cat"] <- relevel(df[, "prev_helminths_cat"], "0%" )
    
    # Prevalence of Pf malaria
    f_hist("prev_pf", df, c(NA, NA))
    df[, "prev_pf_cat"] <- cut(df[, "prev_pf"], c(0, 0.00001, 5, 100), 
      labels = c("0%", "0.1 to 4.9%", ">= 5.0%"), include.lowest = TRUE)
    table(df$prev_pf_cat)
    df[, "prev_pf_cat"] <- relevel(df[, "prev_pf_cat"], "0%" )
    
    # Prevalence of Pv malaria
    f_hist("prev_pv", df, c(NA, NA))
    df[, "prev_pv_cat"] <- cut(df[, "prev_pv"], c(0, 0.00001, 100), 
      labels = c("0%", "> 0%"), include.lowest = TRUE)
    table(df$prev_pv_cat)
    df[, "prev_pv_cat"] <- relevel(df[, "prev_pv_cat"], "0%" )
    
    # Prevalence of schistosomiasis
    f_hist("prev_schisto", df, c(NA, NA))
    df[, "prev_schisto_cat"] <- cut(df[, "prev_schisto"], c(0, 0.00001, 5, 100), 
      labels = c("0%", "0.1 to 4.9%", ">= 5.0%"), include.lowest = TRUE)
    table(df$prev_schisto_cat)
    df[, "prev_schisto_cat"] <- relevel(df[, "prev_schisto_cat"], "0%" )
    
    # Sanitation index
    f_hist("wash_2017", df, c(NA, NA))
    df[, "wash_2017_cat"] <- cut(df[, "wash_2017"], c(0, 25, 50, 75, 1000), 
      labels = c("< 25%", "25 to 49%", "50 to 74%", ">= 75%"), include.lowest = TRUE)
    table(df$wash_2017_cat)
    df[, "wash_2017_cat"] <- relevel(df[, "wash_2017_cat"], "< 25%" )
    
    # Percent of population with increased COVID risk
    f_hist("percent_increased_risk_age_std", df, c(NA, NA))
    df[, "percent_increased_risk_age_std_cat"] <- cut(df[, "percent_increased_risk_age_std"], c(0, 15, 20, 25, 30, 1), 
      labels = c("< 15%", "15 to 19%", "20 to 24%", "25 to 29%", ">= 30%"), include.lowest = TRUE)
    table(df$percent_increased_risk_age_std_cat)
    df[, "percent_increased_risk_age_std_cat"] <- relevel(df[, "percent_increased_risk_age_std_cat"], "< 15%" )
    
    # Human Development Index  
    f_hist("hdi_2018", df, c(NA, NA))
    df[, "hdi_2018_cat"] <- cut(df[, "hdi_2018"], c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 10), 
      labels = c("< 0.5", "0.50 to 0.59", "0.60 to 0.69", "0.70 to 0.79", "0.80 to 0.89", ">= 0.90"), 
      include.lowest = TRUE)
    table(df$hdi_2018_cat)
    df[, "hdi_2018_cat"] <- relevel(df[, "hdi_2018_cat"], "< 0.5" )
    
    # Mean household size
    f_hist("hh_size", df, c(NA, NA))
    df[, "hh_size_cat"] <- cut(df[, "hh_size"], c(0, 3, 4, 5, 100), 
      labels = c("< 3.0", "3.0 to 3.9", "4.0 to 4.9", ">= 5.0"), include.lowest = TRUE)
    table(df$hh_size_cat)
    df[, "hh_size_cat"] <- relevel(df[, "hh_size_cat"], "< 3.0" )
    
    # Mean mobility change  
    f_hist("mean_mobility_change", df, c(NA, NA))
    df[, "mean_mobility_change_cat"] <- cut(df[, "mean_mobility_change"], c(-100, -30, -20, -10, 100), 
      labels = c("-30% or more", "-20 to -29%", "-10 to -19%", "-10% or less"), include.lowest = TRUE)
    table(df$mean_mobility_change_cat)
    df[, "mean_mobility_change_cat"] <- relevel(df[, "mean_mobility_change_cat"], "-30% or more" )
    
    # Population density  
    f_hist("pop_km2_2020", df, c(NA, NA))
    f_hist("pop_km2_2020", df, c(NA, 100))
    f_hist("pop_km2_2020", df, c(NA, 200))
    df[, "pop_km2_2020_cat"] <- cut(df[, "pop_km2_2020"], c(0, 25, 50, 75, 100, 1000000), 
      labels = c("< 25", "25 to 49", "50 to 74", "75 to 99", ">= 100"), include.lowest = TRUE)
    table(df$pop_km2_2020_cat)
    df[, "pop_km2_2020_cat"] <- relevel(df[, "pop_km2_2020_cat"], "< 25" )
    
    # Mean stringency of control measures  
    f_hist("mean_stringency", df, c(NA, NA))
    df[, "mean_stringency_cat"] <- cut(df[, "mean_stringency"], c(0, 40, 50, 60, 100), 
      labels = c("< 40", "40 to 49", "50 to 59", ">= 60"), include.lowest = TRUE)
    table(df$mean_stringency_cat)
    df[, "mean_stringency_cat"] <- relevel(df[, "mean_stringency_cat"], "< 40" )
    
    # Mean testing policy  
    f_hist("mean_testing_policy", df, c(NA, NA))
    df[, "mean_testing_policy_cat"] <- cut(df[, "mean_testing_policy"], c(0, 1, 1.5, 2.0, 10), 
      labels = c("< 1.0", "1.0 to 1.4", "1.5 to 1.9", ">= 2.0"), include.lowest = TRUE)
    table(df$mean_testing_policy_cat)
    df[, "mean_testing_policy_cat"] <- relevel(df[, "mean_testing_policy_cat"], "< 1.0" )
    
    # Testing rate  
    f_hist("test_rate", df, c(NA, NA))
    f_hist("test_rate", df, c(NA, 100))
    f_hist("test_rate", df, c(NA, 500))
    df[, "test_rate_cat"] <- cut(df[, "test_rate"], c(0, 25, 100, 100000), 
      labels = c("< 25 per 1000", "25 to 99 per 1000", ">= 100 per 1000"), include.lowest = TRUE)
    table(df$test_rate_cat)
    df[, "test_rate_cat"] <- relevel(df[, "test_rate_cat"], "< 25 per 1000" )
    

  #...................................      
  ## Identify variables and data types
    # Dataframe of variables
    vars <- data.frame(matrix(NA, ncol = 4, nrow = length(colnames(df)) ))
    colnames(vars) <- c("variable", "variable_long", "role", "type")
    vars[, "variable"] <- colnames(df)
    
    vars[, "variable_long"] <- c("country code", "country", "total population", "mean age", "percent with increased COVID risk",
      "WHO region", "number of cases", "number of deaths", "minimum age category", "maximum age category",
      "number of tests", "prevalence of filaria", "Human Development Index", "prevalence of helminths",
      "mean household size", "percent mobility change", "prevalence of P. falciparum (2-10 yo)", "population density per Km2", 
      "prevalence of P. vivax", "mean reproduction number R_t", "prevalence of schistosomiasis", 
      "mean stringency index", "mean testing policy", "sanitation index", "median age of cases (standardised)", 
      "median age of deaths (standardised)", "crude CFR", "standardised CFR 1", "standardised CFR 2", "ratio of standardised CSR 2",
      "testing rate per 1000", "log crude CFR", "log standardised CFR 1", "mean age", "prevalence of filaria", "prevalence of helminths",
      "prevalence of P. falciparum", "prevalence of P. vivax", "prevalence of schistosomiasis", "sanitation index",
      "percent with increased COVID risk", "Human Development Index", "mean household size", "percent mobility change", 
      "population density per Km2", "mean stringency index", "mean testing policy", "testing rate per 1000")
    
    # General roles variables play
    vars[, "role"] <- "none"
    
    vars[vars$variable %in% c("age_mean", "prev_filaria", "prev_helminths", "prev_pf", 
      "prev_pv", "prev_schisto", "wash_2017"), "role"] <- "exposure"
    vars[vars$variable %in% paste(c("age_mean", "prev_filaria", "prev_helminths", "prev_pf", 
      "prev_pv", "prev_schisto", "wash_2017"), "_cat", sep=""), "role"] <- "exposure"
    
    vars[vars$variable %in% c("percent_increased_risk_age_std", "hdi_2018", "hh_size", "mean_mobility_change", 
      "pop_km2_2020", "mean_stringency", "mean_testing_policy", "test_rate"), "role"] <- "confounder"
    vars[vars$variable %in% paste(c("percent_increased_risk_age_std", "hdi_2018", "hh_size", "mean_mobility_change", 
      "pop_km2_2020", "mean_stringency", "mean_testing_policy", "test_rate"), "_cat", sep=""), "role"] <- "confounder"
    
    vars[vars$variable %in% c("rt_mean", "median_age_deaths", "median_age_cases", 
      "cfr_crude", "cfr_std1", "cfr_std2", "cfr_crude_ln", "cfr_std1_ln"), "role"] <- "outcome"
    
    # Types of variable
    vars[vars$variable %in% c("age_mean", "percent_increased_risk_age_std", "prev_filaria", "hdi_2018", 
      "prev_helminths", "hh_size", "mean_mobility_change", "prev_pf", "pop_km2_2020", "prev_pv", 
      "rt_mean", "prev_schisto", "mean_stringency", "mean_testing_policy", "wash_2017", 
      "median_age_cases", "median_age_deaths", "test_rate", "cfr_crude", "cfr_std1", "cfr_std2", 
      "cfr_crude_ln", "cfr_std1_ln"), "type"] <- "continuous"
    vars[grep("_cat", vars$variable), "type"] <- "categorical"
    vars[vars$variable %in% c("region"), "type"] <- "categorical"

    # Which roles variable play for each model (outcome)
      # R_t
      vars[, "rt_model"] <- NA
      vars[vars$role == "exposure", "rt_model"] <- "exposure"
      vars[vars$variable %in% c("test_rate", "mean_testing_policy", "mean_mobility_change", "hdi_2018", "pop_km2_2020",
        "hh_size", "mean_stringency"), "rt_model"] <- "confounder"
      vars[vars$variable %in% paste(c("test_rate", "mean_testing_policy", "mean_mobility_change", "hdi_2018", "pop_km2_2020",
        "hh_size", "mean_stringency"), "_cat", sep=""), "rt_model"] <- "confounder"

      # median age of cases
      vars[, "age_cases_model"] <- NA
      vars[vars$role == "exposure", "age_cases_model"] <- "exposure"
      vars[vars$variable %in% c("test_rate", "mean_testing_policy", "mean_mobility_change", "hdi_2018", 
        "mean_stringency", "percent_increased_risk_age_std"), "age_cases_model"] <- "confounder"
      vars[vars$variable %in% paste(c("test_rate", "mean_testing_policy", "mean_mobility_change", "hdi_2018",
        "mean_stringency", "percent_increased_risk_age_std"), "_cat", sep=""), "age_cases_model"] <- "confounder"

      # median age of deaths
      vars[, "age_deaths_model"] <- NA
      vars[vars$role == "exposure", "age_deaths_model"] <- "exposure"
      vars[vars$variable %in% c("test_rate", "mean_testing_policy", "mean_mobility_change", "hdi_2018", 
        "mean_stringency", "percent_increased_risk_age_std"), "age_deaths_model"] <- "confounder"
      vars[vars$variable %in% paste(c("test_rate", "mean_testing_policy", "mean_mobility_change", "hdi_2018",
        "mean_stringency", "percent_increased_risk_age_std"), "_cat", sep=""), "age_deaths_model"] <- "confounder"
      
      # crude CFR
      vars[, "cfr_crude_model"] <- NA
      vars[vars$role == "exposure", "cfr_crude_model"] <- "exposure"
      vars[vars$variable %in% c("test_rate", "mean_testing_policy", "hdi_2018", 
        "mean_stringency", "percent_increased_risk_age_std"), "cfr_crude_model"] <- "confounder"
      vars[vars$variable %in% paste(c("test_rate", "mean_testing_policy", "hdi_2018",
        "mean_stringency", "percent_increased_risk_age_std"), "_cat", sep=""), "cfr_crude_model"] <- "confounder"
      
      # age-standardised CFR 1
      vars[, "cfr_std1_model"] <- NA
      vars[vars$role == "exposure", "cfr_std1_model"] <- "exposure"
      vars[vars$variable %in% c("test_rate", "mean_testing_policy", "hdi_2018", 
        "mean_stringency", "percent_increased_risk_age_std"), "cfr_std1_model"] <- "confounder"
      vars[vars$variable %in% paste(c("test_rate", "mean_testing_policy", "hdi_2018",
        "mean_stringency", "percent_increased_risk_age_std"), "_cat", sep=""), "cfr_std1_model"] <- "confounder"
      
      # age-standardised CFR 2
      vars[, "cfr_std2_model"] <- NA
      vars[vars$role == "exposure", "cfr_std2_model"] <- "exposure"
      vars[vars$variable %in% c("test_rate", "mean_testing_policy", "hdi_2018", 
        "mean_stringency", "percent_increased_risk_age_std"), "cfr_std2_model"] <- "confounder"
      vars[vars$variable %in% paste(c("test_rate", "mean_testing_policy", "hdi_2018",
        "mean_stringency", "percent_increased_risk_age_std"), "_cat", sep=""), "cfr_std2_model"] <- "confounder"
      
    # Which variables are plausible effect modifiers for the exposures in each model? (outcome)
    effect_modifiers <- rbind(
          c("rt_model", "age_mean", "test_rate"),
          c("rt_model", "age_mean", "mean_testing_policy"),
          c("age_cases_model", "age_mean", "percent_increased_risk_age_std"),
          c("age_cases_model", "age_mean", "mean_stringency"),
          c("age_cases_model", "age_mean", "mean_mobility_change"),
          c("age_deaths_model", "age_mean", "percent_increased_risk_age_std"),
          c("age_deaths_model", "age_mean", "mean_stringency"),
          c("age_deaths_model", "age_mean", "mean_mobility_change"),
          c("cfr_crude_model", "age_mean", "mean_stringency"),
          c("cfr_crude_model", "age_mean", "mean_mobility_change"),
          c("cfr_crude_model", "age_mean", "percent_increased_risk_age_std"),
          c("cfr_std1_model", "age_mean", "mean_stringency"),
          c("cfr_std1_model", "age_mean", "mean_mobility_change"),
          c("cfr_std1_model", "age_mean", "percent_increased_risk_age_std"),
          c("cfr_std2_model", "age_mean", "mean_stringency"),
          c("cfr_std2_model", "age_mean", "mean_mobility_change"),
          c("cfr_std2_model", "age_mean", "percent_increased_risk_age_std")
      )  
    effect_modifiers <- as.data.frame(effect_modifiers)    
    colnames(effect_modifiers) <- c("model", "exposure", "effect_modifier") 
      
    # Completeness
    vars[, "n_w_data"] <- NA
    vars[, "percent_w_data"] <- NA
    for (i in 1:nrow(vars) ) {
      vars[i, "n_w_data"] <- length(which(complete.cases(df[, vars[i, "variable"]]) ) )
      vars[i, "percent_w_data"] <- round(vars[i, "n_w_data"] * 100 / nrow(df), digits = 1)
    }
    
    # Save variables and data
    write.csv(df, "out_data.csv", row.names = FALSE)
    write.csv(vars, "out_variables_for_analysis.csv", row.names = FALSE)
      
  #...................................      
  ## Impute missing predictor values (above a certain threshold of completeness): to be used for linear regression
  
    # Select continuous versions of independent variables that have at least the desired % completeness
    x1 <- vars[vars$percent_w_data >= threshold_compl & vars$role %in% c("exposure", "confounder"), "variable"]  
    x1 <- x1[! grepl("_cat", x1)]
    x1 <- c("country_iso", x1)
    
    # Perform imputation and merge back into database
    x2 <- mice(df[, x1], m = 20)
    x2 <- complete(x2)
    df_imputed <- merge(df[, vars[vars$role %in% c("none", "outcome"), "variable"]], x2, by = "country_iso", all.x = TRUE)

    # Recategorise independent variables
      # mean age of population
      df_imputed[, "age_mean_cat"] <- cut(df_imputed[, "age_mean"], c(0, 25, 35, 100), 
        labels = c("<25 yo", "25 to 34.9yo", ">= 35yo"), include.lowest = TRUE)
      df_imputed[, "age_mean_cat"] <- relevel(df_imputed[, "age_mean_cat"], "<25 yo" )
          
      # prevalence of filaria
      df_imputed[, "prev_filaria_cat"] <- cut(df_imputed[, "prev_filaria"], c(0, 0.00001, 2.5, 100), 
        labels = c("0%", "0.1 to 2.4%", ">= 2.5%"), include.lowest = TRUE)
      df_imputed[, "prev_filaria_cat"] <- relevel(df_imputed[, "prev_filaria_cat"], "0%" )
        
      # prevalence of helminths
      df_imputed[, "prev_helminths_cat"] <- cut(df_imputed[, "prev_helminths"], c(0, 0.00001, 10, 20, 100), 
        labels = c("0%", "0.1 to 9.9%", "10.0% to 19.9%", ">= 20.0%"), include.lowest = TRUE)
      df_imputed[, "prev_helminths_cat"] <- relevel(df_imputed[, "prev_helminths_cat"], "0%" )
      
      # prevalence of Pf malaria
      df_imputed[, "prev_pf_cat"] <- cut(df_imputed[, "prev_pf"], c(0, 0.00001, 5, 100), 
        labels = c("0%", "0.1 to 4.9%", ">= 5.0%"), include.lowest = TRUE)
      df_imputed[, "prev_pf_cat"] <- relevel(df_imputed[, "prev_pf_cat"], "0%" )
      
      # prevalence of Pv malaria
      df_imputed[, "prev_pv_cat"] <- cut(df_imputed[, "prev_pv"], c(0, 0.00001, 100), 
        labels = c("0%", "> 0%"), include.lowest = TRUE)
      df_imputed[, "prev_pv_cat"] <- relevel(df_imputed[, "prev_pv_cat"], "0%" )
      
      # prevalence of schistosomiasis
      df_imputed[, "prev_schisto_cat"] <- cut(df_imputed[, "prev_schisto"], c(0, 0.00001, 5, 100), 
        labels = c("0%", "0.1 to 4.9%", ">= 5.0%"), include.lowest = TRUE)
      df_imputed[, "prev_schisto_cat"] <- relevel(df_imputed[, "prev_schisto_cat"], "0%" )
      
      # sanitation index
      df_imputed[, "wash_2017_cat"] <- cut(df_imputed[, "wash_2017"], c(0, 25, 50, 75, 1000), 
        labels = c("< 25%", "25 to 49%", "50 to 74%", ">= 75%"), include.lowest = TRUE)
      df_imputed[, "wash_2017_cat"] <- relevel(df_imputed[, "wash_2017_cat"], "< 25%" )
      
      # percent of population with increased COVID risk
      df_imputed[, "percent_increased_risk_age_std_cat"] <- cut(df_imputed[, "percent_increased_risk_age_std"], c(0, 15, 20, 25, 30, 1), 
        labels = c("< 15%", "15 to 19%", "20 to 24%", "25 to 29%", ">= 30%"), include.lowest = TRUE)
      df_imputed[, "percent_increased_risk_age_std_cat"] <- relevel(df_imputed[, "percent_increased_risk_age_std_cat"], "< 15%" )
      
      # Human Development Index  
      df_imputed[, "hdi_2018_cat"] <- cut(df_imputed[, "hdi_2018"], c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 10), 
        labels = c("< 0.5", "0.50 to 0.59", "0.60 to 0.69", "0.70 to 0.79", "0.80 to 0.89", ">= 0.90"), 
        include.lowest = TRUE)
      df_imputed[, "hdi_2018_cat"] <- relevel(df_imputed[, "hdi_2018_cat"], "< 0.5" )
      
      # mean household size
      df_imputed[, "hh_size_cat"] <- cut(df_imputed[, "hh_size"], c(0, 3, 4, 5, 100), 
        labels = c("< 3.0", "3.0 to 3.9", "4.0 to 4.9", ">= 5.0"), include.lowest = TRUE)
      df_imputed[, "hh_size_cat"] <- relevel(df_imputed[, "hh_size_cat"], "< 3.0" )
      
      # mean mobility change  
      df_imputed[, "mean_mobility_change_cat"] <- cut(df_imputed[, "mean_mobility_change"], c(-100, -30, -20, -10, 100), 
        labels = c("-30% or more", "-20 to -29%", "-10 to -19%", "-10% or less"), include.lowest = TRUE)
      df_imputed[, "mean_mobility_change_cat"] <- relevel(df_imputed[, "mean_mobility_change_cat"], "-30% or more" )
      
      # population density  
      df_imputed[, "pop_km2_2020_cat"] <- cut(df_imputed[, "pop_km2_2020"], c(0, 25, 50, 75, 100, 1000000), 
        labels = c("< 25", "25 to 49", "50 to 74", "75 to 99", ">= 100"), include.lowest = TRUE)
      df_imputed[, "pop_km2_2020_cat"] <- relevel(df_imputed[, "pop_km2_2020_cat"], "< 25" )
      
      # mean stringency of control measures  
      df_imputed[, "mean_stringency_cat"] <- cut(df_imputed[, "mean_stringency"], c(0, 40, 50, 60, 100), 
        labels = c("< 40", "40 to 49", "50 to 59", ">= 60"), include.lowest = TRUE)
      df_imputed[, "mean_stringency_cat"] <- relevel(df_imputed[, "mean_stringency_cat"], "< 40" )
      
      # mean testing policy  
      df_imputed[, "mean_testing_policy_cat"] <- cut(df_imputed[, "mean_testing_policy"], c(0, 1, 1.5, 2.0, 10), 
        labels = c("< 1.0", "1.0 to 1.4", "1.5 to 1.9", ">= 2.0"), include.lowest = TRUE)
      df_imputed[, "mean_testing_policy_cat"] <- relevel(df_imputed[, "mean_testing_policy_cat"], "< 1.0" )
      
      # testing rate  
      df_imputed[, "test_rate_cat"] <- cut(df_imputed[, "test_rate"], c(0, 25, 100, 100000), 
        labels = c("< 25 per 1000", "25 to 99 per 1000", ">= 100 per 1000"), include.lowest = TRUE)
      df_imputed[, "test_rate_cat"] <- relevel(df_imputed[, "test_rate_cat"], "< 25 per 1000" )
    
            
#..........................................................................................
### Visualising key patterns in the epidemic outcomes
#..........................................................................................

  #...................................      
  ## Colour coding for WHO regions
    cols <- brewer.pal(length(unique(df$region)), "Set1")
    names(cols) <- sort(unique(df$region))
  
  #...................................      
  ## Visualise the mean R_t by region
      plot_rt <- ggplot(df) +
        geom_boxplot(aes(x = region, y = rt_mean, fill = region), lwd = 1, alpha = 0.5 ) +
        theme_bw() +
        labs(x = "region", y = "mean instantaneous reproduction number (R_t)") +
        scale_fill_manual(values = cols ) +
        scale_y_continuous(breaks = seq(0.0, 1.8, by = 0.2)) +
        geom_hline(yintercept = 1, colour= "red", size = 1, alpha = 0.8)
      plot_rt
      ggsave("rt_mean_by_region.png", height = 15, width = 25, units = "cm", dpi = "print")
         
  #...................................      
  ## visualise age-specific cases and deaths data
      
    # Data availability
    print("number of countries with age-specific cases data:")
    length(unique(cases$country) )
    print("number of countries with age-specific deaths data:")
    length(unique(deaths$country) )
    print("number of countries with age-specific cases and deaths data:")
    length( intersect(unique(cases$country), unique(deaths$country)) )
    
    # Visualise age-specific distributions of cases and deaths, by country
      # cases: country-by-country plots
      plot <- ggplot(cases, aes(x = age_max, y = cases_rate) ) +
        geom_point(colour = palette_cb[3]) +
        geom_line(colour = palette_cb[3], size = 1.5, alpha = 0.5) +
        theme_bw() +
        labs(x = "age (years)", y = "cases per 1000 people") +
        facet_wrap(~country , ncol = 4, scales = "free_y")
      plot
      ggsave("cases_rates_by_country.png", width = 23, height = 30, units = "cm", dpi = "print")    

      # cases: single plot of all countries, colour-coded by region
        # prepare data
        x1 <- cases[, c("region", "country", "age_max", "n_cases")]
        x1 <- x1[order(x1[, "country"], x1[, "age_max"]), ]
          # cumulative percentages of all cases by age
          x2 <- aggregate(x1$n_cases, by = list(country = x1$country), FUN = sum)
          colnames(x2) <- c("country", "total_cases")
          x1[, "cum_cases"] <- do.call(c, tapply(x1[, "n_cases"], x1[, "country"], FUN=cumsum))
          x1 <- merge(x1, x2, by = "country", all.x = TRUE)
          x1[, "cum_prop"] <- x1[, "cum_cases"] / x1[, "total_cases"]
        
        # plot
        plot_cases <- ggplot(x1, aes(x = age_max, y = cum_prop)) +
          geom_line(aes(group = factor(country), colour = region ), size = 1, alpha = 0.3  ) +
          scale_colour_manual(values = cols) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 0.1) ) +
          theme_bw() +
          labs(x = "age (years)", y = "cumulative percent of all-age cases")          
          
        plot_cases
        ggsave("cases_cum_single_plot.png", height = 15, width = 25, units = "cm", dpi = "print")

      
      # deaths: country-by-country plots
      plot <- ggplot(deaths, aes(x = age_max, y = deaths_rate) ) +
        geom_point(colour = palette_cb[7]) +
        geom_line(colour = palette_cb[7], size = 1.5, alpha = 0.5) +
        theme_bw() +
        labs(x = "age (years)", y = "deaths per 1000 people") +
        facet_wrap(~country , ncol = 4, scales = "free_y")
      plot
      ggsave("deaths_rates_by_country.png", width = 23, height = 30, units = "cm", dpi = "print")   
            
      
      # deaths: single plot of all countries, colour-coded by region
        # prepare data
        x1 <- deaths[, c("region", "country", "age_max", "n_deaths")]
        x1 <- x1[order(x1[, "country"], x1[, "age_max"]), ]
          # cumulative percentages of all deaths by age
          x2 <- aggregate(x1$n_deaths, by = list(country = x1$country), FUN = sum)
          colnames(x2) <- c("country", "total_deaths")
          x1[, "cum_deaths"] <- do.call(c, tapply(x1[, "n_deaths"], x1[, "country"], FUN=cumsum))
          x1 <- merge(x1, x2, by = "country", all.x = TRUE)
          x1[, "cum_prop"] <- x1[, "cum_deaths"] / x1[, "total_deaths"]
        
        # plot
        plot_deaths <- ggplot(x1, aes(x = age_max, y = cum_prop)) +
          geom_line(aes(group = factor(country), colour = region ), size = 1, alpha = 0.3  ) +
          scale_colour_manual(values = cols) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 0.1) ) +
          theme_bw() +
          labs(x = "age (years)", y = "cumulative percent of all-age deaths")          
          
        plot_deaths
        ggsave("deaths_cum_single_plot.png", height = 15, width = 25, units = "cm", dpi = "print")

    
  #...................................      
  ## Visualise standardised median ages by region
    
    # Cases
      plot_age_cases <- ggplot(df) +
        geom_boxplot(aes(x = region, y = median_age_cases, fill = region), lwd = 1, alpha = 0.5 ) +
        theme_bw() +
        labs(x = "region", y = "standardised median age of cases") +
        scale_fill_manual(values = cols )
      plot_age_cases
      ggsave("cases_medians_by_region.png", height = 15, width = 25, units = "cm", dpi = "print")

    # Deaths
      plot_age_deaths <- ggplot(df) +
        geom_boxplot(aes(x = region, y = median_age_deaths, fill = region), lwd = 1, alpha = 0.5 ) +
        theme_bw() +
        labs(x = "region", y = "standardised median age of deaths") +
        scale_fill_manual(values = cols )
      plot_age_deaths
      ggsave("deaths_medians_by_region.png", height = 15, width = 25, units = "cm", dpi = "print")
  
  #...................................      
  ## Visualise crude and age-standardised CFRs by region
  
      # Crude CFR
      plot_cfr_crude <- ggplot(df) +
        geom_boxplot(aes(x = region, y = cfr_crude, fill = region), lwd = 1, alpha = 0.5 ) +
        theme_bw() +
        labs(x = "region", y = "crude case-fatality ratio") +
        scale_y_continuous(limits = c(0, 0.075), labels = scales::percent_format(accuracy = 0.1),
          breaks = seq(0, 0.075, by = 0.005) ) +
        scale_fill_manual(values = cols )
      plot_cfr_crude
      ggsave("cfr_crude_by_region.png", height = 15, width = 25, units = "cm", dpi = "print")
  
      # Age-standardised CFR 1
      plot_cfr_std1 <- ggplot(df) +
        geom_boxplot(aes(x = region, y = cfr_std1, fill = region), lwd = 1, alpha = 0.5 ) +
        theme_bw() +
        labs(x = "region", y = "age-standardised case-fatality ratio 1") +
        scale_y_continuous(limits = c(0, 0.075), labels = scales::percent_format(accuracy = 0.1), 
          breaks = seq(0, 0.075, by = 0.005) ) +
        scale_fill_manual(values = cols )
      plot_cfr_std1
      ggsave("cfr_std1_by_region.png", height = 15, width = 25, units = "cm", dpi = "print")
  
      # Age-standardised CFR 2
      plot_cfr_std2 <- ggplot(df) +
        geom_boxplot(aes(x = region, y = cfr_std2, fill = region), lwd = 1, alpha = 0.5 ) +
        theme_bw() +
        labs(x = "region", y = "age-standardised case-fatality ratio 2") +
        scale_y_continuous(limits = c(0, 0.075), labels = scales::percent_format(accuracy = 0.1), 
          breaks = seq(0, 0.075, by = 0.005) ) +
        geom_hline(yintercept = df[df$country_iso == "KOR", "cfr_std2"], colour = "red3", linetype = "longdash",
          size = 1.0, alpha = 0.5) +
        annotate("text", y = 0.014, x = "SEARO", label = "CFR in South Korea", colour = "red3") +
        scale_fill_manual(values = cols )
      plot_cfr_std2
      ggsave("cfr_std2_by_region.png", height = 15, width = 25, units = "cm", dpi = "print")
  

  #...................................      
  ## Create a single panel plot
    plot <- ggarrange(plot_rt + rremove("y.title"), plot_age_cases + rremove("y.title"), plot_cases + rremove("y.title"),
      plot_age_deaths + rremove("y.title"), plot_deaths + rremove("y.title"), plot_cfr_crude + rremove("y.title"), 
      plot_cfr_std1 + rremove("y.title"), plot_cfr_std2 + rremove("y.title"),
      ncol = 2, nrow = 4, labels = c("Reproduction number", "Standardised median age of cases", "Cumulative percent of cases", 
        "Standardised median age of deaths", "Cumulative percent of deaths", "Crude CFR",
        "Age-standardised CFR", "Incidence-standardised CFR"),
      common.legend = TRUE, legend = "right", align = "v", label.y = 0.98, label.x = c(0.14, 0.11, 0.12, 0.1, 0.12, 0.15, 0.13, 0.12 ), 
      hjust = -0.1, font.label = list(size = 10, color = "grey20", face = "bold", family = NULL) )
    plot
    
    ggsave("outcomes_by_region.png", height = 31, width = 25, units = "cm", dpi = "print")
        
 
#..........................................................................................
### Statistical modelling: Univariate analysis
#..........................................................................................

  #...................................      
  ## Visualise collinearity
    # Select variables
    x1 <- vars[vars$role %in% c("exposure", "confounder"), "variable"]
    x1 <- x1[! grepl("_cat", x1)]

    # Scatter plots
    plot <- ggpairs(df[, x1] ) + theme_bw()
    ggsave("collinearity_plots.png", plot, height = 45, width = 45, units = "cm", dpi = "print")    
    
    # Collinearity heatmap
    plot <- ggcorr(df[, x1], method = c("pairwise.complete.obs", "pearson"),
      low = "steelblue", mid = "grey90", high = "darkred", geom = "tile", nbreaks = 5, min_size = 0, max_size = 6, 
      label = TRUE, label_size = 3, label_round = 2, size = 3, hjust = 0.75, layout.exp = 1, legend.position = "off")
    ggsave("collinearity_heatmap.png", plot, height = 25, width = 25, units = "cm", dpi = "print")    
    
    # Collinearity matrix (Pearson correlation coefficients)
    coll_matrix <- round(abs(cor(df[, x1], use = "pairwise.complete.obs", )), 2)
    write.csv(coll_matrix, "out_collinearity_matrix.csv", row.names = TRUE)
    
    # Collinear pairs
    x1 <- data.frame(coll_matrix)
    x2 <- which(x1 >= 0.70, arr.ind = TRUE)
    coll_pairs <- data.frame( cbind(rownames(x1)[x2[, 1]], colnames(x1)[x2[, 2]]) )
    colnames(coll_pairs) <- c("var1", "var2")
    coll_pairs <- subset(coll_pairs, var1 != var2)
    coll_pairs <- t(apply(coll_pairs, 1, sort) )
    coll_pairs <- data.frame(unique(coll_pairs) )
    colnames(coll_pairs) <- c("var1", "var2")
    write.csv(coll_pairs, "out_collinear_vars.csv", row.names = FALSE)    
    
  #...................................      
  ## Reproduction number
  f_univar("rt_mean", "rt_model", vars, df)

  #...................................      
  ## Median age of cases
  f_univar("median_age_cases", "age_cases_model", vars, df)

  #...................................      
  ## Median age of deaths
  f_univar("median_age_deaths", "age_deaths_model", vars, df)
  
  #...................................      
  ## Crude CFR
  f_univar("cfr_crude_ln", "cfr_crude_model", vars, df)
   
  #...................................      
  ## Age-standardised CFR 1
  f_univar("cfr_std1_ln", "cfr_std1_model", vars, df)

  #...................................      
  ## Age-standardised CFR 2
  f_univar("cfr_std2", "cfr_std2_model", vars, df)
  
        
#..........................................................................................
### Statistical modelling: Multivariate analysis - Approach 1: Random Forest Regression
#..........................................................................................
               
  #...................................      
  ## Reproduction number
    # Entire world
      # Random forest with complete cases only
      x1 <- f_rf_prep("rt_mean", "rt_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, NA)
      f_rf("rt_mean", "rt_model", x1, effect_modifiers, "out_rt_mean_rf_complete.csv")

      # Random forest using dataset with missing predictor values imputed
      x1 <- f_rf_prep("rt_mean", "rt_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, NA)
      f_rf("rt_mean", "rt_model", x1, effect_modifiers, "out_rt_mean_rf_impute.csv")

    # AFRO region only
      # Random forest with complete cases only
      x1 <- f_rf_prep("rt_mean", "rt_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, "AFRO")
      f_rf("rt_mean", "rt_model", x1, effect_modifiers, "out_rt_mean_rf_complete_afro.csv")
        
      # Random forest using dataset with missing predictor values imputed
      x1 <- f_rf_prep("rt_mean", "rt_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, "AFRO")
      f_rf("rt_mean", "rt_model", x1, effect_modifiers, "out_rt_mean_rf_impute_afro.csv")

        
  #...................................      
  ## Median age of cases
    # Random forest with complete cases only
    x1 <- f_rf_prep("median_age_cases", "age_cases_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, NA)
    f_rf("median_age_cases", "age_cases_model", x1, effect_modifiers, "out_median_age_cases_rf_complete.csv")
        
    # Random forest using dataset with missing predictor values imputed
    x1 <- f_rf_prep("median_age_cases", "age_cases_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, NA)
    f_rf("median_age_cases", "age_cases_model", x1, effect_modifiers, "out_median_age_cases_rf_impute.csv")


  #...................................      
  ## Median age of deaths
    # Random forest with complete cases only
    x1 <- f_rf_prep("median_age_deaths", "age_deaths_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, NA)
    f_rf("median_age_deaths", "age_deaths_model", x1, effect_modifiers, "out_median_age_deaths_rf_complete.csv")
        
    # Random forest using dataset with missing predictor values imputed
    x1 <- f_rf_prep("median_age_deaths", "age_deaths_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, NA)
    f_rf("median_age_deaths", "age_deaths_model", x1, effect_modifiers, "out_median_age_deaths_rf_impute.csv")


  #...................................      
  ## Crude CFR
    # Entire world
      # Random forest with complete cases only
      x1 <- f_rf_prep("cfr_crude", "cfr_crude_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, NA)
      f_rf("cfr_crude", "cfr_crude_model", x1, effect_modifiers, "out_cfr_crude_rf_complete.csv")
        
      # Random forest using dataset with missing predictor values imputed
      x1 <- f_rf_prep("cfr_crude", "cfr_crude_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, NA)
      f_rf("cfr_crude", "cfr_crude_model", x1, effect_modifiers, "out_cfr_crude_rf_impute.csv")

    # AFRO region only
      # Random forest with complete cases only
      x1 <- f_rf_prep("cfr_crude", "cfr_crude_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, "AFRO")
      f_rf("cfr_crude", "cfr_crude_model", x1, effect_modifiers, "out_cfr_crude_rf_complete_afro.csv")
        
      # Random forest using dataset with missing predictor values imputed
      x1 <- f_rf_prep("cfr_crude", "cfr_crude_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, "AFRO")
      f_rf("cfr_crude", "cfr_crude_model", x1, effect_modifiers, "out_cfr_crude_rf_impute_afro.csv")
  
            
  #...................................      
  ## Age-standardised CFR
    # Random forest with complete cases only
    x1 <- f_rf_prep("cfr_std1", "cfr_std1_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, NA)
    f_rf("cfr_std1", "cfr_std1_model", x1, effect_modifiers, "out_cfr_std1_rf_complete.csv")
        
    # Random forest using dataset with missing predictor values imputed
    x1 <- f_rf_prep("cfr_std1", "cfr_std1_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, NA)
    f_rf("cfr_std1", "cfr_std1_model", x1, effect_modifiers, "out_cfr_std1_rf_impute.csv")

            
  #...................................      
  ## Incidence-standardised CFR
    # Random forest with complete cases only
    x1 <- f_rf_prep("cfr_std2", "cfr_std2_model", vars, df, threshold_compl, c("hdi_2018"), FALSE, NA)
    f_rf("cfr_std2", "cfr_std2_model", x1, effect_modifiers, "out_cfr_std2_rf_complete.csv")
        
    # Random forest using dataset with missing predictor values imputed
    x1 <- f_rf_prep("cfr_std2", "cfr_std2_model", vars, df, threshold_compl, c("hdi_2018"), TRUE, NA)
    f_rf("cfr_std2", "cfr_std2_model", x1, effect_modifiers, "out_cfr_std2_rf_impute.csv")

    
    
#..........................................................................................
### Statistical modelling: Multivariate analysis - Approach 2: Ordinary linear regression (imputed dataset only)
#..........................................................................................

  #...................................      
  ## Reproduction number
    # Select variables for this outcome
      # exposures or confounders
      vars_now <- vars[! is.na(vars$rt_model), "variable"]
      
      # categorical variables when continuous don't seem linearly associated with outcome
      x1 <- c("mean_testing_policy_cat", "test_rate_cat")
      vars_now <- vars_now[! grepl("_cat", vars_now)]
      vars_now <- vars_now[! vars_now %in% gsub("_cat", "", x1) ]
      vars_now <- c(vars_now, x1)
      
      # eliminate collinear variables
      vars_now <- vars_now[! vars_now %in% c("hdi_2018", "hh_size")]
      
      # log or linear outcome
      outcome_now <- "rt_mean"
      
      # interaction terms
      int_terms <- c()
      x1 <- effect_modifiers[effect_modifiers$model == "rt_model", ]
      for (i in 1:nrow(x1)) { int_terms <- c(int_terms, paste(x1[i, "exposure"], x1[i, "effect_modifier"], sep = ":")) }
      
    # Entire world
      # fit model
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed)
      summary(fit)
      f_diag_ols(fit)
        
      # remove influential observations
      x2 <- as.data.frame(influence.measures(fit)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed[-x3, ])
        print(summary(fit))
        f_diag_ols(fit)
      }
      
    
      # reduce model
      vars_red <- vars_now[! vars_now %in% c("mean_mobility_change", "prev_pv", "prev_schisto", "prev_pf", 
        "prev_helminths", "mean_testing_policy_cat")]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed)
      summary(fit_red)
      f_diag_ols(fit_red)
    
      # remove influential observations from reduced model
      x2 <- as.data.frame(influence.measures(fit_red)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed[-x3, ])
        print(summary(fit_red))
        f_diag_ols(fit_red)
      }

      # add interactions one by one
      for (i in int_terms) {
        x2 <- df_imputed
        if (length(x3) > 0) {x2 <- x2[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        
        # retain interaction of age and testing policy
      
      # final reduced fit
      fit_red <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], int_terms[2]), x2)
            
      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute.csv", sep = "") )
  
       
    # AFRO region only
      # fit model
      x1 <- subset(df_imputed, region == "AFRO")
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now],  x1)
      summary(fit)
      f_diag_ols(fit)
        
      # remove influential observations
      x2 <- as.data.frame(influence.measures(fit)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], x1[-x3, ])
        print(summary(fit))
        f_diag_ols(fit)
      }
      
      # reduce model
      vars_red <- vars_now[! vars_now %in% c("mean_mobility_change", "prev_pv", "prev_schisto", "prev_pf", 
        "mean_testing_policy_cat", "prev_helminths", "age_mean", "test_rate_cat")]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], x1)
      summary(fit_red)
      f_diag_ols(fit_red)
    
      # remove influential observations from reduced model
        # [no evidence of influential observations]
      
      # add interactions one by one
      for (i in int_terms) {
        x2 <- x1
        if (length(x3) > 0) {x2 <- x1[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        
        # no interactions improve model
      
      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute_afro.csv", sep = "") )
      
  #...................................      
  ## Median age of cases
    # Select variables for this outcome
      # exposures or confounders
      vars_now <- vars[! is.na(vars$age_cases_model), "variable"]
      
      # categorical variables when continuous don't seem linearly associated with outcome
      x1 <- c("prev_filaria_cat", "prev_helminths_cat", "prev_pf_cat", "prev_pv_cat", "prev_schisto_cat", 
        "mean_mobility_change_cat", "test_rate_cat", "mean_testing_policy_cat")
      vars_now <- vars_now[! grepl("_cat", vars_now)]
      vars_now <- vars_now[! vars_now %in% gsub("_cat", "", x1) ]
      vars_now <- c(vars_now, x1)
      
      # eliminate collinear variables
      vars_now <- vars_now[! vars_now %in% c("hdi_2018", "hh_size")]
      
      # log or linear outcome
      outcome_now <- "median_age_cases"
      
      # interaction terms
      int_terms <- c()
      x1 <- effect_modifiers[effect_modifiers$model == "age_cases_model", ]
      for (i in 1:nrow(x1)) { int_terms <- c(int_terms, paste(x1[i, "exposure"], x1[i, "effect_modifier"], sep = ":")) }
      
    # Model
      # fit model
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed)
      summary(fit)
      f_diag_ols(fit)
        
      # remove influential observations
      x2 <- as.data.frame(influence.measures(fit)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed[-x3, ])
        print(summary(fit))
        f_diag_ols(fit)
      }
      
      # reduce model
      vars_red <- vars_now[! vars_now %in% c("percent_increased_risk_age_std", "prev_pv_cat", "prev_pf_cat", "prev_schisto_cat",
        "prev_filaria_cat", "mean_testing_policy_cat", "prev_helminths_cat", "mean_stringency", "test_rate_cat")]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed)
      summary(fit_red)
      f_diag_ols(fit_red)
    
      # remove influential observations from reduced model
      x2 <- as.data.frame(influence.measures(fit_red)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed[-x3, ])
        print(summary(fit_red))
        f_diag_ols(fit_red)
      }

      # add interactions one by one
      for (i in int_terms) {
        x2 <- df_imputed
        if (length(x3) > 0) {x2 <- x2[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        
        # no interactions improve model
      
      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute.csv", sep = "") )

      
      
  #...................................      
  ## Median age of deaths
    # Select variables for this outcome
      # exposures or confounders
      vars_now <- vars[! is.na(vars$age_deaths_model), "variable"]
      
      # categorical variables when continuous don't seem linearly associated with outcome
      x1 <- c("prev_filaria_cat", "prev_helminths_cat", "prev_pf_cat", "prev_pv_cat", "prev_schisto_cat", 
        "mean_mobility_change_cat", "test_rate_cat", "mean_testing_policy_cat")
      vars_now <- vars_now[! grepl("_cat", vars_now)]
      vars_now <- vars_now[! vars_now %in% gsub("_cat", "", x1) ]
      vars_now <- c(vars_now, x1)
      
      # eliminate collinear variables
      vars_now <- vars_now[! vars_now %in% c("hdi_2018", "hh_size")]
      
      # log or linear outcome
      outcome_now <- "median_age_deaths"
      
      # interaction terms
      int_terms <- c()
      x1 <- effect_modifiers[effect_modifiers$model == "age_deaths_model", ]
      for (i in 1:nrow(x1)) { int_terms <- c(int_terms, paste(x1[i, "exposure"], x1[i, "effect_modifier"], sep = ":")) }
      
    # Model
      # fit model
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed)
      summary(fit)
      f_diag_ols(fit)
        
      # remove influential observations
        # [no evidence of influential observations]

      # reduce model
      vars_red <- vars_now[! vars_now %in% c("mean_mobility_change_cat",
        "mean_testing_policy_cat", "prev_filaria_cat", "prev_pf_cat", "prev_pv_cat", 
        "prev_helminths_cat", "prev_schisto_cat", "test_rate_cat")]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed)
      summary(fit_red)
      f_diag_ols(fit_red)
    
      # remove influential observations from reduced model
      x2 <- as.data.frame(influence.measures(fit_red)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed[-x3, ])
        print(summary(fit_red))
        f_diag_ols(fit_red)
      }
      
      # add interactions one by one
      for (i in int_terms) {
        x2 <- df_imputed
        if (length(x3) > 0) {x2 <- x2[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        
        # no interactions improve model

      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute.csv", sep = "") )

      
  #...................................      
  ## Crude CFR
    # Select variables for this outcome
      # exposures or confounders
      vars_now <- vars[! is.na(vars$cfr_crude_model), "variable"]
      
      # categorical variables when continuous don't seem linearly associated with outcome
      x1 <- c("prev_filaria_cat", "prev_helminths_cat", "prev_pf_cat", "prev_pv_cat", "prev_schisto_cat", 
        "mean_mobility_change_cat", "test_rate_cat", "mean_testing_policy_cat")
      vars_now <- vars_now[! grepl("_cat", vars_now)]
      vars_now <- vars_now[! vars_now %in% gsub("_cat", "", x1) ]
      vars_now <- c(vars_now, x1)
      
      # eliminate collinear variables
      vars_now <- vars_now[! vars_now %in% c("hdi_2018", "hh_size")]
      
      # log or linear outcome
      outcome_now <- "cfr_crude_ln"
      
      # interaction terms
      int_terms <- c()
      x1 <- effect_modifiers[effect_modifiers$model == "cfr_crude_model", ]
      for (i in 1:nrow(x1)) { int_terms <- c(int_terms, paste(x1[i, "exposure"], x1[i, "effect_modifier"], sep = ":")) }
      
    # Entire world
      # fit model
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed)
      summary(fit)
      f_diag_ols(fit)

      # k-folds cross-validation
      f_cv(fit, k_folds, TRUE)
        
      # remove influential observations
      x2 <- as.data.frame(influence.measures(fit)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed[-x3, ])
        print(summary(fit))
        f_diag_ols(fit)
      }
      
      # reduce model
      vars_red <- vars_now[! vars_now %in% c("mean_mobility_change_cat", "prev_pv_cat", "prev_schisto_cat", "prev_pf_cat", 
        "test_rate_cat", "age_mean", "prev_filaria_cat", "mean_stringency", "prev_helminths_cat")]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed)
      summary(fit_red)
      f_diag_ols(fit_red)
    
      # remove influential observations from reduced model
      x2 <- as.data.frame(influence.measures(fit_red)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed[-x3, ])
        summary(fit_red)
        f_diag_ols(fit_red)
      }

      # k-folds cross-validation
      f_cv(fit, k_folds, TRUE)

      # add interactions one by one
      for (i in int_terms) {
        x2 <- df_imputed
        if (length(x3) > 0) {x2 <- x2[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        

      # k-folds cross-validation
      f_cv(fit, k_folds, TRUE)
        # no interactions improve the model

      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute.csv", sep = "") )
  
    # AFRO region only
      # fit model
      x1 <- subset(df_imputed, region == "AFRO")
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now],  x1)
      summary(fit)
      f_diag_ols(fit)

      # k-folds cross-validation
      f_cv(fit, k_folds, TRUE)
        
      # remove influential observations
      x2 <- as.data.frame(influence.measures(fit)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], x1[-x3, ])
        summary(fit)
        f_diag_ols(fit)
      }
      
      # reduce model
      vars_red <- vars_now[! vars_now %in% c("mean_mobility_change_cat", "prev_pv_cat", "prev_schisto_cat", "prev_pf_cat", 
        "test_rate_cat", "age_mean", "prev_filaria_cat", "percent_increased_risk_age_std", "mean_stringency", "prev_helminths_cat" )]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], x1)
      summary(fit_red)
      f_diag_ols(fit_red)

      # k-folds cross-validation
      f_cv(fit, k_folds, TRUE)
    
      # remove influential observations from reduced model
      x2 <- as.data.frame(influence.measures(fit_red)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], x1[-x3, ])
        summary(fit_red)
        f_diag_ols(fit_red)
      }
 
      # add interactions one by one
      for (i in int_terms) {
        x2 <- x1
        if (length(x3) > 0) {x2 <- x2[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        

      # k-folds cross-validation
      f_cv(fit, k_folds, TRUE)
        # no interactions improve model   
        
      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute_afro.csv", sep = "") )

            
  #...................................      
  ## Age-standardised CFR
    # Select variables for this outcome
      # exposures or confounders
      vars_now <- vars[! is.na(vars$cfr_std1_model), "variable"]
      
      # categorical variables when continuous don't seem linearly associated with outcome
      x1 <- c("prev_filaria_cat", "prev_helminths_cat", "prev_pf_cat", "prev_pv_cat", "prev_schisto_cat", 
        "mean_mobility_change_cat", "test_rate_cat", "mean_testing_policy_cat")
      vars_now <- vars_now[! grepl("_cat", vars_now)]
      vars_now <- vars_now[! vars_now %in% gsub("_cat", "", x1) ]
      vars_now <- c(vars_now, x1)
      
      # eliminate collinear variables
      vars_now <- vars_now[! vars_now %in% c("hdi_2018", "hh_size")]
      
      # log or linear outcome
      outcome_now <- "cfr_std1_ln"
      
      # interaction terms
      int_terms <- c()
      x1 <- effect_modifiers[effect_modifiers$model == "cfr_std1_model", ]
      for (i in 1:nrow(x1)) { int_terms <- c(int_terms, paste(x1[i, "exposure"], x1[i, "effect_modifier"], sep = ":")) }
      
    # Entire world
      # fit model
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed)
      summary(fit)
      f_diag_ols(fit)
        
      # remove influential observations
        # [no evidence of influential observations]
      
      # reduce model
      vars_red <- vars_now[! vars_now %in% c("mean_mobility_change_cat", "prev_pv_cat", "prev_schisto_cat", "prev_pf_cat", 
        "test_rate_cat", "prev_filaria_cat", "mean_testing_policy_cat", "age_mean", "prev_helminths_cat",
        "mean_stringency")]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed)
      summary(fit_red)
      f_diag_ols(fit_red)
    
      # remove influential observations from reduced model
        # [no evidence of influential observations]
      
      # add interactions one by one
      for (i in int_terms) {
        x2 <- df_imputed
        if (length(x3) > 0) {x2 <- x2[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        
      
      # final model
      fit_red <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], int_terms[3]), x2)
           
      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute.csv", sep = "") )

            
  #...................................      
  ## Incidence-standardised CFR
    # Select variables for this outcome
      # exposures or confounders
      vars_now <- vars[! is.na(vars$cfr_std2_model), "variable"]
      
      # categorical variables when continuous don't seem linearly associated with outcome
      x1 <- c("prev_filaria_cat", "prev_helminths_cat", "prev_pf_cat", "prev_pv_cat", "prev_schisto_cat", 
        "mean_mobility_change_cat", "test_rate_cat", "mean_testing_policy_cat")
      vars_now <- vars_now[! grepl("_cat", vars_now)]
      vars_now <- vars_now[! vars_now %in% gsub("_cat", "", x1) ]
      vars_now <- c(vars_now, x1)
      
      # eliminate collinear variables
      vars_now <- vars_now[! vars_now %in% c("hdi_2018", "hh_size")]
      
      # log or linear outcome
      outcome_now <- "cfr_std2"
      
      # interaction terms
      int_terms <- c()
      x1 <- effect_modifiers[effect_modifiers$model == "cfr_std2_model", ]
      for (i in 1:nrow(x1)) { int_terms <- c(int_terms, paste(x1[i, "exposure"], x1[i, "effect_modifier"], sep = ":")) }
      
    # Entire world
      # fit model
      fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed)
      summary(fit)
      f_diag_ols(fit)
        
      # remove influential observations
      x2 <- as.data.frame(influence.measures(fit)$infmat)
      x3 <- which(abs(x2$dffit) >= 1)
      if (length(x3) > 0) {
        fit <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_now], df_imputed[-x3, ])
        print(summary(fit))
        f_diag_ols(fit)
      }
      
      # reduce model
      vars_red <- vars_now[! vars_now %in% c("prev_pv_cat", "prev_schisto_cat", "mean_stringency", "prev_pf_cat",
        "age_mean", "prev_filaria_cat", "prev_helminths_cat", "mean_mobility_change_cat", "mean_testing_policy_cat",
        "test_rate_cat")]
      fit_red <- f_lm(outcome_now, colnames(df_imputed)[colnames(df_imputed) %in% vars_red], df_imputed)
      summary(fit_red)
      f_diag_ols(fit_red)
          
      # remove influential observations from reduced model
        # [no evidence of influential observations]

      # add interactions one by one
      for (i in int_terms) {
        x2 <- df_imputed
        if (length(x3) > 0) {x2 <- x2[-3, ]}
        fit_int <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], i), x2)
        print(summary(fit_int))
        f_diag_ols(fit_int)
      }        
      
      # final model
      fit_red <- f_lm(outcome_now, c(colnames(df_imputed)[colnames(df_imputed) %in% vars_red], int_terms[3]), x2)
      
      # write output to file (both full and reduced models)
      f_lm_write(paste("out_", outcome_now, "_lm_impute.csv", sep = "") )


#...................................      
# Age plot

plot_median_age_death <- ggplot(df_imputed, aes(x=median_age_cases, y=median_age_deaths,
                                        color=region)) + 
  geom_point(size=6, shape=15) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15)) +
  labs(y="median age of deaths", x = "median age of cases") +
  theme_bw()
plot_median_age_death
ggsave("median_age_death_by_region.png", height = 15, width = 25, units = "cm", dpi = "print")

#..........................................................................................
### ENDS
#..........................................................................................
    
    
    
