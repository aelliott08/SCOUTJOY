#' Significant Cross-trait OUtliers and Trends in JOint York regression (SCOUTJOY)
#'
#' @description
#' xxx
#'
#' @details
#' xxx
#' 
#' @param BetaOutcome
#' String, giving the column name in `data` for the effect sizes (betas) on the outcome (y-axis trait)
#' 
#' @param BetaExposure 
#' String, giving the column name in `data` for the effect sizes (betas) on the exposure (x-axis trait)
#' 
#' @param SdOutcome,
#' String, giving the column name in `data` for the standard error (se) of the effect sizes on the outcome (y-axis trait)
#' 
#' @param SdExposure, 
#' String, giving the column name in `data` for the standard error (se) of the effect sizes on the exposure (x-axis trait)
#' 
#' @param data
#' Data frame with the specified column names containing effect sizes and standard errors
#' 
#' @param ids 
#' (Optional) string giving the column name in `data` of variant IDs to use for labeling results tables
#' 
#' @param CovIntercept 
#' Default = 0, 
#' 
#' @param ExposureIntercept
#' Default =1, 
#' 
#' @param OutcomeIntercept
#' Default =1, 
#' 
#' @param OutlierNull
#' Default = "estimate"
#' 
#' @param SignifThreshold 
#' Default = 0.05,
#' 
#' @param NullReplicates 
#' Default = NA, 
#' 
#' @param maxOutlierReps
#' Default =10, 
#' 
#' @param seed 
#' Default = NULL, 
#' 
#' @param printProgress 
#' Default = TRUE
#'
#'
#' @return
#' List with elements:
#' \describe{
#'
#'   \item{\code{Global Test}}{xxx}
#'   
#'   \item{\code{xxx}}{xxx}
#'
#' }
#'
#' @export
#'
#' @references
#'
#' xxx
#'
#' @examples
#'
#' x <- stats::rnorm(20)
#' 


scoutjoy <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, ids = NULL, CovIntercept = 0, ExposureIntercept=1, OutcomeIntercept=1, OutlierNull="estimate", SignifThreshold = 0.05, NullReplicates = NA, maxOutlierReps=10, seed = NULL, printProgress = TRUE){

  
	if(!is.null(seed)){
		set.seed(seed)
	}
  
  BonferroniOutliers=TRUE
  IterateOutliers=TRUE
  FixedSlope = NULL
  
  AnalyticLoss = FALSE
  if(is.na(NullReplicates) || NullReplicates < 1){
    AnalyticLoss = TRUE
  }
	
	if((length(BetaExposure)!=1) || (length(SdExposure)!=1) || (length(BetaOutcome)!=1) || (length(SdOutcome)!=1)){
		stop("BetaExposure, SdExposure, BetaOutcome, and SdOutcome must each be a single column name")
	}

  if(class(data)[1] != "data.frame"){
    stop("data must be an object of class data.frame, try converting data to a data.frame \'data = as.data.frame(data)\'")
  }
  
  if((SignifThreshold > 1) || (SignifThreshold <= 0)){
    warning("Significance threshold outside of 0-1; this is probably an error unless manipulating bonferroni correction")
  }
  
	if(!is.null(FixedSlope) && IterateOutliers){
		warning("Cannot iterate with FixedSlope. Ignoring 'IterateOutliers'.")
		IterateOutliers <- FALSE
	}

  if(!is.null(ids)){
    row.names(data) <- as.character(data[,ids])
  }  
  data <- data[, c(BetaExposure, SdExposure,BetaOutcome, SdOutcome)]
  data <- data[rowSums(is.na(data)) == 0, ]
  data[, c(BetaOutcome, BetaExposure)] <- data[, c(BetaOutcome, BetaExposure)] * sign(data[, BetaExposure[1]])

    
  if(AnalyticLoss){
  
  	# 0- Transforming the data + checking number of observations
    colnames(data) <- c("x","sx","y","sy")
    data$re <- rep(CovIntercept/sqrt(ExposureIntercept*OutcomeIntercept), nrow(data))
  
  	if(nrow(data) <= 3)
  		stop("Not enough instrumental variables")
  
  }else{
    data$Weights <- 1/data[, SdOutcome]^2
    
    if(nrow(data) <= length(BetaExposure) + 2)
      stop("Not enough intrumental variables")
    
  	if((1/(NullReplicates+1) > SignifThreshold/nrow(data))){
  		warning(paste0("Outlier test unstable. The Bonferroni significance threshold of ", SignifThreshold/nrow(data), " is not achievable with only ", NullReplicates, " replicates to compute the null distribution. The current minimum is ", 1/(NullReplicates+1), ". Increase NullReplicates."))
  	}
  }

	###
	# Global test
	###	
	if(printProgress){print("Starting global test")}

  if(AnalyticLoss){
    
    mod_all <- yorktools::york(data, intercept=0, gridstartvals=T)
    if(!mod_all$converge$converged){
      stop("york regression failed to converge for global test")
    }
    GlobalTest <- list(RSSobs = mod_all$chi2, Pvalue = mod_all$p.value)
    
  }else{
  	
  	mod_all <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york = TRUE, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)
  	loo <- fit_LOO(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, SdOutcome = SdOutcome, SdExposure = SdExposure, data = data, york = TRUE, yorkResid = TRUE, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)
  	loo_betas <- loo[[2]]
  		
  	# 1- Computing the observed residual sum of squares (RSS)
  	RSSobs <- loo[[1]]
  
  	# 2- Computing the distribution of expected residual sum of squares (RSS)
  	if(printProgress){print("Generating null distribution")}
  	randomData <- replicate(NullReplicates, getRandomData(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, SdOutcome = SdOutcome, SdExposure = SdExposure, loo_betas = loo_betas, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, data = data, yorkResid = TRUE), simplify = FALSE)
  	RSSexp <- sapply(randomData, fit_LOO, BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, SdOutcome = SdOutcome, SdExposure = SdExposure, york = TRUE, yorkResid = TRUE, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope) 
  
  	GlobalTest <- list(RSSobs = RSSobs, Pvalue = (sum(RSSexp[1, ] > RSSobs)+1)/(NullReplicates+1))

  }

	###
	# Outlier test
	###	
	if(printProgress){print("Starting outlier test")}
	
	if(AnalyticLoss){
	
	  OutlierTest <- AnalyticOutlierTest(data, putativeOutliers=NA, OutlierNull=OutlierNull, FixedSlope=FixedSlope)
	  
	  
	}else{
	
  	if(OutlierNull=="estimate"){
  	  OutlierTest <- OutlierTestRep(data, randomData, BetaOutcome, BetaExposure, SdOutcome, SdExposure, loo_betas, york=TRUE, yorkResid=TRUE, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope=FixedSlope)
  	}else if(OutlierNull=="fitted"){
      OutlierTest <- PressoOutlierTestRep(data, randomData, BetaOutcome, BetaExposure, SdOutcome, SdExposure, loo_betas, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept)
  	}
  
	}
  
	# initial outliers
  refOutlier <- which(OutlierTest$Pvalue <= SignifThreshold/nrow(data))
      	
  if(IterateOutliers && is.null(FixedSlope) && (length(refOutlier)>0)){
  		OutlierReps <- matrix(FALSE,2,nrow(data))
  		colnames(OutlierReps) <- rownames(data)
  		OutlierReps[2,refOutlier] <- TRUE
  		OutlierReps <- as.data.frame(OutlierReps)
  		
  		firstOutlierTest <- OutlierTest
  		firstRefOutlier <- refOutlier
  		colnames(firstOutlierTest) <- c("InitialRSSobs","InitialPvalue")
  		
  		nreps <- 1
  		while(!anyDuplicated(OutlierReps) && (nreps < maxOutlierReps)){
  			if(printProgress){print(paste0("Outlier test iteration ",nreps + 1,"..."))}
  
  		  if(AnalyticLoss){
  		    
  		    OutlierTest <- AnalyticOutlierTest(data, putativeOutliers=refOutlier, OutlierNull=OutlierNull, FixedSlope=FixedSlope)
  		    
  		  }else{
  		  
    			nonoutlier_loo_betas <- fit_LOO(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, SdOutcome = SdOutcome, SdExposure = SdExposure, data = data[-refOutlier,], york = TRUE, yorkResid = TRUE, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)[[2]]
    			
    			outlier_betas <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data[-refOutlier,], york = TRUE, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)$slopes
    			
    			iter_loo_betas <- sapply(1:nrow(data),
    				function(i){
    					if(i %in% refOutlier){
    						outlier_betas
    					}else{
    						nonoutlier_loo_betas[which(c(1:nrow(data))[-refOutlier]==i)]
    					}
    				})
    
    			rep_randomData <- replicate(NullReplicates, getRandomData(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, SdOutcome = SdOutcome, SdExposure = SdExposure, loo_betas = iter_loo_betas, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, data = data, yorkResid = TRUE), simplify = FALSE)
  
    			if(OutlierNull=="estimate"){
    			  OutlierTest <- OutlierTestRep(data, rep_randomData, BetaOutcome, BetaExposure, SdOutcome, SdExposure, loo_betas=iter_loo_betas, york=TRUE, yorkResid=TRUE, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope=FixedSlope)
    			}else if(OutlierNull=="fitted"){
    			  OutlierTest <- PressoOutlierTestRep(data, rep_randomData, BetaOutcome, BetaExposure, SdOutcome, SdExposure, loo_betas=iter_loo_betas, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept)
    			}
  		  }
  			  			
  			if(BonferroniOutliers){
  				refOutlier <- which(OutlierTest$Pvalue <= SignifThreshold/nrow(data))
  			}else{
  				refOutlier <- which(OutlierTest$Pvalue <= SignifThreshold)
  			}
  			
  			nreps <- nreps+1
  			OutlierReps[nreps+1,] <- rep(FALSE,nrow(data))
  			OutlierReps[nreps+1,refOutlier] <- TRUE
  			
  		}
  		
  		OutlierReps <- OutlierReps[c(-1),]
  		outlier_converged <- FALSE
  		if(isTRUE(all.equal(OutlierReps[nrow(OutlierReps),],OutlierReps[nrow(OutlierReps)-1,],check.attributes=FALSE))){
  			outlier_converged <- TRUE
  			OutlierReps <- OutlierReps[-nrow(OutlierReps),]
  			if(printProgress){print(paste0("Outlier test iterations converged"))}
  		}else if(nreps == maxOutlierReps){
  			warning(paste0("Outlier test iterations failed to converge in ",nreps," iterations"))
  			OutlierTest$RSSobs <- cbind.data.frame(RSSobs=rep(NA,nrow(data)),Pvalue=rep(NA,nrow(data))) 
  		}else{
  		  warning(paste0("Outlier test iterations detected a loop after ",nreps," iterations"))
  		  OutlierTest$RSSobs <- cbind.data.frame(RSSobs=rep(NA,nrow(data)),Pvalue=rep(NA,nrow(data)))
  		}
  		
  		colnames(OutlierTest) <- c("RSSobs","Pvalue")
  		OutlierTest <- cbind(firstOutlierTest, OutlierTest)
  		rownames(OutlierTest) <- rownames(data)
  		
  		OutlierTestControl <- list(
  			Iterations=OutlierReps,
  			Converged=outlier_converged,
  			NumberIterations=nreps,
  			SignifThreshold=ifelse(BonferroniOutliers,SignifThreshold/nrow(data),SignifThreshold))
  			
	}else{		
		OutlierReps <- matrix(FALSE,1,nrow(data))
		colnames(OutlierReps) <- rownames(data)
		OutlierReps[1,refOutlier] <- TRUE
		OutlierReps <- as.data.frame(OutlierReps)
		
		OutlierTestControl <- list(
			Iterations=OutlierReps,
			Converged=TRUE,
			NumberIterations=1,
			SignifThreshold=ifelse(BonferroniOutliers,SignifThreshold/nrow(data),SignifThreshold))
	}

  
  # reset column names for comparison models
  if(AnalyticLoss){
    colnames(data) <- c(BetaExposure, SdExposure, BetaOutcome, SdOutcome, "re")
    data$Weights <- 1/data[, SdOutcome]^2
    data <- data[,!(colnames(data) %in% "re")]
  }
    
  OutlierTest <- cbind(data[,!(colnames(data) %in% c("Weights","re"))], OutlierTest)

  # get final model without outliers	
	if(BonferroniOutliers){
		refOutlier <- which(OutlierTest$Pvalue <= SignifThreshold/nrow(data))
	}else{
		refOutlier <- which(OutlierTest$Pvalue <= SignifThreshold)
	}

	if(length(refOutlier) > 0){
	  mod_noOutliers <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data[-refOutlier,], york = TRUE, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)
	}



	###
	# Finalize results
	# (formatting and add some slope comparisons)
	###

	# 5- Formatting the results
	if(printProgress){print("Formatting results")}

	methodname <- "York"


	if(is.null(FixedSlope)){
		OriginalMR <- cbind.data.frame(BetaExposure, paste(methodname,"(base)"), "All_variants", Summary_fit(mod_all))

		ivw <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york = FALSE)
		ivw_mr <- cbind.data.frame(BetaExposure, "IVW", "All_variants", Summary_fit(ivw))
		colnames(ivw_mr) <- colnames(OriginalMR)
		OriginalMR <- rbind(OriginalMR, ivw_mr)

		egger <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york = FALSE, MRIntercept=TRUE)
		egger_mr <- cbind.data.frame(BetaExposure, "Egger", "All_variants", Summary_fit(egger))
		colnames(egger_mr) <- colnames(OriginalMR)
		OriginalMR <- rbind(OriginalMR, egger_mr)
		colnames(OriginalMR) <- c("Exposure", "Method", "Data", "Estimate", "SE", "T-stat", "df", "P-value")

		if(exists("mod_noOutliers")){
			if(IterateOutliers && OutlierTestControl$Converged){
				OutlierCorrectedMR <- cbind.data.frame(BetaExposure, methodname, "Exclude_final_outliers", Summary_fit(mod_noOutliers))
				mod_initialoutliers <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data[-firstRefOutlier,], york = TRUE, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)
				initial_mr <- cbind.data.frame(BetaExposure, methodname, "Exclude_initial_outliers", Summary_fit(mod_initialoutliers))
				colnames(initial_mr) <- colnames(OutlierCorrectedMR)
				OutlierCorrectedMR <- rbind(OutlierCorrectedMR, initial_mr)
			}else{
				OutlierCorrectedMR <- cbind.data.frame(BetaExposure, methodname, "Exclude_outliers", Summary_fit(mod_noOutliers))
			}
		}else if (IterateOutliers && !OutlierTestControl$Converged){
			mod_initialoutliers <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data[-firstRefOutlier,], york = TRUE, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)
			OutlierCorrectedMR <- cbind.data.frame(BetaExposure, methodname, "Exclude_initial_outliers", Summary_fit(mod_initialoutliers))
			OutlierCorrectedMR <- cbind.data.frame(BetaExposure, methodname, "Exclude_final_outliers", t(rep(NA, 5)))
		}else{
			OutlierCorrectedMR <- cbind.data.frame(BetaExposure, methodname, "Exclude_outliers", t(rep(NA, 5)))
		}
	
		colnames(OutlierCorrectedMR) <- colnames(OriginalMR)
		MR <- rbind.data.frame(OriginalMR, OutlierCorrectedMR)
		row.names(MR) <- NULL
	}else{
		mod_fixed_all <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york = TRUE, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope=NULL)
	
		SlopeTest <- cbind.data.frame(BetaExposure, paste(methodname,"(base)"), "All_variants", Summary_fit(mod_fixed_all,FixedSlope=FixedSlope))
		names(SlopeTest) <- c("Exposure", "Method", "Data", "Slope", "SE", "df", "Tstat vs 0", "Pvalue vs 0", "Difference vs Fixed", "Tstat vs Fixed", "Pvalue vs Fixed")
		

		mod_fixed_yorkegger <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york = TRUE, MRIntercept=TRUE)
		yegg <- cbind.data.frame(BetaExposure, "York Egger", "All_variants", Summary_fit(mod_fixed_yorkegger,FixedSlope=FixedSlope))
		names(yegg) <- c("Exposure", "Method", "Data", "Slope", "SE", "df", "Tstat vs 0", "Pvalue vs 0", "Difference vs Fixed", "Tstat vs Fixed", "Pvalue vs Fixed")
		SlopeTest <- rbind(SlopeTest,yegg)

		
		mod_fixed_ivw <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york = FALSE)
		ivw <- cbind.data.frame(BetaExposure, "IVW", "All_variants", Summary_fit(mod_fixed_ivw,FixedSlope=FixedSlope))
		names(ivw) <- c("Exposure", "Method", "Data", "Slope", "SE", "df", "Tstat vs 0", "Pvalue vs 0", "Difference vs Fixed", "Tstat vs Fixed", "Pvalue vs Fixed")
		SlopeTest <- rbind(SlopeTest,ivw)
		
		
		mod_fixed_egger <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york = FALSE, MRIntercept=TRUE)
		egger <- cbind.data.frame(BetaExposure, "Egger", "All_variants", Summary_fit(mod_fixed_egger,FixedSlope=FixedSlope))
		names(egger) <- c("Exposure", "Method", "Data", "Slope", "SE", "df", "Tstat vs 0", "Pvalue vs 0", "Difference vs Fixed", "Tstat vs Fixed", "Pvalue vs Fixed")
		SlopeTest <- rbind(SlopeTest,egger)
		
	
		if(length(refOutlier) > 0){
			if(length(refOutlier) < nrow(data)){
				mod_fixed_noOutliers <- RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data[-refOutlier,], york = TRUE, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = NULL)
				SlopeTest_noOut <- cbind.data.frame(BetaExposure, methodname, "Exclude_outliers", Summary_fit(mod_fixed_noOutliers,FixedSlope=FixedSlope))
			}else{
				SlopeTest_noOut <- cbind.data.frame(BetaExposure, methodname, "Exclude_outliers", t(rep(NA,6)))
			}
		}else{
			SlopeTest_noOut <- cbind.data.frame(BetaExposure, methodname, "Exclude_outliers", t(rep(NA,6)))
		}
	
		names(SlopeTest_noOut) <- c("Exposure", "Method", "Data", "Slope", "SE", "df", "Tstat vs 0", "Pvalue vs 0", "Difference vs Fixed", "Tstat vs Fixed", "Pvalue vs Fixed")
		SlopeTest <- rbind(SlopeTest,SlopeTest_noOut)
		row.names(SlopeTest) <- NULL
	}
	
	


	res_tests <- list(`Global Test` = GlobalTest)

	if(IterateOutliers){
		res <- list(`Omnibus Tests` = res_tests, `Slope Estimates` = MR, `Outlier Tests` = OutlierTest, `Outlier Convergence` = OutlierTestControl)
	}else if(!is.null(FixedSlope)){
		res <- list(`Omnibus Tests` = res_tests, `Slope Tests` = SlopeTest, `Outlier Tests` = OutlierTest)
	}else{
		res <- list(`Omnibus Tests` = res_tests, `Slope Estimates` = MR, `Outlier Tests` = OutlierTest)
	}
	

	
	if(printProgress){print("Done!")}
	return(res)
}