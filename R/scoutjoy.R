#' Significant Cross-trait OUtliers and Trends in JOint York regression (SCOUTJOY)
#'
#' @description 
#' Runs SCOUTJOY (Elliott et al., medrxiv) on input effect sizes and standard error to identify a primary regression slope and potential outliers.
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
#' @param SdOutcome
#' String, giving the column name in `data` for the standard error (se) of the effect sizes on the outcome (y-axis trait)
#' 
#' @param SdExposure 
#' String, giving the column name in `data` for the standard error (se) of the effect sizes on the exposure (x-axis trait)
#' 
#' @param data
#' Data frame with the specified column names containing effect sizes and standard errors
#' 
#' @param ids 
#' (Optional) string giving the column name in `data` of variant IDs to use for labeling results tables
#' 
#' @param CovIntercept 
#' Intercept from LD score regression for genetic correlation between the two phenotypes. Default = 0 corresponds to assuming no sample overlap or correlated confounding. 
#' 
#' @param ExposureIntercept
#' Intercept from univariate LD score regression of the exposure phenotype. Used for estimate of correlated error (default = 1).
#' 
#' @param OutcomeIntercept
#' Intercept from univariate LD score regression of the outcome phenotype. Used for estimate of correlated error (default = 1). 
#' 
#' @param OutlierNull
#' String, either "estimate" or "fitted", indicating which null hypothesis to use in outlier testing. Recommended default is "estimate". See Details.
#' 
#' @param SignifThreshold 
#' Significance threshold for testing outlier status,prior to multiple testing correction (default = 0.05).
#' 
#' @param NullReplicates 
#' Number of simulated null replicates to use for hypothesis testing. If NA (default), analytic tests are used rather than simulated replicates. 
#' 
#' @param maxOutlierReps
#' Maximum number of iterations to use when idenfiying stable outlier set (default = 10). 
#' 
#' @param seed 
#' Random seed, for use with simulated null. 
#' 
#' @param printProgress 
#' Logical, whether to print logging messages.
#'
#'
#' @return
#' Structured list with elements:
#' \describe{
#'
#'   \item{\code{Global Test}}{Data frame with results for the test of overall heterogeneity in the York regression. Provides values `RSSobs`, the weighted residual sum of squares observed, and the corresponding `Pvalue` for the test of excess heterogeneity.}
#'   
#'   \item{\code{Slope Estimates}}{Data frame with regression results. Includes York regression results for all variants and excluding outliers, as well as Inverse Variance Weighted (IVW) and Egger regression for comparison.}
#'   
#'   \item{\code{Outlier Tests}}{Data frame containing the input effect sizes and standard errors along with the observed weighted residuals and corresponding P-values for the initial York regression (including all variants) and the final York regression after outlier exclusions. NA values indicate failed convergence. If `ids` were provided, they are present as the row names for this table.}
#'   
#'   \item{\code{Outlier Convergence}}{List with information on `Iterations`, whether the outlier testing `Converged`, the number of iterations performed (`NumberIterations`), and the signficance threshold for identifying outliers after Bonferroni correction for the number of included variants (`SignifTheshold`). `Iterations` is a data frame where each column is a variant and each row is an iteration, with true/false values indicating whether the given variants was a putative outlier in the given iteration. If the procedure converged, the final outlier selection is given by the last row.}
#'
#' }
#'
#' @export
#'
#' @references
#'
#' Elliott A, Walters RK, Pirinen M, Kurki M, Junna N, Goldstein J, Reeve MP, Siirtola H, Lemmelä, Turley P, FinnGen, Palotie A, Daly M, Widén E. (2023). Distinct and shared genetic architectures of gestational diabetes mellitus and Type 2 diabetes mellitus. medRxiv. doi: \url{https://doi.org/10.1101/2023.02.16.23286014}
#' 
#' @examples
#'
#' # data generation settings
#' n <- 50
#' a0 <- 0
#' b <- .75
#' re <- 0.3
#' 
#' # generate underlying data without error
#' xtrue <- rnorm(n,.5,.3)
#' ytrue <- a0 + b*xtrue
#' 
#' # add outliers
#' ytrue[1:3] <- ytrue[1:3] + .35*xtrue[1:3]
#' 
#' # define SEs
#' xsd <- sqrt(runif(n,.05,1))
#' ysd <- sqrt(runif(n,.05,0.5))
#' 
#' # generate data with SEs
#' obs <- t(sapply(1:n, function(i) MASS::mvrnorm(1, mu=c(xtrue[i],ytrue[i]), Sigma = matrix(c(xsd[i]^2, re*xsd[i]*ysd[i], re*xsd[i]*ysd[i], ysd[i]^2),2,2))))
#' dat <- data.frame(x=obs[,1],sx=xsd,y=obs[,2],sy=ysd,re=re)
#'
#' # run scoutjoy
#' fit <- scoutjoy("y", "x", "sy", "sx", dat, CovIntercept = re, ExposureIntercept = 1, OutcomeIntercept = 1)
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
    GlobalTest <- data.frame(RSSobs = mod_all$chi2, Pvalue = mod_all$p.value)
    
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
  
  	GlobalTest <- data.frame(RSSobs = RSSobs, Pvalue = (sum(RSSexp[1, ] > RSSobs)+1)/(NullReplicates+1))

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
	
	

	if(IterateOutliers){
		res <- list(`Global Test` = GlobalTest, `Slope Estimates` = MR, `Outlier Tests` = OutlierTest, `Outlier Convergence` = OutlierTestControl)
	}else if(!is.null(FixedSlope)){
		res <- list(`Global Test` = GlobalTest, `Slope Tests` = SlopeTest, `Outlier Tests` = OutlierTest)
	}else{
		res <- list(`Global Test` = GlobalTest, `Slope Estimates` = MR, `Outlier Tests` = OutlierTest)
	}
	
	class(res) <- "scoutjoy"
	if(printProgress){print("Done!")}
	return(res)
}


#################################

#' Plot SCOUTJOY results
#'
#' @param fit
#' output of \link{scoutjoy}`
#' 
#' @param label
#' string, whether to label "outliers", "all", or "none".
#' 
#' @param colors
#' list of colors to use for plotting base regression, regression excluding outliers, outlier points, and other points, respectively.
#' 
#' @param se_bars
#' logical, whether to include standard error bars
#' 
#' @param legend
#' logical, whether to include a legend
#' 
#' @param legend_names
#' list of names to use for slopes and outliers (only affects legend)
#' 
#' @param se_length
#' length for caps on SE intervals; passed to \link{arrows}
#'
#' @param ...
#' additional parameters to be passed directly to \link{plot()}
#'
#'
#' @export plot.scoutjoy
#' @export
#'

plot.scoutjoy <- function(fit, label="outliers", colors=c("darkred","darkblue","darkorange","black"), se_bars=TRUE, legend=TRUE, legend_names=c("all variants","excluding outliers","outliers"), se_length=0.025, ...){
  
  args <- list(...)
  
  if(!(class(fit)=="scoutjoy") || !is.list(fit) || !hasName(fit,"Slope Estimates") || !hasName(fit,"Outlier Tests") || !hasName(fit,"Outlier Convergence")){
    stop("fit does not appear to be scoutjoy results output.")
  }
    
  if(!label %in% c("outliers", "all", "none")){
    stop("label must be one of: 'outliers', 'all', or 'none'.")
  }
  
  
  x <- fit$`Outlier Tests`[,1]
  y <- fit$`Outlier Tests`[,3]
  sx <- fit$`Outlier Tests`[,2]
  sy <- fit$`Outlier Tests`[,4]
  
  slope_base <- fit$`Slope Estimates`[fit$`Slope Estimates`$Method=="York (base)","Estimate"]
  if("Exclude_final_outliers" %in% fit$`Slope Estimates`$Data){
    slope_outlier <- fit$`Slope Estimates`[fit$`Slope Estimates`$Data=="Exclude_final_outliers","Estimate"]
  }else{
    slope_outlier <- NA
  }
  
  isout <- t(fit$`Outlier Convergence`$Iterations[nrow(fit$`Outlier Convergence`$Iterations),])
  
  if(!("xlab" %in% names(args))){
    args$xlab <- colnames(fit$`Outlier Tests`)[1]
  }
  
  if(!("ylab" %in% names(args))){
    args$ylab <- colnames(fit$`Outlier Tests`)[3]
  }

  if(!("xlim" %in% names(args))){
    args$xlim <- c(min(0,min(x-sx)), 1.1*max(x+sx))
  }
  
  if(!("ylim" %in% names(args))){
    args$ylim <- c(min(0,min(y-sy)), 1.1*max(y+sy))
  }
  
  
  if(!("main" %in% names(args))){
    args$main <- ""
  }
  
  if(!("pch" %in% names(args))){
    args$pch <- 20
  }
  
  if(!("cex" %in% names(args))){
    args$cex <- 1
  }
  
  if(!("lwd" %in% names(args))){
    args$lwd <- 2
  }
  
  if(!("las" %in% names(args))){
    args$las <- 1
  }
  
  
  if(!("bty" %in% names(args))){
    args$bty <- 'l'
  }
  
  if(!("mar" %in% names(args))){
    mar <- c(3.5, 3.5, 1.5, 1.5)
  }else{
    mar <- args$mar
  }
  
  if(!("mgp" %in% names(args))){
    mgp <- c(2.25,1,0)
  }else{
    mgp <- args$mgp
  }
  
  mar_bak <- par()$mar
  mgp_bak <- par()$mgp
  
  par(mar=mar, mgp=mgp)
  
  args$x <- 0
  args$y <- 0
  args$col <- rgb(0,0,0,0)
  do.call(plot, args)
  

  abline(h=0,col="grey20",lwd=0.5)
  abline(v=0,col="grey20",lwd=0.5)
  abline(0,1,col="grey80",lty=2)
  abline(0,-1,col="grey80",lty=2)
  
  
  if(se_bars){    
    arrows(x0=x, x1=x, y0=y-sy, y1=y+sy, length=se_length, angle=90, code=3, col="grey50")
    arrows(x0=x-sx, x1=x+sx, y0=y, y1=y, length=se_length, angle=90, code=3, col="grey50")
  }
  
  points(x[!isout], y[!isout],
         pch=args$pch,
         col=colors[4],
         cex=args$cex)
  
  points(x[isout], y[isout],
         pch=args$pch,
         col=colors[3],
         cex=args$cex)
  

  abline(0, slope_base, col=colors[1], lwd=args$lwd)
  if(!is.na(slope_outlier)){
    abline(0, slope_outlier, col=colors[2], lwd=args$lwd)
  }

  if(legend){
    if(is.na(slope_outlier)){
      legend("topright",
             fill=colors[1],
             border=colors[1],
             cex = 0.9*args$cex,
             bg="white",
             box.col="white",
             legend = legend_names[1])
    }else if(length(legend_names)==2){
      legend("topright", 
             fill=colors[1:2],
             border=colors[1:2],
             cex = 0.9*args$cex, 
             bg="white",
             box.col="white", 
             legend = legend_names)
    }else if(length(legend_names)==3){
      legend("topright", 
             col=colors[1:3],
             lwd=c(rep(args$lwd,2),NA),
             pch=c(NA,NA,args$pch),
             cex = 0.9*args$cex, 
             bg="white",
             box.col="white", 
             legend = legend_names)
    }else if(length(legend_names)==4){
      legend("topright", 
             col=colors[1:4],
             lwd=c(rep(args$lwd,2),NA,NA),
             pch=c(NA,NA,rep(args$pch,2)),
             cex = 0.9*args$cex, 
             bg="white",
             box.col="white", 
             legend = legend_names)
  }
  }
  
  
  if(label=="all"){
    text(x=x, 
         y=y, 
         labels=rownames(fit$`Outlier Tests`), 
         cex=0.8, adj=c(-.4,-.25))
  }else if(label=="outliers"){
    text(x=x[isout], 
         y=y[isout], 
         labels=rownames(fit$`Outlier Tests`)[isout], 
         cex=0.8, adj=c(-.4,-.25))
  }
  
  
  par(mar=mar_bak, mgp=mgp_bak)
  
}

