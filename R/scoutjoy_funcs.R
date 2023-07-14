###
# Functions
# (should probably move these to a separate file)
###


# Makes expected covariance matrix for betas 
# based on reported SDs and the ldsc covariance intercept.
# Works since cov intercept is covariance of Z scores, 
# and intercepts are (unexplained) variances
makeBetaCovariance <- function(SdExposureVal, SdOutcomeVal, CovIntercept, ExposureIntercept, OutcomeIntercept){

	CovOutcomeExposure <- SdOutcomeVal * SdExposureVal * CovIntercept / sqrt(ExposureIntercept*OutcomeIntercept) 

	matrix(c(SdExposureVal^2, CovOutcomeExposure,
			 CovOutcomeExposure, SdOutcomeVal^2),
			 nrow=2)
}


# note: added flexibility here to fit intercept, but is never used currently
# (unclear how to treat in LOO, among other things)
RegressionFit <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, MRIntercept = FALSE, york = FALSE, CovIntercept = 0, ExposureIntercept = 1, OutcomeIntercept = 1, FixedSlope = NULL){
	
	if(!is.null(FixedSlope)){
		out <- list(
			intercept = NA,
			intercept_se = NA,
			slopes = FixedSlope,
			slopes_se = NA,
			df = NA)
		return(out)
	}
	
	if(!york){
		if(!MRIntercept){
			fit <- stats::lm(paste0(BetaOutcome," ~ 0 + ."), data=data[,c(BetaOutcome,BetaExposure)], weights = data$Weights)
			out <- list(
				intercept = NA,
				intercept_se = NA,
				slopes = fit$coefficients[BetaExposure],
				slopes_se = summary(fit)$coefficients[BetaExposure, 'Std. Error'], 
				df = fit$df.residual)
		}else{
			fit <- stats::lm(paste0(BetaOutcome," ~ ."), data=data[,c(BetaOutcome,BetaExposure)], weights = data$Weights)
			out <- list(
				intercept = fit$coefficients['(Intercept)'],
				intercept_se = summary(fit)$coefficients['(Intercept)', 'Std. Error'],
				slopes = fit$coefficients[BetaExposure],
				slopes_se = summary(fit)$coefficients[BetaExposure, 'Std. Error'],
				df = fit$df.residual)
		}
		return(out)
			
	}
	
	else if(york){
		if(length(BetaExposure) == 1){

			dat <- data.frame(
				X=data[,BetaExposure],
				sX=data[,SdExposure],
				y=data[,BetaOutcome],
				sY=data[,SdOutcome],
				rXY=rep(CovIntercept / sqrt(ExposureIntercept * OutcomeIntercept),nrow(data)))
			
			if(!MRIntercept){
				yfit <- yorktools::york(dat,intercept=0)
				ydf <- nrow(dat)-1
			}else{
				yfit <- yorktools::york(dat)
				ydf <- nrow(dat)-2
			}
			if(!yfit$converge$converged){
			  stop("york regression failed to converge in RegressionFit")
			}
			
			out <- list(
				intercept = yfit$a['a'],
				intercept_se = yfit$a['s[a]'],
				slopes = yfit$b['b'],
				slopes_se = yfit$b['s[b]'],
				df = ydf)
				
				names(out$slopes) <- BetaExposure
				names(out$slopes_se) <- BetaExposure
			
			return(out)
			
		}else{
			stop("york only allows 1 exposure. Initial sanity check failed?")
		}
		
	}
}		

Summary_fit <- function(fit,FixedSlope=NULL){
  if(class(fit)=="york"){
    sumstat <- data.frame(Estimate=fit$b[1], SE=fit$b[2], Tstat=fit$b[1]/fit$b[2], df=fit$df, P=2*stats::pt(abs(fit$b[1]/fit$b[2]),df=fit$df,lower=F))
  }else if(is.null(FixedSlope)){
		sumstat <- data.frame(Estimate=fit$slopes, SE=fit$slopes_se, Tstat=fit$slopes/fit$slopes_se, df=fit$df, P=2*stats::pt(abs(fit$slopes/fit$slopes_se),df=fit$df,lower=F))
	}else{
		sumstat <- data.frame(Estimate=fit$slopes, SE=fit$slopes_se, df=fit$df, T0=fit$slopes/fit$slopes_se, P0=2*stats::pt(abs(fit$slopes/fit$slopes_se),df=fit$df,lower=F), Diff=fit$slopes-FixedSlope, TDiff=(fit$slopes-FixedSlope)/fit$slopes_se,  PDiff=2*stats::pt(abs((fit$slopes-FixedSlope)/fit$slopes_se),df=fit$df,lower=F))
	}
	return(sumstat)
}


Summary_fit_intercept <- function(fit,FixedSlope=NULL){
	if(is.null(FixedSlope)){
		sumstat <- data.frame(Slope=fit$slopes, 
			Slope_SE=fit$slopes_se, 
			Slope_Tstat=fit$slopes/fit$slopes_se, 
			df=fit$df, 
			Slope_P=2*stats::pt(abs(fit$slopes/fit$slopes_se),df=fit$df,lower=F), 
			Intercept=fit$intercept, Intercept_SE=fit$intercept_se, 
			Intercept_Tstat=fit$intercept/fit$intercept_se, 
			Intercept_P=2*stats::pt(abs(fit$intercept/fit$intercept_se),df=fit$df,lower=F))
	}else{
		sumstat <- data.frame(Slope=fit$slopes, 
			Slope_SE=fit$slopes_se, 
			Slope_Tstat=fit$slopes/fit$slopes_se, 
			df=fit$df, 
			Slope_P=2*stats::pt(abs(fit$slopes/fit$slopes_se),df=fit$df,lower=F), 
			Intercept=fit$intercept, Intercept_SE=fit$intercept_se, 
			Intercept_Tstat=fit$intercept/fit$intercept_se, 
			Intercept_P=2*stats::pt(abs(fit$intercept/fit$intercept_se),df=fit$df,lower=F))
	}
	return(sumstat)
}


fit_LOO <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, york=FALSE, yorkResid=FALSE, CovIntercept=0, ExposureIntercept = 1, OutcomeIntercept = 1, FixedSlope = NULL){
	

	CausalEstimate_LOO <- sapply(1:nrow(data), function(i){RegressionFit(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data[-i,], york = york, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)$slopes})
	
	Y <- data[,BetaOutcome]
	X <- data[,BetaExposure]
	
	if(yorkResid){
		if(length(BetaExposure) == 1){
			W <- 1/(data[,SdOutcome]^2 + CausalEstimate_LOO^2 * data[,SdExposure]^2 - 2*CausalEstimate_LOO*data[,SdOutcome]*data[,SdExposure]*CovIntercept/sqrt(ExposureIntercept*OutcomeIntercept))
			RSS <- sum(W*(Y-CausalEstimate_LOO*X)^2, na.rm=TRUE)
		}else{
			stop("york only allows 1 exposure. Initial sanity check failed?")
		}
		
	}else{
		if(length(BetaExposure) == 1)  # if theres only 1 E
			RSS <- sum(data$Weights*(Y - CausalEstimate_LOO * X)^2, na.rm = TRUE)  # Y - vector of slope (leaving out snp) *X
		else
			RSS <- sum(data$Weights*(Y - rowSums(t(CausalEstimate_LOO) * X))^2, na.rm = TRUE)
	}
	

	out <- list(RSS, CausalEstimate_LOO)		
	return(out)
}



getRandomData <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, loo_betas, CovIntercept, ExposureIntercept, OutcomeIntercept, data, yorkResid){
	X <- data[,BetaExposure]
	if(yorkResid){
	  if(length(BetaExposure) == 1){
	    rr <- CovIntercept/sqrt(ExposureIntercept*OutcomeIntercept)
	    W <- 1/(data[,SdOutcome]^2 + loo_betas^2 * data[,SdExposure]^2 - 2*loo_betas*data[,SdOutcome]*data[,SdExposure]*rr)
	    predX <- W*( (data[,SdOutcome]^2 - loo_betas*rr*data[,SdOutcome]*data[,SdExposure])*X + (loo_betas*data[,SdExposure]^2 - rr*data[,SdOutcome]*data[,SdExposure])*data[,BetaOutcome])
	    predY <- loo_betas*predX
	  }else{
	    stop("york only allows 1 exposure. Initial sanity check failed?")
	  }	  
	}else{
	  if(length(BetaExposure) == 1){
		  pred <- loo_betas * X
	  }else{
		  pred <- rowSums(t(loo_betas) * X)
	  }
	  predY <- pred
	  predX <- X
	}

	dataRandomBetas <- t(sapply(1:nrow(data), function(i) {MASS::mvrnorm(n=1, mu=c(predX[i], predY[i]), Sigma=makeBetaCovariance(data[i,SdExposure], data[i,SdOutcome], CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept))}))
	
	dataRandom <- data.frame(dataRandomBetas[,1],data[,SdExposure],dataRandomBetas[,2],data[,SdOutcome],data$Weights)
	colnames(dataRandom) <- c(BetaExposure, SdExposure, BetaOutcome, SdOutcome, "Weights")
	return(dataRandom)
}


PressoOutlierTestRep <- function(data, randomData, BetaOutcome, BetaExposure, SdOutcome, SdExposure, loo_betas, CovIntercept=CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept){
  OutlierTest <- do.call("rbind", lapply(1:nrow(data), function(SNV){
    randomSNP <- do.call("rbind", lapply(randomData, function(mat) mat[SNV, ]))

    rr <- CovIntercept/sqrt(ExposureIntercept*OutcomeIntercept)
    W <- 1/(data[SNV,SdOutcome]^2 + loo_betas[SNV]^2 * data[SNV,SdExposure]^2 - 2*loo_betas[SNV]*data[SNV,SdOutcome]*data[SNV,SdExposure]*rr)
    Dif <- data[SNV, BetaOutcome] - data[SNV, BetaExposure] * loo_betas[SNV]
    Exp <- randomSNP[, BetaOutcome] - randomSNP[, BetaExposure] * loo_betas[SNV]

    pval <- (sum(Exp^2 > Dif^2)+1)/(length(randomData)+1)
    pval <- cbind.data.frame(RSSobs = W*Dif^2, Pvalue = pval)
    return(pval)
  }))
  row.names(OutlierTest) <- row.names(data)
  
  return(OutlierTest)
}


OutlierTestRep <- function(data, randomData, BetaOutcome, BetaExposure, SdOutcome, SdExposure, loo_betas, york=FALSE, yorkResid=york, CovIntercept=0, ExposureIntercept = 1, OutcomeIntercept = 1, FixedSlope = NULL){
  
	OutlierTest <- do.call("rbind", lapply(1:nrow(data), function(SNV){
		randomSNP <- do.call("rbind", lapply(randomData, function(mat) mat[SNV, ]))
		randomLOO <- do.call("rbind", lapply(randomData, function(mat) t(RegressionFit(BetaOutcome=BetaOutcome, BetaExposure=BetaExposure, SdOutcome=SdOutcome, SdExposure=SdExposure, data=data[-SNV,], york = york, CovIntercept = CovIntercept, ExposureIntercept=ExposureIntercept, OutcomeIntercept=OutcomeIntercept, FixedSlope = FixedSlope)$slopes)))
		
		if(length(BetaExposure) == 1){
		  Dif <- data[SNV, BetaOutcome] - data[SNV, BetaExposure] * loo_betas[SNV]
		  Exp <- randomSNP[, BetaOutcome] - randomSNP[, BetaExposure] * randomLOO[SNV]
		} else {
		  Dif <- data[SNV, BetaOutcome] - sum(data[SNV, BetaExposure] * loo_betas[, SNV])
		  Exp <- randomSNP[, BetaOutcome] - rowSums(randomSNP[, BetaExposure] * randomLOO[, SNV])
		}
		pval <- (sum(Exp^2 > Dif^2)+1)/(length(randomData)+1)
		
		if(yorkResid){
		  rr <- CovIntercept/sqrt(ExposureIntercept*OutcomeIntercept)
		  W <- 1/(data[SNV,SdOutcome]^2 + loo_betas[SNV]^2 * data[SNV,SdExposure]^2 - 2*loo_betas[SNV]*data[SNV,SdOutcome]*data[SNV,SdExposure]*rr)
		  
		  pval <- cbind.data.frame(RSSobs = W*Dif^2, Pvalue = pval)
		}else{
		  pval <- cbind.data.frame(RSSobs = Dif^2, Pvalue = pval)
		}
		return(pval)
	}))
	row.names(OutlierTest) <- row.names(data)

	return(OutlierTest)
}

AnalyticOutlierTest <- function(data, putativeOutliers, OutlierNull="estimate", FixedSlope=NULL){
  
  if(!is.null(FixedSlope)){

    pred <- yorktools::pred.york(data, intercept=0, slope=FixedSlope)
    
    # RSS <- ((data$x-pred$xpred)^2/data$sx^2 + (data$y-pred$ypred)^2/data$sy^2 - 2*data$re*(data$x-pred$xpred)*(data$y-ypred$pred)/(data$sx*data$sy))/(1-data$re)^2
    RSS <- yorktools::loss_per_obs.york(data, intercept=0, slope=FixedSlope, func="predicted_values")$loss
    
    l0 <- sapply(1:nrow(data), function(a){mvtnorm::dmvnorm(c(data$x[a],data$y[a]), mean=c(pred$xpred[a], pred$ypred[a]), sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2))})
    l1 <- sapply(1:nrow(data), function(a){mvtnorm::dmvnorm(c(data$x[a],data$y[a]), mean=c(data$x[a],data$y[a]), sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2))})
    
    OutlierTest <-data.frame(RSSobs = RSS, Pvalue = stats::pchisq(-2*log(l0/l1),df=1,lower=FALSE) )

        
  }else if(OutlierNull=="fitted"){
    
    OutlierTest <- data.frame(t(sapply(1:nrow(data), function(a){
      mask <- rep(FALSE, nrow(data))
      mask[putativeOutliers] <- TRUE
      mask[a] <- TRUE
      
      ymod <- yorktools::york(data[!mask,], intercept=0, gridstartvals=T)
      if(!ymod$converge$converged){
        stop("york regression failed to converge for global test")
      }else{
        b <- ymod$b[1]
      }

      pred <- yorktools::pred.york(data[a,], intercept=0, slope=b)
      
      # RSS <- ((data$x[a]-pred$predx)^2/data$sx[a]^2 + (data$y[a]-pred$ypred)^2/data$sy[a]^2 - 2*data$re[a]*(data$x[a]-pred$xpred)*(data$y[a]-pred$ypred)/(data$sx[a]*data$sy[a]))/(1-data$re[a])^2
      RSS <- yorktools::loss_per_obs.york(data[a,], intercept=0, slope=b, func="predicted_values")$loss

      l0 <- mvtnorm::dmvnorm(c(data$x[a],data$y[a]), mean=c(pred$xpred, pred$ypred), sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2))
      l1 <- mvtnorm::dmvnorm(c(data$x[a],data$y[a]), mean=c(data$x[a],data$y[a]), sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2))

      data.frame(RSSobs = RSS, Pvalue = stats::pchisq(-2*log(l0/l1),df=1,lower=FALSE) )
    })))

        
  }else if(OutlierNull=="estimate"){
    
    ll_full <- sapply(1:nrow(data), function(a){
      mask <- rep(FALSE, nrow(data))
      mask[putativeOutliers] <- TRUE
      mask[a] <- FALSE
      ymod <- yorktools::york(data[!mask,], intercept=0, gridstartvals=T)
      if(!ymod$converge$converged){
        stop(paste0("york regression failed to converge for outlier test a=",a))
      }else{
        return(ymod$chi2)
      }
    })
    ll_loo <- sapply(1:nrow(data), function(a){
      mask <- rep(FALSE, nrow(data))
      mask[putativeOutliers] <- TRUE
      mask[a] <- TRUE
      ymod <- yorktools::york(data[!mask,], intercept=0, gridstartvals=T)
      if(!ymod$converge$converged){
        stop(paste0("york regression failed to converge for outlier test a=",a))
      }else{
        return(ymod$chi2)
      }
    })
    
    OutlierTest <- data.frame(RSSobs = ll_full-ll_loo, Pvalue = stats::pchisq(ll_full-ll_loo, df=1, lower=F))
    
  }else{
    stop("Unsupported OutlierNull")
  }
  
  row.names(OutlierTest) <- row.names(data)
  return(OutlierTest)
  
}




