qc.measure<-function(target, background, type = c("zpfactor", "robust_zpfactor", "ssmd", "robust_ssmd", 
                                      "cnr", "gcnr", "pgssmd", "pgssmd_kde", "pgssmd_ksmooth", "pgssmd_plugin", "gssmd", 
                                      "gssmd_subsample", "pgssmd_subsample")) {
  # Author: Seongyong Park (psychemistz@gmail.com). 
  # Released under MIT License. Retrieved from Github, 1-14-25. https://github.com/psychemistz/gssmd 
  # Reference: S. Park, S. Khan, M. Moinuddin and U. M. Al-Saggaf doi: 10.1109/BIBM49941.2020.9313582.

  nsample <- length(target)
  nbins <- round(1+log2(nsample))
  mu_t <- mean(target)
  mu_b <- mean(background)
  med_t <- median(target)
  med_b <- median(background)
  v_t <- var(target)
  v_b <- var(background)
  sd_t <- sd(target)
  sd_b <- sd(background)
  mad_t <- mad(target)
  mad_b <- mad(background)
  
  if(type == "zpfactor") {
    zpfactor <- 1 - ((3 * (sqrt(v_t) + sqrt(v_b)))/abs(mu_t - mu_b))
    return(zpfactor)
  } else if(type == "robust_zpfactor") {
    robust_zpfactor <- 1 - ((3 * (sqrt(mad_t) + sqrt(mad_b)))/abs(med_t - med_b))
    return(robust_zpfactor)
  } else if(type == "ssmd") {
    ssmd <- (mu_t-mu_b)/sqrt(v_t + v_b)
    return(ssmd)
  } else if(type == "robust_ssmd") {
    robust_ssmd <- (med_t-med_b)/sqrt(mad_t^2 + mad_b^2)
    return(robust_ssmd)
  } else if(type == "cnr") {
    cnr <- abs(mu_t-mu_b)/sqrt(v_t + v_b)
    return(cnr)
  } else if(type == "gcnr") {
    x <- range(min(target, background), max(target, background))
    x <- seq(x[1], x[2], length.out=nbins)
    pdf_t <- hist(target, x, plot=FALSE)
    pdf_b <- hist(background, x, plot=FALSE)
    OVL <- sum(apply((rbind(pdf_t$counts/sum(pdf_t$counts), pdf_b$counts/sum(pdf_b$counts))),2,min))
    GCNR <- 1-OVL
    return(GCNR)
  } else if (type == "pgssmd"){
    require(fitdistrplus)
    x <- range(min(target, background), max(target, background))
    x <- seq(x[1], x[2], length.out=nbins)
    fit_b <- fitdist(background, "norm")
    fit_t <- fitdist(target, "norm")
    pd_b <- pnorm(x, fit_b$estimate[1], fit_b$estimate[2])
    pd_t <- pnorm(x, fit_t$estimate[1], fit_t$estimate[2])
    OVL <- sum(apply(cbind(diff(pd_b), diff(pd_t)), 1, min))
    PCNR <- 1-OVL
    pGSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * PCNR
    return(pGSSMD)
  } else if (type == "pgssmd_kde"){
    require(sfsmisc)
    lower <- min(c(target, background)) - 1 
    upper <- max(c(target, background)) + 1
    density_background <- density(target, from=lower, to=upper)
    density_target <- density(background, from=lower, to=upper)    
    OVL <- integrate.xy(density_background$x, pmin(density_background$y, density_target$y))
    PCNR <- 1-OVL
    pGSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * PCNR
    return(pGSSMD)
  } else if (type == "pgssmd_ksmooth"){
    require(KernSmooth)
    require(sfsmisc)
    x <- range(min(target, background), max(target, background))
    density_background <- KernSmooth::bkde(background, range.x = x)
    density_target <- KernSmooth::bkde(target, range.x = x)
    OVL <- integrate.xy(density_background$x, pmin(density_background$y, density_target$y))
    PCNR <- 1-OVL
    pGSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * PCNR
    return(pGSSMD)
  } else if (type == "pgssmd_plugin"){
    min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
      f1 <- dnorm(x, mean=mu1, sd=sd1)
      f2 <- dnorm(x, mean=mu2, sd=sd2)
      pmin(f1, f2)
    }
    OVL = integrate(min.f1f2, -Inf, Inf, mu1=mu_t, mu2=mu_b, sd1=sd_t, sd2=sd_t)$value
    PCNR <- 1-OVL
    pGSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * PCNR
    return(pGSSMD)
  } else if (type == "gssmd"){
    x <- range(min(target, background), max(target, background))
    x <- seq(x[1], x[2], length.out=nbins)
    pdf_t <- hist(target, x, plot=FALSE)
    pdf_b <- hist(background, x, plot=FALSE)
    OVL <- sum(apply((rbind(pdf_t$counts/sum(pdf_t$counts), pdf_b$counts/sum(pdf_b$counts))),2,min))
    GSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * (1-OVL)
    return(GSSMD)
  } else if (type == "gssmd_subsample"){
    gssmd_ss <- function(target, background){
      tg <- sample(target, length(target) * 0.9)
      bg <- sample(background, length(background) * 0.9)
      mu_t <- mean(tg)
      mu_b <- mean(bg)
      x <- range(min(tg, bg), max(tg, bg))
      x <- seq(x[1], x[2], length.out=nbins)
      pdf_t <- hist(tg, x, plot=FALSE)
      pdf_b <- hist(bg, x, plot=FALSE)
      OVL <- sum(apply((rbind(pdf_t$counts/sum(pdf_t$counts), pdf_b$counts/sum(pdf_b$counts))),2,min))
      temp_GSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * (1-OVL)
      return(temp_GSSMD)
    }
    GSSMD <- lapply(1:1000, function(x) gssmd_ss(target, background))
    sd <- sd(unlist(GSSMD))
    GSSMD <- mean(unlist(GSSMD))    
    return(list(GSSMD=GSSMD, sd=sd))
  } else if (type == "pgssmd_subsample"){
    pgssmd_ss <- function(target, background){
      tg <- sample(target, length(target) * 0.9)
      bg <- sample(background, length(background) * 0.9)
      mu_t <- mean(tg)
      mu_b <- mean(bg)
      x <- range(min(tg, bg), max(tg, bg))
      x <- seq(x[1], x[2], length.out=nbins)
      fit_b <- fitdist(background, "norm")
      fit_t <- fitdist(target, "norm")
      pd_b <- pnorm(x, fit_b$estimate[1], fit_b$estimate[2])
      pd_t <- pnorm(x, fit_t$estimate[1], fit_t$estimate[2])
      OVL <- sum(apply(cbind(diff(pd_b), diff(pd_t)), 1, min))
      PCNR <- 1-OVL
      temp_GSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * PCNR
      return(temp_GSSMD)
    }
    pGSSMD <- lapply(1:1000, function(x) pgssmd_ss(target, background))
    sd <- sd(unlist(pGSSMD))
    pGSSMD <- mean(unlist(pGSSMD))    
    return(list(GSSMD=pGSSMD, sd=sd))
  }else {
    warning("Type is incorrect")
  }
}
