# ---
# title: Modeling (Global CPR)
# author: Bryan Alpecho
# date: 2025-
# output: model summary and predictions
# ---

library(tidyverse)
library(glmmTMB)
library(GGally)
library(MuMIn)
library(DHARMa)
library(marginaleffects)
library(reprex)


# Import data
  df <- read_rds("data_input/global/df_complete_11062025.rds")
  df_withGroup <- read_rds("data_input/global/df_complete_withGroup_11062025.rds")
  
  #mutate Chla_sqrt_withCeiling
  df <- df %>% 
    mutate(Chla_sqrt_withCeilingAt3 = case_when(
      Chla >= 3 ~ 3,
      .default = Chla
    )) %>%  
    mutate(Chla_sqrt_withCeilingAt5 = case_when(
      Chla >= 5 ~ 5,
      .default = Chla
    )) 
  
  df_withGroup <- df_withGroup %>% 
    mutate(Chla_sqrt_withCeilingAt3 = case_when(
      Chla >= 3 ~ 3,
      .default = Chla
    )) %>%  
    mutate(Chla_sqrt_withCeilingAt5 = case_when(
      Chla >= 5 ~ 5,
      .default = Chla
    )) 
  
  #to remove NA values
  df_noNA <- df %>% filter(!is.na(Chla)) 
  df_noNA_withGroup <- df_withGroup %>% filter(!is.na(Chla)) 
  df_1997 <- df_withGroup %>% filter(date(SampleTime_UTC) >= "1997-09-01") %>% 
    filter(!is.na(Chla)) #should be the same # of objects as df_withGroup_noNA
  
## Exploratory data analysis
  #Look at density plots
    #Filter-feeder, Carnivore, Omnivore
    par(mfrow = c(1,3))
    hist((df$Carni_sum/df$Omni_sum))
    hist(df$FF_sum, xlab = "Filter feeder abundance", main = "", prob = TRUE, 
         cex.lab = 0.9, col = "grey", 
         border = "black") # Heavily LEFT skewed
    pu <- par("usr")
    text(pu[2], pu[4], "(A)", adj = c(3.5,1), cex = 1.5)
    
    hist(df$Carni_sum, xlab = "Carnivore abundance", main = "", prob = TRUE, 
         cex.lab = 0.9, col = "grey", 
         border = "black") # Heavily LEFT skewed
    pu <- par("usr")
    text(pu[2], pu[4], "(B)", adj = c(3.5,1), cex = 1.5)
    
    hist(df$Omni_sum, xlab = "Omnivore abundance", main = "", prob = TRUE, 
         cex.lab = 0.9, col = "grey", 
         border = "black") # Heavily RIGHT skewed
    pu <- par("usr")
    text(pu[2], pu[4], "(C)", adj = c(4,1), cex = 1.5)
    
    #Filter feeder and Non-Filter feeder
    par(mfrow = c(1,2))
    hist(df$FF_sum, xlab = "Filter feeder abundance", main = "", prob = TRUE, 
         cex.lab = 0.9, col = "grey", 
         border = "black") # Heavily LEFT skewed
    pu <- par("usr")
    text(pu[2], pu[4], "(A)", adj = c(3.5,1), cex = 1.5)
    
    hist(df$NonFF_sum, xlab = "Non-Filter Feeder abundance", main = "", prob = TRUE, 
         cex.lab = 0.9, col = "grey", 
         border = "black") # 
    pu <- par("usr")
    text(pu[2], pu[4], "(B)", adj = c(4,1), cex = 1.5)


  # Look at spearman correlations and distributions
    ggpairs_plot <- GGally::ggpairs(df, columns = c("Latitude", "Chla"), progress = FALSE, 
                                    upper = list(continuous = GGally::wrap("cor", method = "spearman", stars = FALSE))
    )
    ggpairs_plot

  #Density Plots - Feeding Group proportions ~ Chlorophyll-a
    par(mfrow = c(1,2))
    plot(density((df$Omni_sum/df$Carni_sum), bw = 0.37), #Omnivores
         col = "darkgreen", main = "", xlab = "Omnivores", 
         lwd = 1.5, cex.axis=1.2, cex.lab = 1.3)
    lines(density((df$Carni_sum/df$Omni_sum), bw = 0.37), #Carnivores
          lwd = 1.5, col = "darkred")
    lines(density((df$FF_sum/df$NonFF_sum), bw = 0.37), #Filter-feeders
          lwd = 1.5, col = "darkblue")
    legend(7, 0.3, legend=c("Omnivores", "Carnivores", "Filter-feeders"),
           col=c("darkgreen", "darkred","darkblue"), lty = c(1,1), lwd = c(1.5, 1.5))
    
    with(df, {
      par(mfrow = c(1,1))
      
      # Chl-a
      hist(Chla, xlab = "Chl-a", main = "", prob = TRUE, 
           cex.lab = 0.9, col = "grey", 
           border = "black") # Heavily right skewed
      pu <- par("usr")
      text(pu[2], pu[4], "(A)", adj = c(1.5,1), cex = 1.5)
      
    })


# Run models 
  #zero-one-inflated beta regression model
  
  #regular beta reg on transformed data
    #SVT - Smithson & Verkuilen (2006) transformation
  Filter_mdl_beta <- glmmTMB(RFF_SVT ~ Chla_sqrt + Survey +
                          (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                          data = df, family = beta_family(link = "logit"))
  
  Omni_mdl_beta <- glmmTMB(ROC_SVT ~ Chla_sqrt + Survey +
                        (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                        data = df, family = beta_family(link = "logit"))
  
  Carni_mdl_beta <- glmmTMB(RCO_SVT ~ Chla_sqrt + Survey +
                         (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                       data = df, family = beta_family(link = "logit"))
  
  #zero-inflated beta regression model 
    #SVT transformation for 1s only; non-transformed [0,1) 
  Filter_mdl_zib <- glmmTMB(RFF_SVT_zib ~ Chla_sqrt_withCeilingAt3 + Survey + (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                        ziformula = ~1,
                        data = df, family = beta_family(link = "logit"))
  
  Omni_mdl_zib <- glmmTMB(ROC_SVT_zib ~ Chla_sqrt + Survey + (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                      ziformula = ~1,
                      data = df, family = beta_family(link = "logit"))
  
  Carni_mdl_zib <- glmmTMB(RCO_SVT_zib ~ Chla_sqrt + Survey + (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                       ziformula = ~1,
                       data = df, family = beta_family(link = "logit"))
  
  #with Group; regular beta
    Omni_withGroup_mdl_beta <- glmmTMB(ROC_SVT ~ Chla_sqrt * Group + Survey +
                               (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                             data = df_withGroup, family = beta_family(link = "logit"))
    
    Carni_withGroup_mdl_beta <- glmmTMB(RCO_SVT ~ Chla_sqrt * Group + Survey +
                                (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                              data = df_withGroup, family = beta_family(link = "logit"))
  #with Group; zero-inflated beta
    Omni_withGroup_mdl_zib <- glmmTMB(ROC_SVT_zib ~ Chla_sqrt_withCeilingAt3 * Group + Survey + (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                            ziformula = ~1,
                            data = df_withGroup, family = beta_family(link = "logit"))
    
    Carni_withGroup_mdl_zib <- glmmTMB(RCO_SVT_zib ~ Chla_sqrt_withCeilingAt3 * Group + Survey + (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst), 
                             ziformula = ~1,
                             data = df_withGroup, family = beta_family(link = "logit"))
    
  
  #sqrt-chla
    Filter_mdl <- glmmTMB(cbind(FF_sum, NonFF_sum) ~ Chla_sqrt + Survey +
                            (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst),
                          data = df, family = binomial)
    
    Omni_mdl <- glmmTMB(cbind(Omni_sum, Carni_sum) ~ Chla_sqrt + Survey +
                          (1 +Tow_Days | Survey: Tow_No) + (1 | Longhurst),
                        data = df, family = binomial)
    
    Carni_mdl <- glmmTMB(cbind(Carni_sum, Omni_sum) ~ Chla_sqrt + Survey +
                           (1 + Tow_Days | Tow_No) + (1 | Longhurst),
                         data = df, family = binomial)
  
  #log-transformed chla; no ceiling
    Filter_mdl_log <- glmmTMB(cbind(FF_sum, NonFF_sum) ~ Chla_eightday_log + Survey +
                            (1 + Tow_Days | Survey: Tow_No) + (1 | Longhurst),
                          data = df, family = binomial)
    
    Omni_mdl_log <- glmmTMB(cbind(Omni_sum, Carni_sum) ~ Chla_eightday_log + Survey +
                          (1 +Tow_Days | Survey: Tow_No) + (1 | Longhurst),
                        data = df_noNA, family = binomial)
    
    Carni_mdl_log <- glmmTMB(cbind(Carni_sum, Omni_sum) ~ Chla_eightday_log + Survey +
                           (1 + Tow_Days | Tow_No) + (1 | Longhurst),
                         data = df, family = binomial)
  
  #with Ceiling
    #withTowDays
  
    Filter_mdl_withCeiling <- glmmTMB(cbind(FF_sum, NonFF_sum) ~ Chla_withCeilingAt5 + Survey +
                            (1 + Tow_Days | Tow_No) + (1 | Longhurst),
                          data = df, family = binomial)
    
    Omni_mdl_withCeiling <- glmmTMB(cbind(Omni_sum, Carni_sum) ~ Chla_withCeilingAt5 + Survey +
                          (1 +Tow_Days | Tow_No) + (1 | Longhurst),
                        data = df, family = binomial)
    
    Carni_mdl_withCeiling <- glmmTMB(cbind(Carni_sum, Omni_sum) ~ Chla_withCeilingAt5 + Survey +
                           (1 + Tow_Days | Tow_No) + (1 | Longhurst),
                         data = df, family = binomial)
    

#Visualize the predictions    
  #Filter-feeders proportion
    pop_preds_FF<- predictions(Filter_mdl_zib, 
                               newdata = datagrid(Survey = unique(df_noNA$Survey), 
                                                  Chla_sqrt = seq(floor(min(df_noNA$Chla_sqrt)),ceiling(max(df_noNA$Chla_sqrt)),0.05)), 
                               re.form = NA) # Zeros out random effects
    
    ff_plot <- ggplot(data = pop_preds_FF) + 
      geom_point(data = df_noNA, 
                 aes(x = Chla_sqrt, y = RFF_SVT_zib, colour = Survey), alpha = 0.1) +
      geom_ribbon(aes(x = Chla_sqrt, y = estimate, 
                      ymin = conf.low, ymax = conf.high, 
                      colour = Survey, fill = Survey), alpha = 0.3, colour = NA) + 
      geom_line(aes(x = Chla_sqrt, y = estimate, colour = Survey)) + 
      ylab("Filter-feeders proportion")
    
    ggsave("Output/predictions/global/Filter-feeder-prop-zib-12062025.png", plot = ff_plot, 
           width = 8, height = 5, dpi = 300)
    
  #Omnivore proportion
    pop_preds_omni <- predictions(Omni_withGroup_mdl_zib, 
                                  newdata = datagrid(Chla_sqrt = seq(floor(min(df_noNA_withGroup$Chla_sqrt)),ceiling(max(df_noNA_withGroup$Chla_sqrt)),0.05),
                                                     Group = unique(df_noNA_withGroup$Group)), 
                                  re.form = NA) # Zeros out random effects
    
    omni_plot <- ggplot(data = pop_preds_omni) + 
      geom_point(data = df_noNA_withGroup, 
                 aes(x = Chla_sqrt, y = ROC_SVT_zib, colour = Group), alpha = 0.1) +
      geom_ribbon(aes(x = Chla_sqrt, y = estimate, 
                      ymin = conf.low, ymax = conf.high, 
                      colour = Group, fill = Survey), alpha = 0.3, colour = NA) + 
      geom_line(aes(x = Chla_sqrt, y = estimate, colour = Group)) + 
      ylab("Omnivores proportion")
    
    ggsave("Output/predictions/global/Omnivore-prop-withGroup-zib-11062025.png", plot = omni_plot, 
           width = 8, height = 5, dpi = 300)
    
  #Carnivore proportion
    #Visualize the predictions
    pop_preds_carni <- predictions(Carni_withGroup_mdl_zib, 
                                   newdata = datagrid(Chla_sqrt = seq(floor(min(df_noNA_withGroup$Chla_sqrt)),ceiling(max(df_noNA_withGroup$Chla_sqrt)),0.05),
                                                      Group = unique(df_noNA_withGroup$Group)),  
                                   re.form = NA) # Zeros out random effects
    
    carni_plot <- ggplot(data = pop_preds_carni) + 
      geom_point(data = df_noNA_withGroup, 
                 aes(x = Chla_sqrt, y = RCO_SVT_zib, colour = Group), alpha = 0.1) +
      geom_ribbon(aes(x = Chla_sqrt, y = estimate, 
                      ymin = conf.low, ymax = conf.high, 
                      colour = Group, fill = Group), alpha = 0.3, colour = NA) + 
      geom_line(aes(x = Chla_sqrt, y = estimate, colour = Group)) + 
      ylab("Carnivores proportion")
    
    ggsave("Output/predictions/global/Carnivore-prop-withGroup-zib-11062025.png", plot = carni_plot,
           width = 8, height = 5, dpi = 300)

#Assess model fit
    summary(Filter_mdl)
    summary(Carni_mdl)
    summary(Omni_mdl)
    summary(Filter_mdl)$coefficients
    summary(Omni_mdl)$coefficients
    summary(Carni_mdl)$coefficients
    MuMIn::r.squaredGLMM(Filter_mdl)
    MuMIn::r.squaredGLMM(Omni_mdl)
    MuMIn::r.squaredGLMM(Carni_mdl)
    
    #
    summary(Filter_mdl_beta)
    summary(Carni_mdl_beta)
    summary(Omni_mdl_beta)
    MuMIn::r.squaredGLMM(Filter_mdl_beta)
    MuMIn::r.squaredGLMM(Omni_mdl_beta)
    MuMIn::r.squaredGLMM(Carni_mdl_beta)
    #
    Filter_mdl_zib_coeff <- summary(Filter_mdl_zib)$coefficients$cond %>% as.data.frame()
    Carni_mdl_zib_summary <- summary(Carni_mdl_zib)$coefficients$cond %>% as.data.frame()
    Omni_mdl_zib_summary <- summary(Omni_mdl_zib)$coefficients$cond %>% as.data.frame()
    Filter_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib)
    Omni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib)
    Carni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib)
    #
    summary(Carni_withGroup_mdl_beta)
    summary(Omni_withGroup_mdl_beta)
    MuMIn::r.squaredGLMM(Omni_withGroup_mdl_beta)
    MuMIn::r.squaredGLMM(Carni_withGroup_mdl_beta)
    #
    Carni_withGroup_mdl_zib_coeff <- summary(Carni_withGroup_mdl_zib)$coefficients$cond %>% as.data.frame()
    Omni_withGroup_mdl_zib_coeff <- summary(Omni_withGroup_mdl_zib)$coefficients$cond %>% as.data.frame()
    Carni_withGroup_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Omni_withGroup_mdl_zib)
    Omni_withGroup_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Carni_withGroup_mdl_zib)
    #
    
    ## COMPARISON OF MODELS
    tibble_models <- tibble(
      mdl = c("Filter_mdl_zib", "Omni_withGroup_mdl_zib", "Carni_withGroup_mdl_zib"),
      Estimate = c(Filter_mdl_zib_coeff$Estimate[2],  Omni_withGroup_mdl_zib_coeff$Estimate[2], Carni_withGroup_mdl_zib_coeff$Estimate[2]),
      StdError = c(Filter_mdl_zib_coeff$`Std. Error`[2],  Omni_withGroup_mdl_zib_coeff$`Std. Error`[2], Carni_withGroup_mdl_zib_coeff$`Std. Error`[2]),
      z_value = c(Filter_mdl_zib_coeff$`z value`[2], Omni_withGroup_mdl_zib_coeff$`z value`[2], Carni_withGroup_mdl_zib_coeff$`z value`[2]),
      p_value = c(Filter_mdl_zib_coeff$`Pr(>|z|)`[2], Omni_withGroup_mdl_zib_coeff$`Pr(>|z|)`[2], Carni_withGroup_mdl_zib_coeff$`Pr(>|z|)`[2]))
      PseudoR2_marginal = c(Filter_mdl_zib_R2[1], Omni_withGroup_mdl_zib_R2[1], Carni_withGroup_mdl_zib_R2[1]),
      PseudoR2_conditional = c(Filter_mdl_zib_R2[2], Omni_withGroup_mdl_zib_R2[2], Carni_withGroup_mdl_zib_R2[2]))
      
    
#Distribution of random intercepts
    RE <- names(Filter_mdl$sdr$par.random) == "b" # Select  random intercepts for each Longhurst province (this slot has fixed parameters too)
    Mn <- Filter_mdl$sdr$par.random[RE] 
    
    Sd <- sqrt(Filter_mdl$sdr$diag.cov.random)[RE] # Standard deviations of each random intercept
    
    df <- data.frame(Rows = 1:length(Mn),
                     Survey = levels(Filter_mdl$frame$Survey),
                     Mn = Mn, 
                     CIlow =  Mn - 1.96*Sd, # Calculate Confidence Intervals
                     CIhigh = Mn + 1.96*Sd) %>% 
      arrange(Mn) %>% # Arrange in ascending order
      mutate(Survey = factor(Survey, levels = unique(Survey))) # Order factors as they are in df, not alphabetically
    
    ggplot(data = df) +
      geom_errorbarh(aes(xmin = CIlow, xmax = CIhigh, y = Longhurst)) +
      geom_point(aes(x = Mn, y = 1:25))

##plot by Survey
    library(marginaleffects)
    
    plot_predictions(Filter_mdl, 
                     newdata = df, 
                     type = "response", 
                     x = "Chla", 
                     by = "Survey") +
      labs(x = "Chlorophyll-a", y = "Probability of Filter-feeders") +
      theme_bw()
    
    #Filter feeders
    df$pred_fltr_mdl <- predict(Filter_mdl, type = "response") 
    #Error in `$<-`:
    #! Assigned data `predict(Filter_mdl, type = "response")` must be compatible with existing data.
    ggplot(df) +
      aes(x = Chla, y = FF_sum, color = Survey)+
      geom_point(size = 0.2, alpha = 0.3) +
      geom_line(aes(y = pred_fltr_mdl), size = 1) + 
      facet_grid(.~Survey) + 
      theme_bw()
    
    #Carnivores
    df$pred_carni_mdl <- predict(Carni_mdl, type = "response")
    ggplot(df) +
      aes(x = Chla, y = Carni_sum, color = Survey)+
      geom_point(size = 0.2, alpha = 0.3) +
      geom_line(aes(y = pred_carni_mdl), size = 1) + 
      facet_grid(.~Survey) + 
      theme_bw()
    
    #Omnivores
    df$pred_omni_mdl <- predict(Omni_mdl, type = "response")
    ggplot(df) +
      aes(x = Chla, y = Omni_sum, color = Survey)+
      geom_point(size = 0.2, alpha = 0.3) +
      geom_line(aes(y = pred_omni_mdl), size = 1) + 
      facet_grid(.~Survey) + 
      theme_bw()

#Assumptions

    dharma = function(model){
      simulationOutput <- simulateResiduals(fittedModel = model, plot = FALSE)
      plotQQunif(simulationOutput, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
      plotResiduals(simulationOutput, quantreg = FALSE) 
    }
    
    dharma(Filter_mdl)
    dharma(Omni_mdl)
    dharma(Carni_mdl)

  # Fitted vs residual plot
  plot(Filter_mdl, type = c("p","smooth"), col.line=1)
  
  # Scale-location plot
  plot(Filter_mdl,
       sqrt(abs(resid(.)))~fitted(.),
       type = c("p","smooth"), col.line=1)
  
  ## QQplot
  qqnorm(resid(Filter_mdl, type = "deviance"), main = "")
  
  # Residuals vs leverage
  plot(Filter_mdl, rstudent(.) ~ hatvalues(.))  
  
  summary(Filter_mdl)
  
  #Carni_model
  # Fitted vs residual plot
  plot(Carni_mdl, type = c("p","smooth"), col.line=1)
  
  # Scale-location plot
  plot(Carni_mdl,
       sqrt(abs(resid(.)))~fitted(.),
       type = c("p","smooth"), col.line=1)
  
  ## QQplot
  qqnorm(resid(Carni_mdl, type = "deviance"), main = "")
  
  # Residuals vs leverage
  plot(Carni_mdl, rstudent(.) ~ hatvalues(.))  
  
  summary(Carni_mdl)
  
  #Omni_model
  # Fitted vs residual plot
  plot(Omni_mdl, type = c("p","smooth"), col.line=1)
  
  # Scale-location plot
  plot(Omni_mdl,
       sqrt(abs(resid(.)))~fitted(.),
       type = c("p","smooth"), col.line=1)
  
  ## QQplot
  qqnorm(resid(Omni_mdl, type = "deviance"), main = "")
  
  # Residuals vs leverage
  plot(Omni_mdl, rstudent(.) ~ hatvalues(.))  
  
  summary(Omni_mdl)
  
  #NonFilter_model
  # Fitted vs residual plot
  plot(Filter_mdl, type = c("p","smooth"), col.line=1)
  
  # Scale-location plot
  plot(Filter_mdl,
       sqrt(abs(resid(.)))~fitted(.),
       type = c("p","smooth"), col.line=1)
  
  ## QQplot
  qqnorm(resid(Filter_mdl, type = "deviance"), main = "")
  
  # Residuals vs leverage
  plot(Filter_mdl, rstudent(.) ~ hatvalues(.))  
  
  summary(Filter_mdl)

