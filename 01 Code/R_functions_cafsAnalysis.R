## Convenience functions used by the code 
# Author: Chris McGraw, chris.mcgraw@gmail.com, christopher.mcgraw@northwestern.edu
# Version 1. 1-14-2025
# For release with manuscript, McGraw et al EJP 2025.
# Reference (Biorxiv) doi: https://doi.org/10.1101/2024.08.01.606228

fct_trainModel <- function(trainData, testData, modelFormula,
                           cv_number=10,
                           positiveClass="X1", metric="accuracy",
                           preprocess=TRUE, weighted=TRUE)
{
  spikes_forCV <- trainData;
  
  spikes_forCV$Conditions_names <- as.factor(spikes_forCV$Conditions_code)
  levels(spikes_forCV$Conditions_names) <- make.names(levels(spikes_forCV$Conditions_names))
  
  testData$Conditions_names <- as.factor(testData$Conditions_code)
  levels(testData$Conditions_names) <- make.names(levels(testData$Conditions_names))
  
  # folds<-createFolds(spikes_forCV$Conditions, k=5)    #1-8-25
  
  ### calculate weights
  Y= as.matrix(as.numeric(as.character(spikes_forCV$Conditions_code)))   # I think this works -- need to take the "name", not level
  # xtabs(~Conditions_code,data=spikes_forCV)
  
  fraction_0 <- rep(1 - sum(Y == 0) / nrow(Y), sum(Y == 0))
  fraction_1 <- rep(1 - sum(Y == 1) / nrow(Y), sum(Y == 1))
  # assign that value to a "weights" vector
  weights <- numeric(nrow(Y))
  if (weighted == TRUE) {
    weights[Y == 0] <- fraction_0
    weights[Y == 1] <- fraction_1
  } else {
    weights <- rep(1, nrow(Y))
  }
  
  glmnetGrid <- expand.grid(alpha=c(0, .5, 1), lambda=c(.1, 1, 10))
  
  if (metric == "accuracy" | metric =="") {
    train_control <- trainControl( method="cv",
                                   number=cv_number,
                                   classProbs = TRUE,
                                   #summaryFunction = f1,
                                   savePredictions = TRUE)
    
    if (preprocess) {
      thisModel<- train(eval(modelFormula),
                        data=spikes_forCV,
                        trControl=train_control,
                        method="glmnet",
                        family="binomial",
                        tuneGrid=glmnetGrid,
                        metric="Accuracy",
                        weights=weights,
                        preProcess=c("center","scale")
      )
    } else {   # no preprocess
      thisModel<- train(eval(modelFormula),
                        data=spikes_forCV,
                        trControl=train_control,
                        method="glmnet",
                        family="binomial",
                        tuneGrid=glmnetGrid,
                        metric="Accuracy",
                        weights=weights
      )
    }
  } else if (metric == "F1") {
    train_control <- trainControl( method="cv",
                                   number=cv_number,
                                   classProbs = TRUE,
                                   summaryFunction = f1,
                                   savePredictions = TRUE)
    
    glmnetGrid <- expand.grid(alpha=c(0, .5, 1), lambda=c(.1, 1, 10))
    
    if (preprocess) {
      thisModel<- train(eval(modelFormula),
                        data=spikes_forCV,
                        trControl=train_control,
                        method="glmnet",
                        family="binomial",
                        tuneGrid=glmnetGrid,
                        metric="F1",
                        weights=weights,
                        preProcess=c("center","scale")
      )
    } else {  # no preprocess
      thisModel<- train(eval(modelFormula),
                        data=spikes_forCV,
                        trControl=train_control,
                        method="glmnet",
                        family="binomial",
                        tuneGrid=glmnetGrid,
                        metric="F1",
                        weights=weights
      )
      
    }
    
  } else {
    
    # no metric specified 
    warning("Invalid 'metric' requested.")
  }
  
  # output results to console
  thisModel
  print(modelFormula)
  print(summary(thisModel))
  # effects of elastic net on factors
  print("Parameter tuning")
  print(coef(thisModel$finalModel, thisModel$bestTune$lambda))
  
  spikes_forCV$Class_predicted_CV<-predict(thisModel, newdata=spikes_forCV,type = "raw")
  spikes_forCV$Class_predicted_CV_num<-predict(thisModel, newdata=spikes_forCV,type = "prob")[[2]]   # just select X1 condition
  testData$Class_predicted_CV<-predict(thisModel, newdata=testData,type = "raw")
  testData$Class_predicted_CV_num<-predict(thisModel, newdata=testData,type = "prob")[[2]] 
  
  ## training data
  print("Model performance on training data")
  print(confusionMatrix(reference=spikes_forCV$Conditions_names,
                        data=spikes_forCV$Class_predicted_CV,
                        mode = "prec_recall", positive=positiveClass))  # needs (obs, predict) both as
  
  confmat_train<-confusionMatrix(reference=spikes_forCV$Conditions_names,
                                 data=spikes_forCV$Class_predicted_CV,positive=positiveClass)  # needs (obs, predict) both as factors
  
  print(confmat_train$byClass)
  #print(confmat_train$table)
  
  dat<-data.frame(obs = spikes_forCV$Conditions_names, pred=spikes_forCV$Class_predicted_CV)
  print(prSummary(dat,lev=levels(dat$obs)))
  
  ## test data
  print("Model performance on hold-out test data")
  print(confusionMatrix(reference=testData$Conditions_names,
                        data=testData$Class_predicted_CV,
                        mode = "prec_recall", positive=positiveClass))  # needs (obs, predict) both as
  
  confmat_test<-confusionMatrix(reference=testData$Conditions_names,
                                data=testData$Class_predicted_CV, positive=positiveClass)  # needs (obs, predict) both as
  print(confmat_test$byClass)
  #print(confmat_test$table)
  dat<-data.frame(obs = testData$Conditions_names, pred=testData$Class_predicted_CV)
  print(prSummary(dat,lev=levels(dat$obs)))
  
  return(thisModel)
  #return(c(thisModel,spikes_forCV,testData))
}

fct_saveImagesFromMLEval <- function(res, exptTitle, subtitle, imgSubDir,img_h,img_w)
  #V2 for publication
{
  ifelse(!dir.exists(file.path(imgSubDir)), dir.create(file.path(imgSubDir)), FALSE)
  for (j in 1:4) {
    thisPlot<-res[[j]]
    saveimages(thisPlot,paste0(exptTitle,subtitle," fig"),file.path(imgSubDir),h=img_h,w=img_w)    
  }
}

fct_makeTableFromMLEval <- function(res,filename, exptTitle, subtitle, imgSubDir,item=6)
  #V2 for publication
  {
  #item=6;# corresponding to list of data.frames with results for each model
  # ifelse(!dir.exists(file.path("images", imgSubDir)), dir.create(file.path("images", imgSubDir)), FALSE)
  ifelse(!dir.exists(file.path(imgSubDir)), dir.create(file.path(imgSubDir)), FALSE)
  
  t<-res[[item]]  
  col_names<-names(t)
  t_l =length(t)
  
  theseDF<-bind_cols(t[1:t_l])
  indx <- seq(from=1, to=length(colnames(theseDF))-1, by= 2)
  colnames(theseDF)[indx]<-col_names
  
  print(xtable(theseDF, type="latex"),file=paste0(file.path(imgSubDir),"/","table_",exptTitle,subtitle," ",filename,".tex"))
  print(xtable(theseDF),type="html",file=paste0(file.path(imgSubDir),"/","table_",exptTitle,subtitle," ",filename,".html"))
  write.xlsx(as.data.frame(theseDF), file=paste0(file.path(imgSubDir),"/","table_",exptTitle,subtitle," ",filename,".xlsx"),
             col.names=TRUE, row.names=TRUE, append=FALSE,showNA=FALSE)
  
  # ugly output
  #library(stargazer)
  #print(stargazer(theseDF,type="html"),file=paste0(file.path("images", imgSubDir),"/","table.html"))
  
  # library(knitr)
  #  kable(theseDF)
  
}

#deprecated -- did not handle fish with no events properly.
fct_perFishAvg <- function(spikes5,cafs3,thisModel){
  # take model, spikes, cafs data
  
  #thisModel <- model
  spikes_f1 <- droplevels(filter(spikes5, 
                                 totalCentroidSize_mode_mm2>0.05,
                                 #Genotype2=="wt" | Genotype2=="hom" | Genotype2=="het",
                                 #MaxIntensity_F_F0_centroid >0.1,
                                 (Alive_after == "y" | Alive_after == "yes")))
  
  # apply model to classify spikes
  spikes_f1$classPredicted3 <- predict(thisModel, newdata=spikes_f1,type = "raw")
  
  t1<-round(prop.table(ftable(xtabs(~ Genotype + Conditions + expt_date+ classPredicted3, data=spikes_f1)),1),3)*100
  print(t1)
  
  t2<-round(prop.table(ftable(xtabs(~ Genotype_combined + Conditions + expt_date+ classPredicted3, data=spikes_f1)),1),3)*100
  
  print(t2[2])
  
  theseSpikes <- spikes_f1;
  
  
  require(dplyr)
  
  spikes_summary_byFish <- theseSpikes %>% group_by(animalID,classPredicted3) %>% dplyr::summarise(
    count= n(),
    rate_perMin = n()/30, #30min in duration
    totalTime_min = sum(duration_sec)/60,
    AUC_F_centroid_rate_perMin = sum(AUC_F_centroid)/30,  # 30min in duration
    mean_MaxChange_Intensity_F_F0_centroid = mean(MaxIntensity_F_F0_centroid-MinIntensity_F_F0_centroid),
    MaxIntensity_F_F0_centroid_avg = mean(MaxIntensity_F_F0_centroid),
    distance_xy_mm_avg = mean(distance_xy_mm),
    maxVeloc_xy_mm_sec_avg = mean(maxVeloc_xy_mm_sec),
    meanVeloc_xy_mm_sec_avg = mean(meanVeloc_xy_mm_sec),
    duration_sec_avg = mean(duration_sec),
    timeSpentMoving_xy_sec_prct_avg = mean(timeSpentMoving_xy_sec_prct),
    totalRevolutions_perEvent_avg = mean(totalRevolutions_perEvent),
    SNR_F_F0_centroid_perEvent_avg = mean(SNR_F_F0_centroid),
    MaxIntensity_F_centroid_avg = mean(MaxIntensity_F_centroid),
    # mean_centroidSize_prctChg = mean(baseline_centroidSize_change_prct)
  )
  
  
  spikes_summary_byFish_all <- theseSpikes %>% group_by(animalID) %>% dplyr::summarise(
    count= n(),
    rate_perMin = n()/30, #30min in duration
    totalTime_min = sum(duration_sec)/60,
    AUC_F_centroid_rate_perMin = sum(AUC_F_centroid)/30,  # 30min in duration
    mean_MaxChange_Intensity_F_F0_centroid = mean(MaxIntensity_F_F0_centroid-MinIntensity_F_F0_centroid),
    MaxIntensity_F_F0_centroid_avg = mean(MaxIntensity_F_F0_centroid),
    distance_xy_mm_avg = mean(distance_xy_mm),
    maxVeloc_xy_mm_sec_avg = mean(maxVeloc_xy_mm_sec),
    meanVeloc_xy_mm_sec_avg = mean(meanVeloc_xy_mm_sec),
    duration_sec_avg = mean(duration_sec),
    timeSpentMoving_xy_sec_prct_avg = mean(timeSpentMoving_xy_sec_prct),
    totalRevolutions_perEvent_avg = mean(totalRevolutions_perEvent),
    SNR_F_F0_centroid_perEvent_avg = mean(SNR_F_F0_centroid),
    MaxIntensity_F_centroid_avg = mean(MaxIntensity_F_centroid),
    # mean_centroidSize_prctChg = mean(baseline_centroidSize_change_prct)
  )
  
  #spikes_summary_byFish_all$classPredicted3 <- 2
  spikes_summary_byFish_all$classPredicted3 <- 'X2'
  
  spikes_summary_byFish2 <- bind_rows(spikes_summary_byFish,spikes_summary_byFish_all)
  spikes_summary_byFish <- spikes_summary_byFish2  
  spikes_summary_byFish$classPredicted3 <- as.factor(spikes_summary_byFish$classPredicted3)
  
  cafs_clean <- cafs3 %>%
    filter(totalCentroidSize_mode_mm2 > 0.05,
           #PTZadded =="y",
           (Alive_after == "y" | Alive_after == "yes"))
  
  
  cafs_remerged_class0 <- strict_left_join(cafs_clean, filter(spikes_summary_byFish, classPredicted3=='X0'), by=c("animalID"))
  levels(spikes_summary_byFish$classPredicted3)
  
  # start and end of the new columns from spikes_summary_byFish
  startCol <-length(cafs_clean)+2
  # endCol <- startCol+length(spikes_summary_byFish)-4   # don't redefine threshold
  endCol <- startCol+length(spikes_summary_byFish)-4   # don't redefine threshold
  # lastCol <- startCol+length(spikes_summary_byFish)-3   # don't redefine threshold
  
  thisCafs <- cafs_remerged_class0[,startCol:endCol]
  thisCafs[is.na(thisCafs)]<-0
  #thisCafs$thresh_maxIntensity_F_F0_ma_centroid<-thisThreshold
  
  # print(endCol-startCol+1)
  # print(colnames(thisCafs))
  
  cafs_remerged_class0[,startCol:(endCol)]<- thisCafs
  cafs_remerged_class0$classPredicted3 <- 'X0'   # will relabel fish that never had events to X0
  
  #fish_toExclude2 <- filter(cafs_remerged_class0, rate_perMin == 0, PTZ==40);
  #cafs_remerged_class0 <- filter(cafs_remerged_class0, 
  #                          !is.element(animalID_date,fish_toExclude2$animalID_date))
  
  ### events from X1 #############################
  cafs_remerged_class1 <- strict_left_join(cafs_clean, filter(spikes_summary_byFish, classPredicted3=='X1'), by=c("animalID"))
  
  thisCafs <- cafs_remerged_class1[,startCol:endCol]
  thisCafs[is.na(thisCafs)]<-0
  
  cafs_remerged_class1[,startCol:(endCol)]<- thisCafs
  cafs_remerged_class1$classPredicted3 <- 'X1'   # will relabel fish that never had events to X0
  
  ### all events w/o classification  #####################
  cafs_remerged_class2 <- strict_left_join(cafs_clean, filter(spikes_summary_byFish, classPredicted3=='X2'), by=c("animalID"))
  
  thisCafs <- cafs_remerged_class2[,startCol:endCol]
  thisCafs[is.na(thisCafs)]<-0
  thisCafs$thresh_maxIntensity_F_F0_ma_centroid<-thisThreshold
  
  cafs_remerged_class2[,startCol:(endCol+1)]<- thisCafs
  cafs_remerged_class2$classPredicted3 <- 'X2'   # will relabel fish that never had events to X2
  
  cafs_merged_all <- droplevels(bind_rows(cafs_remerged_class0,cafs_remerged_class1,cafs_remerged_class2))
  
  # need to fix rate_perMin 
  
  # saveRDS(cafs_merged_all,file=paste0(dataDirec,"cafs_merged_all.RData"))
  
  theseGroups <- c("X0","X1","X2")
  theseGroups_labels <- c("classifiedX0","classifiedX1","unclassified")
  
  cafs_forXLS <- cafs_merged_all
  # cafs_forXLS$Genotype2 <- cafs_forXLS$Genotype
  # cafs_forXLS$Genotype_combined <- cafs_forXLS$Genotype
  return(cafs_forXLS)
  
}

# properly handle fish that have no events
fct_perFishAvg_NA <- function(spikes5,cafs3,thisModel){
  # take model, spikes, cafs data
  
  # spikes5<-spikes5_TGB_CT
  # cafs3<-cafs3_TGB_CT
  # thisModel <- model.glmnet.mf
  
  spikes_f1 <- droplevels(filter(spikes5, 
                                 totalCentroidSize_mode_mm2>0.05,
                                 #Genotype2=="wt" | Genotype2=="hom" | Genotype2=="het",
                                 #MaxIntensity_F_F0_centroid >0.1,
                                 (Alive_after == "y" | Alive_after == "yes")))
  
  # apply model to classify spikes
  spikes_f1$classPredicted3 <- predict(thisModel, newdata=spikes_f1,type = "raw")
  
  t1<-round(prop.table(ftable(xtabs(~ Genotype + Conditions + expt_date+ classPredicted3, data=spikes_f1)),1),3)*100
  print(t1)
  
  t2<-round(prop.table(ftable(xtabs(~ Genotype_combined + Conditions + expt_date+ classPredicted3, data=spikes_f1)),1),3)*100
  
  print(t2[2])
  
  theseSpikes <- spikes_f1;
  
  require(dplyr)
  
  spikes_summary_byFish <- theseSpikes %>% group_by(animalID,classPredicted3) %>% dplyr::summarise(
    count= n(),
    rate_perMin = n()/30, #30min in duration
    totalTime_min = sum(duration_sec)/60,
    AUC_F_centroid_rate_perMin = sum(AUC_F_centroid)/30,  # 30min in duration
    mean_MaxChange_Intensity_F_F0_centroid = mean(MaxIntensity_F_F0_centroid-MinIntensity_F_F0_centroid),
    MaxIntensity_F_F0_centroid_avg = mean(MaxIntensity_F_F0_centroid),
    distance_xy_mm_avg = mean(distance_xy_mm),
    maxVeloc_xy_mm_sec_avg = mean(maxVeloc_xy_mm_sec),
    meanVeloc_xy_mm_sec_avg = mean(meanVeloc_xy_mm_sec),
    duration_sec_avg = mean(duration_sec),
    timeSpentMoving_xy_sec_prct_avg = mean(timeSpentMoving_xy_sec_prct),
    totalRevolutions_perEvent_avg = mean(totalRevolutions_perEvent),
    SNR_F_F0_centroid_perEvent_avg = mean(SNR_F_F0_centroid),
    MaxIntensity_F_centroid_avg = mean(MaxIntensity_F_centroid),
    MinIntensity_F_centroid_avg = mean(MinIntensity_F_centroid),
    totalDistance_events = sum(distance_xy_mm)
    # mean_centroidSize_prctChg = mean(baseline_centroidSize_change_prct)
  )
  
  spikes_summary_byFish_all <- theseSpikes %>% group_by(animalID) %>% dplyr::summarise(
    count= n(),
    rate_perMin = n()/30, #30min in duration
    totalTime_min = sum(duration_sec)/60,
    AUC_F_centroid_rate_perMin = sum(AUC_F_centroid)/30,  # 30min in duration
    mean_MaxChange_Intensity_F_F0_centroid = mean(MaxIntensity_F_F0_centroid-MinIntensity_F_F0_centroid),
    MaxIntensity_F_F0_centroid_avg = mean(MaxIntensity_F_F0_centroid),
    distance_xy_mm_avg = mean(distance_xy_mm),
    maxVeloc_xy_mm_sec_avg = mean(maxVeloc_xy_mm_sec),
    meanVeloc_xy_mm_sec_avg = mean(meanVeloc_xy_mm_sec),
    duration_sec_avg = mean(duration_sec),
    timeSpentMoving_xy_sec_prct_avg = mean(timeSpentMoving_xy_sec_prct),
    totalRevolutions_perEvent_avg = mean(totalRevolutions_perEvent),
    SNR_F_F0_centroid_perEvent_avg = mean(SNR_F_F0_centroid),
    MaxIntensity_F_centroid_avg = mean(MaxIntensity_F_centroid),
    MinIntensity_F_centroid_avg = mean(MinIntensity_F_centroid),
    totalDistance_events = sum(distance_xy_mm)
    # rate_perMin_n = rate_perMin*30/duration
    # mean_centroidSize_prctChg = mean(baseline_centroidSize_change_prct)
  )
  
  #spikes_summary_byFish_all$classPredicted3 <- 2
  spikes_summary_byFish_all$classPredicted3 <- 'X2'
  
  spikes_summary_byFish2 <- bind_rows(spikes_summary_byFish,spikes_summary_byFish_all)
  spikes_summary_byFish <- spikes_summary_byFish2  
  spikes_summary_byFish$classPredicted3 <- as.factor(spikes_summary_byFish$classPredicted3)
  
  cafs_clean <- cafs3 %>%
    filter(totalCentroidSize_mode_mm2 > 0.05,
           #PTZadded =="y",
           (Alive_after == "y" | Alive_after == "yes")) %>%
    droplevels()%>%ungroup()
  
  ##### for X0 
  
  # spikes_summary_X0 <- filter(spikes_summary_byFish, classPredicted3=='X0')    # the summary does not have NAs
  
  cafs_remerged_class0 <- strict_left_join(cafs_clean, filter(spikes_summary_byFish, classPredicted3=='X0'), by=c("animalID"))
  levels(spikes_summary_byFish$classPredicted3)
  
  # start and end of the new columns from spikes_summary_byFish
  startCol <-length(cafs_clean)+2
  # endCol <- startCol+length(spikes_summary_byFish)-4   # don't redefine threshold
  endCol <- startCol+length(spikes_summary_byFish)-4   # don't redefine threshold
  # lastCol <- startCol+length(spikes_summary_byFish)-3   # don't redefine threshold
  
  thisCafs <- cafs_remerged_class0[,startCol:endCol]
  # thisCafs[is.na(thisCafs)]<-0      ### 3-23-24
  thisCafs$count[is.na(thisCafs$count)]<-0
  thisCafs$rate_perMin[is.na(thisCafs$rate_perMin)]<-0
  thisCafs$totalTime_min[is.na(thisCafs$totalTime_min)]<-0

  cafs_remerged_class0[,startCol:(endCol)]<- thisCafs
  cafs_remerged_class0$classPredicted3 <- 'X0'   # will relabel fish that never had events to X0
  
  #fish_toExclude2 <- filter(cafs_remerged_class0, rate_perMin == 0, PTZ==40);
  #cafs_remerged_class0 <- filter(cafs_remerged_class0, 
  #                          !is.element(animalID_date,fish_toExclude2$animalID_date))
  
  ### events from X1 #############################
  cafs_remerged_class1 <- strict_left_join(cafs_clean, filter(spikes_summary_byFish, classPredicted3=='X1'), by=c("animalID"))
  
  thisCafs <- cafs_remerged_class1[,startCol:endCol]
  thisCafs$count[is.na(thisCafs$count)]<-0
  thisCafs$rate_perMin[is.na(thisCafs$rate_perMin)]<-0
  thisCafs$totalTime_min[is.na(thisCafs$totalTime_min)]<-0
  
  cafs_remerged_class1[,startCol:(endCol)]<- thisCafs
  cafs_remerged_class1$classPredicted3 <- 'X1'   # will relabel fish that never had events to X0
  
  ### all events w/o classification  #####################
  cafs_remerged_class2 <- strict_left_join(cafs_clean, filter(spikes_summary_byFish, classPredicted3=='X2'), by=c("animalID"))
  
  thisCafs <- cafs_remerged_class2[,startCol:endCol]
  thisCafs$count[is.na(thisCafs$count)]<-0
  thisCafs$rate_perMin[is.na(thisCafs$rate_perMin)]<-0
  thisCafs$totalTime_min[is.na(thisCafs$totalTime_min)]<-0
  
  cafs_remerged_class2[,startCol:(endCol)]<- thisCafs
  cafs_remerged_class2$classPredicted3 <- 'X2'   # will relabel fish that never had events to X2
  
  cafs_merged_all <- droplevels(bind_rows(cafs_remerged_class0,cafs_remerged_class1,cafs_remerged_class2))
  
  # need to fix rate_perMin 
  
  # saveRDS(cafs_merged_all,file=paste0(dataDirec,"cafs_merged_all.RData"))
  
  theseGroups <- c("X0","X1","X2")
  theseGroups_labels <- c("classifiedX0","classifiedX1","unclassified")
  
  cafs_forXLS <- cafs_merged_all
  # cafs_forXLS$Genotype2 <- cafs_forXLS$Genotype
  # cafs_forXLS$Genotype_combined <- cafs_forXLS$Genotype
  return(cafs_forXLS)
  
}

fct_outputFigures<-function(theseCafs,theseCafs_labels,theseEventtype,thisFilter, 
                            thisX, thisGroup,thisFacet, thisFile, thisDirec,figTitle0, 
                            theseMeasures1,theseMeasures2,thisStatRef,
                            img_h,img_w,
                            keepXAxis="n",rotateXLabels="n",groupLines="y") {
  require(xlsx)
  # start for loop
  for (k in 1:length(theseCafs)) {
    
    thisData <-theseCafs[[k]]
    
    for (j in 1:length(theseEventtype)){
      
      theseEventtype[j]
      thisFilter_events <- quos(classPredicted3==theseEventtype[j])
      
      figTitle=paste0(figTitle0, theseEventtype[j],"\n",
                      theseCafs_labels[k])
      
      thisData_f <- thisData %>% filter(!!!thisFilter_events) %>% droplevels() %>% ungroup() %>%suppressMessages()
      thisData_f <- thisData_f %>% filter(!!!thisFilter) %>% droplevels() %>% ungroup() %>% suppressMessages()
      
      p3_scn1lab1 <- suppressMessages(fct_plotGroupAvgs4(thisX,theseMeasures1,thisFilter,thisGroup,thisFacet,
                                        thisData_f,groupLines=groupLines,
                                        keepXAxis = keepXAxis,rotateXLabels = rotateXLabels, stats="y",
                                        stats_ref_group =thisStatRef,generateSummaryData="y"))
      
      p3_scn1lab2 <- suppressMessages(fct_plotGroupAvgs4(thisX,theseMeasures2,thisFilter,thisGroup,thisFacet,
                                        thisData_f,groupLines=groupLines,
                                        keepXAxis = keepXAxis,rotateXLabels = rotateXLabels, stats="y",
                                        stats_ref_group = thisStatRef,generateSummaryData="y"))
      
      p3_scn1lab_plots <- p3_scn1lab1[[1]]
      p3_scn1lab_data <- p3_scn1lab1[[2]]
      
      p4 <- p3_scn1lab_plots + plot_annotation(figTitle,paste(as.character(thisFilter),collapse=" "))
      
      suppressMessages(saveimages(p4,thisFile,thisDirec,h=img_h,w=img_w))
      suppressMessages(savedata(p3_scn1lab_data, thisFile, thisDirec, h=6, w=4, figTitle=figTitle))
      
      p3_scn1lab_plots <- p3_scn1lab2[[1]]
      p3_scn1lab_data <- p3_scn1lab2[[2]]
      
      p4 <- p3_scn1lab_plots + plot_annotation(figTitle,paste(as.character(thisFilter),collapse=" "))
      
      suppressMessages(saveimages(p4,thisFile,thisDirec,h=img_h,w=img_w))
      suppressMessages(savedata(p3_scn1lab_data, thisFile, thisDirec, h=6, w=4, figTitle=figTitle))
      
    }
  }
}

ggsummarytable_cm<- function (data, x, y, digits = 0, size = 3, color = "black", palette = NULL, facet.by = NULL, labeller = "label_value", 
                              position = "identity", ggtheme = theme_pubr(), ...) 
{
  require(ggpubr)
  require(tidyverse)
  
  if (missing(ggtheme) & !is.null(facet.by)) {
    ggtheme <- theme_pubr(border = TRUE)
  }
  if (is.null(names(y))) 
    names(y) <- y
  df <- as.data.frame(data)
  
  # df$x <- df[[x]]
  df$x <- df[[x]]
  # df$x <- df[[paste0(",!!!x,")]]     # changing
  # print(df$x)
  if (color %in% colnames(df)) {
    if (missing(position)) 
      position <- position_dodge(0.8)
    group <- color
  }
  else {
    group <- 1
  }
  # df <- df %>% mutate_if(is.double, round, digits) %>% unite(col = "label", 
  # !!!syms(y), sep = "\n") %>% mutate(y = paste(names(y), 
  # collapse = "\n"))
  df <- df %>% mutate_if(is.double, round, digits) %>% 
    unite(col = "label", !!!(y), sep = "\n") %>% 
    mutate(y = paste(as.character(y),collapse = "\n"))
  
  # print(df$label)
  # print(df$y)
  
  # p <- ggplot(data, aes(x, y)) + geom_exec(geom_text, data = df, 
  #   label = "label", size = size, color = color, group = group, 
  #   position = position)
  # p <- ggpar(p, ggtheme = ggtheme, palette = palette, xlab = x, 
  #   ...)
  
  p <- ggplot(data, aes(x, y)) + geom_exec(geom_text, data = df,
                                           label = "label", size = size, color = color, group = group,
                                           position = position)
  
  p <- ggpar(p, ggtheme = ggtheme, palette = palette, xlab = x,
             ...)
  
  if (!is.null(facet.by)) 
    p <- facet(p, facet.by = facet.by, labeller = labeller, 
               ...)
  p + rremove("ylab")
}

# allowing more complex faceting
plot_PTZresponse_qq2 <- function(theseCafs,x1, y1, thisFilter, thisGroup, thisFacet,title="",groupLines="y")
{ x1<-enquo(x1)
x1_n<-quo_name(x1)
y1<-enquo(y1)
y1_n<-quo_name(y1)

#thisFacet<-enquo(thisFacet)
thisGroup<-enquo(thisGroup)

if (groupLines=="y") {
  theseCafs %>%
    filter(!!!thisFilter) %>%
    ggplot(aes(x= !!x1, y=	!!y1)) +
    geom_boxplot(aes(color=!!thisGroup),na.rm = TRUE)+
    #geom_line(aes(color=Genotype2, group=Genotype2),span=10)+
    # geom_smooth(aes(color=!!thisGroup, group=!!thisGroup),span=10,method='loess', se=FALSE)+
    stat_summary(fun.y="mean", geom="line", linewidth=1.5,aes(group=!!thisGroup, color=!!thisGroup),na.rm = TRUE)+
    #facet_wrap(.~Genotype_combined,labeller="label_both")+
    facet_wrap(vars(!!!thisFacet),dir="h",nrow=1,labeller="label_both")+
    stat_summary(fun.data = give.n, geom = "text", position = position_stack(),na.rm = TRUE)+
    ggtitle(paste0(exptTitle, {{title}})) %>%
    suppressMessages()
} else if (groupLines=="n"){
  theseCafs %>%
    filter(!!!thisFilter) %>%
    ggplot(aes(x= !!x1, y=	!!y1)) +
    geom_boxplot(aes(color=!!thisGroup),na.rm = TRUE)+
    #geom_line(aes(color=Genotype2, group=Genotype2),span=10)+
    # geom_smooth(aes(color=!!thisGroup, group=!!thisGroup),span=10,method='loess', se=FALSE)+
    # stat_summary(fun.y="mean", geom="line", linewidth=1.5,aes(group=!!thisGroup, color=!!thisGroup))+
    #facet_wrap(.~Genotype_combined,labeller="label_both")+
    facet_wrap(vars(!!!thisFacet),dir="h",nrow=1,labeller="label_both")+
    stat_summary(fun.data = give.n, geom = "text", position = position_stack(),na.rm = TRUE)+
    ggtitle(paste0(exptTitle, {{title}}))%>%
    suppressMessages()
} else {
  print("error")
  } 
# return(theseCafs)
}

# updated to allow stats, and output summary data 
fct_plotGroupAvgs4 <- function(thisX, theseMeasures, thisFilter, thisGroup, thisFacet, thisCafs, groupLines="y", keepXAxis="n",rotateXLabels="n",
                               stats="n", stats_ref_group =NaN, stats_method = NaN, stats_paired=NaN, 
                               generateSummaryData = "n") {
  
  thesePlots = results = vector("list", length(theseMeasures));
  
  require(patchwork)
  require(ggpubr)
  
  for (j in 1:length(theseMeasures)){
    
    thisMeasure_q <- theseMeasures[[j]]
    
    p<-plot_PTZresponse_qq2(thisCafs, !!thisX, !!thisMeasure_q, thisFilter, !!thisGroup, thisFacet,'test',groupLines=groupLines)
    
    if (j==length(theseMeasures)){
      ## for the last plot, resurrect the legend
      
      if (keepXAxis=="n") {
        p  <- p+theme_hideFacetLegendXLabel(hideLegend = "n")
      }
      else {
        p <- p+theme_hideFacetLegend(hideLegend = "n")
      }
    }
    else {
      ## not the last plot, keep the legend hidden
      
      if (keepXAxis=="n") {
        p  <- p+theme_hideFacetLegendXLabel()
      } else {
        p <- p+theme_hideFacetLegend()
      }
      
    }
    
    if (stats=="y")  {
      
      if (is.nan(stats_ref_group)){
        stats_ref_group=c("0") # default = 0, probably will crash eventually
      }
      
      if (is.nan(stats_method)) {
        stats_method = "wilcox.test"
      }
      
      if (is.nan(stats_paired)) {
        stats_paired = FALSE;
      }
      
      p <- p + stat_compare_means(label="p.signif",
                                  ref.group=stats_ref_group,  
                                  method=stats_method, 
                                  paired = stats_paired, na.rm = TRUE)
    }
    
    if (rotateXLabels=="y") {
      p <- p + theme(axis.text.x=element_text(angle=90, hjust=1))
    }
    
    thesePlots[[j]] <-p
    
    if (generateSummaryData=="y"){
      
      thisData<-thisCafs
      
      summary.stats <- thisData %>%
        filter(!!!thisFilter) %>% droplevels() %>%
        group_by(!!!thisX,!!!thisGroup,!!!thisFacet) %>%
        select(!!!thisMeasure_q) %>%
        get_summary_stats()
      
      # return grouped data (can be used for Prism)
      summary.data <- thisData %>%
        filter(!!!thisFilter) %>% droplevels() %>%
        group_by(!!!thisX,!!!thisGroup,!!!thisFacet) %>%
        select(!!!thisMeasure_q, animalID_date)   # hard-coded animalID_date -- may break
      
      thisData_filtered <- thisData %>% filter(!!!thisFilter) %>% droplevels() %>%
        group_by(!!!thisX,!!!thisGroup,!!!thisFacet) #
      
      # return stats from grouped data
      # summary.tests <- compare_means(
      summary.tests <- compare_means(
        # as.formula(paste0("'", as_label(thisMeasure_q), "~", as_label(thisGroup), "'")),
        # as.formula(paste0("'", as_label(thisMeasure_q), "~", as_label(thisGroup), "+", as_label(thisFacet), "'")),
        as.formula(paste0("'", as_label(thisMeasure_q), "~", as_label(thisX), "'")),
        data=thisData_filtered,
        # group.by=c("Genotype2"),
        method=stats_method, 
        paired = stats_paired,
        p.adjust.method="fdr")
      
      theseVars <- quos(n, median,mad,iqr,mean,se)
      
      summary.plot <- ggsummarytable_cm(
        summary.stats, x = thisGroup, 
        # y = c("n","median"),
        y=theseVars,
        digits=3, 
        ggtheme = theme_bw()
      ) +
        clean_table_theme()    
      
      p4_stat <- p 
      
      p_summary<-ggarrange(
        p4_stat, 
        summary.plot,
        ncol = 1, align = "v",
        heights = c(0.80, 0.20)
      )
      
      # summary files 
      #  OUTPUT
      #   summary.data,  the full data set, grouped by facet/group etc.
      #   summary.stats, descriptive statistics for the data set, grouped by facet/group etc.
      #   summary.tests, the result of Wilcoxon test, w FDR adjustment
      #   p_summary,     the desired plot with ggsummary table added
      results[[j]]<-list(p_summary, summary.data,summary.stats,summary.tests)
    
  }
  
    }
  
  thisFinalPlot<-wrap_plots(thesePlots, ncol=length(thesePlots))    # this is the object returned 
  
  list(thisFinalPlot, results)

}



### fct_bootstrapSSMD2
# Performs bootstrap resampling of rate_perMin data from cafs data to determine ecdf of RSSMD values for treatment vs vehicle -- relative to an off-plate control/vehicle/background set.
#INPUT
#  thisData, filtered to be two groups of interest
#  thisX, an expr indicating column from thisData as measure of interest
#  target, an expr 
#  background, an expr
# colors, a list of colors 
# bootSampleSize
# replacement
# reps
#OUTPUT
# 
fct_bootstrapSSMD2 <- function (thisData, thisX, thisGroup, target, bkgd, theseColors=c("black","red"), bootSampleSize, replacement=TRUE, reps=1000,thisMetric="robust_ssmd",outputGraph=FALSE,thisTitle="") {
  # target_df <- droplevels(filter(thisData, classPredicted3 == "X0", PTZ==15, Treatment =="03_VPA_long"))
  
  require(moderndive)
  
  target_df <- thisData %>% filter(!!!target) %>% droplevels()
  background_df <- thisData %>% filter(!!!bkgd) %>% droplevels()
  # bootSampleSize = 6
  # replacement=TRUE
  # reps=1000
  set.seed(1234) 
  
  target_rs <- target_df %>%
    select(!!thisX,!!thisGroup) %>%
    rep_sample_n(size=bootSampleSize,replace=replacement,reps=reps)
  
  set.seed(2345) 
  
  background_rs <- background_df  %>%
    select(!!thisX, !!thisGroup) %>%
    rep_sample_n(size=bootSampleSize,replace=replacement,reps=reps)
  
  # # modify to allow background to remain the same group, no bootstrapping
  # background_rs <- background_df  %>% ungroup() %>%
  #   select(!!thisX, !!thisGroup) 
  # 
  # background_size <-nrow(background_rs)
  # 
  # background_rs<-background_rs %>%
  #   rep_sample_n(size=background_size,replace=replacement,reps=reps)
  
  rs_dt<- bind_rows(target_rs, background_rs)
  
  rssmd <- list()
  for (i in 1:length(unique(rs_dt$replicate))) {
    
    thisDf <- dplyr::filter(rs_dt, replicate==i)
    # rssmd<- c(rssmd, fct_qcmeasure_helper(thisDf,target,bkgd,!!thisX,"robust_ssmd"))
    rssmd<- c(rssmd, fct_qcmeasure_helper(thisDf,target,bkgd,!!thisX,thisMetric))
  }
  rssmd_df <-do.call(rbind.data.frame,rssmd)
  colnames(rssmd_df)<-c('rssmd')
  # plot(ecdf(rssmd_df$rssmd),xlim=c(-4,1),main='ecdf of RSSMD, TGB vs DMSO, 1000 resamples, N=6 per group')
  # abline(v=-1)
  # abline(h=0.9)
  # quantile(rssmd_df$rssmd, c(0.5,0.8,0.9,0.95))
  
  set.seed(3456) 
  
  ### make null distribution
  target_rs_null <- background_df  %>%   # take from null distribution
    select(!!thisX, !!thisGroup) %>%
    rep_sample_n(size=bootSampleSize,replace=replacement,reps=reps)
  
  # pretend it came from the txt group
  targetName<- strsplit(quo_name(target[[1]]),split="\"")[[1]][2]   # wow, what a pain in the ass
  eval(parse(text=paste0("levels(target_rs_null$",quo_name(thisGroup),')<-\"',targetName,"\"")))
  
  rs_dt<- bind_rows(target_rs_null, background_rs)
  
  rssmd_null <- list()
  for (i in 1:length(unique(rs_dt$replicate))) {
    
    thisDf <- dplyr::filter(rs_dt, replicate==i)
    # rssmd_null<- c(rssmd_null, fct_qcmeasure_helper(thisDf,target,bkgd,!!thisX,"robust_ssmd"))
    rssmd_null<- c(rssmd_null, fct_qcmeasure_helper(thisDf,target,bkgd,!!thisX,thisMetric))
  }
  rssmd_null_df <-do.call(rbind.data.frame,rssmd_null)
  colnames(rssmd_null_df)<-c('rssmd')
  
  print(quantile(rssmd_df$rssmd, c(0.05,0.8,0.9,0.95)))
  print(quantile(rssmd_null_df$rssmd, c(0.05,0.8,0.9,0.95)))
  
  r_1 <-ecdf(rssmd_df$rssmd)
  r_2 <- ecdf(rssmd_null_df$rssmd)
  
  rssmd_threshold = quantile(rssmd_null_df$rssmd,0.05)[[1]];
  print(paste0("RSSMD threshold to maintain FPR5% = ", rssmd_threshold))
  TPR_at_rssmd_t = r_1(rssmd_threshold);
  print(paste0("TPR at RSSMDt = ", TPR_at_rssmd_t))
  
  theseNames<-eval(parse(text=paste0("levels(rs_dt$",quo_name(thisGroup),")")))  # wow what a headache
  
  color1<-theseColors[[1]]
  color2<-theseColors[[2]]
  # xlim1 <- quantile(rssmd_df$rssmd,0.005)[[1]];
  xlim2 <- quantile(rssmd_df$rssmd,0.05)[[1]];
  
  xlim1 <- min(quantile(rssmd_df$rssmd,0.005)[[1]],quantile(rssmd_null_df$rssmd,0.005)[[1]]);
  # xlim2 <- max(quantile(rssmd_df$rssmd,0.05)[[1]],quantile(rssmd_null_df$rssmd,0.05)[[1]]);
  
  # dev.new(width=6, height=4,unit="in")
  # p1<-recordPlot()
  if (outputGraph==TRUE) {
    plot(ecdf(rssmd_df$rssmd),ylim=c(0,1),xlim=c(xlim1,1.5),
         main=paste0(thisTitle, ", ", paste0(theseNames,collapse=" ")),
         sub=paste0('Resamples= ', reps, ' N= ', bootSampleSize, ' RSSMD =', round(rssmd_threshold,2), ' --> TPR=', TPR_at_rssmd_t),
         cex.sub=0.55, lwd=5,
         cex.main=1,col=color1,xlab="Simulated RSSMD",ylab="ECDF (0-1)")
    lines(ecdf(rssmd_null_df$rssmd),col=color2,lwd=5)
    abline(v=rssmd_threshold, col = "gray60",lty=3)
    text(x=rssmd_threshold, y=0.1, adj=1,paste0("RSSMD=",round(rssmd_threshold,2)))   # SSMD threshold line
    abline(h=0.05, col = "gray60",lty=2)
    abline(h=TPR_at_rssmd_t, col = "gray80",lty=2)  # TPR line
    text(x=xlim2, y=round(TPR_at_rssmd_t,digits=2), paste0("TPR=",TPR_at_rssmd_t))
    # legend("left", inset=.02, legend= theseNames, title=paste0("Resample N=", bootSampleSize), fill=theseColors)
    legend(x=xlim1,y=0.4, inset=.02, legend= theseNames, title=paste0("Resample N=", bootSampleSize), fill=theseColors)
  }
  
  dataObj = list(
    rssmd_df=rssmd_df,
    rssmd_null_df=rssmd_null_df,
    r_1 = r_1,
    r_2 = r_2,
    theseNames = theseNames,
    # thisPlot=p1,
    thisMetric = thisMetric,
    bootSampleSize=bootSampleSize,
    rssmd_threshold=rssmd_threshold,
    TPR_at_rssmd_t=TPR_at_rssmd_t
  )
  
  return(dataObj)
  
  # ggplot(rssmd_df,aes(rssmd)) + 
  #   stat_ecdf(geom = "point",color="green",data=rssmd_df) +
  #   stat_ecdf(geom = "point",color="red",data=rssmd_null_df) +
  #   theme(legend.position = "right")
  
}

# for generating more complete dataset
groupwise_test2<- function(thisData, thisMeasure, thisGroup) {

  #n,p,n_g,chisq.signif
  #N = length(p)
  
  b<-thisData[which(is.element(names(thisData),thisGroup))]
  cM<-thisData[which(is.element(names(thisData),thisMeasure))]
  
  data<-cbind(b,cM)
  
  n <- unique(b)[[1]]  # names of groups
  N <- length(n)  # total # of groups
  
  value1 =value2=value3 = comp = group1 = group2 =group1n = group2n =  c()
  
  ## Compute critical values.
  for (i in 1:(N-1))
  { for (j in (i+1):N)
  {
    group1 = c(group1,as.character(n[i]))
    
    d1  = data[which(is.element(data[,1],n[i])),]
    group1n = c(group1n, length(d1[,1]))
    
    group2 = c(group2,as.character(n[j]))
    d2 = data[which(is.element(data[,1],n[j])),]
    group2n = c(group2n,length(d2[,1]))
    
    comp = c(comp,paste("g",i,"_g",j,sep=""))
    
    value1 <- c(value1,qc.measure(d1[,2],d2[,2],type="ssmd"))
    value2 <- c(value2,qc.measure(d1[,2],d2[,2],type="robust_ssmd"))

  }
  }
  
  resultTable<-as.data.frame(cbind(value1,value2))
  resultTable$comp <- comp
  resultTable$group1 <- group1
  resultTable$group1n <- group1n
  resultTable$group2 <- group2
  resultTable$group2n <- group2n
  
  # switch g1 for g2 to calculate opposite value
  for (i in 1:(N-1))
  { for (j in (i+1):N)
  {
    group1 = c(group1,as.character(n[j]))
    
    d1  = data[which(is.element(data[,1],n[j])),]
    group1n = c(group1n, length(d1[,1]))
    
    group2 = c(group2,as.character(n[i]))
    d2 = data[which(is.element(data[,1],n[i])),]
    group2n = c(group2n,length(d2[,1]))
    
    comp = c(comp,paste("g",j,"_g",i,sep=""))
    
    value1 <- c(value1,qc.measure(d1[,2],d2[,2],type="ssmd"))
    value2 <- c(value2,qc.measure(d1[,2],d2[,2],type="robust_ssmd"))
  
  }
  }
  
  resultTable2<-as.data.frame(cbind(value1,value2))
  resultTable2$comp <- comp
  resultTable2$group1 <- group1
  resultTable2$group1n <- group1n
  resultTable2$group2 <- group2
  resultTable2$group2n <- group2n
  
  resultTable<-rbind(resultTable,resultTable2)
  names(resultTable)[1]<-"ssmd"
  names(resultTable)[2]<-"robustssmd"

  return(resultTable)
}

function(dt_forfct,target,bkgd,measure,typeoftest)
{
  target_dt <- dt_forfct %>% filter(!!!target) %>% data.frame() %>% select({{measure}}) # have to do this to make dplyr output 1 column
  #target_dt <- dt %>% filter(!!!target) %>% data.frame() %>% select(median_iGABA_mBFP)
  bkgd_dt <- dt_forfct %>% filter(!!!bkgd)%>% data.frame() %>% select({{measure}})
  
  #return(bkgd_dt)
  z<-qc.measure(target_dt[,],bkgd_dt[,], type={{typeoftest}})
  #z<-qc.measure(target_dt,bkgd_dt, type="robust_ssmd")
  
  return(z)
}

savedata <- function(d,thisFile, thisDirec, figTitle="", w=12,h=8,dpi=300){
  
  # take a results object, list of lists, and output to XLS files and images in defined directory 
  # INPUT
  # d, the results object generated from fct_plotGroupAvgs4
  #     structure is list of N-length, 
  #     each element is a list of 4-length corresponding to N measures
  #       each element of the 4-length list is: 
  #         1. a ggplot2 graphic with p-values and descriptive stats added
  #         2. summary.data,   the grouped/faceted data used to generate the plot
  #         3. summary.stats,  descriptive stats from the grouped/facted data
  #         3. summary.tests,  statistical testing for groupwise differences among levels of group
  #
  # thisFile,  the prefix to append to each file  #
  # thisDirec, the images subdirectory to save
  #
  # w, width for individual plot
  # h, height for individual plot
  #
  
  require(stringr)
  require(openxlsx)
  require(patchwork)
  
  ifelse(!dir.exists(file.path(thisDirec)), dir.create(file.path(thisDirec)), FALSE)
  
  thisFile<-stringr::str_replace_all(thisFile,'[ ]','_')
  # filename = paste0(thisDirec,"/",next_file(thisFile,'png',thisDirec))
  filename = paste0(thisDirec,"/",next_file_2(thisFile,'png',thisDirec))
  
  rootName<-paste0(substr(filename, 1, (nchar(filename)-4)),"_")
  
  for (i in 1:length(d)){
    thisMeasure_result <- d[[i]]
    
    # the plot
    p <- thisMeasure_result[[1]] + plot_annotation(figTitle)
    ggsave(plot=p, filename = paste0(rootName,i,".png"), width=w, height=h, dpi=dpi)
    ggsave(plot=p, filename = paste0(rootName,i,".svg"), width=w, height=h, dpi=dpi)
    try(ggsave(plot=p, filename = paste0(rootName,i,".eps"), width=w, height=h, dpi=dpi))
    # save(p, file=paste0(rootName,i,'.RData',thisDirec), )
    saveRDS(p, file=paste0(rootName,i,'.rds'))
    
    # the data    -- started throwing errors 
    # write.xlsx(list("data" = as.data.frame(thisMeasure_result[[2]]),
    #              "stats" = as.data.frame(thisMeasure_result[[3]]),
    #              "tests" = as.data.frame(thisMeasure_result[[4]])),
    #            file=paste0(rootName,i,"_data.xlsx"),
    #            col.names=TRUE, row.names=TRUE, append=TRUE,showNA=TRUE)
    # 
    
    write.xlsx2(list(data = as.data.frame(thisMeasure_result[[2]])), sheetName="data",
                file = paste0(rootName, i, "_data.xlsx"), col.names = TRUE, row.names = TRUE, append = FALSE, showNA = TRUE)
    
    write.xlsx2(list(stats = as.data.frame(thisMeasure_result[[3]])), sheetName="stats",
                file = paste0(rootName, i, "_data.xlsx"), col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE)
    
    write.xlsx2(list(tests = as.data.frame(thisMeasure_result[[4]])), sheetName="tests",
                file = paste0(rootName, i, "_data.xlsx"), col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE)
    
    
  }
  
}

# saveimages 
# v2. saves RDS, take a figTitle, improved logic re sequential numbering
saveimages <- function(p,thisFile2='', thisDirec2='', w=12,h=8,dpi=300,figTitle='')
{
  require(svglite)
  require(stringr)
  require(patchwork)
  thisFile2<-stringr::str_replace_all(thisFile2,'[ ]','_')
  
  if (figTitle!='') {
    p <- p + plot_annotation(figTitle)
  }
  
  filename <- next_file(thisFile2,'png',thisDirec2)
  rootName<-paste0(substr(filename, 1, (nchar(filename)-4)))
  print(rootName)
  
  # ggsave(plot=p, filename = paste0(thisDirec,"/",next_file(thisFile,'png',thisDirec)), width=w, height=h, dpi=dpi)
  # ggsave(plot=p, filename = paste0(thisDirec,"/",next_file(thisFile,'svg',thisDirec)), width=w, height=h, dpi=dpi)
  # #eps sometimes fails if Courier not available
  # try(ggsave(plot=p, filename = paste0(thisDirec,"/",next_file(thisFile,'eps',thisDirec)), width=w, height=h, dpi=dpi))
  # # save as RDS
  # saveRDS(p, file=paste0(thisDirec,"/",next_file(thisFile,'rds',thisDirec)))
  
  ifelse(!dir.exists(file.path(thisDirec2)), dir.create(file.path(thisDirec2)), FALSE)
  
  ggsave(plot=p, filename = paste0(thisDirec2,"/",rootName,'.png'), width=w, height=h, dpi=dpi)
  ggsave(plot=p, filename = paste0(thisDirec2,"/",rootName,'.svg'), width=w, height=h, dpi=dpi)
  #eps sometimes fails if Courier not available
  try(ggsave(plot=p, filename = paste0(thisDirec2,"/",rootName,'.eps'), width=w, height=h, dpi=dpi))
  # save as RDS
  saveRDS(p, file=paste0(thisDirec2,"/",rootName,'.rds'))
  
  }

next_file = function(basename = 'myfile', fileext = 'png', filepath = '.'){
  
  old.fnames = grep(paste0(basename,' \\d+\\.', fileext,'$'), 
                    list.files(filepath), value = T)
  lastnum = gsub(paste0(basename,' (\\d+)\\.', fileext,'$'), '\\1', old.fnames)
  if (!length(lastnum)) { 
    lastnum = 1 
  } else {
    lastnum = sort(as.integer(lastnum),T)[1] + 1L 
  }
  return(paste0(basename, ' ', sprintf('%03i', lastnum), '.', fileext))
}

# next_file_2, derived from next_file -- accounts for specific pattern at end of filename
next_file_2 = function(basename = 'myfile', fileext = 'png', filepath = '.', prefixpattern = ' \\d{3}'){
  
  # basename="FDSS_ML_scn1lab_wstats_fig"
  # fileext='png'
  # filepath = thisDirec
  # prefixpattern = ' \\d{3}'
  
  # Build the regex pattern dynamically based on arguments
  pattern = paste0(basename, prefixpattern, '_\\d+\\.', fileext, '$')
  #extract_pattern = paste0(basename, '_', prefixpattern, '\\.', fileext, '$')
  
  # List files and find those matching the pattern
  old.fnames = grep(pattern, list.files(filepath), value = TRUE)
  
  # Extract the numeric part from the filenames (the part after 'fig ')
  # lastnum = gsub(paste0(basename, ' (\\d{3})_\\d+\\.', fileext, '$'), '\\1', old.fnames)
  lastnum = gsub(paste0(basename, '(', prefixpattern, ')', '_\\d+\\.', fileext, '$'), '\\1', old.fnames)
  
  # Determine the next file number
  if (!length(lastnum)) {
    lastnum = 1 
  } else {
    lastnum = max(as.integer(lastnum)) + 1L 
  }
  
  # Generate the new filename
  # return(paste0(basename, '_fig ', sprintf('%03i', lastnum), '.', fileext))
  return(paste0(basename, " ", sprintf('%03i', lastnum), '.', fileext))
}

check_data_types <- function(...){
  # Get the list of data frames
  dfs <- list(...)
  
  # Create a list to store the unique column names
  all_cols <- unique(unlist(lapply(dfs, names)))
  
  # Initialize a list to store mismatched column information
  type_mismatches <- list()
  
  # Loop through each column
  for(col in all_cols){
    # Get the data types of this column in each data frame where it exists
    types <- sapply(dfs, function(df) if(col %in% names(df)) class(df[[col]]) else NA)
    
    # Check if there is more than one unique type (excluding NA)
    unique_types <- unique(types[!is.na(types)])
    if(length(unique_types) > 1){
      type_mismatches[[col]] <- unique_types
    }
  }
  
  # Return the list of mismatches
  return(type_mismatches)
}
theme_hideFacetLegendXLabel <- function(hideLegend="y",hideXAxis="y") {
  require(grid)
  require(ggthemes)
  
  theme_Publication()  %+replace% 
    
    if (hideLegend=="y") {
      theme(legend.position = "none", plot.title = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x=element_blank(), #remove x axis labels
            axis.ticks.x=element_blank()) #remove x axis ticks
      
    } else {
      theme(legend.position = "right",  legend.direction ="vertical",
            strip.background = element_blank(), 
            plot.title = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x=element_blank(), #remove x axis labels
            axis.ticks.x=element_blank()) #remove x axis ticks
    }
}

theme_hideFacetLegend <- function(hideLegend="y",hideXAxis="y") {
  require(grid)
  require(ggthemes)
  
  theme_Publication()  %+replace% 
    
    if (hideLegend=="y") {
      theme(legend.position = "none", plot.title = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank())
      
    } else {
      theme(legend.position = "right",  legend.direction ="vertical",
            strip.background = element_blank(), 
            plot.title = element_blank(),
            strip.text.x = element_blank())
      
    }
}

strict_left_join <- function(x, y, by = NULL, ...){
  by <- common_by(by, x, y)
  if(any(duplicated(y[by$y]))) {
    stop("Duplicate values in foreign key")
  } else left_join(x, y, by = by, ...)
}

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# function for mean labels
mean.n <- function(x){
  return(c(y = median(x)*0.97, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}

theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

geom_text_repel2 <- function(...) {
  layer <- ggrepel::geom_text_repel(...)
  layer$ggrepel <- TRUE
  class(layer) <- c("ggrepel", class(layer))
  return(layer)
}

ggplot_add.ggrepel <- function(object, plot, object_name) {
  if (any(do.call(c, lapply(plot$layer, function(x) x$ggrepel)))) {
    warning(
      "There is more than one ggrepel layers. ",
      "This may cause overlap of labels"
    )
  }
  # Optionally, one may modify `object` here.
  NextMethod("ggplot_add")
}

fct_qcmeasure_helper<-function(dt_forfct,target,bkgd,measure,typeoftest)
{
  target_dt <- dt_forfct %>% filter(!!!target) %>% data.frame() %>% select({{measure}}) # have to do this to make dplyr output 1 column
  bkgd_dt <- dt_forfct %>% filter(!!!bkgd)%>% data.frame() %>% select({{measure}})
  
  z<-qc.measure(target_dt[,],bkgd_dt[,], type={{typeoftest}})

    return(z)
}
