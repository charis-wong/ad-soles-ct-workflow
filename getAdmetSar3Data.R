#for drugs with SMILES available, get data from admetSAR3 http://lmmd.ecust.edu.cn/admetsar3/predict.php

library(dplyr)

# 1. write compatible files for admetsar (txt files containing smiles only, in batches of 1000)-----
allDrugs <- read.csv("output/2024-08-04masterDrugList.csv")
drugsWithSmiles <- allDrugs %>% filter(SMILES != "NA")
drugsSMILES<- drugsWithSmiles$SMILES
batches <-split(drugsSMILES, ceiling(seq_along(drugsSMILES)/1000))
for(i in 1:length(batches)){
  write.table(batches[i], 
              paste0("data/admetSAR3/input/SMILES_batch", i, ".txt"), 
              quote = FALSE,
              row.names = FALSE, 
              col.names = FALSE)
  
  
}

# 2. get, upload, and compile predictions from admetSAR3------ 
# Go to http://lmmd.ecust.edu.cn/admetsar3/predict.php, Select 'Predict' tab, Select 'ADMET Screening', enter txt file in text box, click 'submit' then get results and download txt files. 
#upload to juniper under 'data/admetSAR3/predictions/' folder

predFolderPath <- "data/admetSAR3/predictions/"
pred <- list.files(predFolderPath)



for(i in 1:length(pred)){
  
  predData <- read.delim(paste0(predFolderPath, pred[i]))
  if(!exists("predictionOutput")){ predictionOutput <- predData} else {predictionOutput <- rbind(predictionOutput, predData)}
  
}


# 3. write predictionOutputs

write.csv(predictionOutput, "data/admetSAR3Predictions.csv", row.names = FALSE)


# # 4. combine masterDrugList with relevant data (BBB + ro5) from admetSAR3 - COMBINE in generateMasterDrugList.R
# 
# allDrugs <- allDrugs %>% 
#   # select(-admetSAR_BBB, -admetSAR_p_CNSPenetrance, -BBB)%>%
#   # rename(B3DB_BBB = BBB1)%>%
#   unique()
# 
# admetData <- predictionOutput%>%
#   select(SMILES, MW, HBA, HBD, SlogP, admetSAR3_BBB = BBB)%>%unique()%>%
#   mutate(ro5_MW = ifelse(MW<500, 0, 1),
#          ro5_HBA = ifelse(HBA<10, 0, 1), 
#          ro5_HBD = ifelse(HBD<5, 0, 1),
#          ro5_logP = ifelse(SlogP<5, 0, 1))
# 
# admetData <- admetData%>%
#   mutate(across(starts_with('ro5'), as.numeric))%>%
#   mutate(admetSAR_ro5Violations = ro5_MW+ro5_HBA + ro5_HBD+ ro5_logP)
# 
# 
# allDrugs1 <- left_join(allDrugs, admetData, by = "SMILES")
# 
# allDrugs1 <- allDrugs1%>%unique()
# 
# write.csv(allDrugs1, paste0("output/", Sys.Date(), "masterDrugList.csv"), row.names = FALSE)
