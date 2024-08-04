#Helper functions for dictionary validation

require(tidyr)
require(AutoAnnotation)
require(dplyr)
#---- MetaMap helpers ----
# Process MetaMap results from downloaded Metamap data 
# MetaMap data are batch processed data from ASCII TiAb data 
# updated (SyRF) and historcial data downloaded from metaMap are saved in folder "method2MetaMap"
# Convert metaMap json file to data frame
ReadInMetaMapData <- function(dataSource, dataFilenames){
  metaTypeDrug <- "antb|crbs|chem|clnd|dora|elii|food|fngs|hops|horm|inch|mosq|nnon|nusq|orch|phsu|phpr|sbst|vita"
  metaTypeDisease <- "dsyn"
  if(dataSource == "syrf"){
    diseaseFile <- dataFilenames$SyRFDiseaseMetaMapDF
    drugFile <- dataFilenames$SyRFDrugMetaMapDF
  } else {
    diseaseFile <- dataFilenames$HistoricalDiseaseMetaMapDF
    drugFile <- dataFilenames$HistoricalDrugMetaMapDF
  }
  if(!file.exists(diseaseFile) | !file.exists(drugFile)){
    if(dataSource == "syrf"){
      method2JsonResults <- fromJSON("data/historicalMetaMapJsonOutput.txt")
    }else{
      method2JsonResults <- fromJSON("data/updatedMetaMapJsonOutput.txt")
    }
    utterances <- method2JsonResults$AllDocuments$Document$Utterances
    for(i in seq_along(utterances)){
      utterance <- utterances[[i]]
      df <- unnest(as.data.frame(utterance))[,c("PMID","PhraseText","Mappings"  # ,   "SyntaxUnits" ,"Candidates" "UttStartPos","UttLength","PhraseStartPos","PhraseLength", "UttText", "UttNum","UttSection", 
      )]
      df$nMappings <- sapply(df$Mappings, length)
      df <- df[df$nMappings > 0, ]
      if(nrow(df) == 0) next
      df <- unnest(unnest(df))[,c("PMID","Sources", "PhraseText" , "CandidateCUI","CandidateMatched" ,"CandidatePreferred","MatchedWords","SemTypes" # , "IsHead" , "UttSection" ,"UttNum","IsOverMatch" ,"ConceptPIs"   ,"Status"    ,"Negated"    ,"MatchMaps","MappingScore"  ,"CandidateScore"
      )]
      
      df$MatchedWords <- sapply(df$MatchedWords, paste, collapse = " " )
      uniqueDf <- df %>%
        group_by(PMID, PhraseText, CandidateCUI, CandidatePreferred) %>%
        summarise(MatchedWords = first(MatchedWords), SemTypes = list(unique(unlist(SemTypes)))) %>%
        group_by(PMID,  CandidateCUI, CandidatePreferred) %>%
        summarise(frequency = n(), SemTypes = paste(unique(unlist(SemTypes)), collapse = "|")
                  ,MatchedWords = paste(unique(unlist(MatchedWords)), collapse = ","))
      
      disease <-  uniqueDf %>%
        filter(grepl(metaTypeDisease, SemTypes)) %>%
        group_by(PMID) %>%
        mutate(maxFrequency = max(frequency), ID = as.character(PMID)) %>%
        filter(frequency == maxFrequency)
      
      drug <-  uniqueDf %>%
        filter(grepl(metaTypeDrug, SemTypes)) %>%
        group_by(PMID) %>%
        mutate(maxFrequency = max(frequency) , ID = as.character(PMID)) %>%
        filter(frequency == max(frequency))
      
      if(i == 1 | !exists("diseases")) {
        diseases <- disease
      } else   {
        diseases <- rbind(diseases, disease)  
      }
      
      if(i == 1 | !exists("drugs")) {
        drugs <- drug
      } else   {
        drugs <- rbind(drugs, drug)  
      }
    }
    
    write.csv(diseases, diseaseFile, row.names = F)
    write.csv(drugs, drugFile, row.names = F)
  } else{
    diseases <- read.csv(diseaseFile, stringsAsFactors = F)
    drugs <- read.csv(drugFile, stringsAsFactors = F)
  }
  return(list(diseasesMetaMap = diseases, drugsMetaMap = drugs))
}
#---- SyRF helpers ----
ReadStudiesFromSyRF <- function(syrfConnection, projectIdStr, recentUpdateDate = "2016-04-01T00:00:00.000Z", onlyIncluded = F){
  CalculateDecisionForIOE <- function(myStudies, myProject){
    # include: 1, exclude: 0, unknown: 99
    return( 
      ifelse(myStudies$ScreeningInfo$AgreementMeasure$NumberScreened < myProject$AgreementThreshold$NumberScreened |
               myStudies$ScreeningInfo$AgreementMeasure$AbsoluteAgreementRatio < myProject$AgreementThreshold$AbsoluteAgreementRatio
             , 99
             , ifelse(myStudies$ScreeningInfo$Inclusion > 0.5, 1, 0))
    )
  }
  
  myProject <- GetProject(syrfConnection$conProject, projectIdStr, raw = F)
  myStudies <- GetStudiesForProject(syrfConnection$conStudy, projectIdStr, raw = T,recentUpdateDate)
  myStudies$Label <- CalculateDecisionForIOE(myStudies, myProject)
  
  if(onlyIncluded) myStudies <- myStudies[myStudies$Label == 1,]
  
  myStudies$CleanTitle <- gsub("[[:punct:]]","", myStudies$Title,perl=T)
  myStudies$CleanAbstract <- gsub("[[:punct:]]","", myStudies$Abstract,perl=T)
  
  return(myStudies)
}
ReadInIncludedStudiesFromSyRF <- function(syrfConnection, projectName, projects){
  return(ReadStudiesFromSyRF(syrfConnection, projectName, projects, onlyIncluded = T))
}
KeptColumnsForInlcudedStudies <- function(includedStudies){
  keptColumns <- c( "Title","Abstract","Authors" ,"PdfRelativePath", "Year","DOI",  "idStr"   ,"SystematicSearchIdStr" ,"Journal" )
  
  validColumns <- intersect(keptColumns, names(includedStudies))
  
  return( includedStudies[, validColumns])
}
LoadFullData <- function(){
  read.csv(file = outputFilenames$CombinedStudiesHumanDiseaseDrugForValidation, stringsAsFactors = F)
}
LoadUpdatedData <- function(){
  allStudies <- LoadFullData()
  allStudies[which(is.na(as.numeric(allStudies$HistoricalID))),]
}
LoadHistoricalData <- function(){
  allStudies <- LoadFullData()
  allStudies[which(!is.na(as.numeric(allStudies$HistoricalID))),]
  # allStudies[which(allStudies$HistoricalID != ""),]
}
LoadMSTypeData <- function(){
  read.csv(file = outputFilenames$SyRFMSTypeAnnotation, stringsAsFactors = F)
}

#---- Annotation Extraction helpers ----
ExtractAnnotation  <- function(myStudies , dictionaryName, ignoreCase = T, textSearchingHeaders = c( "CleanTitle" ,"CleanAbstract"), idColumn ="idStr"){
  results <- CountTermsInStudies(searchingData = as.data.frame(myStudies)
                                 , dictionary = dictionaryName
                                 , dictionaryNameHeader = "Id"
                                 , dictionaryRegexHeader = "Regex"
                                 , textSearchingHeaders = textSearchingHeaders
                                 , ignoreCase = ignoreCase
  )
  resultsWide <- data.frame(sapply(results[, -1,drop=F], function(x) as.numeric(as.character(x))))
  resultsWide[resultsWide==0] <- NA
  resultsWide[,"Id"] <- myStudies[,idColumn]
  dictionary <- GetData(dictionaryName) %>%
    mutate(Id = as.character(Id))
  resultsLong <- resultsWide %>%
    gather(key = "Var", value="Frequency", -"Id", na.rm = T)%>%
    mutate(Var = gsub("X","",Var)) %>%
    left_join(dictionary, by = c("Var"="Id"))
  
  names(resultsLong)[names(resultsLong) =="Id"] <- idColumn
  return(resultsLong)
}
TagFrequentAnnotation  <- function(resultsLong, idColumn ="idStr", minimumIncludeFrequency = 10 , maximumExcludeFrequency = 1){
  names(resultsLong)[names(resultsLong) == idColumn] <- "Id"
  
  resultsLongFiltered <- resultsLong %>%
    group_by(Id) %>%
    mutate(
      Count = n(),
      MaxFrequency = max(Frequency),
      Max = Frequency == MaxFrequency,
      Vars = paste(Var, collapse = ","),
      VarFrequency = paste(Frequency, collapse = ","),
      MaxVar = paste(Var[which(Frequency == MaxFrequency)], collapse = ","),
      MaxVarFrequency = paste(Frequency[which(Frequency == MaxFrequency)], collapse = ","),
      Selected = (Frequency >= minimumIncludeFrequency) | ((MaxFrequency - Frequency) < maximumExcludeFrequency)
    )
  names(resultsLongFiltered)[names(resultsLongFiltered) == "Id"] <- idColumn
  # write.csv(resultsLongFiltered, "3_drugBank/resultsLongFiltered.csv")
  return(resultsLongFiltered)
}

GroupSelectedAnnotation <- function(resultsLongFiltered, idColumn ="idStr", varnames = c("Terms","Frequencies"),groupAnnotation=F){
  names(resultsLongFiltered)[names(resultsLongFiltered) == idColumn] <- "Id"
  
  if(groupAnnotation){
    print("groupAnnotation")
    groupedData <- as.data.frame(resultsLongFiltered) %>%
      filter(Selected) %>%
      group_by(Id) %>%
      summarise(
        Terms = paste(Name, collapse = "||"),
        Frequencies = paste(Frequency, collapse = "||")
      )
  } else{
    groupedData <- as.data.frame(resultsLongFiltered) %>%
      filter(Selected) %>%
      select(Id, Terms = Name, Frequencies = Frequency)
  }
  
  names(groupedData) <- c(idColumn,  varnames)
  
  return(groupedData) 
}

AnnotateStudiesWithDictionary <- function(resultsLong, idColumn = "idStr", varnames = c("Terms",    "Frequencies"), minimumIncludeFrequency = 9, maximumExcludeFrequency = 1, groupAnnotation = F){
  resultsLongFiltered <- TagFrequentAnnotation(resultsLong, idColumn = idColumn, minimumIncludeFrequency = minimumIncludeFrequency , maximumExcludeFrequency = maximumExcludeFrequency)
  groupedData <- GroupSelectedAnnotation(resultsLongFiltered, idColumn = idColumn, varnames = varnames, groupAnnotation=groupAnnotation)
  
  annotatedData <- unique(resultsLongFiltered[,idColumn ,drop=F]) %>%
    full_join(groupedData)
  
  return(annotatedData)
}
ExtractDisease <- function(myStudies, dictionaryName, idColumn = "idStr", textSearchingHeaders = c("Title","Abstract"), minimumIncludeFrequency = 9,  maximumExcludeFrequency = 1, varnames = c("Disease", "DiseaseFrequency"), groupAnnotation=groupAnnotation, ignoreCase=T){
  resultsLong <- ExtractAnnotation(myStudies, dictionaryName, idColumn = idColumn, ignoreCase = ignoreCase, textSearchingHeaders = textSearchingHeaders)
  
  myDiseases <- AnnotateStudiesWithDictionary(resultsLong, idColumn = idColumn, varnames = varnames, minimumIncludeFrequency = 9, maximumExcludeFrequency = 1, groupAnnotation=groupAnnotation)
  return(myDiseases)
}
ExtractDrug <- function(myStudies, dictionaryName, idColumn = "idStr", textSearchingHeaders = c("Title","Abstract"), minimumIncludeFrequency = 9,  maximumExcludeFrequency = 1, varnames = c("Drug", "DrugFrequency"), groupAnnotation=groupAnnotation, ignoreCase=T){
  resultsLong <- ExtractAnnotation(myStudies, dictionaryName, idColumn = idColumn, ignoreCase = ignoreCase, textSearchingHeaders = textSearchingHeaders)
  
  myDrugs <- AnnotateStudiesWithDictionary(resultsLong, idColumn = idColumn, minimumIncludeFrequency = minimumIncludeFrequency, maximumExcludeFrequency = maximumExcludeFrequency , varnames = varnames, groupAnnotation=groupAnnotation)
  return(myDrugs)
}
ExtractDrugDisease <- function(myStudies, diseaseDictionary, drugDictionary, idColumn = "idStr", minimumIncludeFrequency = 9,  maximumExcludeFrequency = 1, groupAnnotation =F ){  
  myDiseases <- ExtractDisease(myStudies, dictionaryName= diseaseDictionary, idColumn = idColumn, minimumIncludeFrequency = minimumIncludeFrequency,  maximumExcludeFrequency = maximumExcludeFrequency,groupAnnotation=groupAnnotation)
  
  myDrugs <- ExtractDrug(myStudies, dictionaryName=drugDictionary, idColumn = idColumn, minimumIncludeFrequency = minimumIncludeFrequency,  maximumExcludeFrequency = maximumExcludeFrequency,groupAnnotation=groupAnnotation)
  
  myStudiesDiseaseDrugAnnotated <- myStudies %>%
    full_join(myDiseases) %>%
    full_join(myDrugs)
  
  return(myStudiesDiseaseDrugAnnotated)
}

GridSearch <- function(myStudies, dictionaryName, textSearchingHeadersHuman, textSearchingHeadersDic, idColumn, minimumIncludeFrequencies , maximumExcludeFrequenciess, outputFolder, na.correct = FALSE, ignoreCase = FALSE, candidates = NULL)
{
  resultsLongHuman <- ExtractAnnotation(myStudies, dictionaryName, idColumn = idColumn, ignoreCase = ignoreCase, textSearchingHeaders = textSearchingHeadersHuman)
  resultsLongDic <- ExtractAnnotation(myStudies, dictionaryName, idColumn = idColumn, ignoreCase = ignoreCase, textSearchingHeaders = textSearchingHeadersDic)
  # run grid search to identify the best parameters 
  performances <- data.frame(matrix(nrow=0,ncol=8))
  for(minimumIncludeFrequency in minimumIncludeFrequencies ){
    for(maximumExcludeFrequency in maximumExcludeFrequencies )  {
      print(paste0("minimumIncludeFrequenc: ", minimumIncludeFrequency))
      print(paste0("maximumExcludeFrequency: ", maximumExcludeFrequency))
      
      humanAnnotation <- AnnotateStudiesWithDictionary(resultsLongHuman, idColumn = idColumn, varnames = c("HumanTerm","HumanFrequency"), minimumIncludeFrequency = 1, maximumExcludeFrequency = 10, groupAnnotation = T)
      dicAnnotation <- AnnotateStudiesWithDictionary(resultsLongDic, idColumn = idColumn, varnames = c("DicTerm","DicFrequency"), minimumIncludeFrequency = minimumIncludeFrequency, maximumExcludeFrequency = maximumExcludeFrequency, groupAnnotation = T)
      
      annotatedFull <- myStudies[, c(idColumn,"Title" ,"Abstract","Diseases", "Drugs")] %>%
        left_join(humanAnnotation) %>%
        left_join(dicAnnotation) 
      annotatedFull$OriginalHumanTerm <- annotatedFull$HumanTerm
      annotatedFull$OriginalDicTerm <-annotatedFull$DicTerm
      
      if(!is.null(candidates)) {
        annotatedFull$HumanTerm <- sapply(annotatedFull$OriginalHumanTerm, function(term){
          termList <- unlist(strsplit(term, "||", fixed = T))
          finalList <- intersect(termList, candidates)
          return(paste(finalList, collapse = "||") )
        })
        annotatedFull$DicTerm <- sapply(annotatedFull$OriginalDicTerm, function(term){
          termList <- unlist(strsplit(term, "||", fixed = T))
          finalList <- intersect(termList, candidates)
          return(paste(finalList, collapse = "||") )
        })
      }
      annotatedFull <- annotatedFull %>%
        mutate(correct1 = HumanTerm == DicTerm,
               correct2 = na.correct & is.na(HumanTerm) & is.na(DicTerm),
               correct = correct1 | correct2)
      
      # filter(!is.na(HumanTerm)) 
      nStudies <- nrow(annotatedFull)
      
      performance <- c(minimumIncludeFrequency= minimumIncludeFrequency, maximumExcludeFrequency = maximumExcludeFrequency, rate = sum(annotatedFull$correct, na.rm = T)/nStudies, filename = paste0(outputFolder, "/",minimumIncludeFrequency, "_",maximumExcludeFrequency, ".csv",sep=""))
      
      write.csv(annotatedFull, file = performance["filename"], row.names = F, na = "")
      print(performance)
      performances <- rbind(performances, performance)
    }
  }
  
  names(performances) <- names(performance)
  print(paste0("best performance is:", max(performances$rate)))
  write.csv(performances, file = paste(outputFolder,"performances.csv",sep="/"))
  results <- read.csv(file=performances$filename[which(performances$rate == max(performances$rate))[1]], stringsAsFactors = F)
  outputResults <- myStudies[,c("HistoricalID", "PublicationID", "OldIdStr")] %>% right_join(results) 
  write.csv(unique(outputResults),file= paste0(outputFolder, "/forValidatin.csv"), row.names = F)
  
  return(outputResults)
}
