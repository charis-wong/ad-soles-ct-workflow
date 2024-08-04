createClinicalProgressSummary <- function(clinicalOutputCrossTable, clinicalProgressSummary, ad_longlist){
  
  
  nDrugMeetLogic = nrow(DrugsMeetLogicTable)
  nPublicationsMeetLogic <- sum(DrugsMeetLogicTable$nPub)
  nCoreDrugs <- nrow(ad_longlist)
  CoreDrugPubs <- clinicalStudyList%>%filter(Drug %in% ad_longlist$Drug) 
  nCoreDrugPublications <- length(unique(CoreDrugPubs$idStr))
  nSingleAnnotated <- CoreDrugPubs%>%filter(nReviews == 1)%>%select(idStr)%>%unique()%>%nrow()
  nDualAnnotated <- CoreDrugPubs%>%filter(nReviews == 2)%>%select(idStr)%>%unique()%>%nrow()
  nReconciled <- CoreDrugPubs %>% filter(reconcile == TRUE)%>%select(idStr)%>%unique()%>%nrow()
  nUniquePublications = clinicalProgressSummary$nUniquePublications
  nIncludedPublications = clinicalProgressSummary$nIncludedPublications
  
  clinicalProgressSummary <- data.frame(
    ReviewType = 'Clinical', 
    nUniquePublications = nUniquePublications, 
    nIncludedPublications = nIncludedPublications, 
    nDrugMeetLogic = nDrugMeetLogic,
    nPublicationsMeetLogic = nPublicationsMeetLogic,
    nCoreDrugs = nCoreDrugs,
    nCoreDrugPublications = nCoreDrugPublications,
    nSingleAnnotated = nSingleAnnotated,
    nDualAnnotated = nDualAnnotated,
    nReconciled = nReconciled, 
    lastUpdate = Sys.Date()
  )
  
  return(clinicalProgressSummary)
}

createAnimalProgressSummary <- function(ad_soles_con, con, ad_longlist){
  unique_citations <- tbl(ad_soles_con, "unique_citations")
  unique_citations_smaller <- unique_citations %>%
    select(date, uid, title, journal, year, doi, uid, url, author, abstract, keywords)
  
  unique_citations_doi <- unique_citations_smaller %>%
    select(doi, uid, year)
  
  unique_citations_doi_df <- unique_citations_doi %>%
    collect()
  
  nUniquePublications <- length(unique_citations_doi_df$uid)
  
  included <- tbl(ad_soles_con, "study_classification")  %>%
    select(uid, decision) %>%
    left_join(unique_citations_doi, by = "uid") %>%
    select(-year) %>%
    filter(decision == "include") %>%
    # anti_join(conf_abstracts, by = "doi") %>%
    select(-doi)
  
  included_with_metadata <- unique_citations_smaller %>%
    inner_join(included, by = "uid") %>%
    collect()
  
  nIncludedPublications <- length(included_with_metadata$uid)
  
  ad_intervention_citations <- dbReadTable(con, "ad_invivo_citations")
  
  nDrugMeetLogic <-  nrow(DrugsMeetLogicTable)
  nPublicationsMeetLogic <- ad_intervention_citations %>% filter(intervention %in% drugMeetsLogic)%>%select(uid)%>%unique()%>%nrow()
  nCoreDrugPublications <- ad_intervention_citations %>% filter(intervention %in% ad_longlist$Drug)%>%select(uid)%>%unique()%>%nrow()
  animalProgressSummary <- data.frame(
    ReviewType = 'In vivo', 
    nUniquePublications = nUniquePublications, 
    nIncludedPublications = nIncludedPublications, 
    nDrugMeetLogic = nDrugMeetLogic,
    nPublicationsMeetLogic = nPublicationsMeetLogic,
    nCoreDrugs = nCoreDrugs,
    nCoreDrugPublications = nCoreDrugPublications,
    nSingleAnnotated = 0,
    nDualAnnotated = 0,
    nReconciled = 0, 
    lastUpdate = Sys.Date()
  )
  
  return(animalProgressSummary)
}