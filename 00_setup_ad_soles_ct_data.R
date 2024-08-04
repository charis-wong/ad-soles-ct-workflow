#script to set up AD-SOLES-CT data

library(dplyr)
library(fst)
library(DBI)
library(RMySQL)
library(RPostgres)
library(RSQLite)
library(googlesheets4)

source("configure.R")
source("R/functions.R")

con <- DBI::dbConnect(
  RMySQL::MySQL(),
  dbname = Sys.getenv("relisyr_dbname"),
  host = Sys.getenv("relisyr_host"),
  port = 3306,
  user = Sys.getenv("relisyr_user"),
  password = Sys.getenv("relisyr_password")
)


#  longlist tab----
ad_longlist <- read_sheet(Sys.getenv("relisyr_mnd_gsheet"), "entityOfInterest")%>%filter(Type == "drugOfInterestAD")%>%select(Drug = Item)

ad_longlist <- data.frame(Drug = ad_longlist, DateAdded = "2024-08-01")

write_sheet(ad_longlist, Sys.getenv("relisyr_ad_gsheet"), "longlist")

# progressSummary tab----
### clinical----
#add clinical n corpus, unique publications and included publications from progress summary from relisyr-mnd
relisyrmndProgressSummary <- read_sheet(Sys.getenv("relisyr_mnd_gsheet"), "ProgressSummary")
clinicalProgressSummary <- relisyrmndProgressSummary[1,]

#add ad logic and longlist pubs from db

reviewerSession <- googlesheets4::read_sheet(Sys.getenv("relisyr_mnd_gsheet"), sheet = "reviewerSession")
reviewStage <- reviewerSession%>%
  filter(StageIdStr == "2c400348-d871-4055-9c68-2bc529ac9ccc")%>%
  group_by(StudyIdStr)%>%
  summarise(nReviews = length(InvestigatorIdStr))%>%
  rename(idStr = StudyIdStr)

reconciledPapers <-reviewerSession%>%
  filter(StageIdStr == "2c400348-d871-4055-9c68-2bc529ac9ccc")%>%
  filter(InvestigatorIdStr %in% reconcilerIdStrs)%>%
  group_by(StudyIdStr)%>%
  unique()%>%
  select(idStr = StudyIdStr)

clinicalStudyList <- RMySQL::dbReadTable(con, "ReLiSyRClinicalStudies") 
clinicalStudyList <- left_join(clinicalStudyList, reviewStage, by = "idStr")
clinicalStudyList <- clinicalStudyList%>%mutate(nReviews = ifelse(is.na(nReviews), 0, nReviews),
                                                reconcile = ifelse(idStr%in%reconciledPapers$idStr, TRUE, FALSE))

entityOfInterest <- googlesheets4::read_sheet(Sys.getenv("relisyr_mnd_gsheet"), sheet="entityOfInterest")

diseaseOfInterest <- entityOfInterest[entityOfInterest$Type == "diseaseOfInterest", ]$Item

clinicalStudyList$Disease <- factor(clinicalStudyList$Disease, levels = diseaseOfInterest)
clinicalOutputCrossTable <- as.data.frame.matrix(table(clinicalStudyList[,c("Drug","Disease")]))

clinicalOutputCrossTable$select1  <- F
clinicalOutputCrossTable$score1 <- rowSums(clinicalOutputCrossTable[, "AD", drop = F])
clinicalOutputCrossTable$score2 <- rowSums(clinicalOutputCrossTable[, setdiff(diseaseOfInterest, "AD"), drop=F] > 0)
clinicalOutputCrossTable$select1 <- clinicalOutputCrossTable$select | (clinicalOutputCrossTable$score1 > 0 | clinicalOutputCrossTable$score2 >= 2)

DrugsMeetLogicTable <- clinicalOutputCrossTable%>%filter(select1 == TRUE)%>%
  mutate(nPub = AD+PD+HD+MND+FTD+MS)

drugMeetsLogic <- rownames(DrugsMeetLogicTable)
nDrugMeetLogic = nrow(DrugsMeetLogicTable)
nCoreDrugs <- nrow(ad_longlist)

clinicalProgressSummary <- createClinicalProgressSummary(clinicalOutputCrossTable, clinicalProgressSummary, ad_longlist)
## animal progress summary----
ad_soles_con<- pool::dbPool(RPostgres::Postgres(),
                            dbname = Sys.getenv("ad_soles_dbname"),
                            host = Sys.getenv("ad_soles_host"),
                            port = 5432,
                            user = Sys.getenv("ad_soles_user"),
                            password = Sys.getenv("ad_soles_password"))




animalProgressSummary <- createAnimalProgressSummary(ad_soles_con, con, ad_longlist)

progressSummary <- rbind(clinicalProgressSummary, animalProgressSummary)

# write results
write_sheet(progressSummary, Sys.getenv("relisyr_ad_gsheet"), "progressSummary")


# annotated publications----
relisyrmndAnnotatedPublications <- read_sheet(Sys.getenv("relisyr_mnd_gsheet"), "PublicationList2")
adlonglistAnnotatedPublications <- relisyrmndAnnotatedPublications%>% filter(Drug %in% tolower(ad_longlist$Drug))
write_sheet(adlonglistAnnotatedPublications, Sys.getenv("relisyr_ad_gsheet"), "PublicationList")

# drug summary----
relisyrmndDrugSummary<- read_sheet(Sys.getenv("relisyr_mnd_gsheet"), "DrugSummary2")
adLonglistDrugSummary <- relisyrmndDrugSummary %>% filter(Drug %in% tolower(ad_longlist$Drug))
write_sheet(adLonglistDrugSummary, Sys.getenv("relisyr_ad_gsheet"), "DrugSummary")
