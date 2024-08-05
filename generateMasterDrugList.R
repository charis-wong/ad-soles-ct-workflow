
source('R/dictionaryValidationHelpers.R')
library(AutoAnnotation)
library(DBI)
library(dplyr)
library(googlesheets4)
library(purrr)
library(readr)
library(readxl)
library(RMySQL)
library(rvest)
library(stringr)
library(tibble)
library(tidyr)
library(xml2)

drugBankDictionary <- read_sheet(Sys.getenv("relisyr_ad_gsheet"), "drugDictionary")
drugBankVocabulary <- read.csv("data/drugbank vocabulary.csv")
allDrugs <- left_join(drugBankDictionary, drugBankVocabulary, by=c("Name" = "Common.name"))
allDrugs$CAS <-as.character(allDrugs$CAS)


#get chembl data----
allDrugs <- allDrugs%>% mutate(row = row_number())

for (i in 1:nrow(allDrugs)){
  drugi<-allDrugs[i, ]
  drug<- drugi$Name%>%tolower()
  drug <-gsub(" ", "%20", drug)
  row<-i
  url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule?pref_name__iexact=", drug)
  
  data <- read_html(url)
  molecules <- data %>%html_node("molecules")%>%html_children()
  if(length(molecules)==0){
    url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule?molecule_synonyms__molecule_synonym__iexact=", drug)
    data <- read_html(url)
    molecules <- data %>%html_node("molecules")%>%html_children()
  }
  if(length(molecules)>0){
    for(n in 1:length(molecules)){
      data<-molecules[n]
      name <- data %>%html_nodes("pref_name")%>%html_text()%>%paste(collapse="; ")
      chemblIDs <- data %>%html_nodes("molecule_chembl_id")%>%html_text()%>%unique()%>%paste(collapse="; ")
      synonyms <- data %>%html_nodes("molecule_synonym")%>%html_text()%>%unique()%>%paste(collapse="; ")
      SMILES <- data %>% html_node("canonical_smiles")%>%html_text()%>%paste(collapse="; ")  
      mol_formula<-data%>%html_node("full_molformula")%>%html_text()%>%paste(collapse="; ")
      mol_weight<-data%>%html_node("full_mwt")%>%html_text()%>%paste(collapse="; ")
      mol_inchi<- data%>%html_node("standard_inchi")%>%html_text()%>%paste(collapse="; ")
      mol_type<-data%>%html_node("molecule_type")%>%html_text()%>%paste(collapse="; ")
      natural_product <- data %>%html_node("natural_product")%>%html_text()%>%paste(collapse="; ")
      natural_product <- ifelse(natural_product =="", FALSE, TRUE)%>%paste(collapse="; ")
      oral <- data %>%html_node("oral")%>%html_text()
      oral <- ifelse(oral =="", FALSE, TRUE)%>%paste(collapse="; ")
      class <- data%>%html_node("usan_stem_definition")%>%html_text()%>%paste(collapse="; ")
      ro5Violations <- data%>%html_node("num_ro5_violations")%>%html_text()
      # ro5 <- ifelse(ro5Violations == "", TRUE, FALSE)%>%paste(collapse="; ")
      ro3_pass <- data %>%html_node("ro3_pass")%>%html_text()
      ro3 <- ifelse(ro3_pass == "Y", TRUE, FALSE)%>%paste(collapse="; ")
      max_phase <- data %>% html_node("max_phase")%>%html_text()%>%paste(collapse="; ")
      availability_type <- data %>%html_node("availability_type")%>%html_text()
      df<- data.frame(row, name, chemblIDs, synonyms, SMILES, mol_formula, mol_weight, mol_type, natural_product, class, ro5Violations, ro3, max_phase, availability_type, oral)
      
    }
    
    if(!exists("outputDF")) outputDF <- df else outputDF<-rbind(outputDF, df)
  }}


outputDF <- outputDF%>%
  rename(chemblName = name)

allDrugs <-left_join(allDrugs, outputDF,  by = 'row')


#get BBB data from B3DB------
# https://www.nature.com/articles/s41597-021-01069-5
B3DB <- read_tsv("data/B3DB_classification.tsv")

myDrug <- as.data.frame(unique(B3DB$compound_name))
names(myDrug) <- "Name"
myDrug$id <- seq_along(myDrug$Name)

minimumIncludeFrequencies <- 6
maximumExcludeFrequencies <- 1

results1 <- ExtractDrug(myDrug, dictionaryName =  as.data.frame(drugBankDictionary), idColumn = "id", textSearchingHeaders = c("Name"), minimumIncludeFrequency = 6,  maximumExcludeFrequency = 1, varnames = c("Drug", "DrugFrequency"), groupAnnotation=F, ignoreCase=T)


B3DBRegex <- merge(myDrug, results1)
B3DBRegex <- B3DBRegex%>%
  select(compound_name = Name, 
         Name = Drug)

B3DB_data <- B3DB%>%select(compound_name, `BBB+/BBB-`)
B3DBRegexData <- right_join(B3DB_data, B3DBRegex, by = "compound_name", relationship = "many-to-many") %>%select(Name, compound_name, B3DB = `BBB+/BBB-`)

allDrugsB3DB_BBB <- left_join(allDrugs1, B3DBRegexData, by="Name", relationship = "many-to-many")

allDrugsBBB <- allDrugsB3DB_BBB%>%
  group_by(row)%>%
  summarise(BBB1 = ifelse(length(unique(B3DB)) == 1, first(B3DB), "BBB+"))


allDrugs <- left_join(allDrugs, allDrugsBBB, by = "row")

allDrugs <- allDrugs %>% 
  mutate(BBB1 = case_match(BBB1, 
                                    "BBB+" ~ TRUE,
                                    "BBB-" ~ FALSE
                                    ))


#get admetsar data from previous
ican_mnd_druglist <- read.csv("data/2024-03-13ICAN-MNDallDrugsList.csv")
admetSarBBB <- ican_mnd_druglist%>%select(Name, 
                                          admetSAR_BBB = CNSPenetrance, 
                                          admetSAR_p_CNSPenetrance = p_CNSPenetrance)
allDrugs <- left_join(allDrugs, admetSarBBB, by = "Name", relationship = "many-to-many")
naBBB <- allDrugs %>% filter(is.na(BBB1) & is.na(admetSAR_BBB))

#combine both sources, if discrepancy, bring forward BBB=TRUE leaning towards being overinclusive
allDrugs1 <- allDrugs%>%
  filter(!row %in% naBBB$row) %>%
  mutate(BBB = ifelse(
    is.na(BBB1), admetSAR_BBB, ifelse(
      is.na(admetSAR_BBB), BBB1, ifelse(   
        BBB1 == "BBB+"|BBB1=="unclear"|admetSAR_BBB == TRUE, TRUE, FALSE)
    ))
  )%>%
  select(row, BBB)%>%
  unique()

allDrugs <- left_join(allDrugs, allDrugs1, by = "row")

#get BNF Snomed data-----
#from NHS business authority: https://www.nhsbsa.nhs.uk/prescription-data/understanding-our-data/bnf-snomed-mapping
BNF <- read_xlsx("data/BNF Snomed Mapping data 20240722.xlsx", "Sheet1")

# download BNF to dmd mapping (source: https://raw.githubusercontent.com/ebmdatalab/bnf-code-to-dmd/master/data/bnf_to_dmd.csv;
#info: https://www.bennett.ox.ac.uk/blog/2023/11/bnf-to-dictionary-of-medicines-and-devices-dm-d-map-now-available/)
bnf_to_dmd <- read.csv("data/bnf_to_dmd.csv")

bnf_to_vtm_nm <- bnf_to_dmd %>% select(bnf_code, vtm_nm)


BNF <- left_join(BNF, bnf_to_vtm_nm, by = c("BNF Code" = "bnf_code"))

BNF1 <- BNF%>%
  mutate(generic = ifelse(str_sub(`BNF Code`, 10,11) == "AA", T, F)) %>%
  group_by(vtm_nm)%>%
  summarise(generic = ifelse(sum(generic)>0, T, F ))

bnf_to_vtm_nm1 <- data.frame(vtm_nm = setdiff(bnf_to_vtm_nm$vtm_nm, BNF1$vtm_nm),
                             generic = NA)

bnfDrugs <- rbind(BNF1, bnf_to_vtm_nm1)
                             
bnfDrugs <- bnfDrugs%>%filter(vtm_nm != "")

myBnfDrug <- as.data.frame(unique(bnfDrugs$vtm_nm))
names(myBnfDrug) <- "Name"
myBnfDrug$id <- seq_along(myBnfDrug$Name)


results2 <- ExtractDrug(myBnfDrug, dictionaryName =  as.data.frame(drugBankDictionary), idColumn = "id", textSearchingHeaders = c("Name"), minimumIncludeFrequency = 6,  maximumExcludeFrequency = 1, varnames = c("Drug", "DrugFrequency"), groupAnnotation=F, ignoreCase=T)

BNFRegex <- merge(myBnfDrug, results2)

bnfDrugs<- left_join(bnfDrugs, select(BNFRegex, c(Name, Drug)), by = c("vtm_nm" = "Name"))

allDrugs <- allDrugs %>%
  mutate(BNF = ifelse(Name %in% bnfDrugs$Drug, TRUE, FALSE),
         BNFgeneric = ifelse(Name %in% filter(bnfDrugs, generic == TRUE)$Drug, TRUE, FALSE))

# write output----
write.csv(allDrugs, paste0("output/", 
                           Sys.Date(),
                           "masterDrugList.csv"), row.names = FALSE)
