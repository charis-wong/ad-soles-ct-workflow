library(dplyr)
library(fst)
library(DBI)
library(RMySQL)
library(RPostgres)
library(RSQLite)
source("getADSOLESData/configure.R")

##generate list of studies from ad-soles tagged for intervention for relisyr-ad

con <- pool::dbPool(RPostgres::Postgres(),
                    dbname = Sys.getenv("ad_soles_dbname"),
                    host = Sys.getenv("ad_soles_host"),
                    port = 5432,
                    user = Sys.getenv("ad_soles_user"),
                    password = Sys.getenv("ad_soles_password"))



unique_citations <- tbl(con, "unique_citations")
unique_citations_smaller <- unique_citations %>%
  select(date, uid, title, journal, year, doi, uid, url, author, abstract, keywords)

unique_citations_doi <- unique_citations_smaller %>%
  select(doi, uid, year)

unique_citations_doi_df <- unique_citations_doi %>%
  collect()

nUniquePublications <- length(unique_citations_doi_df$uid)

included <- tbl(con, "study_classification")  %>%
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

included_with_metadata$year <- as.numeric(included_with_metadata$year)


pico_dictionary<- tbl(con, "pico_dictionary")
pico_dictionary <- pico_dictionary %>%select(id, name)%>%collect()

pico_ontology <- tbl(con, "pico_ontology")



interventionlist <- pico_ontology%>%filter(type == "intervention") %>% filter(regex_id != "9999993")%>% select(regex_id)%>% collect()
interventionlist <- left_join(interventionlist, pico_dictionary, by = c("regex_id" = "id"))
interventionlist <- interventionlist%>%unique()

pico_tag <-tbl(con, "pico_tag")
pico_tag <- pico_tag%>%filter(regex_id%in% !!interventionlist$regex_id,
                              method == "tiabkw_regex",
                              uid %in% !!included_with_metadata$uid)%>%select(-method, -frequency)%>%collect()


intervention_tag <- left_join(pico_tag, interventionlist, by = "regex_id")



intervention_citations<- unique_citations%>%filter(uid %in% local(intervention_tag$uid))
intervention_citations <- intervention_citations%>%
  select(uid, 
         title,
         author,
         journal, 
         abstract,
         year,
         doi,
         secondarytitle,
         author_affiliation,
         ptype,
         keywords,
         doi, 
         url
         )%>% collect()


intervention_citations <- left_join(intervention_citations, intervention_tag, by = "uid")

intervention_citations <- intervention_citations%>%select(-regex_id)%>%rename(intervention = name)

intervention_citations <- intervention_citations%>%dplyr::mutate_if(is.character, stringi::stri_enc_toascii)

relisyrCon <- DBI::dbConnect(
  RMySQL::MySQL(),
  dbname = Sys.getenv("relisyr_dbname"),
  host = Sys.getenv("relisyr_host"),
  port = 3306,
  user = Sys.getenv("relisyr_user"),
  password = Sys.getenv("relisyr_password")
)

DBI::dbWriteTable(relisyrCon, "ad_invivo_citations", intervention_citations, overwrite = TRUE)

#check upload
# x <- dbReadTable(relisyrCon, "ad_invivo_citations")


dbDisconnect(relisyrCon)


