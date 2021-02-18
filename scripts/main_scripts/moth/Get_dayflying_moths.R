
####   Getting day flying moths   ####

library(tidyverse)
library(readr)

spdf <- read_csv("data/edited_insect_data/moth/moth_no_duplicates.csv")

# day fliers as based on the rmarkdown document and 
# the expert decision makers
dfm <- c("Jordanita globulariae",
         "Adscita statices",
         "Adscita geryon",
         "Zygaena purpuralis",
         "Zygaena loti",
         "Zygaena exulans",
         "Zygaena viciae",
         "Zygaena filipendulae",
         "Zygaena lonicerae",
         "Zygaena lonicerae/trifolii",
         "Zygaena trifolii",
         # "Lasiocampa quercus",
         "Hemaris tityus",
         "Hemaris fuciformis",
         "Idaea muricata",
         "Scotopteryx bipunctaria",
         "Epirrhoe tristata",
         "Minoa murinata",
         "Rheumaptera hastata",
         "Odezia atrata",
         "Perizoma albulata",
         "Eupithecia pygmaeta",
         "Archiearis parthenias",
         "Archiearis notha",
         "Macaria carbonaria",
         "Macaria brunneata",
         "Chiasmia clathrata",
         "Pseudopanthera macularia",
         "Lycia lapponaria",
         "Lycia zonaria",
         "Ematurga atomaria",
         "Glacies coracina",
         "Siona lineata",
         "Perconia strigillaria",
         "Orgyia antiqua",
         "Orgyia recens",
         "Parasemia plantaginis",
         "Tyria jacobaeae",
         "Phytometra viridaria",
         "Euclidia glyphica",
         "Euclidia mi",
         "Tyta luctuosa",
         "Panemeria tenebrata",
         "Heliothis viriplaca",
         "Heliothis maritima",
         # "Eremobia ochroleuca",
         "Photedes captiuncula",
         # "Cerapteryx graminis",
         "Anarta melanopa",
         "Anarta myrtilli",
         "Anarta cordigera")

dfm <- dfm[!duplicated(dfm)]

dfm_df <- spdf[spdf$sp_n %in% dfm,]
getwd()
# write.csv(dfm_df, file = "data/edited_insect_data/moth/DayFlyingMoths.csv")