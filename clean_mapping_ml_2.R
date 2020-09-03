
keepImportantColumns <- function(mapping, metabol = F){
  mapping <- as.data.frame(mapping)
  #keep <- colnames(ps@sam_data)[c(1, 4, 11, 12, 17, 22, 23, 38, 39, 40, 41, 42, 43, 48, 49, 53, 54, 105, 106, 107,
  #                      112:118, 138, 140)]
  keep_names <- c('Biospecimen.Barcode', 'Biospecimen.Name', 'Host.Name',"Within.study.sampling.date..Biospecimen.",
                  "Family.group.ID..Biospecimen.",
                  "Pet.dog..Biospecimen.","Pet.cat..Biospecimen.", 
                  "Minimum.time.since.antibiotics..Biospecimen.",
                  "Stool.frequency..Biospecimen.","Born.by.caesarian.section..Biospecimen.",
                  "Racial.group..Biospecimen.","Specific.food.allergy..Biospecimen.","GI.issues.this.week..M3...Biospecimen." ,
                  "Other.GI.symptoms..M3...Biospecimen.","Biological.sex..Biospecimen.",
                  "Gluten.allergy..Biospecimen.","Non.celiac.gluten.sensitivity..Biospecimen.",
                  "Whole.grain..consumption.frequency...Biospecimen.","Fermented.vegetable..consumption.frequency...Biospecimen.",
                  "Dairy..consumption.frequency...Biospecimen.","Fruit..consumption.frequency...Biospecimen.",
                  "Meals.prepared.at.home..consumption.frequency...Biospecimen.","Ready.to.eat.meals..consumption.frequency...Biospecimen.",
                  "Meat..consumption.frequency...Biospecimen.","Olive.oil.used.in.cooking..M3...Biospecimen.",
                  "Seafood..consumption.frequency...Biospecimen.","Sweetened.drink..consumption.frequency...Biospecimen.",
                  "Vegetable..consumption.frequency...Biospecimen.","Restaurant.prepared.meals..consumption.frequency...Biospecimen.",
                  "Sugary.food..consumption.frequency...Biospecimen.",
                  "Diet.during.infancy..Biospecimen.","Lactose.intolerance..Biospecimen.","Age..months...Biospecimen.",
                  "Multivitamin..Biospecimen.","Born.prematurely..Biospecimen.","Probiotic..consumption.frequency...Biospecimen.",
                  "Dietary.restrictions..M3...Biospecimen.","Environmental.tobacco.smoke.exposure..Biospecimen.",
                  "Dietary.supplement..Biospecimen.","Vitamin.B.complex.supplement..consumption.frequency...Biospecimen.",
                  "Vitamin.D..consumption.frequency...Biospecimen.",
                  "Diarrhea..Biospecimen.","Constipation..Biospecimen.",
                  "Starchy.food..consumption.frequency...longitudinal...Biospecimen.",
                  "Meats.and.seafood..consumption.frequency...longitudinal...Biospecimen.",
                  "Bread..consumption.frequency...longitudinal...Biospecimen.","Dairy..consumption.frequency...longitudinal...Biospecimen.",
                  "Dietary.fat.and.oil..consumption.frequency...longitudinal...Biospecimen.",
                  "Vegetable..consumption.frequency...longitudinal...Biospecimen.",
                  "Fruit..consumption.frequency...longitudinal...Biospecimen.",
                  "Natural.father.current.age..months...Biospecimen.", "Natural.mother.current.age..months...Biospecimen.",
                  "Recently.ill..Biospecimen.",
                  "Mobile.Autism.Risk.Assessment.Score..Biospecimen.", "Biospecimen.Date.Collected")
  
  if(metabol){
    #annoyingly, the metabolomics mapping file has colnames with different format
    keep_names <- gsub("\\.$", "",  gsub("\\.\\.Biospecimen", "", keep_names))
  }
  
  
  renames <- c("biospecimen_id", "biospecimen_name", "host_name", "timepoint",
               "familyID", 
               "dog", "cat", "min_time_antibiotics", 
               "stool_freq", "csection", 
               "racial_group","specific_food_allergy", "GI_issues_this_week", 
               "other_GI_symptoms", "sex", 
               "gluten_allergy", "nonceliac_sensitivity",
               "whole_grain", "fermented_vegetables", "dairy", "fruit",
               "meal_home_prep", "meal_ready_eat", "meat", "olive_oil", "seafood", "sweetened_drink", "vegetable",
               "restaurant", "sugary_food", 
               "infant_diet", "lactose_intolerance", "age",
               "multivitamin", "prematurely_born", "probiotic", "dietary_restriction", "env_tobacco", 
               "dietary_supplement", "vitamin_B", "vitamin_D", 
               "diarrhea","constipation",
               "starchy_food", "meats_and_seafood", "bread", "dairy_freq", "fat_oil_freq", "vegetable_freq", 
               "fruit_freq", "father_age", "mother_age", "recently_ill", "MARA", "date")
  mapping <- mapping[, keep_names]
  phenotypes <- sapply(strsplit(as.character(mapping$Host.Name), "_"), `[`,2)
  colnames(mapping) <- renames
  mapping$phenotype <- as.character(phenotypes)
  return(mapping)
}

makeFieldsNumeric <- function(map){
  handleNAs <- function(vec){
    vec[vec == ""] <- "NA"
    vec[is.na(vec)] <- "NA"
    return(vec)
  }
  
  breast_fed <- c("089_A","054_N", "158_N" )
  map <- map[!(as.character(map$host_name) %in% breast_fed), ] #Child is breastfed 
  
  map$sex <- ifelse(map$sex == "Male", 0, 1)
  map$stool_freq <- handleNAs(as.character(map$stool_freq))
  map$stool_freq[as.character(map$stool_freq) == "Less than 1"] = 0
  map$stool_freq[as.character(map$stool_freq) == "5 or more"] = 5
  map$dairy_freq[map$dairy_freq == 5] <- "3-4 meals per week"
  #map$LR2[map$LR2 == "1 (ASD)"] = 1
  #map$LR2[map$LR2 == "0 (non-ASD"] = 0
  
  mother_milk_dict <- list("Mother's milk diet" = 1, "Mixed infant formula and mother's milk diet" = 2, "Infant formula diet" = 3, "NA" = NA)
  map$infant_diet <- handleNAs(map$infant_diet)
  map$infant_diet <- unlist(mother_milk_dict[map$infant_diet])
  
  freq_dict_2 <- list("Never" = 0, "Rarely" = 1, "Occasionally" = 2, "Regularly" = 3, "Weekly" = 4,
                      "Several time weekly" = 5, "Several times weekly" = 5, "Daily" = 6, "NA" = NA)
  dict_2_items <- c("whole_grain", "fermented_vegetables", "dairy", "fruit",
                    "meal_home_prep", "meal_ready_eat", "meat", "olive_oil", "seafood", "sweetened_drink", "vegetable",
                    "restaurant", "sugary_food", "probiotic", "vitamin_B", "vitamin_D")
  for(item in dict_2_items){
    print(item)
    tmp <- rep(NA, nrow(map))
    freqs <- handleNAs(map[,item])
    numeric_rep <- unlist(freq_dict_2[freqs])
    print(paste("Numeric rep length: ", length(numeric_rep)))
    print(sum(!is.na(freqs)))
    tmp[!is.na(freqs)] <- as.numeric(numeric_rep)  
    map[ , item] <- tmp
  }
  
  
  freq_dict_1 <- list("Never or less than once per week" = 0, "3-4 meals per week" = 1, "5" = 2, "7-10 meals per week" = 3, "Almost every meal" = 4, "NA" = NA)
  dict_1_items <- c("starchy_food", "meats_and_seafood", "bread", "dairy_freq", "fat_oil_freq", "vegetable_freq", 
                    "fruit_freq")
  for(item in dict_1_items){
    print(item)
    tmp <- rep(NA, nrow(map))
    freqs <- handleNAs(map[ , item])
    numeric_rep <- unlist(freq_dict_1[freqs])
    print(paste("Numeric rep length: ", length(numeric_rep)))
    print(sum(!is.na(freqs)))
    tmp[!is.na(freqs)] <- as.numeric(numeric_rep)  
    map[ , item] <- tmp
  }
  map <- map[!duplicated(map$biospecimen_id), ]
  rownames(map) <- map$biospecimen_name
  
  return(map)
  
}





#These people have really terrible metadata, and are not part of the 16S data seemingly
#unique(map$host_name[!(map$host_name %in% map_old$host_name)])
#[1] "047_A" "047_N" "016_A" "010_N" "010_A" "046_N" "046_A" "098_A"
#[9] "098_N" "191_A" "191_N" "029_A" "029_N" "041_A" "041_N" "049_N"
#[17] "049_A" "072_A" "072_N" "074_A" "074_N" "106_A" "106_N" "124_A"
#[25] "124_N" "127_A" "127_N" "119_A" "119_N" "149_A" "158_A" "149_N"
#[33] "158_N" "165_N" "165_A" "141_A" "141_N" "170_A" "170_N" "133_A"
#[41] "133_N" "132_A" "132_N"


