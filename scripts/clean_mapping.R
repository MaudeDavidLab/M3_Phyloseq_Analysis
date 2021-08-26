
makeFieldsNumeric <- function(map){
  handleNAs <- function(vec){
    vec[vec == ""] <- "NA"
    vec[is.na(vec)] <- "NA"
    return(vec)
  }
  
  map$Stool.frequency <- handleNAs(as.character(map$Stool.frequency))
  map$Stool.frequency[as.character(map$Stool.frequency) == "Less than 1"] = 0
  map$Stool.frequency[as.character(map$Stool.frequency) == "5 or more"] = 5
  map$Dairy..consumption.frequency...longitudinal.[map$Dairy..consumption.frequency...longitudinal. == 5] <- "3-4 meals per week"
  #map$LR2[map$LR2 == "1 (ASD)"] = 1
  #map$LR2[map$LR2 == "0 (non-ASD"] = 0
  
  
  freq_dict_2 <- list("Never" = "0", "Rarely" = "1", "Occasionally" = "2", "Regularly" = "3", "Weekly" = "4", "weekly" = "4",
                      "Several time weekly" = "5", "Several times weekly" = "5", "Daily" = "6", "NA" = "NA")
  dict_2_items <- c("Whole.grain..consumption.frequency.", "Fermented.vegetable..consumption.frequency.", "Dairy..consumption.frequency.","Meat..consumption.frequency." , "Fruit..consumption.frequency.", "Meals.prepared.at.home..consumption.frequency.",   "Ready.to.eat.meals..consumption.frequency.", "Red.meat..consumption.frequency.", "Olive.oil.used.in.cooking..M3.", "Seafood..consumption.frequency.",   "Sweetened.drink..consumption.frequency.", "Vegetable..consumption.frequency.",
                    "Restaurant.prepared.meals..consumption.frequency.", "Sugary.food..consumption.frequency.", "Probiotic..consumption.frequency.", "Vitamin.B.complex.supplement..consumption.frequency.", "Vitamin.D..consumption.frequency.")
  for(item in dict_2_items){
    print(item)
    #tmp <- rep(NA, nrow(map))
    #freqs <- handleNAs(map[,item])
    #numeric_rep <- unlist(freq_dict_2[freqs])
    tmp <- map
    tmp[,item] <- as.character(tmp[,item])
    for (n in names(freq_dict_2)) {
      tmp[,item][which(tmp[,item] == n )] <- rep(freq_dict_2[[n]], length(tmp[,item][which(tmp[,item] == n )]))
    }
    map[,item] <- as.numeric(tmp[,item])
    #numeric_rep <-as.numeric(numeric_rep)
    #print(paste("Numeric rep length: ", length(numeric_rep)))
    #print(sum(!is.na(freqs)))
    #tmp[!is.na(freqs)] <- as.numeric(numeric_rep)  
    #map[ , item] <- tmp
  }
  
  freq_dict_1 <- list("Never or less than once per week" = "0", "3-4 meals per week" = "1", "5" = "2", "7-10 meals per week" = "3", "Almost every meal" = "4", "NA" = "NA")
  dict_1_items <- c("Starchy.food..consumption.frequency...longitudinal.", "Meats.and.seafood..consumption.frequency...longitudinal.", "Bread..consumption.frequency...longitudinal.", "Dairy..consumption.frequency...longitudinal.", "Dietary.fat.and.oil..consumption.frequency...longitudinal.", "Vegetable..consumption.frequency...longitudinal.", 
                    "Fruit..consumption.frequency...longitudinal.")
  for(item in dict_1_items){
    print(item)
    #tmp <- rep(NA, nrow(map))
    #freqs <- handleNAs(map[ , item])
    #numeric_rep <- unlist(freq_dict_1[freqs])
    #numeric_rep <- as.numeric(numeric_rep)
    #print(paste("Numeric rep length: ", length(numeric_rep)))
    #print(sum(!is.na(freqs)))
    #tmp[!is.na(freqs)] <- as.numeric(numeric_rep)  
    #map[ , item] <- tmp
    
    tmp <- map
    tmp[,item] <- as.character(tmp[,item])
    for (n in names(freq_dict_1)) {
      tmp[,item][which(tmp[,item] == n )] <- rep(freq_dict_1[[n]], length(tmp[,item][which(tmp[,item] == n )]))
    }
    map[,item] <- as.numeric(tmp[,item])
  }
  
  #may add more, but these variable only apply to phenotype for autism
  freq_dict_2 <- list("Able to speak fluently" = 4,"Phrase speech"=3, "Single word speech"=2, "Little to no speech" = 1,"Able to have conversation" = 4, "Limited conversation ability" = 3, "Difficulty with conversation" = 2, "Cannot have a conversation" = 1, "Understands about half of words" = 1, "Understands few or no words"= 0, "Understands many words" = 2,  "Understands most words"= 3, "Understands nearly all words" = 4 , "Never" = 1, "Rarely" = 2, "Sometimes" = 3 , "Regularly" = 4, "No opportunity to play with other children" = NA, "Consistent eye contact" = 4,  "Little or no eye contact" = 1, "Some eye contact" = 3, "Uncertain about behavioral development" = NA, "Developmentally delayed (autism)" = 1, "Developmentally delayed (not clearly autism)" = 2, "Some developmental delays" = 3, "Met all developmental milestones" = 4 , "Does not imitate others" =1, "Imitates others when prompted" = 2, "Imitates actions or gestures of others" = 4 ,"NA" = 4)
  dict_2_items <- c("Language.ability.and.use", "Conversation.ability", "Understands.speech", "Plays.imaginatively.when.alone", "Plays.imaginatively.with.others",  "Plays.in.a.group.with.others", "Eye.contact.finding", "Childhood.behavioral.development.finding", "Picks.up.objects.to.show.to.others", "Imitation.behavior")
  for(item in dict_2_items){
    print(item)
    #tmp <- rep(NA, nrow(map))
    #freqs <- handleNAs(map[ , item])
    #numeric_rep <- unlist(freq_dict_2[freqs])
    #print(paste("Numeric rep length: ", length(numeric_rep)))
    #print(sum(!is.na(freqs)))
    #tmp[!is.na(freqs)] <- as.numeric(numeric_rep)  
    #tmp[is.na(tmp)] <- 4
    #map[ , item] <- tmp
    
    tmp <- map
    tmp[,item] <- as.character(tmp[,item])
    for (n in names(freq_dict_2)) {
      tmp[,item][which(tmp[,item] == n )] <- rep(freq_dict_2[[n]], length(tmp[,item][which(tmp[,item] == n )]))
    }
    map[,item] <- as.numeric(tmp[,item])
  }
  
  #may add more, but these variable only apply to phenotype for autism
  freq_dict_3 <- list( "Never" = 1, "Sometimes" = 2, "Regularly" = 3,  "NA" = 1, "Constant sleep difficulties" = 3,"Some sleep difficulties" = 2, "Healthy sleep pattern" = 1, "Highly sensitive to typical sounds" = 3,  "Sensitive to typical sounds" = 2, "Not bothered by typical sounds" = 1, "No self-injurious behavior" = 1,"Mild self-harming behavior" = 2, "Dangerous or frequent self-harming behavior" = 3,  "No issues" = 1, "Continuous" = 3, "No elevated anxiety" = 1, "Somewhat elevated anxiety" = 2, "Elevated anxiety" = 3)
  dict_3_items <- c("Repetitive.motion", "Sleep.pattern.finding", "Response.to.typical.sounds", "Self.injurious.behavior.finding", "Gastrointestinal.problems..M3.", "Recent.anxiety..caretaker.reported.")
  for(item in dict_3_items){
    print(item)
    #tmp <- rep(NA, nrow(map))
    #freqs <- handleNAs(map[ , item])
    #numeric_rep <- unlist(freq_dict_3[freqs])
    #print(paste("Numeric rep length: ", length(numeric_rep)))
    #print(sum(!is.na(freqs)))
    #tmp[!is.na(freqs)] <- as.numeric(numeric_rep)
    #tmp[is.na(tmp)] <- 1
    #map[ , item] <- tmp
    
    tmp <- map
    tmp[,item] <- as.character(tmp[,item])
    for (n in names(freq_dict_3)) {
      tmp[,item][which(tmp[,item] == n )] <- rep(freq_dict_3[[n]], length(tmp[,item][which(tmp[,item] == n )]))
    }
    map[,item] <- as.numeric(tmp[,item])
  }
  
  
  map <- map[!duplicated(map$Biospecimen.Barcode), ]
  rownames(map) <- map$Biospecimen.Barcode
  map$Stool.frequency <- as.numeric(map$Stool.frequency)
  return(map)
  
}

dict_1_items <- c("Starchy.food..consumption.frequency...longitudinal.", "Meats.and.seafood..consumption.frequency...longitudinal.", "Bread..consumption.frequency...longitudinal.", "Dairy..consumption.frequency...longitudinal.", "Dietary.fat.and.oil..consumption.frequency...longitudinal.", "Vegetable..consumption.frequency...longitudinal.", 
                  "Fruit..consumption.frequency...longitudinal.")
dict_2_items <- c("Language.ability.and.use", "Conversation.ability", "Understands.speech")


clean_mapping_file <- function(map){
  map <- makeFieldsNumeric(map)
  
  
  map_levels<-sapply(map, levels)
  map_levelscount<-sapply(map_levels, length)
  mapnotfac <- names(map_levelscount[which(map_levelscount >= 18)])
  
  for (i in mapnotfac){
    map[,i]<-as.character(map[,i])
  }
  
  #Round years
  map$Age..years. <-round(map$Age..years.)
  
  #create scores for social, learning, and anxiety
  social <- c("Language.ability.and.use", "Conversation.ability", "Plays.imaginatively.when.alone", "Plays.imaginatively.with.others", "Plays.in.a.group.with.others", "Eye.contact.finding", "Picks.up.objects.to.show.to.others", "Imitation.behavior")
  
  learning <- c("Understands.speech", "Childhood.behavioral.development.finding")
  
  anxiety <- c("Repetitive.motion", "Sleep.pattern.finding", "Response.to.typical.sounds", "Self.injurious.behavior.finding", "Recent.anxiety..caretaker.reported.")
  
  map$social_score_sum<-rowSums(map[, social])
  map$learning_score_sum<-rowSums(map[, learning])
  map$anxiety_score_sum <-rowSums(map[, anxiety])
  
  map$social_score_avg<-rowMeans(map[, social])
  map$learning_score_avg<-rowMeans(map[, learning])
  map$anxiety_score_avg <-rowMeans(map[, anxiety])
  
  
  
  #add season during collection
  map$Season <- gsub("-01-", "Winter", map$Biospecimen.Date.Collected)
  wint <- c("-02-", "-12-")
  for (i in wint){
    map$Season <- gsub(i, "Winter", map$Season)
  }
  
  spring <- c("-03-", "-04-", "-05-")
  for (i in spring){
    map$Season <- gsub(i, "Spring", map$Season)
  }
  
  summer <- c("-06-", "-07-", "-08-")
  for (i in summer){
    map$Season <- gsub(i, "Summer", map$Season)
  }
  
  fall<- c("-09-", "-10-", "-11-")
  for (i in fall){
    map$Season <- gsub(i, "Fall", map$Season)
  }
  
  map$Season <- gsub("2018","", map$Season)
  map$Season <- gsub("2017","", map$Season)
  
  for (i in c("1","2","3", "4", "5", "6", "7", "8", "9", "0")){
    map$Season <- gsub(i, "", map$Season)
  }
  map$Season <- as.factor(map$Season)
  
  map$Annual.household.income_rank <- gsub("150001-200000", 5, map$Annual.household.income)
  
  map$Annual.household.income_rank <- gsub("More than 200000", 6, map$Annual.household.income_rank)
  
  map$Annual.household.income_rank <- gsub("Less than 20000", 1, map$Annual.household.income_rank)
  
  map$Annual.household.income_rank <- gsub("20001-40000", 2, map$Annual.household.income_rank)
  
  map$Annual.household.income_rank <- gsub("40001-80000", 3, map$Annual.household.income_rank)
  
  map$Annual.household.income_rank <- gsub("80001-150001", 4, map$Annual.household.income_rank)
  
  map$Annual.household.income_rank <- as.numeric(map$Annual.household.income_rank)
  
  
  return(map)
}
