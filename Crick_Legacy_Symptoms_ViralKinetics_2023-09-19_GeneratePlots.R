#--------------------------------------------------------------------#
#
# Crick/UCLH Legacy Cohort Analysis:
#
# COVID-19 in non-hospitalised adults caused by either SARS-CoV-2 sub-variants 
# Omicron BA.1, BA.2, BA.5 or Delta associates with similar illness duration, 
# symptom severity and viral kinetics, irrespective of vaccination history
#
#
# 24 October 2022
# Updated 19 September 2023
#
# Creator, author:
# Ed Carr
# edward.carr@crick.ac.uk
#
# Author:
# Hermaleigh Townsley
# hermaleigh.townsley@crick.ac.uk
#
# Author:
# Timothy Russell
# Timothy.Russell@lshtm.ac.uk
#
# Author:
# Joshua Gahir
# joshua.gahir@crick.ac.uk
#
# Author:
# David LV Bauer
# david.bauer@crick.ac.uk
#
# Francis Crick Institute, London, UK
#
# This work is licensed under a CC BY 4.0 license and is free to use with attribution.
#
#--------------------------------------------------------------------#

library(tidyverse)
library(lubridate)
library(ggpubr)
library(rstatix)
library(epitools)
library(pheatmap)
library(ggplotify)
#library(systemfonts)
library(khroma)
library(patchwork)

#--------------------------------------------------------------------#
# Set working directory to the folder this script is in:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Start in peace
rm(list = ls())

#--------------------------------------------------------------------#
## Load pre-processed legacy data ####
#--------------------------------------------------------------------#

rdat.file <- list.files(
  path = "./",
  pattern = "Crick_Legacy_Symptoms_ViralKinetics_2023-09-19_PUBLIC.RData",
  full.names=T)

load(file = rdat.file)

#--------------------------------------------------------------------#
## Functions
#--------------------------------------------------------------------#
  
# Function to format clustering matrix
cluster_mat <- function(data,VOC=NULL){
  
  if(!is.null(VOC)){
    data <- data %>%
      filter(str_detect(episode_variant_summarised,{{VOC}}))
  }
  
  mat <- data %>%
    filter(!Sx_severity == "asymptomatic") %>%
    select(Sx_severity,functionalDosesAtTimeOfInfection, episode_variant_summarised,
           fever, lossSmell, cold, fatigue, diarrhoea,
           ache, cough, shrtbrth) %>%
    ## remove rows where all are NA ##
    janitor::remove_empty(., which = "rows") %>%
    ## remaining NA are zero ##
    mutate(across(c(fever, lossSmell, cold, fatigue, diarrhoea,
                    ache, cough, shrtbrth), ~ as.numeric(str_replace_na(.x, "0")))) %>%
    mutate(across(c(fever, lossSmell, cold, fatigue, diarrhoea,
                    ache, cough, shrtbrth), ~ ifelse(.x>0, 1, 0))) %>%
    arrange(episode_variant_summarised,Sx_severity)
  
  
  row_anno <- mat %>% select(functionalDosesAtTimeOfInfection, episode_variant_summarised, Sx_severity)
  row_anno <- data.frame(row_anno)
  rownames(row_anno) <- paste0("ID",rownames(row_anno)) ## this is important
  colnames(row_anno) <- c("Doses", "VOC","Severity")
  
  row_anno$Doses <- factor(row_anno$Doses)
  row_anno$VOC <- factor(row_anno$VOC)
  row_anno$Severity <- factor(row_anno$Severity)
  
  mat <- mat %>% select(!c(functionalDosesAtTimeOfInfection, episode_variant_summarised, Sx_severity))
  
  colnames(mat) <- c("Fever", "Anosmia", "Coryza", "Fatigue", "Diarrhoea", 'Myalgia', "Cough", "Dyspnoea")
  rownames(mat) <- paste0("ID",rownames(mat)) ## this is important
  
  return(list(mat = mat,
              row_anno = row_anno))
  
  
}

#--------------------------------------------------------------------#
## Plot helpers and presets ####
#--------------------------------------------------------------------#

no <- 5
weak <- 10
compl <- 5120

strainOrder <- c("wildtype", "D614G", "Alpha", "Beta", "Delta", "Omicron BA.1", "Omicron BA.2", 
                 "Omicron BA.4/5")
vocOrder <- c("\nwildtype", "\nD614G", "\nAlpha", "\nBeta", "\nDelta", "\nOmicron BA.1", "\nOmicron BA.2",
              "\nOmicron BA.4/5")

pd <- position_dodge(width=0.4)

## Adapted from khroma::muted()
manual.pal <- 
  c("#CC6677",
             "#332288",
             "#DDCC77",
             "#117733",
             "#88CCEE",
             "#882255",
             "#44AA99",
             "black")
             
#--------------------------------------------------------------------#
## PLOT: Fig 2A - How many doses received at the time of their infection? ####
#--------------------------------------------------------------------#

cts <- dosesAtTimeofInfectionAfterD2_14d %>%
  group_by(dosesAtTimeofInfection, episode_variant_summarised) %>%
  tally()

panel.DaysSinceDose <- 
  dosesAtTimeofInfectionAfterD2_14d %>%
  left_join(cts) %>%
  mutate(dose_VOC_n=
           paste0(
             episode_variant_summarised,
             "\n(n=",n, ")")) %>%
  ## remove any groups of VOC/doses with less than 5 ##
  filter(n>4) %>%
  ## force empty levels preceding Delta to give standard colours ##
  mutate(episode_variant_summarised =
           factor(episode_variant_summarised,
                  levels = c("Ancestral", "D614G", "Alpha", "Beta",
                             levels(factor(episode_variant_summarised))))) %>%
  mutate(dose_VOC_n=factor(dose_VOC_n)) %>%
  mutate(dose = paste("infection\nafter", dosesAtTimeofInfection, "doses")) %>%
  ggplot(aes(x = dose_VOC_n,
             y = daysSincedose, 
             col = episode_variant_summarised,
             group=elig_study_id)) + 
  scale_color_manual(values = manual.pal, drop=F) +
  geom_boxplot(aes(group=NULL), outlier.shape = NA) + 
  ylim(0, NA) + 
  geom_point(shape=20, alpha=0.3, position=pd, size=2) + 
  scale_x_discrete(name="") + 
  facet_grid(.~dose, scales = "free_x", space = "free") + 
  labs(y="Time since dose (days)") + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5))

panel.DaysSinceDose

#--------------------------------------------------------------------#
## PLOT: Fig 2B - Symptom severity ####
#--------------------------------------------------------------------#

panel.SxSeverity <- dosesAtTimeofInfectionAfterD2_14d %>%
  filter(! str_detect(VOC, "unkno")) %>%
  mutate(Sx_severity=str_replace_na(Sx_severity, replacement="no diary")) %>%
  mutate(Sx_severity=factor(Sx_severity, 
                            levels = c(
                              "asymptomatic",
                              "grade I",
                              "grade II",
                              "grade III","no diary"
                            ), 
                            exclude = "NULL")) %>%
  mutate(VOC = gsub(VOC, pattern="\n", replacement="")) %>%
  mutate(VOC = str_replace_all(VOC, " ", "\n")) %>%
  mutate(VOC=factor(VOC, levels = c("Delta", "Omicron\nBA.1", "Omicron\nBA.2", 
                                    "Omicron\nBA.4/5"))) %>%
  group_by(VOC,functionalDosesAtTimeOfInfection) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n>4) %>%
  ggplot(aes(x=functionalDosesAtTimeOfInfection,fill=Sx_severity)) +
  geom_bar(position=position_fill(reverse = T),col = "black",size=0.05) + 
  scale_fill_manual(values=rev(c("white",grey.colors(n=4,start=0.1,end=0.9)))) + 
  guides(fill=guide_legend(ncol=3)) +
  facet_grid(.~VOC, drop=T, scales="free_x", space = "free_x") + 
  theme(strip.text.x = element_text(angle = 90)) +
  labs(y="Proportion of\ninfection episodes") +
  # theme_bw() +
  labs(fill = "Symptom\nseverity", x = "Doses at time of infection") +
  theme(axis.title.x = element_text(vjust=+2))
panel.SxSeverity

## Chisq test for discussion ##
panel.SxSeverity$data %>%
  filter(str_detect(VOC, "Om")) %>%
  filter(functionalDosesAtTimeOfInfection==3) %>%
  filter(Sx_severity != "no diary") %>%
  group_by(VOC,Sx_severity) %>%
  tally() %>%
  group_by(VOC) %>%
  mutate(total_case_for_VOC = sum(n)) %>%
  filter(str_detect(Sx_severity, "asym")) %>%
  ungroup() %>%
  select(n,total_case_for_VOC) %>%
  mutate(n_not_asym = total_case_for_VOC-n) %>%
  select(n,n_not_asym) %>%
  t(.) %>%
  chisq_test()

# Chisq test with OR - comparing asymptomatic between VOC
dosesAtTimeofInfectionAfterD2_14d %>%
  
  mutate(VOC = episode_variant_summarised) %>% 
  
  mutate(Sx_severity=str_replace_na(Sx_severity, replacement="no diary")) %>%
  filter(Sx_severity != "no diary") %>%
  
  # Count by VOC / severity 
  group_by(VOC,Sx_severity) %>%
  tally() %>%
  group_by(VOC) %>%
  mutate(n = as.numeric(n),
         tot = as.numeric(sum(n))) %>% 
  mutate(pct = 100* n/tot) %>%
  
  # Filter to asym only 
  filter(str_detect(Sx_severity, "asym")) %>%
  ungroup() %>%
  
  # Get difference between total an asymptomatic 
  mutate(n_not = tot-n) %>% 
  select(VOC, Sx_severity, n, n_not) %>%
  
  # Make table 
  pivot_wider(id_cols = VOC, names_from = Sx_severity, values_from = c(n,n_not),values_fill = 0) %>%
  tibble::column_to_rownames('VOC') %>% 
  select(where(is.numeric))->
  stats

# BA.1 compared to BA.2 
stats %>% 
  # Filter to comparison of interest
  filter(grepl("BA.1$|BA.2$",rownames(.))) %>% 
  as.matrix(.) %>% 
  epitools::oddsratio.wald(x=., rev="rows", correction = F,verbose=T)


# BA.5 compared to BA.1 
stats %>% 
  # Filter to comparison of interest
  filter(grepl("BA.1$|BA.4/5$",rownames(.))) %>% 
  as.matrix(.) %>% 
  epitools::oddsratio.wald(x=.,correction = F,verbose=T)

#--------------------------------------------------------------------#
## PLOT: Fig 2C -  Duration of episodes ####
#--------------------------------------------------------------------#

### Note: The episode_length has been redacted. 

# cts <- dosesAtTimeofInfectionAfterD2_14d %>%
#   group_by(dosesAtTimeofInfection,episode_variant_summarised) %>%
#   tally()
#
# panel.DaysEpisodeLength <- 
#   
#   dosesAtTimeofInfectionAfterD2_14d %>%
#   left_join(cts) %>%
#   
#   ## remove any groups of VOC/doses with less than 5 ##
#   filter(n>4) %>%
#   
#   filter(Sx_severity != "asymptomatic") %>% 
#   
#   mutate(dose_VOC_n=
#            paste0(episode_variant_summarised,"\nafter ",
#                   dosesAtTimeofInfection, " doses")) %>%
#   mutate(dose_VOC_n=factor(dose_VOC_n)) %>%
#   ## force empty levels preceding Delta to give standard colours ##
#   mutate(episode_variant_summarised =
#            factor(episode_variant_summarised,
#                   levels = c("Ancestral", "D614G", "Alpha", "Beta",
#                              levels(factor(episode_variant_summarised))))) %>%
#   ggplot(aes(x = dose_VOC_n,
#              y = episode_length, col = episode_variant_summarised,
#              group = elig_study_id)) + 
#   scale_color_manual(values = manual.pal, drop=F) + 
#   geom_boxplot(aes(group=NULL), outlier.shape = NA) + 
#   geom_point(shape=20, alpha=0.3, position=pd, size=2) + 
#   scale_x_discrete(name="") + 
#   labs( y="\nDuration of symptoms (days)") + 
#   theme_bw() + 
#   theme(legend.position = "none", axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5)) 
# 
# panel.DaysEpisodeLength
# 
# ## KW test
# dosesAtTimeofInfectionAfterD2_14d %>%
#   mutate(dose_VOC_n=
#            paste0(episode_variant_summarised,"\nafter ",
#                   dosesAtTimeofInfection, " doses")) %>%
#   group_by(dose_VOC_n) %>% 
#   filter(n()>4) %>% 
#   ungroup() %>% 
#   kruskal.test(episode_length~dose_VOC_n, data =.)


#--------------------------------------------------------------------#
## PLOT: Fig 2E - part 1 - Symptom heatmap ####
#--------------------------------------------------------------------#

## extras: anosmia coryza fatigue fever chest pain iritis nausea
SxDistro <- 
  dosesAtTimeofInfectionAfterD2_14d %>% 
  ## we only want presence / absence
  ## absence is indicated by NA, presence = # of days duration
  ## a few participants wrote free text that was obviously fever (eg night sweats)
  ## So here we add +1 to relevant columns.
  mutate(cold = case_when(str_detect(parsed_description, "coryza") ~ as.numeric(str_replace_na(cold, "0")) + 1,
                          ! str_detect(parsed_description, "coryza") ~ cold,
                          is.na(parsed_description) ~ cold)) %>%
  mutate(lossSmell = case_when(str_detect(parsed_description, "anosmia") ~ as.numeric(str_replace_na(lossSmell, "0")) + 1,
                               ! str_detect(parsed_description, "anosmia") ~ lossSmell,
                               is.na(parsed_description) ~ lossSmell)) %>%
  mutate(fatigue = case_when(str_detect(parsed_description, "fatigue") ~ as.numeric(str_replace_na(fatigue, "0")) + 1,
                             ! str_detect(parsed_description, "fatigue") ~ fatigue,
                             is.na(parsed_description) ~ fatigue)) %>%
  mutate(fever = case_when(str_detect(parsed_description, "fever") ~ as.numeric(str_replace_na(fever, "0")) + 1,
                           ! str_detect(parsed_description, "fever") ~ fever,
                           is.na(parsed_description) ~ fever)) 

panel.SxDistro <- SxDistro %>%
  group_by(functionalDosesAtTimeOfInfection,
           episode_variant_summarised) %>%
  summarise(total=n(),
            Coryza_prop = sum(!is.na(cold))/total,
            Fatigue_prop = sum(!is.na(fatigue))/total,
            Myalgia_prop = sum(!is.na(ache))/total,
            Anosmia_prop = sum(!is.na(lossSmell))/total,
            Cough_prop = sum(!is.na(cough))/total,
            Fever_prop = sum(!is.na(fever))/total,
            `Dyspnoea_prop` = sum(!is.na(shrtbrth))/total,
            Diarrhoea_prop = sum(!is.na(diarrhoea))/total) %>%
  ## make sure we report summary stats for voc:doses combinations of 5 or more ##
  filter(total > 4) %>%
  ungroup() %>%
  pivot_longer(contains("_prop"),
               names_sep = "_",
               names_to = c("sx", ".value")) %>%
  mutate(Symptom = factor(sx,
                          levels = rev(c("Anosmia",
                            "Cough", "Fever",
                            "Coryza", "Fatigue",
                            "Myalgia",
                            "Dyspnoea",
                            "Diarrhoea")))) %>%
  rename(Doses=functionalDosesAtTimeOfInfection) %>%
  rename(`Participants\nReporting`=prop) %>%
  mutate(dose_VOC_n=
           paste0(
             episode_variant_summarised,
             "\n(n=",total, ")")) %>%
  mutate(dose_VOC_n = str_replace_all(dose_VOC_n, 
                                      c("Omicron-BA.1"= "Omicron\nBA.1",
                                        "Omicron-BA.2"= "Omicron\nBA.2",
                                        "Omicron-BA.4/5"= "Omicron\nBA.4/5"))) %>%
  mutate(dose_VOC_n=factor(dose_VOC_n)) %>%
  
  mutate(dose = paste("infection\nafter", Doses, "doses")) %>%
  filter(! (episode_variant_summarised == "Omicron BA.2" & Doses == 2)) %>%
  ggplot(aes(y=Symptom, x=Doses,
             fill=`Participants\nReporting`)) + geom_tile() + 
  geom_text(aes(label=100*round(`Participants\nReporting`,2)), col="white", 
            size=8/ggplot2::.pt) + 
  scale_fill_viridis_c(limits=c(0,1),
                       breaks = c(0,0.25,0.50, 0.75, 1),
                       labels = c("0%", "25%","50%", "75%", "100%")) + 
  facet_grid(~dose_VOC_n, scales="free_x", space="free") + 
  theme_bw() +
  labs(x="Doses at time of infection")

nhs.labels <- ifelse(str_detect(levels(panel.SxDistro$data$Symptom),"Cough|Fever|Anosmia"), 
                     rgb(red = 23,green = 72,blue = 169,max=255),"black")
 
panel.SxDistro + theme(axis.text.y  =element_text(colour = nhs.labels),  panel.grid = element_blank()) +
   theme(legend.position = "bottom")

#--------------------------------------------------------------------#
## PLOT: Fig 2D - Number of symptoms reported, amongst symptomatic ####
#--------------------------------------------------------------------#

SxDistroNumberOfSx <- SxDistro %>%
  filter(! str_detect(Sx_severity, "asym")) %>%
  group_by(functionalDosesAtTimeOfInfection,
           episode_variant_summarised) %>%
  summarise(total=n(),
            AggregateSxReported = 
              
              paste(cold,
                    fatigue,
                    ache,
                    lossSmell,
                    cough,
                    fever,
                    shrtbrth,
                    diarrhoea
              )) %>%
  mutate(numberOfSxReported =
           8 - str_count(AggregateSxReported, "NA"))

panel.NumbOfSx <- SxDistroNumberOfSx %>%
  filter(total >4) %>%
  
  ## force empty levels preceding Delta to give standard colours ##
  mutate(episode_variant_summarised =
           factor(episode_variant_summarised,
                  levels = c("Ancestral", "D614G", "Alpha", "Beta",
                             "Delta", "Omicron-BA.1", "Omicron-BA.2",
                             "Omicron-BA.4/5"))) %>%
  mutate(clean_x_labels = paste0(episode_variant_summarised,
                                 "\nafter ",
                                 functionalDosesAtTimeOfInfection,
                                 " doses")) %>%
  
  ggplot(
    aes(x = clean_x_labels,
        y = numberOfSxReported,
        col = episode_variant_summarised)) +
  
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  geom_point(position = position_jitter(0.2), alpha=0.3, shape=20,size=2, show.legend = F) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_manual(values = manual.pal, drop=F) + 
  labs(x = NULL, 
       #"Number of doses at time of infection per VOC",
       y = "\nNumber of symptoms reported")


# Perform tests
stat.test <- compare_means(
  numberOfSxReported ~ clean_x_labels, data = panel.NumbOfSx$data,
  method = "wilcox.test") %>% 
  mutate(p.lab = ifelse(after_stat(as.numeric(p))<=0.001,paste0("< 0.001"),
                        format(round(after_stat(as.numeric(p)),3), nsmall=2)))

# Select comparisons for plotting   
stat.test<-stat.test %>% 
  filter(
    (group1 == "Delta\nafter 2 doses" & group2 == "Omicron-BA.1\nafter 3 doses")|
      (group1 == "Omicron-BA.1\nafter 2 doses" & group2 == "Omicron-BA.1\nafter 3 doses")|
      (group1 == "Omicron-BA.1\nafter 3 doses" & group2 == "Omicron-BA.2\nafter 3 doses")|
      (group1 == "Omicron-BA.1\nafter 3 doses" & group2 == "Omicron-BA.4/5\nafter 3 doses"))

# Add to plot 
panel.NumbOfSx <- 
  panel.NumbOfSx +
  stat_pvalue_manual(stat.test,tip.length = 0.01,
                     label = "p.lab",
                     y.position = c(8,9,10,11),
                     size = 8/ggplot2::.pt)

panel.NumbOfSx


#--------------------------------------------------------------------#
## PLOT: Fig 2E - part 2 - Symptom heatmap signif tests ####
#--------------------------------------------------------------------#

SxDistroFreq <- SxDistro %>%
  group_by(functionalDosesAtTimeOfInfection,
           episode_variant_summarised) %>%
  summarise(total=n(),
            Coryza = sum(!is.na(cold)),
            NotCoryza = total-sum(!is.na(cold)),
            Fatigue = sum(!is.na(fatigue)),
            NotFatigue = total-sum(!is.na(fatigue)),
            Myalgia = sum(!is.na(ache)),
            NotMyalgia = total-sum(!is.na(ache)),
            Anosmia = sum(!is.na(lossSmell)),
            NotAnosmia = total-sum(!is.na(lossSmell)),
            Cough = sum(!is.na(cough)),
            NotCough = total-sum(!is.na(cough)),
            Fever = sum(!is.na(fever)),
            NotFever = total-sum(!is.na(fever)),
            Dyspnoea = sum(!is.na(shrtbrth)),
            NotDyspnoea = total-sum(!is.na(shrtbrth)),
            Diarrhoea = sum(!is.na(diarrhoea)),
            NotDiarrhoea = total-sum(!is.na(diarrhoea))) %>% 
  mutate(functionalDosesAtTimeOfInfection = paste0(functionalDosesAtTimeOfInfection)) %>% 
  mutate(dose_variant = paste0(functionalDosesAtTimeOfInfection,"_",episode_variant_summarised)) %>% 
  tibble::column_to_rownames("dose_variant")

## Chisq test with OR ## 
# For each VOC and symptom, this tests the odds of experiencing the
# symptom in that group, compared to the Omicron BA.4/5 group (3 doses)
# This generates the OR, CI, and a table (object name: SxDistroFreqOR)
# Retaining the SxDistroFreqPvalues object used in plotting. 
SxDistroFreqPvalues <- matrix(nrow = nrow(SxDistroFreq), ncol = 8)
SxDistroFreqOR <- list()
for(j in 1:8) {
  for(i in 1:nrow(SxDistroFreq)) {
    start <- (2 * j)+2
    end <- (2* j)+3
    
    # This takes the frequency table, makes a 2x2, and applies a 
    zz<-SxDistroFreq[c(i,nrow(SxDistroFreq)), c(start,end)] %>% 
      as.matrix(.) %>% 
      t(.) %>% 
      epitools::oddsratio.wald(x=.,rev = "both",correction = T,verbose=T)
    
    label = paste0(colnames(zz$x)[2],"_",rownames(zz$x)[2])
    
    SxDistroFreqOR[[label]] <- data.frame(VOC = colnames(zz$x)[2],
                                          referece = "d3+Omicron-BA.4/5",
                                          Symptom =  rownames(zz$x)[2],
                                          n = zz$data[2,2],
                                          total = zz$data[3,2],
                                          prop = zz$p.exposed[2,2],
                                          ref_n = zz$data[2,1],
                                          ref_total = zz$data[3,1],
                                          ref_prop = zz$p.exposed[2,1],
                                          p.value = zz$p.value[2,3],
                                          OR = zz$measure[2,1],
                                          Lower = as.numeric(zz$measure[2,2]),
                                          Upper = as.numeric(zz$measure[2,3]))
    
    cat(zz$p.value[2,3], sep = "\n\n")
    SxDistroFreqPvalues[i,j] <- zz$p.value[2,3]
    rm(zz)
    rm(label)
  }
}


# Format for plot 
colnames(SxDistroFreqPvalues) <- colnames(SxDistroFreq)[ !grepl(colnames(SxDistroFreq), pattern = "Not")] %>%
  .[!grepl(., pattern="functional|summar|tot")]

SxDistroFreqPvalues <- as_tibble(SxDistroFreqPvalues) %>%
  mutate(episode_variant_summarised=SxDistroFreq$episode_variant_summarised) %>%
  mutate(functionalDosesAtTimeOfInfection=SxDistroFreq$functionalDosesAtTimeOfInfection) %>%
  mutate(total = SxDistroFreq$total) %>%
  mutate(grps = paste(episode_variant_summarised, functionalDosesAtTimeOfInfection))

panel.SxDistroPs <- SxDistroFreqPvalues %>%
  filter(total >4) %>%
  select(!total) %>%
  pivot_longer(cols = ! c(grps, episode_variant_summarised, functionalDosesAtTimeOfInfection), names_to = "symptom") %>%
  mutate(value = value) %>%
  mutate(label = case_when(value < 0.001 ~ "***",
                           value < 0.01 ~ "**",
                           value < 0.05 ~ "*",
                           (value == 1 & 
                              str_detect(grps, "Omicron-BA.5") & 
                              str_detect(grps, "3$")) ~ "ref")) %>%
  mutate(Symptom = factor(symptom,
                          levels = rev(c(
                            "Anosmia",
                            "Cough", "Fever",
                            "Coryza", "Fatigue",
                            "Myalgia",
                            "Dyspnoea",
                            "Diarrhoea")))) %>%
  mutate(episode_variant_summarised = str_replace_all(episode_variant_summarised, 
                                                      c("Omicron-BA.1"= "BA.1",
                                                        "Omicron-BA.2"= "BA.2",
                                                        "Omicron-BA.4"= "BA.4",
                                                        "Omicron-BA.4/5"= "BA.4/5",
                                                        "Omicron-BA.5"= "BA.5"))) %>%
  mutate(episode_variant_summarised = factor(episode_variant_summarised,
                                             levels = c("Delta",
                                                        "BA.1",
                                                        "BA.2",
                                                        "BA.4",
                                                        "BA.4/5",
                                                        "BA.5"))) %>%
  ggplot(aes(x=functionalDosesAtTimeOfInfection, y=Symptom, fill=-log10(value), label = label)) +
  geom_tile() +
  geom_text(col="white", size=8/ggplot2::.pt) +
  scale_fill_viridis_c(option="E", name = "-log<sub>10</sub>(<i>P</i>)") +
  facet_grid(.~episode_variant_summarised, space="free_x", scales = "free_x") + 
  theme(        legend.title = ggtext::element_markdown()) +
  xlab("Doses at time of infection")


panel.SxDistro + theme(axis.text.y  =element_text(colour = nhs.labels),
                       panel.grid = element_blank()) +
  theme(legend.position = "bottom") + panel.SxDistroPs




#--------------------------------------------------------------------#
## PLOT: Fig 3 - symptom clustering
#--------------------------------------------------------------------#

# Annotation colours 
ann_cols <- list(Doses = c("2" ="white",
                           "3"= "black"),
                 VOC = c(Delta=colour("muted")(9)[[5]],
                         `Omicron-BA.1`=colour("muted")(9)[[6]], 
                         `Omicron-BA.2`=colour("muted")(9)[[7]],
                         `Omicron-BA.4/5`= colour("muted")(9)[[9]]), #`unknown` = "grey"),
                 Severity = c("grade I"="grey63",
                              "grade II"="grey37",
                              "grade III"="grey10"))


## All VOC 
out <- cluster_mat(SxDistro)
p1 <- as.ggplot(pheatmap(t(out$mat), 
                         border_color = NA,
                         color = c("grey90", color("muted")(9)[[1]]),
                         clustering_distance_rows="binary",
                         clustering_distance_cols="binary",
                         cluster_rows = T, 
                         cluster_cols = T,
                         cellheight = 14,
                         #cellwidth = 1,
                         cutree_rows = c(2),
                         #cutree_cols = 2,
                         legend = F,
                         annotation_legend = F,
                         show_rownames = T,
                         show_colnames = F,
                         annotation_col = out$row_anno,
                         annotation_colors = ann_cols
))

## Delta
out <- cluster_mat(SxDistro,VOC = "Delta")
p2 <- as.ggplot(pheatmap(t(out$mat), #
                         border_color = NA,
                         color = c("grey90", color("muted")(9)[[1]]),
                         clustering_distance_rows="binary",
                         clustering_distance_cols="binary",
                         cluster_rows = T, 
                         cluster_cols = F,
                         cellheight = 14,
                         #cellwidth = 1,
                         cutree_rows = 2,
                         #cutree_cols = 2,
                         legend = F,
                         annotation_legend = F,
                         show_rownames = T,
                         show_colnames = F,
                         annotation_col = out$row_anno,
                         annotation_colors = ann_cols
))

## Ba.1
out <- cluster_mat(SxDistro,VOC = "BA.1")
p3 <- as.ggplot(pheatmap(t(out$mat), 
                         border_color = NA,
                         color = c("grey90", color("muted")(9)[[1]]),
                         clustering_distance_rows="binary",
                         clustering_distance_cols="binary",
                         cluster_rows = T, 
                         cluster_cols = F,
                         cellheight = 14,
                         #cellwidth = 1,
                         cutree_rows = 2,
                         #cutree_cols = 2,
                         legend = F,
                         annotation_legend = F,
                         show_rownames = T,
                         show_colnames = F,
                         annotation_col = out$row_anno,
                         annotation_colors = ann_cols
))

## Ba.2
out <- cluster_mat(SxDistro,VOC = "BA.2")
p4 <- as.ggplot(pheatmap(t(out$mat), 
                         border_color = NA,
                         color = c("grey90", color("muted")(9)[[1]]),
                         clustering_distance_rows="binary",
                         clustering_distance_cols="binary",
                         cluster_rows = T, 
                         cluster_cols = F,
                         cellheight = 14,
                         #cellwidth = 1,
                         cutree_rows = 2,
                         #cutree_cols = 2,
                         legend = F,
                         annotation_legend = F,
                         show_rownames = T,
                         show_colnames = F,
                         annotation_col = out$row_anno,
                         annotation_colors = ann_cols
))

## Ba.4/5
out <- cluster_mat(SxDistro,VOC = "BA.4/5")
p5 <- as.ggplot(pheatmap(t(out$mat), 
                         border_color = NA,
                         color = c("grey90", color("muted")(9)[[1]]),
                         clustering_distance_rows="binary",
                         clustering_distance_cols="binary",
                         cluster_rows = T, 
                         cluster_cols = F,
                         cellheight = 14,
                         #cellwidth = 1,
                         cutree_rows = 2,
                         #cutree_cols = 2,
                         legend = F,
                         annotation_legend = F,
                         show_rownames = T,
                         show_colnames = F,
                         annotation_col = out$row_anno,
                         annotation_colors = ann_cols
))



# A - 9 w x 6 h
# B:E - 6 w x 4 h 

lay <- 
  "AAAAAAAAA
   AAAAAAAAA
   BBBBFFFFF
   BBBBFFFFF
   CCCCFFFFF
   CCCCFFFFF
   DDDDFFFFF
   DDDDFFFFF
   EEEEFFFFF
   EEEEFFFFF"


p<-p1+p2+p3+p4+p5+plot_layout(design = lay,guides = "collect") +
  plot_annotation(tag_levels = c('A')) &
  theme(# axis.text.y = element_text(size=8),
    legend.margin = margin(0,0,0,0, unit="cm"),
    plot.tag = element_text(face="bold", size=18),
    plot.tag.position = c(0.02, 0.95)
  )
print(p)



#--------------------------------------------------------------------#
## PLOT: Fig 4A - CT values over time since symptom onset 
#--------------------------------------------------------------------#


panel.traj <- episodeCts %>%
  
  filter(dose >=2) %>%
  filter(dose <4 ) %>%
  filter(AllPCRs_swab_type == 'dry') %>%
  ## remove any less than 5 cases variants ##
  group_by(dose, episode_variant_summarised) %>%
  filter(nlevels(factor(elig_study_id)) > 4) %>%
  ungroup() %>%
  ##
  ## remove any less than 20 point curves ##
  group_by(dose, episode_variant_summarised) %>%
  filter(n() > 20) %>%
  ungroup() %>%
  ##
  filter(! is.na(episode_variant_summarised)) %>%
  mutate(no_vaccines = paste("infection\nafter", dose, "doses")) %>%
  ## force empty levels preceding Delta to give standard colours ##
  mutate(episode_variant_summarised =
           factor(episode_variant_summarised,
                  levels = c("Ancestral", "D614G", "Alpha", "Beta",
                             "Delta", "Omicron-BA.1", "Omicron-BA.2",
                             "Omicron-BA.4/5"))) %>%
  ggplot(aes(x = days_since_episode_start, y = AllPCRs_tacpath_orf1ab, col = episode_variant_summarised)) +
  facet_grid(~episode_variant_summarised+no_vaccines, drop = T) +
  scale_y_reverse() +
  geom_smooth(show.legend = F, se=F,colour="black") +
  geom_point(shape=20, alpha=0.3,size=2,show.legend=F,colour="black") +
  theme_bw()+
  labs(x = "Time since symptom onset (days)", y = "Ct")

panel.traj


#--------------------------------------------------------------------#
## PLOT: Fig 4B - peak ct values by variant 
#--------------------------------------------------------------------#

peak.Cts <- panel.traj$data %>%
  filter(days_since_episode_start> 1 & days_since_episode_start<6) %>%
  group_by(dose, episode_variant_summarised, elig_study_id) %>%
  arrange(AllPCRs_tacpath_orf1ab, .by_group = T) %>%
  slice_tail(n=1) %>%
  ungroup()

## Plot infection after 2 doses 

panel.peaksCts2 <- peak.Cts %>%
  filter(str_detect(no_vaccines, "2")) %>%
  ## force empty levels preceding Delta to give standard colours ##
  mutate(episode_variant_summarised =
           factor(episode_variant_summarised,
                  levels = c("Ancestral", "D614G", "Alpha", "Beta",
                             "Delta", "Omicron-BA.1", "Omicron-BA.2",
                             "Omicron-BA.4/5"))) %>%
  ggplot(aes(x=episode_variant_summarised, y = AllPCRs_tacpath_orf1ab,
             col=episode_variant_summarised)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(0.2), shape=20, alpha=0.3, size=2) +
  facet_grid(.~no_vaccines) + 
  #scale_y_reverse(limits=c(40,3)) +
  labs(y = "Ct", x = "") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_color_manual(values = manual.pal, drop=F) + 
  guides(colour="none")

# Perform tests
stat.test <- compare_means(AllPCRs_tacpath_orf1ab ~ episode_variant_summarised,
                           data = panel.peaksCts2$data,
                           method = "wilcox.test") %>% 
  mutate(p.lab = ifelse(after_stat(as.numeric(p))<=0.001,paste0("< 0.001"),
                        format(round(after_stat(as.numeric(p)),3), nsmall=2)),
         no_vaccines = unique(panel.peaksCts2$data$no_vaccines))


# Add to plot 
y.pos = c(-16)
panel.peaksCts2<-panel.peaksCts2 +
  stat_pvalue_manual(stat.test,
                     tip.length = 0.01,
                     label = "p.lab",
                     y.position = y.pos,
                     size = 8/ggplot2::.pt)+
  scale_y_reverse()+
  coord_cartesian(ylim = c(40, 3))

## Plot infection after 3 doses 

panel.peaksCts3 <- peak.Cts %>%
  filter(str_detect(no_vaccines, "3")) %>%
  ## force empty levels preceding Delta to give standard colours ##
  mutate(episode_variant_summarised =
           factor(episode_variant_summarised,
                  levels = c("Ancestral", "D614G", "Alpha", "Beta",
                             "Delta", "Omicron-BA.1", "Omicron-BA.2",
                             "Omicron-BA.4/5"))) %>%
  ggplot(aes(x=episode_variant_summarised, y = AllPCRs_tacpath_orf1ab,
             col=episode_variant_summarised)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(0.2), shape=20, alpha=0.3, size=2) +
  facet_grid(.~no_vaccines) + 
  labs(y = "Ct", x = "") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_color_manual(values = manual.pal, drop=F) + 
  guides(colour="none") 

# Get medians
panel.peaksCts3$data %>%
  group_by(episode_variant_summarised,no_vaccines) %>%
  summarise(median_ct = median(AllPCRs_tacpath_orf1ab,na.rm=T))


# Perform tests
stat.test <- compare_means(AllPCRs_tacpath_orf1ab ~ episode_variant_summarised,
                           data = panel.peaksCts3$data,
                           method = "wilcox.test") %>% 
  mutate(p.lab = ifelse(after_stat(as.numeric(p))<=0.001,paste0("< 0.001"),
                        format(round(after_stat(as.numeric(p)),3), nsmall=2)),
         no_vaccines = unique(panel.peaksCts3$data$no_vaccines))


# Add to plot 
y.pos = c(-12,-9,-5)
panel.peaksCts3<-panel.peaksCts3 +
  stat_pvalue_manual(stat.test,
                     tip.length = 0.01,
                     label = "p.lab",
                     y.position = y.pos,
                     size = 8/ggplot2::.pt)+
  scale_y_reverse()+
  coord_cartesian(ylim = c(40, 3))



l <- "aabbb"
panel.peaksCts2+panel.peaksCts3 +plot_layout(design = l)

#--------------------------------------------------------------------#
## PLOT: Supplemental Fig 1 - Peak CT by infection and fever status
#--------------------------------------------------------------------#

SxDistroAll <-
  dosesAtTimeofInfectionAfterD2 %>%
  mutate(cold = case_when(str_detect(parsed_description, "coryza") ~ as.numeric(str_replace_na(cold, "0")) + 1,
                          ! str_detect(parsed_description, "coryza") ~ cold,
                          is.na(parsed_description) ~ cold)) %>%
  mutate(lossSmell = case_when(str_detect(parsed_description, "anosmia") ~ as.numeric(str_replace_na(lossSmell, "0")) + 1,
                               ! str_detect(parsed_description, "anosmia") ~ lossSmell,
                               is.na(parsed_description) ~ lossSmell)) %>%
  mutate(fatigue = case_when(str_detect(parsed_description, "fatigue") ~ as.numeric(str_replace_na(fatigue, "0")) + 1,
                             ! str_detect(parsed_description, "fatigue") ~ fatigue,
                             is.na(parsed_description) ~ fatigue)) %>%
  mutate(fever = case_when(str_detect(parsed_description, "fever") ~ as.numeric(str_replace_na(fever, "0")) + 1,
                           ! str_detect(parsed_description, "fever") ~ fever,
                           is.na(parsed_description) ~ fever))


FeverPlusSxDistroAll <- SxDistroAll %>%
  group_by(functionalDosesAtTimeOfInfection,
           episode_variant_summarised) %>%
  mutate(Plus = coalesce(lossSmell,cold,cough,ache,
                         shrtbrth,fatigue, diarrhoea, oth)) %>%
  mutate(Plus = ifelse(is.na(Plus)==T, "without other", "with other")) %>%
  mutate(Fever= ifelse(is.na(fever)==T, "afebrile", "febrile")) %>%
  mutate(FeverPlus = paste(Fever, Plus, sep=", ")) %>%
  mutate(FeverPlus = ifelse(is.na(Sx_severity)==T, NA, FeverPlus)) %>%
  mutate(FeverPlus = str_replace_na(FeverPlus, "no diary")) %>%
  mutate(FeverPlus = factor(FeverPlus,
                            levels = c("febrile, with other",
                                       "febrile, without other",
                                       "afebrile, with other",
                                       "afebrile, without other",
                                       "no diary")))



peakCts <- peak.Cts %>%
  filter(elig_study_id %in% FeverPlusSxDistroAll$elig_study_id) %>%
  ## Take the categories from feverplus
  left_join(FeverPlusSxDistroAll %>%
              group_by(elig_study_id) %>%
              arrange(episode_start, .by_group = T) %>%
              slice_tail() %>%
              ungroup() %>%
              select(elig_study_id, episode_number, FeverPlus)) %>%
  filter(str_detect(no_vaccines, "3"))


peakFeverPlusCts <- peakCts %>%
  filter(episode_variant_summarised != "Omicron-BA.4/5") %>% 
  filter(FeverPlus %in% c("febrile, with other", "afebrile, with other")) %>%
  ggplot(aes(x=FeverPlus, y=AllPCRs_tacpath_orf1ab,col=FeverPlus)) +
  facet_grid(.~ episode_variant_summarised + no_vaccines) +
  scale_y_reverse(name= "Ct") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape=20, size=2,alpha = 0.3, position=position_jitter(width = 0.3)) +
  scale_colour_devon(reverse=F, discrete=T, drop=F) +
  guides(colour="none") +
  stat_compare_means(method = "wilcox", paired = F, comparisons = list(c(1,2)),label.y = -11) +
  theme(axis.text.x = element_text(angle=90, hjust=1))

# Get medians
peakFeverPlusCts$data %>% 
  group_by(episode_variant_summarised,FeverPlus,no_vaccines) %>% 
  summarise(median = median(AllPCRs_tacpath_orf1ab))

peakFeverPlusCts

#--------------------------------------------------------------------#
## PLOT: Fig 4C - CT values over time since symptom onset by fever status
#--------------------------------------------------------------------#

panel.FeverPlusCtTrajectories <- episodeCts %>%
  filter(elig_study_id %in% FeverPlusSxDistroAll$elig_study_id) %>%
  ## Take the categories from feverplus
  left_join(FeverPlusSxDistroAll %>%
              group_by(elig_study_id) %>%
              arrange(episode_start, .by_group = T) %>%
              slice_tail() %>%
              ungroup() %>%
              select(elig_study_id, episode_number, FeverPlus)) %>%
  mutate(no_vaccines = paste("infection\nafter", dose, "doses")) %>%
  filter(FeverPlus %in% c("febrile, with other", "afebrile, with other")) %>%
  filter(str_detect(episode_variant_summarised, "BA.1|BA.2")) %>%
  filter(str_detect(no_vaccines, "3")) %>%
  mutate(FeverPlus = gsub(", with other",", with other\nsymptoms",FeverPlus)) %>%
  ggplot(aes(x=days_since_episode_start, y=AllPCRs_tacpath_orf1ab,col=FeverPlus)) +
  facet_grid(.~episode_variant_summarised+no_vaccines) +
  scale_y_reverse(limits = c(40,10), name= "Ct") +
  xlim(-2,20) +
  geom_point(shape=20, size=2,alpha = 0.3, position=position_jitter(width = 0.3)) +
  geom_smooth(se = FALSE, span=1) +
  scale_colour_manual(values = c("#788BD9","#2B194C"))+
  labs(x = "Time since symptom onset (days)") +
  theme_bw()+
  theme(legend.position = "bottom")#axis.text.x = element_text(angle=90, hjust=1),

panel.FeverPlusCtTrajectories


#--------------------------------------------------------------------#
## PLOT: Fig 4D - peak ct values by variant and severity 
#--------------------------------------------------------------------#

## i.e. does asking +ve individuals who feel too unwell to conduct their day-to-day lives to stay at home,
## change the RNA loads of individuals at large in the community.

peakSeverityCts.stats <- peakCts %>%
  mutate(Severity_group = case_when(str_detect(Sx_severity, "mod|seve") ~ "II or III",
                                    str_detect(Sx_severity, "mild") ~ "I",
                                    str_detect(Sx_severity, "asym") ~ "asym")) %>%
  mutate(Severity_group=factor(Severity_group)) %>%
  filter(!str_detect(episode_variant_summarised, "BA.5")) %>%
  filter(!str_detect(Severity_group, "asym")) %>%
  group_by(episode_variant_summarised) %>%
  wilcox_test(AllPCRs_tacpath_orf1ab ~ Severity_group) %>%
  add_xy_position(scales = "free") %>%
  mutate(episode_variant_summarised = factor(episode_variant_summarised, levels = c("Omicron BA.1", "Omicron BA.2")))


peakSeverityCts <- peakCts %>%
  mutate(Severity_group = case_when(str_detect(Sx_severity, "mod|seve") ~ "II or III",
                                    str_detect(Sx_severity, "mild") ~ "I",
                                    str_detect(Sx_severity, "asym") ~ "asym")) %>%
  mutate(Severity_group=factor(Severity_group)) %>%
  filter(!str_detect(episode_variant_summarised, "BA.4/5")) %>%
  filter(!str_detect(Severity_group, "asym")) %>%
  ## force empty levels preceding Delta to give standard colours ##
  mutate(episode_variant_summarised =
           factor(episode_variant_summarised,
                  levels = c("Ancestral", "D614G", "Alpha", "Beta",
                             "Delta", "Omicron-BA.1", "Omicron-BA.2",
                             "Omicron-BA.4/5"))) %>%
  ggplot(aes(x=Severity_group, y=AllPCRs_tacpath_orf1ab , col=episode_variant_summarised)) +
  facet_grid(.~episode_variant_summarised+no_vaccines, drop=T) +
  scale_color_manual(values = manual.pal, drop=F) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape=20, size=2,alpha = 0.3, position=position_jitter(width = 0.3)) +
  labs(x = "Severity Grade") +
  guides(color = guide_none())+
  theme_bw()

# Get medians
peakSeverityCts$data %>% 
  group_by(episode_variant_summarised,no_vaccines,Severity_group) %>% 
  summarise(median = median(AllPCRs_tacpath_orf1ab ,na.rm=T))


# BA.1
stat.test.1 <- compare_means(AllPCRs_tacpath_orf1ab  ~ Severity_group,
                             data = peakSeverityCts$data %>% filter(
                               no_vaccines == "infection\nafter 3 doses",
                               episode_variant_summarised == unique(peakSeverityCts$data$episode_variant_summarised)[1]),
                             method = "wilcox.test") %>% 
  mutate(p.lab = ifelse(after_stat(as.numeric(p))<=0.001,paste0("< 0.001"),
                        format(round(after_stat(as.numeric(p)),3), nsmall=2)),
         no_vaccines = unique(peakSeverityCts$data$no_vaccines),
         episode_variant_summarised = unique(peakSeverityCts$data$episode_variant_summarised)[1])
# BA.2
stat.test.2 <- compare_means(AllPCRs_tacpath_orf1ab ~ Severity_group,
                             data = peakSeverityCts$data %>% filter(
                               no_vaccines == "infection\nafter 3 doses",
                               episode_variant_summarised == unique(peakSeverityCts$data$episode_variant_summarised)[2]),
                             method = "wilcox.test") %>% 
  mutate(p.lab = ifelse(after_stat(as.numeric(p))<=0.001,paste0("< 0.001"),
                        format(round(after_stat(as.numeric(p)),3), nsmall=2)),
         no_vaccines = unique(peakSeverityCts$data$no_vaccines),
         episode_variant_summarised = unique(peakSeverityCts$data$episode_variant_summarised)[2])
stat.test <- bind_rows(stat.test.1,stat.test.2)


# Add to plot 
y.pos = c(-10)
peakSeverityCts <- 
  peakSeverityCts +
  stat_pvalue_manual(stat.test,
                     tip.length = 0.01,
                     label = "p.lab",
                     y.position = y.pos,
                     size = 8/ggplot2::.pt)+
  scale_y_reverse(name="Ct")+
  coord_cartesian(ylim = c(40, 3))

peakSeverityCts

#--------------------------------------------------------------------#
## PLOT: Fig 4E - linear fit peak ct values by time since dose 
#--------------------------------------------------------------------#

ct.supp <- ct.supp.dat %>% 
  mutate(episode_variant_summarised =
           factor(episode_variant_summarised,
                  levels = c("Delta", "Omicron-BA.1", "Omicron-BA.2",
                             "Omicron-BA.4/5"))) %>%
  ggplot( aes(x=daysSincedose, y=AllPCRs_tacpath_orf1ab,col=episode_variant_summarised)) +
  facet_grid(~no_vaccines, scales="free_x") +
  scale_y_reverse(limits = c(45,10), name= "Ct") +
  geom_smooth(method="lm", #aes(col=NULL),
              #col="black",
              show.legend = T, se = F) +
  geom_point(shape=20, alpha = 0.3, size=2) +
  scale_color_manual(values = c("#88CCEE","#882255","#44AA99","black" ), drop=T) + 
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  labs(x = "Time since last dose (days)", y="Ct",colour="VOC")+
  theme_bw()+
  theme(legend.position = "bottom")

ct.supp


# Linear fits 
stats <- ct.supp$data %>% 
  group_by(no_vaccines,episode_variant_summarised) %>%
  nest() %>% 
  mutate(model = map(data, ~lm(AllPCRs_tacpath_orf1ab ~ daysSincedose, data = .x) %>% tidy)) %>%
  unnest(model) %>% 
  filter(term == 'daysSincedose') %>% 
  mutate(lower = format(round(estimate-1.96*std.error,2),nsmall=2),
         upper = format(round(estimate+1.96*std.error,2),nsmall=2),
         CI = paste0("(",lower," - ",upper,")"),
         estimate = format(round(estimate,2),nsmall=2),
         p.value = format(round(p.value,3),nsmall=3),
         label = paste0("β ",estimate," ",CI,"\nP = ",p.value))


# A data frame with labels for each facet
f_labels <- data.frame(no_vaccines = stats %>% filter(p.value<0.05) %>% pull(no_vaccines),
                       episode_variant_summarised = stats %>% filter(p.value<0.05) %>% pull(episode_variant_summarised),
                       label = stats %>% filter(p.value<0.05) %>% pull(label),
                       daysSincedose = 100,
                       AllPCRs_tacpath_orf1ab =11)
ct.supp<-ct.supp +
  geom_text(aes(label = label),show.legend = F, data = f_labels,size=3)

ct.supp


#--------------------------------------------------------------------#
## END #### 
#--------------------------------------------------------------------#