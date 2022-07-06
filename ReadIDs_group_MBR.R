

################## User Input Section start ####################################

setwd("X:/Share/Xiaofeng/6-6")   #where you save the file, use "/" instead of "\"




################## User Input Section  ends#####################################

library(tidyverse)
library(dplyr)
library(ggplot2)

#functions
read_file_name = function(PDrunName){
  file_assign = read_tsv(paste0(PDrunName, "_InputFiles.txt")) %>%
    separate(`File Name`, into = c("drive","year","month","project", "rawfile","raw"), sep = "[\\\\\\.]" )
  return(file_assign)
}

read_group_info = function(PDrunName){
  df_protein_name = read_tsv(paste0(PDrunName, "_Proteins.txt")) 
  group_names = colnames(df_protein_name)[which(str_detect(colnames(df_protein_name), "Found in Sample:"))]
  group_names = str_replace_all(group_names, c("Found in Sample: .+] "="", ": Sample"=""))
  group_names = as.data.frame(group_names) %>%
    separate(group_names, into = c("FileID", "Cell", "Gradient"), sep = ", ")
  return(group_names)
}

read_protein = function(PDrunName){
  df_protein_name = read_tsv(paste0(PDrunName, "_Proteins.txt")) 
  df_protein_name = df_protein_name %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Master == "IsMasterProtein") %>%
    filter(`Protein FDR Confidence: Combined` == "High") %>%
    select(Accession, str_subset(colnames(df_protein_name), "Found in Sample: "))  
  colnames(df_protein_name) = str_replace_all(colnames(df_protein_name), ": Sample.+", "")  #delete ": Sample.+"
  colnames(df_protein_name) = str_replace_all(colnames(df_protein_name), "Found in Sample: .+] ", "")  #delete "Found in Sample: "
  return(df_protein_name)
}

read_peptides = function(PDrunName){
  df_peptide_name = read_tsv(paste0(PDrunName, "_PeptideGroups.txt"))
  df_peptide_name = df_peptide_name %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Confidence == "High") %>%
    select(str_subset(colnames(df_peptide_name), "Found in Sample: "))    
  colnames(df_peptide_name) = str_replace_all(colnames(df_peptide_name), ": Sample.+", "")  #delete ": Sample.+"
  colnames(df_peptide_name) = str_replace_all(colnames(df_peptide_name), "Found in Sample: .+] ", "")  #delete "Found in Sample: "
  return(df_peptide_name)
}

read_PSMs = function(PDrunName){
  df_PSM = read_tsv(paste0(PDrunName, "_PSMs.txt"))
  df_PSM = df_PSM %>%
    filter(Contaminant == "FALSE")%>%
    filter(Confidence == "High")
  return(df_PSM)
}

read_MSMS = function(PDrunName){
  df_MSMS = read_tsv(paste0(PDrunName, "_MSMSSpectrumInfo.txt"))
  return(df_MSMS)
}

#######


PDrunName = "6-6"
InputFiles = read_file_name(PDrunName)
proteins = read_protein(PDrunName)
peptides = read_peptides(PDrunName)
PSMs = read_PSMs(PDrunName)
MSMS = read_MSMS(PDrunName)
group_info = read_group_info(PDrunName)

FileIDs = unique(InputFiles$`File ID`)
Identification = data.frame()
i = 1
for (file in FileIDs){
  Identification[i,"FileID"] = file
  Identification[i,"rawfile"] = pull(InputFiles[which(InputFiles$`File ID`== file), "rawfile"])
  Identification[i,"proteinID_MSMS"] = sum(proteins[,which(colnames(proteins)== file)] == "High")
  Identification[i,"proteinID_MBR"] = sum(proteins[,which(colnames(proteins)== file)] == "Peak Found")
  Identification[i,"proteinID_Total"] = Identification[i,"proteinID_MSMS"] + Identification[i,"proteinID_MBR"]
  Identification[i,"peptideID_MSMS"] = sum(peptides[,which(colnames(peptides)== file)] == "High")
  Identification[i,"peptideID_MBR"] = sum(peptides[,which(colnames(peptides)== file)] == "Peak Found")
  Identification[i,"peptideID_Total"] = Identification[i,"peptideID_MSMS"]+Identification[i,"peptideID_MBR"]
  Identification[i,"PSMs"] = sum(PSMs[,"File ID"] == file)
  Identification[i,"MSMS"] = sum(MSMS[,"File ID"] == file)
  i=i+1
}
Identification = left_join(Identification, group_info, by = "FileID")

write_csv(Identification, "madi_id.csv")