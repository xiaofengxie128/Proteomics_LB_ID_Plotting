

################## User Input Section start ####################################
# type in your full inputfile path or leave it blank to choose from the pop 
#out selector,the file path need to be "\\", not single "\"

inputfile_path =""
#e.g., "D:\\R-Programming_temp\\6-6\\6-6_InputFiles.txt"

if (inputfile_path ==""){
  
  
  inputfile_path= file.choose()
}

# 
# # output_path ="D:\\R-Programming_temp\\6-6\\result.csv", set it manually or
# #leave it blank so it generated based on input_file_path
# 
Output_path= ""

if (Output_path ==""){

  Output_path = file.path(dirname(inputfile_path),
                          paste0("Result_",format(Sys.time(),
                                                  "%Y_%m_%d_%H_%M_%OS"),
                                 ".csv"))
}



################## User Input Section  ends#####################################




library(tidyverse)
library(dplyr)
library(ggplot2)

#functions
read_file_name = function(PDrunName){
  file_assign = read_tsv(PDrunName) %>%
    separate(`File Name`, into = c("drive","year","month","project",
                                   "rawfile","raw"), sep = "[\\\\\\.]" )
  return(file_assign)
}

read_group_info = function(PDrunName){
  df_protein_name = read_tsv(str_replace(PDrunName,"_InputFiles.txt",
                                         "_Proteins.txt"))
  group_names = colnames(df_protein_name)[
    which(str_detect(colnames(df_protein_name), "Found in Sample:"))]
  group_names = str_replace_all(
    group_names, c("Found in Sample: .+] "="", ": Sample"=""))
  group_names = as.data.frame(group_names) %>%
    separate(group_names, into = c("FileID", "Cell", "Gradient"), sep = ", ")
  return(group_names)
}

read_protein = function(PDrunName){
  df_protein_name  = read_tsv(str_replace(
    PDrunName,"_InputFiles.txt", "_Proteins.txt")) 
  df_protein_name = df_protein_name %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Master == "IsMasterProtein") %>%
    filter(`Protein FDR Confidence: Combined` == "High") %>%
    select(Accession, str_subset(colnames(
      df_protein_name), "Found in Sample: "))  
  colnames(df_protein_name) = str_replace_all(
    colnames(df_protein_name), ": Sample.+", "")  
  #delete ": Sample.+" for assigned sample name
  colnames(df_protein_name) = str_replace_all(
    colnames(df_protein_name), ": Sample", "")  
  #delete ": Sample" for no assigned sample name
  
  colnames(df_protein_name) = str_replace_all(
    colnames(df_protein_name), "Found in Sample: .+] ", "")  
  #delete "Found in Sample: "
  return(df_protein_name)
}

read_peptides = function(PDrunName){
  df_peptide_name = read_tsv(str_replace(
    PDrunName,"_InputFiles.txt", "_PeptideGroups.txt")) 
  df_peptide_name = df_peptide_name %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Confidence == "High") %>%
    select(str_subset(colnames(df_peptide_name), "Found in Sample: "))    
  colnames(df_peptide_name) = str_replace_all(
    colnames(df_peptide_name), ": Sample.+", "")  
  #delete ": Sample.+" for assigned sample name
  colnames(df_peptide_name) = str_replace_all(
    colnames(df_peptide_name), ": Sample", "")  
  #delete ": Sample.+"  for no assigned sample name
  
  colnames(df_peptide_name) = str_replace_all(
    colnames(df_peptide_name), "Found in Sample: .+] ", "")  
  #delete "Found in Sample: "
  return(df_peptide_name)
}

read_PSMs = function(PDrunName){
  df_PSM = read_tsv(str_replace(PDrunName,"_InputFiles.txt", "_PSMs.txt")) 
  df_PSM = df_PSM %>%
    filter(Contaminant == "FALSE")%>%
    filter(Confidence == "High")
  return(df_PSM)
}

read_MSMS = function(PDrunName){
  df_MSMS = read_tsv(str_replace(
    PDrunName,"_InputFiles.txt", "_MSMSSpectrumInfo.txt"))
  return(df_MSMS)
}

#######


InputFiles = read_file_name(inputfile_path)
proteins = read_protein(inputfile_path)
peptides = read_peptides(inputfile_path)
PSMs = read_PSMs(inputfile_path)
MSMS = read_MSMS(inputfile_path)
group_info = read_group_info(inputfile_path)

FileIDs = unique(InputFiles$`File ID`)
Identification = data.frame()
i = 1
for (file in FileIDs){
  Identification[i,"FileID"] = file
  Identification[i,"rawfile"] = pull(
    InputFiles[which(InputFiles$`File ID`== file), "rawfile"])
  Identification[i,"proteinID_MSMS"] = sum(
    proteins[,which(colnames(proteins)== file)] == "High")
  Identification[i,"proteinID_MBR"] = sum(
    proteins[,which(colnames(proteins)== file)] == "Peak Found")
  Identification[i,"proteinID_Total"] = Identification[
    i,"proteinID_MSMS"] + Identification[i,"proteinID_MBR"]
  Identification[i,"peptideID_MSMS"] = sum(
    peptides[,which(colnames(peptides)== file)] == "High")
  Identification[i,"peptideID_MBR"] = sum(
    peptides[,which(colnames(peptides)== file)] == "Peak Found")
  Identification[i,"peptideID_Total"] = Identification[
    i,"peptideID_MSMS"]+Identification[i,"peptideID_MBR"]
  Identification[i,"PSMs"] = sum(PSMs[,"File ID"] == file)
  Identification[i,"MSMS"] = sum(MSMS[,"File ID"] == file)
  i=i+1
}
Identification = left_join(Identification, group_info, by = "FileID")

write_csv(Identification, Output_path)