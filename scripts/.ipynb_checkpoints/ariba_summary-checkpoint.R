library(dplyr)
#library(ggplot)


### to load ariba report
load_ariba = function(ariba_report_path) {
  d_full = read.table(ariba_report_path,header=TRUE,sep="\t",comment.char = "",check.names = FALSE) %>%
    dplyr::rename(ariba_ref_name = "#ariba_ref_name")
  return(d_full)
}

### to load all ariba reports from folder
load_ariba_from_folder = function(input_folder) {
  samples = c()
  d_full = data.frame()
  samples = list.files(input_folder)
  for (sample in samples) {
    print(sample)
    ariba_report_path = file.path(input_folder,sample,"report.tsv")
    if (file.exists(ariba_report_path)) {
      d = load_ariba(ariba_report_path)
      if (nrow(d) > 0) {
        d$sample_name = sample
        d_full = rbind(d_full,d)
        
      }
    }
  }
  return(list(ariba_data = d_full,sample_names = samples))
}

### load SNPs from ariba reference
load_SNPS_from_ariba_ref = function(ariba_ref_path) {
  d = read.table(ariba_ref_path,sep="\t",header=FALSE)
  SNPs = paste0(d$V1,"::",d$V4)
  return(SNPs)
}


#d_full = load_ariba_from_folder("/Volumes/data/MPV/LAEH/reads/sweden_carrier_and_PJI_20221223_pointmutations_2/")
args_full =commandArgs(trailingOnly = FALSE)

args =commandArgs(trailingOnly = TRUE)

script_dir_arg = args_full[grep('^--file=',args_full)]
script_dir = dirname(substr(script_dir_arg,8,nchar(script_dir_arg)))
ariba_metadata_path = file.path(script_dir,"..","resources","ariba_refs",args[3],"01.filter.check_metadata.tsv")

reference_SNPS = load_SNPS_from_ariba_ref(ariba_metadata_path)

d_list = load_ariba_from_folder(args[1])

d = d_list$ariba_data %>% filter(!var_description == ".",!ref_nt == ctg_nt)
d$Variant = paste0(d$ref_name,"::",d$known_var_change)
d$Variant = factor(d$Variant,levels = reference_SNPS)
d$AMR_class = unlist(lapply(d$var_description, function(x) strsplit(x,"::")[[1]][6]))
d$sample_name = factor(d$sample_name,levels = d_list$sample_names)

var_table = d %>% select(sample_name,Variant) %>% table() %>% as.data.frame.matrix()
var_table$sample_name = rownames(var_table)
var_table = var_table %>% relocate(sample_name)

write.table(var_table, file = args[2], sep="\t",row.names = FALSE,quote = FALSE)
