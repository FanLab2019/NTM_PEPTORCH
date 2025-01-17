# install the needed packages before run
library(dplyr)
library(jsonlite)
library(tidyr)
# set your own working directionary
setwd('C:\\Users\\TulaneProteomicsCore\\Downloads') 

##############################################################
##### Section 1. Setting all necessary functions #############
##############  functions zone start #########################

# 1.1 mapping the peptide to their taxons from Unipept API
unip<-function(demo){
  data1<-fromJSON(paste0('https://api.unipept.ugent.be/api/v1/pept2taxa.json?input[]=',demo[1,1],'&extra=true&equate_il=true')) # loop for input[]=
  for (i in 2:nrow(demo)) {
    temp<-fromJSON(paste0('https://api.unipept.ugent.be/api/v1/pept2taxa.json?input[]=',demo[i,1],'&extra=true&equate_il=true'))
    data1<-rbind(data1, temp)
    Sys.sleep(0.05) # set the sleeping time for accessing the Unipept API continuously
  }
  data_taxa<- data1 %>% filter(grepl('Myco',taxon_name)) #%>% filter(taxon_rank=='species')
  return(data_taxa)
}

# 1.2 converting the long table to the peptide-taxon matrix

refer<-read.csv('partic_wgs_humanhost.csv',header = TRUE) # The data described the human host infected taxons

dotdata<-function(temp){
  referlst <- dplyr::pull(refer, TaxonID)
  temp<-temp %>% filter(taxon_id %in% referlst)
  myco_tax<-temp[c('peptide','taxon_name')]
  myco_tax['count']<-1
  myco_tax<-na.omit(unique(myco_tax))
  myco_new<-tidyr::spread(myco_tax,taxon_name,count) # get summary table for taxon_name
  myco_new[is.na(myco_new)] = 0  # replace NA to '0'
  rownames(myco_new)<-myco_new[,1]
  res2<-subset(myco_new,select=-c(peptide))
  return(res2)
}

# 1.3 generating table contained unique peptide combinations for specific taxons
findjobk<-function(temp,numb){
  # temp<-res2 from dotdata function
  # numb<-3 # the maximum number of combinations we want to summarize from the unique combination, we used '3' in the study
  
  # ignoring 'I' and 'L' in the main pipeline, regarding the I to L from all the peptides.
  temp$peptide<-rownames(temp)
  temp$peptide<-gsub('I','L',temp$peptide)
  temp <- temp[!duplicated(temp),]
  temp <- temp[ , !(names(temp) %in% "peptide")]
  
  resk_outer<-data.frame(matrix(ncol = 3)) # initialize new table for identification
  colnames(resk_outer)<-c('taxon', 'identifier', 'n_marker')
  for (i in 1:numb){
    rlst<-c(1:nrow(temp))
    tempr<-combn(rlst,i)  # generating the combination order number of the peptide set
    resk_total<-data.frame(matrix(ncol = 3)) # initialize new table for identification
    colnames(resk_total)<-c('taxon', 'identifier', 'n_marker')
    for (k in 1:ncol(tempr)){
      tempk<-temp[tempr[,k],]
      # return result small table
      resk<-data.frame(matrix(ncol = 3)) # initialize new table for identification
      colnames(resk)<-c('taxon', 'identifier', 'n_marker')
      # tempk is the checking table
      if (ncol(tempk[colSums(tempk)==i])==1){ 
        res4<-tempk[colSums(tempk)==i]   # res4 can help us to locate the coordinate 
        resk["taxon"]<-colnames(res4) # species or subspecies for identification
        resk["identifier"]<-paste(rownames(res4),collapse = '_') # peptide combination for identification
        resk["n_marker"]<-i # the peptide number in the combination
        resk_total<-rbind(resk_total,resk)
        rm(resk)
      } else { next }
      print(paste0('combination ',nrow(res4),' ',rownames(res4))) 
    }
    resk_total<-na.omit(unique(resk_total))
    resk_outer<-rbind(resk_outer, resk_total)
    if (nrow(resk_total)>0){
      for (r in 1:nrow(resk_total)){
        x<-unique(strsplit(as.character(resk_total["identifier"][r,]),'_')[[1]])
        y<-resk_total['taxon'][r,]
        unirow<-which(rownames(temp) %in% x)
        temp[x,y]<-1 # for future calculation optimizing purpose
        print(paste0('processing ',temp[x,y]))
      }
    }
  }
  return(resk_outer)
}  

##############  functions zone finished ###################
###########################################################


##############################################################
##### Section 2. Managing the input data   ###################
##############################################################

setwd('C:\\Users\\TulaneProteomicsCore\\Downloads')
# The table contained two main columns: sample/patient, peptide
filename<-'Abscessus_3_subspecies_db_PEAKS'

tb_total<-read.csv(paste0(filename,'_DB search psm.csv'))
#tb_res1<-read.csv('human_samples_slim_DB.csv')

# filter the data by FDR  %>% filter(X.10lgP > FDR001)
tb_res1 <- tb_total %>% filter(X.10lgP >35 & Area >0)  %>% dplyr::select(c("Peptide", "Source.File")) %>% drop_na()
tb_res1$Peptide<-gsub("\\+57.02|\\+15.99|\\(|\\)", "", tb_res1$Peptide) # remove some modification strings from the table
tb_res1 <- distinct(tb_res1) %>% drop_na() # remove the duplicate of peptides by modification strings in previous step
names(tb_res1) <- c('peptide', 'patient')
tbname<-unique(tb_res1['patient']) 
n <- nrow(tbname) # identify how many samples from the table

# count the peptide numbers before and after the logp and area criteria
before<-unique(tb_total[c('Peptide','Source.File')])
res_before<-data.frame(before %>% group_by(Source.File)%>% summarise(count=n()))
after<-unique(tb_res1[c('peptide','patient')])
res_after<-data.frame(after %>% group_by(patient)%>% summarise(count=n()))
PEAKS_file_count<-merge(res_before,res_after,by.x='Source.File',by.y = 'patient',all = TRUE)
write.csv(PEAKS_file_count,paste0(filename,'_PEAKS_peptide_count.csv'))
# extra counting finished


tb_total<-data.frame()

for (i in 1:n){
  ### bugs when the input peptides had no unipept records, fixed it later
  ### issue when the feedback results containing long taxonomy list
  tbind<-tb_res1 %>% filter(patient==tbname[i,]) %>% dplyr::select('peptide')
  # peptide number before the PEPTORCH processing
  pepn_before<-nrow(tbind)
  tb_ind<-unip(tbind)
  tb_ind['patient']<-tbname[i,]
  tb_ind['pep_num_PEAKS']<-pepn_before
  tb_total<-rbind(tb_total, tb_ind)
  print(paste0('processing..', tbname[i,]))
  
}
tb_total<-unique(tb_total)
# export the unipept searching results
write.csv(tb_total, paste0(filename,'_large_DB_searching_unipept1.csv'))
# read the report into the analysis again
# filename <- 'tb_paper'
tb_total<-read.csv(paste0(filename,'_large_DB_searching_unipept1.csv'),header = TRUE)


##############################################################
##### Section 3. Processing the species level analysis  ######
##############################################################

comb_total<-data.frame() # initializing the output table

tbname<-unique(tb_total['patient']) # unique patient name list

# To get the main NTM species with high prevalence, we can merge certain species within those species and do manual operations to analyze the details
own_lst_new<-c('Mycobacterium gastri','Mycobacterium asiaticum',
           'Mycobacterium persicum', 'Mycobacterium attenuatum', 'Mycobacterium pseudokansasii',
           'Mycobacterium sp. MAC_080597_8934','Mycobacterium lacus','Mycobacterium simulans',
           'Mycobacterium sp. ACS4054', 'Mycobacterium celatum','Mycobacterium nebraskense','Mycobacterium paraffinicum',
           'Mycobacterium sp. MAC_011194_8550', 'Mycobacterium sp. JS623','Mycobacterium paraense',
           'Mycobacterium sp. MOTT36Y') # filter(!taxon_name %in% )

for (i in 1:nrow(tbname)){ #
  res1 <- tb_total %>%  filter(patient==tbname[i,1])  %>%  filter(taxon_rank == 'species')%>% filter(grepl('Mycobac',taxon_name))  
  # filter out the main NTM species with high prevalence
  res1<-res1 %>% filter(!taxon_name %in% own_lst_new)
  if(nrow(res1)>1){
    tryCatch({
      res2 <- dotdata(res1) 
      # generate new folder named 'dotdata' for save the peptide-taxon matrix
      write.csv(res2,paste0('dotdata/comb_',tbname[i,1],'_pep_taxon_species.csv')) # generate the transitional peptide-taxon matrix table
      # when the peptide number was huge, we can adjust the number after 'res2' to decrease the running time.
      res3<-findjobk(res2,3) 
      res3['patient_id']<-tbname[i,1]
      res3['pep_num_w_UniPept']<-nrow(res2)
      # generate new folder named 'comb' for save the unique peptide or combination information
      write.csv(res3,paste0('comb/comb_',tbname[i,1],'_3comb_species.csv'))      
      comb_total<-rbind(comb_total,res3)
      # remove the res3 to release the RAM and avoid the wrong table output
      rm(res3)
    }, error=function(e){})
  } else {print(paste0('need to manually check the data of ',tbname[i,1]))}
}
comb_total<-na.omit(comb_total)

write.csv(comb_total,paste0(filename,'_p35_DB_species_analysis.csv'))

###############################################################################################
##### Section 4. Calculating the subspecies level's scores from previous outputs ##############
###############################################################################################
res<-read.csv(paste0(filename,'_p35_DB_species_analysis.csv'), header = TRUE)
res<-na.omit(res)
p_freq<-data.frame(res %>% group_by(patient_id,taxon)%>%summarize(Freq=n()))
total_freq<-data.frame(res %>% group_by(patient_id)%>%summarize(Freq=n()))
s_freq<-merge(p_freq, total_freq, by='patient_id')
s_freq['score']<-round(s_freq['Freq.x']/s_freq['Freq.y'],digits = 4)
s_freq<-s_freq[order(s_freq[,'patient_id'], -s_freq[,'score']), ]
res <- res %>% select(c('patient_id','pep_num_wo_UniPept', 'pep_num_w_UniPept')) %>% distinct()
s_freq <- merge(s_freq, res, by.x = 'patient_id', by.y = 'patient_id', all.x=TRUE)
write.csv(s_freq, paste0(filename,'_DB_taxon_score.csv'))



###########################################################################
##### Section 5. Processing the subspecies or lower levels analysis  ######
###########################################################################

setwd('C:\\Users\\TulaneProteomicsCore\\Downloads')
# The table contained two main columns: sample/patient, peptide
filename<-'MTB12345_subspecies' # change to your own data names
tb_total<-read.csv(paste0(filename,'_DB search psm.csv'))

tb_res1 <- tb_total %>% filter(X.10lgP >35 & Area >0)  %>% dplyr::select(c("Peptide", "Source.File")) %>% drop_na()
tb_res1$Peptide<-gsub("\\+57.02|\\+15.99|\\(|\\)", "", tb_res1$Peptide) # remove some modification strings from the table
tb_res1 <- distinct(tb_res1) %>% drop_na() # remove the duplicate of peptides by modification strings in previous step
names(tb_res1) <- c('peptide', 'patient')
tbname<-unique(tb_res1['patient']) 
n <- nrow(tbname) # identify how many samples from the table

# count the peptide numbers before and after the logp and area criteria
before<-unique(tb_total[c('Peptide','Source.File')])
res_before<-data.frame(before %>% group_by(Source.File)%>% summarise(count=n()))
after<-unique(tb_res1[c('peptide','patient')])
res_after<-data.frame(after %>% group_by(patient)%>% summarise(count=n()))
PEAKS_file_count<-merge(res_before,res_after,by.x='Source.File',by.y = 'patient',all = TRUE)
write.csv(PEAKS_file_count,paste0(filename,'_PEAKS_peptide_count.csv'))
# extra counting finished

tb_total<-data.frame()

for (i in 1:n){
  tbind<-tb_res1 %>% filter(patient==tbname[i,]) %>% dplyr::select('peptide')
  # peptide number before the PEPTORCH processing
  pepn_before<-nrow(tbind)
  tb_ind<-unip(tbind)
  tb_ind['patient']<-tbname[i,]
  tb_ind['pep_num_PEAKS']<-pepn_before
  tb_total<-rbind(tb_total, tb_ind)
  print(paste0('processing..', tbname[i,]))
  
}
tb_total<-unique(tb_total)
# export the unipept searching results
write.csv(tb_total, paste0(filename,'_large_DB_searching_unipept1.csv'))

tb_total<-read.csv(paste0(filename,'_large_DB_searching_unipept1.csv'),header = TRUE)
comb_total_sub<-data.frame() # initializing the output table
tbname<-unique(tb_total['patient']) # unique patient name list
for (i in 1:nrow(tbname)){ 
  # certain main species name: "Mycobacteroides abscessus", "Mycobacterium tuberculosis", 'Mycobacterium avium','Mycobacterium intracellulare'
  # replace the main species names into the next 'grepl' function to identify specific subspecies within the identified species in the former step
  res1 <- tb_total %>%filter(patient==tbname[i,1]) %>% filter(grepl("Mycobacterium tuberculosis",taxon_name))%>% filter(taxon_rank %in% c('subspecies','strain','no rank'))
  if(nrow(res1)>1){
    tryCatch({
      res2 <- dotdata(res1) 
      # generate new folder named 'dotdata' for save the peptide-taxon matrix
      write.csv(res2,paste0('dotdata/comb_sub_',tbname[i,1],'_pep_taxon_subspecies.csv')) # generate the transitional peptide-taxon matrix table
      res3<-findjobk(res2,3)
      res3['patient_id']<-tbname[i,1]
      res3['pep_num_w_UniPept']<-nrow(res2)
      res3['pep_num_wo_UniPept']<-res1['pep_num_PEAKS'][1,1]
      write.csv(res3,paste0('comb/comb_sub_',tbname[i,1],'_3comb_subspecies.csv'))
      comb_total_sub<-rbind(comb_total_sub,res3)
      rm(res3)
    }, error=function(e){})
  } else {print(paste0('need to manually check the data of ',tbname[i,1]))}
}

comb_total_sub<-na.omit(comb_total_sub)
write.csv(comb_total_sub,paste0(filename,'_DB_subspecies_analysis.csv'))

###############################################################################################
##### Section 6. Calculating the subspecies level's scores from previous outputs ##############
###############################################################################################

res<-read.csv(paste0(filename,'_DB_subspecies_analysis.csv'), header = TRUE)
res<-na.omit(res)
p_freq<-data.frame(res %>% group_by(patient_id,taxon)%>%summarize(Freq=n()))
total_freq<-data.frame(res %>% group_by(patient_id)%>%summarize(Freq=n()))
s_freq<-merge(p_freq, total_freq, by='patient_id')
s_freq['score']<-round(s_freq['Freq.x']/s_freq['Freq.y'],digits = 4)
s_freq<-s_freq[order(s_freq[,'patient_id'], -s_freq[,'score']), ]
res <- res %>% select(c('patient_id','pep_num_wo_UniPept', 'pep_num_w_UniPept')) %>% distinct()
s_freq <- merge(s_freq, res, by.x = 'patient_id', by.y = 'patient_id', all.x=TRUE)
write.csv(s_freq, paste0(filename,'_DB_subspecies_score.csv'))

