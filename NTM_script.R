# install the needed packages before run
library(dplyr)
library(jsonlite)
library(tidyr)
setwd('/Users/duranbao/Desktop/projects/NTM_project/daily_updating/1001_mtb_subspecies_analysis')

##### Section 1. Setting all necessary functions #############
##############  functions zone start #########################

# 1.1 mapping the peptide to their taxons from Unipept API
unip<-function(demo){
  data1<-fromJSON(paste0('https://api.unipept.ugent.be/api/v1/pept2taxa.json?input[]=',demo[1,1],'&extra=true&equate_il=true')) # loop for input[]=
  for (i in 2:nrow(demo)) {
    temp<-fromJSON(paste0('https://api.unipept.ugent.be/api/v1/pept2taxa.json?input[]=',demo[i,1],'&extra=true&equate_il=true'))
    data1<-rbind(data1, temp)
    Sys.sleep(0.5) # set the sleeping time for accessing the Unipept API continuously
  }
  data_taxa<- data1 %>% filter(grepl('Mycobac',taxon_name)) #%>% filter(taxon_rank=='species')
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



##### Section 2. Managing the input data   ###################

# The table contained two main columns: sample/patient, peptide



##### Section 3. Processing the species level analysis  ######
### here I used the 'tuberculosis' as demo
### tb_total was the table contained samples/patient and peptides

comb_total<-data.frame() # initializing the output table

tbname<-unique(tb_total['patient']) # unique patient name list

for (i in 1:nrow(tbname)){ #
  res1 <- tb_total %>%  filter(patient==tbname[i,1])  %>%  filter(taxon_rank == 'species')%>% filter(grepl('Mycobacterium tuberculosis',taxon_name))  
  if(nrow(res1)>1){
    tryCatch({
      res2 <- dotdata(res1) 
      # generate new folder named 'dotdata' for save the peptide-taxon matrix
      write.csv(res2,paste0('dotdata/comb_',tbname[i,1],'_pep_taxon_species.csv')) # generate the transitional peptide-taxon matrix table
      res3<-findjobk(res2,3)
      res3['patient_id']<-tbname[i,1]
      res3['input_pep_num']<-nrow(res2)
      write.csv(res3,paste0('comb/comb_',tbname[i,1],'_3comb_species.csv'))
      comb_total<-rbind(comb_total,res3)
      rm(res3)
    }, error=function(e){})
  } else {print(paste0('need to manually check the data of ',tbname[i,1]))}
}

write.csv(comb_total,'species_analysis.csv')


##### Section 4. Processing the subspecies or lower levels analysis  ######

comb_total<-data.frame() # initializing the output table

tbname<-unique(tb_total['patient']) # unique patient name list

for (i in 1:nrow(tbname)){ 
  res1 <- tb_total %>%filter(patient==tbname[i,1]) %>% filter(grepl('Mycobacterium tuberculosis',taxon_name)%>% filter(taxon_rank %in% c('subspecies','strain','no rank')))
  if(nrow(res1)>1){
    tryCatch({
      res2 <- dotdata(res1) 
      # generate new folder named 'dotdata' for save the peptide-taxon matrix
      write.csv(res2,paste0('dotdata/comb_',tbname[i,1],'_pep_taxon_subspecies.csv')) # generate the transitional peptide-taxon matrix table
      res3<-findjobk(res2,3)
      res3['patient_id']<-tbname[i,1]
      res3['input_pep_num']<-nrow(res2)
      write.csv(res3,paste0('comb/comb_',tbname[i,1],'_3comb_subspecies.csv'))
      comb_total<-rbind(comb_total,res3)
      rm(res3)
    }, error=function(e){})
  } else {print(paste0('need to manually check the data of ',tbname[i,1]))}
}

write.csv(comb_total,'subspecies_analysis.csv')


##### Section 5. Calculating the scores from previous outputs ##############

res<-read.csv('species_analysis.csv', header = TRUE)
res<-na.omit(res)
p_freq<-data.frame(res %>% group_by(patient_id,taxon)%>%summarize(Freq=n()))
total_freq<-data.frame(res %>% group_by(patient_id)%>%summarize(Freq=n()))
s_freq<-merge(p_freq, total_freq, by='patient_id')
s_freq['score']<-round(s_freq['Freq.x']/s_freq['Freq.y'],digits = 4)
s_freq<-s_freq[order(s_freq[,'patient_id'], -s_freq[,'score']), ]

write.csv(s_freq, 'taxon_score.csv')








