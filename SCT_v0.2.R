###
### 
###
library(bcp)
library(ChIPpeakAnno)
library(plyr)
library(ComplexHeatmap)
library(arules)
library(stats)

setwd("C:/Users/guida/Desktop/stratification_tcga/drugdb")
tabDrugs<-read.delim(file="interactions.tsv")
colnames(tabDrugs)[1]<-"genes"

print("Step1: Read the parameters")

##First Step: Read the configuration files always located in the path of software 
setwd("C:/Users/guida/Desktop/stratification_tcga")
conf_file<-read.delim(file="configuration_file2.txt",stringsAsFactors=F,header=T)

#1) Read TF file: read the reference files with the list of file of TFs
reference_TF<-read.table(paste("./","sampleinfo_db.txt",sep=""),sep="\t")

fileTF<-reference_TF[reference_TF[,1]==conf_file[1,2],]
fileTF_parse<-paste(fileTF[,1],fileTF[,2],sep="")

# 2) GMID consider only the start end coordinates
input_TF<-read.table(paste("./db/",fileTF_parse,sep=""),sep="\t",stringsAsFactors=F)
input_TF_fr<-input_TF[,c(1:3)]

# 3) read the gene expression data 
#row genes #column patients
file_ge<-conf_file[conf_file[,"PARAMETERS"]=="path_GE",2]
input_GE<-read.table(file=file_ge,sep="\t",stringsAsFactors=FALSE,header=T,quote=NULL,fill=F) #read empty values with 0 
input_GE[is.na(input_GE)]<-0

# extract the gene expression data for the current TF
input_GE_tf<-input_GE[input_GE[,"Hugo_Symbol"]== conf_file[1,2],]
ncol_input_GE<-ncol(input_GE)

# 4) read the copy-number alteration data
file_CNV<-conf_file[conf_file[,"PARAMETERS"]=="path_CNV",2]
input_CNV<-read.table(file=file_CNV,sep="\t",stringsAsFactors=F,header=T,fill=T)

input_CNV_tf<-input_CNV[input_CNV[,1]==as.character(fileTF[1,1]),]


# 5) read the mutation data
file_MUTATION<-conf_file[conf_file[,"PARAMETERS"]=="path_MUTATION",2]
input_MUTATION<-read.table(file=file_MUTATION,sep="\t",stringsAsFactors=F,header=T,fill=T)

# read the mutation data of current tf
input_MUTATION_tf<-input_MUTATION[input_MUTATION[,1]==as.character(fileTF[1,1]),]

# 6) read the methylation data
file_METH<-conf_file[conf_file[,"PARAMETERS"]=="path_METHYLATION",2]
input_METH<-read.table(file=file_METH,sep="\t",stringsAsFactors=F,header=T,fill=T)

#read the methylation data of the current TF
input_METH_tf<-input_METH[input_METH[,1]==as.character(fileTF[1,1]),]

# 7) read the clinical data
file_CLINICAL<-conf_file[conf_file[,"PARAMETERS"]=="path_CLINICAL",2]
input_CLINICAL<-read.table(file=file_CLINICAL,sep="\t",stringsAsFactors=F,header=T,fill=T)

# 8) Read parameter mad_behavior 
mad_behavior<-conf_file[conf_file[,"PARAMETERS"]=="mad_thr_estimation",2]
mad_thr<-conf_file[conf_file[,"PARAMETERS"]=="mad_thr",2]

# 9) Read parameter perc_patient 
parameter_perc_patient<-as.numeric(conf_file[conf_file[,"PARAMETERS"]=="perc_patients",2])

# 10) Read the type of mutations to used 
parameter_type_mutation<-as.character(conf_file[conf_file[,"PARAMETERS"]=="type_mutation",2])

# 11) Read the path of the output directory
output_dir<-as.character(conf_file[conf_file[,"PARAMETERS"]=="output_dir",2])

# 12) Read the string of the output file
output_file<-as.character(conf_file[conf_file[,"PARAMETERS"]=="output_file",2])

# 13) Read the distance for the annotation
distance_max<-as.numeric(conf_file[conf_file[,"PARAMETERS"]=="distance_max",2])

# 14) categorization parameter
parameter_discr<-as.character(conf_file[conf_file[,"PARAMETERS"]=="param_discr",2])

#remove the chr string from the columns of chromosome, only chromosome wihtout chr are allowed
input_TF_fr[,1]<-gsub(input_TF_fr[,1],pattern="chr",replacement="")
colnames(input_TF_fr)<-c("seqnames","start","end")
peaks <- toGRanges(input_TF_fr, format="broadPeak")

print("Step2: Association of epigenetic components with the genes")

#load annotation data format: chr, start, end, gene symbol
annotation_file<-as.character(conf_file[conf_file[,"PARAMETERS"]=="annotation_library",2])
input_annotation<-read.table(annotation_file,sep="\t",header=T,stringsAsFactors=F)
colnames(input_annotation)<-c("seqnames","start","end","gene_name")
input_annotation[,1]<-gsub(input_annotation[,1],pattern="chr",replacement="")

rownames(input_annotation)<-paste(input_annotation[,4],seq(1,nrow(input_annotation)),sep="_")

annotationGR<-makeGRangesFromDataFrame(input_annotation)

anno <- annotatePeakInBatch(peaks, AnnotationData=annotationGR,output="nearestLocation",
                            maxgap=distance_max, ignore.strand=TRUE) #time consuming with huge number of peaks

annotation_results_current_dist<-as.data.frame(anno)

subsetAnnoDF<-annotation_results_current_dist[,c("feature","distancetoFeature")]
subsetAnnoDF[,1]<-sapply(strsplit((subsetAnnoDF$feature),split="_"),"[[",1)

#first column genes #second column distance
subsetAnnoDF_unique<-aggregate(distancetoFeature ~ feature, data=subsetAnnoDF, FUN=mean)

####
####extract the gene expression values of the annotated genes: obviously only some genes are available.
####

input_GE_selectedgenes<-input_GE[input_GE[,1]%in%subsetAnnoDF_unique[,1],]

#remove possible NA values
input_GE_selectedgenes[is.na(input_GE_selectedgenes)]<-0

###
### The genes derived in the first step are first selected by expression value. Filter out the genes with low expression
###

  #estimate variance of selected genes
  sdevestimation<-apply(input_GE_selectedgenes[,3:ncol(input_GE_selectedgenes)],1,sd)
  #estimate mean of selected genes
  mean_genes<-apply(input_GE_selectedgenes[,3:ncol(input_GE_selectedgenes)],1,mean)

  #create a new.data.frame with the coefficient of variation
  input_GE_selectedgenes<-data.frame(input_GE_selectedgenes,mean=mean_genes,mad =apply(input_GE_selectedgenes[,3:ncol(input_GE_selectedgenes)],1,mad))
  #remove na
  input_GE_selectedgenes_filter<-input_GE_selectedgenes[which(!is.na(input_GE_selectedgenes$mad)),]
  #sort by cv
  input_GE_selectedgenes2<-input_GE_selectedgenes_filter[order(input_GE_selectedgenes_filter$mad,decreasing=TRUE),]
  #select all genes with a cv of variation greater than 0 

if(isTRUE(mad_behavior=="AUTO")){
  
#find the flex point of the distribution of CV
bcp_results<-bcp(input_GE_selectedgenes2$mad,mcmc = 1000)

#find where is the minimum break points (this correspond with the flex point)
min_bp<-min(bcp_results$blocks)
#find the corresponding values of c.v. in min_bp
min_cv_automatic_thresholds<-input_GE_selectedgenes2$mad[min_bp]

#i selected all genes with a thresholds less than or greater to the min cv defined after break points changes identification
input_GE_selectedgenes2_mad_select_afterBCP<-input_GE_selectedgenes2[which(input_GE_selectedgenes2$mad>=min_cv_automatic_thresholds),]

} else {
  
  input_GE_selectedgenes2_mad_select_afterBCP<-input_GE_selectedgenes2[which(input_GE_selectedgenes2$mad>=as.numeric(mad_thr)),]
  
}


input_GE_afterCVBCP_filtering<-input_GE_selectedgenes2_mad_select_afterBCP

check_ge_for_patients<-input_GE_afterCVBCP_filtering #old code and then deprecated

###### Finish filtering

###
### Select genes filtered before from CNV, METHYLATION, MUTATION DATA
###

selected_genes_after_GE_filtering<-check_ge_for_patients[,c("Hugo_Symbol","Entrez_Gene_Id")]

input_CNV_selected<-input_CNV[which(input_CNV[,1] %in% selected_genes_after_GE_filtering[,1]), ]
input_METH_selected<-input_METH[which(input_METH[,1] %in% selected_genes_after_GE_filtering[,1]), ]

sub_input_MUTATION<-input_MUTATION[which(input_MUTATION[,1] %in% selected_genes_after_GE_filtering[,1]),]

#for mutation a further level of filtering is require, in fact we want select only particular mutations defined in the conf files 
defined_mutation<-paste(unlist(strsplit(parameter_type_mutation,split=";")),collapse="|")
index_defined_mutation<-grep(x=sub_input_MUTATION$Variant_Classification,pattern=defined_mutation)
input_MUTATION_selected<-sub_input_MUTATION[index_defined_mutation,]

print("Step3: Select the common ID patients between all patients")

##
## Control
##

#
#Find the common patients between the experiment
#

#define a function to unify the names of patients. The function it is based on TCGA format to define the IDpatients.

parseID<-function(x){
  
 x2<- gsub(x,pattern="-",replacement=".") #replace all - with . 
 x3<- substr(x2,0,15)
 
return(x3)
 
}

pts_exp_ge<-parseID(colnames(check_ge_for_patients))
pts_exp_cnv<-parseID(colnames(input_CNV_selected[,-c(1,2)]))
pts_meth_cnv<-parseID(colnames(input_METH_selected))
pts_mut<-parseID(input_MUTATION_selected$Tumor_Sample_Barcode) #tumor sample barcode it is the mandatory colnames

list_pts_experiments<-list(ge=pts_exp_ge,cnv=pts_exp_cnv,meth=pts_meth_cnv,mut=pts_mut)

tts_for_reducing<-list()

for(tts in 1:length(list_pts_experiments)){

    tts_current<-list_pts_experiments[[tts]]
      
    if(length(tts_current)>1){
      
      tts_for_reducing[[tts]]<- tts_current
        
    } else {next}
    
}

common_patient_GE_CNV_METH_MUT<-Reduce(intersect, tts_for_reducing)


###
### Start the analysis
###

print("Step4: Run the analysis")

ALL_samples_UNIQUE<-common_patient_GE_CNV_METH_MUT

for(asu in 1:length(ALL_samples_UNIQUE)){
  
  se_patient_selection<-ALL_samples_UNIQUE[asu]
  
  GE_current_patient<-check_ge_for_patients[,c("Hugo_Symbol","Entrez_Gene_Id",se_patient_selection)]
  CNV_current_patient<-input_CNV_selected[,c("Hugo_Symbol","Entrez_Gene_Id",se_patient_selection)]
  METH_current_patient<-input_METH_selected[,c("Hugo_Symbol","Entrez_Gene_Id",se_patient_selection)]
  MUT_current_patient<-input_MUTATION_selected[which(input_MUTATION_selected$Tumor_Sample_Barcode==se_patient_selection),]
  
  list_DF<-c("GE_current_patient","CNV_current_patient","METH_current_patient","MUT_current_patient")
  
  #check the presence of empty data.frame
  res_nrow<-NULL
  
  for(dfe in list_DF){
  
  current_df<-get(dfe)
  
  res_nrow<-c(res_nrow,dim(current_df)[1])
  
  }
  
  #check2: control which samples does not have genes
  nogenesindataset<-which(res_nrow==0)
  
  if(length(nogenesindataset)==0){
  
    list_DF_clean<-list_DF
    
    DF_notpresent<-"AlldataAvailable"
    
  } else{
    
  list_DF_clean<-list_DF[-nogenesindataset]
  
  DF_notpresent<-list_DF[nogenesindataset]
  
  }
  
  list_df_patients<-list()
  
  for(ldfc in 1:length(list_DF_clean)){
  
  list_df_patients[[ldfc]]<-get(list_DF_clean[ldfc])
    
  }
  
  ###
  ### merge the different experiments for the same patient
  ###
  merge_experiment_patient_df<-Reduce(function(...) merge(...,by="Hugo_Symbol",all=T),list_df_patients)
  
  ##
  ## Start to parse the object merge_experiment_patient_df
  ##
  
  #find the index of columns with the string of patients: Select Data of Copy Number. -> Mutation does not have a column with the name of sample
  index_mepd<-grep(colnames(merge_experiment_patient_df),pattern=se_patient_selection)
  
  #create a new data.frame with hugo symbol and entrez and the values of experiments
  dfPatientForAnalysis<-cbind(merge_experiment_patient_df[1:2],merge_experiment_patient_df[index_mepd],DF_notpresent=rep(0,nrow(merge_experiment_patient_df)))

  colnames(dfPatientForAnalysis)<-c("Hugo_Symbol","Entrez",list_DF_clean[c(1:3)])
  
  #Change the last column, the experiment without data is always in the last column, see previously line of code
  colnames(dfPatientForAnalysis)[ncol(dfPatientForAnalysis)]<-DF_notpresent
  
  print("Step5: Categorization of the genes's properties for the current patient")
  
  #remove NA values 
  dfPatientForAnalysis[is.na(dfPatientForAnalysis)]<-0
  
  rownames(dfPatientForAnalysis)<-paste(dfPatientForAnalysis[,1],seq(1:nrow(dfPatientForAnalysis)))
  
  #### Parse data of transcriptional factors 

  if(length(intersect(colnames(input_GE_tf),se_patient_selection)==1)){
    
    ge_TF_current_patient<-as.numeric(input_GE_tf[se_patient_selection]) 
    
  } else {ge_TF_current_patient<-0}
  
  
  if(length(intersect(colnames(input_CNV_tf),se_patient_selection)==1)){
    
    cnv_TF_current_patient<-as.numeric(input_CNV_tf[se_patient_selection])
    
  } else {cnv_TF_current_patient<-0}
  
  
  if(length(intersect(colnames(input_METH_tf),se_patient_selection)==1)){
    
    meth_TF_current_patient<-as.numeric(input_METH_tf[se_patient_selection])
    
  } else {meth_TF_current_patient <-0 }
  
  
  if(length(intersect(input_MUTATION_tf$Tumor_Sample_Barcode,se_patient_selection)==1)){
    
    mutation_TF_current_patient<-input_MUTATION_tf[input_MUTATION_tf$Tumor_Sample_Barcode==se_patient_selection,]
    index_defined_mutation<-grep(x=mutation_TF_current_patient$Variant_Classification,pattern=defined_mutation)
    mutation_TF_current_patient_variant<-1
    
  } else {mutation_TF_current_patient_variant <-0 }
  
  #extract the other data of experiment
  TF_ge_rep<-rep(as.numeric(ge_TF_current_patient),length(dfPatientForAnalysis[,"GE_current_patient"]))
  TF_CNV_rep<-rep(as.numeric(cnv_TF_current_patient),length(dfPatientForAnalysis[,"CNV_current_patient"]))
  TF_METH_rep<-rep(as.numeric(meth_TF_current_patient),length(dfPatientForAnalysis[,"METH_current_patient"]))

  
  ###
  ### Step1: Identify if the gene expression, cnv, methylation of TFs predicted expression of genes of interest
  ###
  parameter_discr_unlist<-as.numeric(unlist(strsplit(parameter_discr,split=";")))
  ge_d<-parameter_discr_unlist[1]
  cnv_d<-parameter_discr_unlist[2]
  meth_d<-parameter_discr_unlist[3]
  
  ##
  ## Step 1.1: categorize the genes associated with the expression or not of the tf using a fold-change values cut-off
  ##
  FC_GE_TF<-abs(dfPatientForAnalysis[,"GE_current_patient"]-ge_TF_current_patient)
  FC_GE_TF_categorization<-rep(0,length(FC_GE_TF))
  
  #where the fold-change is less than to a thresholds then those genes are related with the transcriptional factors
  FC_GE_TF_categorization[which(FC_GE_TF<=ge_d)]<-1 #absolute values
  FC_GE_TF_categorization[which(FC_GE_TF>ge_d)]<-0 
  
  ###
  ### Discreterization of ghe gene expression values according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC151169/
  ### works for that median-centered at the levels of genes and z-score
  ###
  
  genes_overexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
  names(genes_overexpressed)<-as.character(dfPatientForAnalysis[,1])#give the names to the vector
  rows_to_selected<-round((length(dfPatientForAnalysis[,1])*5)/100)
  #select the over expressed genes in the current_patient
  genes.overexpressed.strings<-dfPatientForAnalysis[order(dfPatientForAnalysis[,"GE_current_patient"],decreasing=T),][1:rows_to_selected,1]
  genes_overexpressed[which(names(genes_overexpressed) %in% genes.overexpressed.strings)]<-1
  
  
  genes_underexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
  names(genes_underexpressed)<-as.character(dfPatientForAnalysis[,1])#give the names to the vector
  rows_to_selected<-round((length(dfPatientForAnalysis[,1])*5)/100)
  #select the under expressed genes in the current_patient
  genes.underexpressed.strings<-dfPatientForAnalysis[order(dfPatientForAnalysis[,"GE_current_patient"]),][1:rows_to_selected,1]
  genes_underexpressed[which(names(genes_underexpressed) %in% genes.underexpressed.strings)]<-1
  
  
  ##
  ## Step 1.2: categorize the copy-number alteration
  ##
  
  #here i have a different problem, if one gene has a copy number greater than 
  # 
  # #first: test in which genes the copy-number is alterated respect with a copy-number thr
   CNV_TF_categorization<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
  # 
  # tcnv<-which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_d |  dfPatientForAnalysis[,"CNV_current_patient"]<= -cnv_d)
  # 
  # if(length(tcnv)==0){
  #   
  # tcnv<-CNV_TF_categorization #all genes are not alterated by the copy number variation 
  # 
  # } else {
  #   
  # #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
  # CNV_TF_categorization[which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_d)]<-1 #if true that the copy number is greater than of a threshold there is amplification
  # CNV_TF_categorization[which(dfPatientForAnalysis[,"CNV_current_patient"] <= -cnv_d)]<-1 #if true that the copy number is less than of a threshold there is depletion
  # 
  # }
  # 
  ###
  ### Step3: Find Gain-amplication
  ###

  CNV_TF_gain<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
  
  tcnv<-which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_d |  dfPatientForAnalysis[,"CNV_current_patient"]<= -cnv_d)
  
  if(length(tcnv)==0){
    
    tcnv<-CNV_TF_gain #all genes are not alterated by the copy number variation 
    
  } else {
    
    #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
    CNV_TF_gain[which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_d)]<-1
    CNV_TF_gain[which(dfPatientForAnalysis[,"CNV_current_patient"] <= -cnv_d)]<-0
    
  }
  
  ###
  ### Step3: Find Depletion
  ###
  
  CNV_TF_depletion<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
  
  tcnv<-which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_d |  dfPatientForAnalysis[,"CNV_current_patient"]<= -cnv_d)
  
  if(length(tcnv)==0){
    
    tcnv<-CNV_TF_depletion #all genes are not alterated by the copy number variation 
    
  } else {
    
    #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
    CNV_TF_depletion[which(dfPatientForAnalysis[,"CNV_current_patient"] <= -cnv_d)]<-1
    CNV_TF_depletion[which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_d)]<-0
    
  }
  
  #now i can test in which case the copy-number is greater or less than the copy number of TFs
  #0 The CNV of the genes is greater than the copy-number variation of TF
  #1 The CNV of the genes is less than the copy-number variation of the TF, interest this case because allow to identify TFs that are alterated
  #and that regulates the genes
  
if(cnv_TF_current_patient>=cnv_d | cnv_TF_current_patient<=-cnv_d){
    #create a new object to fill
    CNV_TF_categorization_TF<-CNV_TF_categorization
    CNV_TF_categorization_TF[which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_TF_current_patient)]<-0
    CNV_TF_categorization_TF[which(dfPatientForAnalysis[,"CNV_current_patient"] < cnv_TF_current_patient)]<-1
    
  } else {
    
    CNV_TF_categorization_TF<-CNV_TF_categorization
    
    }
  
# ##
# ## Step 1.3: categorize the methylation
# ##
# 
   METH_TF_categorization<-rep(0,length(dfPatientForAnalysis[,"METH_current_patient"]))
# 
# tmeth<-which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d |  dfPatientForAnalysis[,"METH_current_patient"]< meth_d)
# 
# if(length(tmeth)==0){
#   
#   tmeth<-METH_TF_categorization #all genes are not alterated by the copy number variation 
#   
# } else {
#   
#   #if is false the condition tmeth ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
#   METH_TF_categorization[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-1 #check hypermethylation
#   METH_TF_categorization[which(dfPatientForAnalysis[,"METH_current_patient"] <= -meth_d)]<-1 #check hypomethylation
#   
# }

###
### Step 1.4: categorize the hyper-hypomethylation
###

METH_TF_hyper<-rep(0,length(dfPatientForAnalysis[,"METH_current_patient"]))

tmeth<-which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d |  dfPatientForAnalysis[,"METH_current_patient"]< meth_d)

if(length(tmeth)==0){
  
  tmeth<-METH_TF_hyper #all genes are not alterated by the copy number variation 
  
} else {
  
  METH_TF_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-1
  METH_TF_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] < meth_d)]<-0
  
}

METH_TF_hypo<-rep(0,length(dfPatientForAnalysis[,"METH_current_patient"]))

tmeth<-which(dfPatientForAnalysis[,"METH_current_patient"] >=meth_d |dfPatientForAnalysis[,"METH_current_patient"] < meth_d)

if(length(tmeth)==0){
  
  tmeth<-METH_TF_hyper #all genes are not alterated by the copy number variation 
  
} else {
  
  METH_TF_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] <  meth_d)]<-1
  METH_TF_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-0
  
}


#now i can test in which case the methylation is greater or less than the methylation of TFs
#0 The METH of the genes is greater than the methylation variation of TF
#1 The METH of the genes is less than the methylation variation of the TF, interest this case because allow to identify TFs that are alterated
#and that regulates the genes

if(meth_TF_current_patient>=meth_TF_current_patient | meth_TF_current_patient< meth_TF_current_patient){
  
  METH_TF_categorization_TF<-METH_TF_categorization
  
  METH_TF_categorization_TF[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_TF_current_patient)]<-0
  METH_TF_categorization_TF[which(dfPatientForAnalysis[,"METH_current_patient"] < meth_TF_current_patient)]<-1
  
} else {
  
  
  METH_TF_categorization_TF<-METH_TF_categorization
  
}


##
## Step 1.4: categorize the MUTATION
##

MUT_TF_categorization<-rep(0,nrow(dfPatientForAnalysis))
umutationallpatient<-(unique(MUT_current_patient[,1]))
#find which genes are mutated 
index_MUT_genes_allPatients<-which(dfPatientForAnalysis[,1] %in% umutationallpatient)
MUT_TF_categorization[index_MUT_genes_allPatients]<-1 #1 is mutated the gene

#create a data.frame with the update data
dfPatientForAnalysis_GAC<-cbind(dfPatientForAnalysis,
                                
                                GexpTF=rep(ge_TF_current_patient,nrow(dfPatientForAnalysis)),
                                FC_GE_TF=FC_GE_TF_categorization,
                                Genes_overexpressed=genes_overexpressed,
                                Genes_underexpressed=genes_underexpressed,
                                
                                CNV_TF=rep(cnv_TF_current_patient,nrow(dfPatientForAnalysis)),
                                CNV_EC_gain=if(cnv_TF_current_patient>=cnv_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                CNV_EC_depletion=if(cnv_TF_current_patient<=-cnv_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                  
                                CNV_gain=CNV_TF_gain,
                                CNV_depletion=CNV_TF_depletion,
                                CNV_TF_categorization_TF=CNV_TF_categorization_TF,
                                
                                METH_TF=rep(meth_TF_current_patient,nrow(dfPatientForAnalysis)),
                                METH_EC_hyper=if(meth_TF_current_patient>=meth_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                METH_EC_hypo=if(meth_TF_current_patient<meth_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                  
                                METH_hyper=METH_TF_hyper,
                                METH_hypo=METH_TF_hypo,
                                METH_TF_categorization_TF=METH_TF_categorization_TF,
                                
                                MUT_genes=MUT_TF_categorization,
                                MUT_TF=rep(mutation_TF_current_patient_variant,nrow(dfPatientForAnalysis)))

#columns that describe relation between genes and TF considering other data
col_relTF<-c("Hugo_Symbol","FC_GE_TF",
             "Genes_overexpressed","Genes_underexpressed",
             "CNV_EC_gain","CNV_EC_depletion",
             "CNV_gain","CNV_depletion",
             "CNV_TF_categorization_TF",
             "METH_EC_hyper","METH_EC_hypo",
             "METH_hyper","METH_hypo",
             "METH_TF_categorization_TF","MUT_genes","MUT_TF")

dfPatientForAnalysis_GAC_rel_TF<-dfPatientForAnalysis_GAC[,col_relTF]
###################################
rownames(dfPatientForAnalysis_GAC_rel_TF)<-dfPatientForAnalysis_GAC_rel_TF[,1] #26/10/2017

#do a control, if the properties of the genes are always equal to 0 it you can remove these columns
#this step it is important to reduce the computational cost.
resSumControl<-apply(dfPatientForAnalysis_GAC_rel_TF[,2:ncol(dfPatientForAnalysis_GAC_rel_TF)],2,sum)
names.good.properties<-names(resSumControl[which(resSumControl!=0)])

input_for_apriori<-cbind(Hugo_Symbol=dfPatientForAnalysis_GAC_rel_TF[,1],dfPatientForAnalysis_GAC_rel_TF[,names.good.properties])

input_for_apriori2 <- data.frame(sapply(input_for_apriori,as.factor))

###################################
  

print("Step 6: Run apriori")

#a confidence of 0.05 is a good solution
# support: a numeric value for the minimal support of an item set (default: 0.1)
rules = apriori(input_for_apriori2, parameter=list(support=0.0001, confidence=0.90,minlen=3,maxlen=ncol(input_for_apriori2),
                                      target = "closed frequent itemset"))

# support: a numeric value for the minimal support of an item set (default: 0.1)
# rules = apriori(trans , parameter=list(support=0.80, confidence=0.80,
#                                        target = "closed frequent itemset"))

ruledf = DATAFRAME(rules)
subrules.sort<-ruledf[order(ruledf$support,decreasing=T),]

subrules.sort2<-subrules.sort[grep(x=subrules.sort$items,pattern="Hugo_Symbol"),]

###
### Filtering of the rules
###
print("Step 7: Run filtering")

list_rules<-unique(unlist(lapply(X=strsplit(as.character(subrules.sort2[,1]),split=","),FUN=function(X){paste(X[2:length(X)],collapse=",")})))

genes_in_rules<-sapply(strsplit(sapply(strsplit(as.character(subrules.sort2[,1]),split=","),"[[",1),split="="),"[[",2)#

if(length(genes_in_rules)==nrow(input_for_apriori2)){ #if the number of genes detected in rules is the same of the initial rules do the follow processing
  

# Create a list in which the number of element of the list correspond with the number of rule, the content of each elements are the genes with a given rule

      association_rules_genes<-sapply(list_rules,FUN=function(x){
      
          gsub(sapply(strsplit(unlist(strsplit(grep(pattern=x,x=subrules.sort2[,1],value=T),split="Hugo_Symbol=")),split=","),"[[",1),pattern="\\{",replacement="")
        
        })

#parse the previously results in a data.frame
ARG<-list(1:length(association_rules_genes))

      for(l in 1:length(association_rules_genes)){
        
        current_rule<-names(association_rules_genes[l])
        number_of_genes<-length(association_rules_genes[[l]])
        data_arg_module<-data.frame(rule=rep(current_rule,number_of_genes),genes=association_rules_genes[[l]])
        ARG[[l]]<-data_arg_module
      }
      
      ARG2<-do.call(rbind,ARG)
      ARG2clean<-ARG2[ARG2[,2]!="",]
      ARG2clean[,1]<-as.character(ARG2clean[,1])
      
      #assign for each rule a number (module) that identify same rules with one id
      association_rule_module<-data.frame(rule=names(association_rules_genes),Groups_Apriori=1:length(association_rules_genes))
      association_rule_module[,1]<-as.character(association_rule_module[,1])
      
      ARM<-merge(x=association_rule_module,y=ARG2clean,by="rule")

} else { #statement if, code to decide how to treat the genes with multiple rules: slow to run
    
ugir<-unique(genes_in_rules)
    
rules_identified_from_multiple_genes<-sapply(X=ugir,FUN=function(X){ #start function sapply
    
    idx<-grep(subrules.sort2[,1],pattern=X)
  
    subrules.sort.gene<-subrules.sort2[idx,]
    
    if(nrow(subrules.sort.gene)==1){ # if for the current gene there is one rule
    
      rules<-as.character(subrules.sort.gene[,1]) # gmid find the rule and mantain it
    
    } else { # if for one genes are available more rules
      
      maxsupport<-max(subrules.sort.gene[,2]) #check the maximum values of support
      minsupport<-max(subrules.sort.gene[,2]) #check the minimum values of support
      
      if(minsupport!=maxsupport){
      
        rules<-as.character(subrules.sort.gene[which.max(subrules.sort.gene[,2]),1]) #get the rules with higher support value
      
      } else {
        
      }
      nostrings_rules<-unlist(lapply(X=strsplit(as.character(subrules.sort.gene[,1]),split=","),FUN=function(X){ #get the rules with the maximum sum of properties. 
      nostrings_rules<-gsub(X,pattern="[^0-9]",replacement="")                                                   #this is useful because GMID prefer to consider properties with 1 than 0
      sum(as.numeric(nostrings_rules[!nostrings_rules==""]))                                                     #at biological level i prefer features that explain a phenoma.
      }))
      
      idxrule<-which.max(nostrings_rules)
      rules<-as.character(subrules.sort.gene[idxrule,1])
    } #close internal if statement for genes with many rules
  
    } #close function sapply
    
    ) #close sapply

  #repeat the same processing of before
  rules_identified_from_multiple_genes_df<- data.frame(genes=names(rules_identified_from_multiple_genes),rules=rules_identified_from_multiple_genes)
  
  list_rules<-unique(unlist(lapply(X=strsplit(as.character(rules_identified_from_multiple_genes),split=","),FUN=function(X){paste(X[2:length(X)],collapse=",")})))
  
  association_rules_genes<-sapply(list_rules,FUN=function(x){
    
    gsub(sapply(strsplit(unlist(strsplit(grep(pattern=x,x=as.character(rules_identified_from_multiple_genes_df[,2]),value=T),split="Hugo_Symbol=")),split=","),"[[",1),pattern="\\{",replacement="")
    
  })
  
  ARG<-list(1:length(association_rules_genes))
  
  for(l in 1:length(association_rules_genes)){
    
    current_rule<-names(association_rules_genes[l])
    number_of_genes<-length(association_rules_genes[[l]])
    data_arg_module<-data.frame(rule=rep(current_rule,number_of_genes),genes=association_rules_genes[[l]])
    ARG[[l]]<-data_arg_module
  }
  
  ARG2<-do.call(rbind,ARG)
  ARG2clean<-ARG2[ARG2[,2]!="",]
  ARG2clean[,1]<-as.character(ARG2clean[,1])
  
  #assign for each rule a number (module) that identify same rules with one id
  association_rule_module<-data.frame(rule=names(association_rules_genes),Groups_Apriori=1:length(association_rules_genes))
  association_rule_module[,1]<-as.character(association_rule_module[,1])
  
  ARM<-merge(x=association_rule_module,y=ARG2clean,by="rule")
  
} #close statemenent if for one genes are avialable multiple rules

mergeGAC_COM<-merge(dfPatientForAnalysis_GAC,ARM,by.x="Hugo_Symbol",by.y="genes")
mergeGAC_COM_sort<-mergeGAC_COM[order(mergeGAC_COM$Groups_Apriori,decreasing=F),]


###
### Analysis 2: Use k-means for clustering of gene-expression, copy-number variation, methylation data.
###

  mergeGAC_COM_KGE<-kmeans(x=mergeGAC_COM_sort[,"GE_current_patient"],4)
  mergeGAC_COM_KCNV<-kmeans(x=mergeGAC_COM_sort[,"CNV_current_patient"],4)
  mergeGAC_COM_KMETH<-kmeans(x=mergeGAC_COM_sort[,"METH_current_patient"],4)
  
  mergeGAC_COM_res_K<-cbind(mergeGAC_COM_sort,clusterGE=mergeGAC_COM_KGE$cluster,clusterCNV=mergeGAC_COM_KCNV$cluster,clusterMETH=mergeGAC_COM_KMETH$cluster)

  #obtain the list of clusters
  clusters_list_for_2<-unique(mergeGAC_COM_res_K$Groups_Apriori)
  
  mergeGAC_COM_res_K_2<-mergeGAC_COM_res_K
    
  ###
  ### Merge the results of analysis with the drug database.
  ###

  #search in which rows are present the genes in the table of drugs-genes interactions
  #two columns in the tables of database must be present: genes and drug_primary_name
  indexSubDrug<-grep(x=tabDrugs$genes,pattern=paste(mergeGAC_COM_res_K_2[,1],collapse="|"))
  subtabDrugs<-tabDrugs[indexSubDrug,]
  
  #i use the symbol "#" to concatenate the strings, because when i will count the number of drugs for gene if are present "," inside the name of
  #drugs i will obtain a mistake number of genes (es. 1-2ethyl,diol,benze)
  collapseDrugTable<-aggregate(drug_primary_name ~ genes, data = subtabDrugs, paste,collapse = "#")
  
  mergeGAC_COM_res_K_2_drugs<-merge(mergeGAC_COM_res_K_2,collapseDrugTable,by.x="Hugo_Symbol",by.y="genes",all.x=T)
 
  countDrugsFunc<-function(x){
  #drugs name are repeated for the same genes, i used unique to manage this issue: the reason is that different database have the same drugs.
  ld<-as.numeric(length(unique(unlist(strsplit(x,split="#")))))
  
  return(ld)
  }
  resCountDrugs<-as.numeric(sapply(mergeGAC_COM_res_K_2_drugs$drug_primary_name,countDrugsFunc))
  resCountDrugs[which(is.na(mergeGAC_COM_res_K_2_drugs$drug_primary_name))]<-0
  
  #test the druggability of a modules
  mergeGAC_COM_res_K_2_drugs<-cbind(mergeGAC_COM_res_K_2_drugs,Count_Drugs_For_Gene=resCountDrugs)
 
  ### estimate the druggability of the modules
  
  print("Step 8: Estimate the scores")
  
  TOTAL_score_module<-NULL
  TOTAL_score_module_drugs<-NULL
  mergeGAC_COM_res_K_2_drugs$Groups_Apriori<-as.numeric(mergeGAC_COM_res_K_2_drugs$Groups_Apriori)
  
  mergeGAC_COM_res_K_2_drugs<-mergeGAC_COM_res_K_2_drugs[order(mergeGAC_COM_res_K_2_drugs$Groups_Apriori),]
  
  uniqGA<-as.numeric(unique(mergeGAC_COM_res_K_2_drugs$Groups_Apriori))

  for(mga in uniqGA){
   
   print(mga)
  
   smgcrkl<- mergeGAC_COM_res_K_2_drugs[mergeGAC_COM_res_K_2_drugs$Groups_Apriori==mga,]
   #check the presenc of at least one column with 1 for a gene
   
   #check the number of modules alterated: the first column is hugo symbol, the second is the fold-change between the expression (is not an alteration)
   #of tf and target genes. I remove these columns because are not useful in the categorization of alterated and not alterated genes
   alterated_genes_in_module<-apply(smgcrkl[,col_relTF][,-c(1:2)],1,FUN=function(x){ifelse(sum(x)==0,"Not_altered","Alterated")})
   
  #count the number of not alterated genes in modules
  check_alteration<-cbind(smgcrkl,alterated_genes_in_module)
  
  nrow_module<-nrow(smgcrkl[,col_relTF][,-c(1:2)])
  
  ncol_module<-ncol(smgcrkl[,col_relTF][,-c(1:2)])
  
  #total size module
  total_cell<-nrow_module*ncol_module
  
  number_not_alteration_modules<-length(which(smgcrkl[,col_relTF][,-c(1:2)]==0))
  
  number_alteration_modules<-length(which(smgcrkl[,col_relTF][,-c(1:2)]==1))
  
  #estimate LNAM
  ratio_na_inmodule_with_totsize<-number_not_alteration_modules/total_cell
  #estimate LAM
  ratio_alt_inmodule_with_totsize<-number_alteration_modules/total_cell
  
  #estimate DELTA-A
  
  #estimate the difference between the number of cell without alteration and with alteration in modules
  #positive values indicate that the modules is again integrate otherwise not. [range-1,1]
  
  deltamodulesinalt<-ratio_na_inmodule_with_totsize-ratio_alt_inmodule_with_totsize
  scorestatusmodule<-as.numeric(rep(deltamodulesinalt,nrow_module))
  
  TOTAL_score_module<-c(TOTAL_score_module,scorestatusmodule)
  
  #estimate rdg
  genes_with_drugs<-length(which(smgcrkl[,"Count_Drugs_For_Gene"]>=1))/nrow_module
  
  #estimate NRDG
  genes_without_drugs<-length(which(smgcrkl[,"Count_Drugs_For_Gene"]==0))/nrow_module
  
  #estimate DELTA-D
  deltamodulegenesdrugs<-genes_with_drugs-genes_without_drugs
  deltamodulegenesdrugs <-rep(deltamodulegenesdrugs,nrow_module)
  
  TOTAL_score_module_drugs<-c(TOTAL_score_module_drugs,deltamodulegenesdrugs)
    
  }
  
  #estimate SAD
  mergeGAC_COM_res_K_2_drugs <-cbind(cbind(cbind(mergeGAC_COM_res_K_2_drugs,scorescorestatusmodule=TOTAL_score_module),TOTAL_score_module_drugs),combinedscore=TOTAL_score_module_drugs+TOTAL_score_module)
  
  assign(as.character(paste(se_patient_selection,".analysisGMIEC",sep="")),mergeGAC_COM_res_K_2_drugs) 
  
}

setwd(output_dir)

#1) save the .RData with all results
res_analysis_each_patient<-grep(ls(),pattern=".analysisGMIEC",value=T,fixed=T)
save(list=res_analysis_each_patient,file=paste(output_file,".analysis_single_patient.RData",sep=""))
save.image(file=paste(output_file,".ALL_analysis.RData",sep=""))

#save annotation results
write.table(selected_genes_after_GE_filtering,file=paste(paste("annotation_results_GMIEC",output_file,sep="."),".txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F) 

###
### In this step i want create the main output of the analysis the MATRIX with all  data
###

###extract number of modules in total
res_analysis_each_patient<-grep(ls(),pattern=".analysisGMIEC",value=T,fixed=T)

GPM_TOT<-NULL

for(bpm in res_analysis_each_patient){
  
  gpm<-unique(get(bpm)[,"Groups_Apriori"])
  GPM_TOT<-c(GPM_TOT,gpm)
}

#create the colnames to save the results
patient_id<-"patient_id"
total_genes_patients<-"total_genes_patients"
number_modules<-"number_modules"
total_modules<-unique(GPM_TOT)

number_genes_for_module_size<-paste("#genes_module",total_modules,sep="")
genes_in_each_module<-paste("genes_in_",paste("module",total_modules,sep=""),sep="")
number_drugs_for_module_size<-paste("#drugs_module",total_modules,sep="")
drugs_in_each_module<-paste("drugs_in_",paste("module",total_modules,sep=""),sep="")
scores_module_unique<-paste("score_alteration_module(deltaa)",total_modules,sep="")
scores_module_drugs_unique<-paste("score_alteration_drugs(deltad)",total_modules,sep="")
scores_module_sad_unique<-paste("score_sad",total_modules,sep="")

ALL_colnames<-c(patient_id,total_genes_patients,number_modules
                ,number_genes_for_module_size,genes_in_each_module,
                number_drugs_for_module_size,
                drugs_in_each_module,
                scores_module_unique,
                scores_module_drugs_unique,
                scores_module_sad_unique
)

MATRIX_RESULTS_ALL<-data.frame(matrix(,ncol=length(ALL_colnames)))
colnames(MATRIX_RESULTS_ALL)<-ALL_colnames

for(bpm in 1:length(res_analysis_each_patient)){
  
  data_current_patient<-get(res_analysis_each_patient[bpm])
  patient_id<- res_analysis_each_patient[bpm]
  names(patient_id)<-"patient_id"
  
  genes_patients<-data_current_patient[,1] 
  total_genes_patients<-length(data_current_patient[,2])
  names(total_genes_patients)<-"total_genes_patients"
  
  number_modules<-length(unique(data_current_patient[,"Groups_Apriori"]))
  names(number_modules)<-"number_modules"
  
  #extract the number of genes for modules and the genes inside each modules
  number_genes_for_module<-aggregate(Hugo_Symbol ~ Groups_Apriori, data =  data_current_patient[,c("Groups_Apriori","Hugo_Symbol")], length)
  
  number_genes_for_module_size<-number_genes_for_module[,2]
  names(number_genes_for_module_size)<-paste("#genes_module",number_genes_for_module[,1],sep="")
  
  genes_in_each_module<-aggregate(Hugo_Symbol ~ Groups_Apriori, data =  data_current_patient[,c("Groups_Apriori","Hugo_Symbol")], paste,collapse=",")
  genes_in_each_module<-genes_in_each_module[,2]
  names(genes_in_each_module)<-paste("genes_in_",paste("module",number_genes_for_module[,1],sep=""),sep="")
  
  #extract the number of drugs for modules and the genes inside each modules
  number_drugs_for_module<-aggregate(drug_primary_name ~ Groups_Apriori, data =  data_current_patient[,c("Groups_Apriori","drug_primary_name")], paste,collapse=",")
  number_drugs_for_module[,2]<-apply(number_drugs_for_module,1,FUN=function(x){length(unlist(strsplit(x[2],split="#")))})
  
  number_drugs_for_module_size<-number_drugs_for_module[,2]
  names(number_drugs_for_module_size)<-paste("#drugs_module",number_drugs_for_module[,1],sep="")
  
  drugs_in_each_module<-aggregate(drug_primary_name ~ Groups_Apriori, data =  data_current_patient[,c("Groups_Apriori","drug_primary_name")], paste,collapse="@") #use this character to merge drugs because it is important in the next step.
  drugs_in_each_module<-drugs_in_each_module[,2]
  names(drugs_in_each_module)<-paste("drugs_in_",paste("module",number_drugs_for_module[,1],sep=""),sep="")
  
  #extract the alteration scores for each module
  scores_module<-aggregate(scorescorestatusmodule ~ Groups_Apriori, data =  data_current_patient[,c("Groups_Apriori","scorescorestatusmodule")], unique)
  scores_module_unique<-scores_module[,2]
  names(scores_module_unique)<-paste("score_alteration_module(deltaa)",scores_module[,1],sep="")
  
  #extract the drugs scores for each module
  scores_module_drugs<-aggregate(TOTAL_score_module_drugs ~ Groups_Apriori, data =  data_current_patient[,c("Groups_Apriori","TOTAL_score_module_drugs")], unique)
  scores_module_drugs_unique<-scores_module_drugs[,2]
  names(scores_module_drugs_unique)<-paste("score_alteration_drugs(deltad)",scores_module_drugs[,1],sep="")
  
  #extract the sad scores for each module
  scores_module_sad<-aggregate(combinedscore ~ Groups_Apriori, data =  data_current_patient[,c("Groups_Apriori","combinedscore")], unique)
  scores_module_sad_unique<-scores_module_sad[,2]
  names(scores_module_sad_unique)<-paste("score_sad",scores_module_sad[,1],sep="")
  
  #extract rule for each module 
  rule_in_each_module<-aggregate(Groups_Apriori ~ rule, data =  data_current_patient[,c("Groups_Apriori","rule")], unique)
  rule_save<-rule_in_each_module[,1]
  names(rule_save)<-paste(paste("rule_module",rule_in_each_module[,2],sep=""))
  
  row_for_each_patient<-c(patient_id,
                          total_genes_patients,
                          number_modules,
                          number_genes_for_module_size,
                          genes_in_each_module,
                          number_drugs_for_module_size,
                          drugs_in_each_module,
                          scores_module_unique,
                          scores_module_drugs_unique,
                          scores_module_sad_unique,
                          rule_save
  )
  row_for_each_patient_t<-t(data.frame(row_for_each_patient))
  colnames(row_for_each_patient_t)<-names(row_for_each_patient)
  
  dfrow<-data.frame(row_for_each_patient_t,stringsAsFactors=F)
  colnames(dfrow)<-as.character(names(row_for_each_patient))
  # MATRIX_RESULTS_ALL[[bpm]]
  
  MATRIX_RESULTS_ALL<-rbind.fill(MATRIX_RESULTS_ALL,dfrow)
}

MATRIX_RESULTS_ALL$patient_id<-gsub(MATRIX_RESULTS_ALL$patient_id,pattern=".analysisGMIEC",replacement="")
  
MATRIX_RESULTS_ALL_CLINICAL<-merge(MATRIX_RESULTS_ALL[-1,],input_CLINICAL,by.x="patient_id",by.y="SAMPLE_ID")

#2) save the main results of analysis
write.table(t(MATRIX_RESULTS_ALL_CLINICAL[-1,]),file="Analysis_GMIEC_main_results.txt",sep="\t",row.names=T,col.names=F,quote=F) # the first row is always empty

### warning about the previously code: in the first step the total number of modules estimated during the analysis is saved.
### Data are combined by the same columns. Genes that are in the same modules (e.g. modules1) between two patients can have different 
### genomics featurest explore the data


#3) define a function to save the results of the modules
extractModules<-function(MATRIX_RESULTS_ALL_CLINICAL,selection=c("which.max","which.min","custom",thrsad=NULL)){
  
  selectColumns<-grep(grep(colnames(MATRIX_RESULTS_ALL_CLINICAL),pattern="#",invert=T,value=T),pattern="score_alteration",invert=T,value=T)
  filterMRAC<-MATRIX_RESULTS_ALL_CLINICAL[-1,selectColumns]
  idexSADcolumn<-grep(colnames(filterMRAC),pattern="sad") #index of SAD columns
  idXcutValue<-apply(filterMRAC[,idexSADcolumn],1,selection)
  
  ###
  ### define some controls
  ###
  
  if(selection=="which.max"){
    output<-"maxSAD"
    newcolnames<-paste(c("patientID","genes_in_module_with_","drugs_in_module_with_","score_of_module","rule_of_module"),output,sep="")
    
  }
  
  if(selection=="which.min"){
    
    output<-"minSAD"
    newcolnames<-paste(c("patientID","genes_in_module_with_","drugs_in_module_with_","score_of_module","rule_of_module"),output,sep="")
  
  }
  
  if(selection=="custom"){
    
            output1<-paste("custom_max",thrsad,sep="")
            output2<-paste("custom_min",-thrsad,sep="") 
            newcolnames1<-paste(c("patientID","genes_in_module_with_","drugs_in_module_with_","score_of_module","rule_of_module"),output1,sep="")
            newcolnames2<-paste(c("patientID","genes_in_module_with_","drugs_in_module_with_","score_of_module","rule_of_module"),output2,sep="")
  }
  
  
  simple.output.results<-data.frame()
  
  if(selection!="custom"){
    
  for(i in 1:nrow(filterMRAC[,idexSADcolumn])){
    
    cutValue.index<-idXcutValue[i] #the number of cutimum value is the same of the number of rows used for their estimation
    cutValue.extract<-filterMRAC[i,idexSADcolumn][cutValue.index]
    names.cutValues<-gsub(pattern="score_sad",names(cutValue.extract),replacement="")
    
    stringtogrep<-paste(c("genes_in_module","drugs_in_module","score_sad","rule_module"),names.cutValues,sep="")
    
    newcol1<- colnames(filterMRAC)[which(colnames(filterMRAC)%in% stringtogrep)]
    idx_start_clinical<-tail(grep(colnames(filterMRAC),pattern="rule"),1)+1 # see the index in which start the clinical data, i known that come after the rules
    newcol2<-colnames(filterMRAC)[idx_start_clinical:ncol(filterMRAC)]
    totalCol<-c(colnames(filterMRAC)[1],newcol1,newcol2)
    
    subRow<-cbind(patientID=rownames(filterMRAC)[i],data.frame(filterMRAC[i,totalCol]),n.module.cut.sad=names.cutValues,stringsAsFactors=F)
    #the first four columns are always the same
    colnames(subRow)[2:6]<-newcolnames
    
    simple.output.results<-rbind(simple.output.results,subRow)
    
  }
  
  simple.output.results<-cbind(simple.output.results[,-1],type.of.rule=as.numeric(as.factor(simple.output.results$rule_of_module)))
  
  write.table(simple.output.results,file=paste(paste("Analysis_GMIEC_simplified_results",output,sep=""),".txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F) # the first row is always empty
  
  } else {
    

    #define the output for output custom mode
    selectColumns<-grep(grep(colnames(MATRIX_RESULTS_ALL_CLINICAL),pattern="#",invert=T,value=T),pattern="score_alteration",invert=T,value=T)
    filterMRAC<-MATRIX_RESULTS_ALL_CLINICAL[-1,selectColumns]
    idexSADcolumn<-grep(colnames(filterMRAC),pattern="sad") #index of SAD columns

    idXMinValue<-apply(X=filterMRAC[,idexSADcolumn],1, FUN=function(X){
      if(length(which(as.numeric(X) < -thrsad))!=0) 
      {which(as.numeric(X) < -thrsad)} else {
        which.min(as.numeric(X))
      }
    }
    )
    
    idXMaxValue<-apply(X=filterMRAC[,idexSADcolumn],1, FUN=function(X){
      if(length(which(as.numeric(X) > thrsad))!=0)
      {which(as.numeric(X) > thrsad)} else {
        which.max(as.numeric(X))
      }
    }
    )
    
    list_parameter_for_custom_analysis<-list()
    list_parameter_for_custom_analysis[[1]]<-c(-thrsad,thrsad)
    list_parameter_for_custom_analysis[[2]]<-list(idXMinValue,idXMaxValue)
    list_parameter_for_custom_analysis[[3]]<-data.frame(newcolnames1,newcolnames2)
    list_parameter_for_custom_analysis[[4]]<-data.frame(output2,output1)
    
  for(ts in 1:length(c(-thrsad,thrsad))){
    
    simple.output.results<-data.frame()
    
    idXcutValue<-list_parameter_for_custom_analysis[[2]][ts]
    
    for(i in 1:nrow(filterMRAC[,idexSADcolumn])){
      
      cutValue.index<-unlist(idXcutValue[[1]][i]) #the number of cutimum value is the same of the number of rows used for their estimation
      cutValue.extract<-filterMRAC[i,idexSADcolumn][cutValue.index]
      names.cutValues<-gsub(pattern="score_sad",names(cutValue.extract),replacement="")
      
      stringtogrep<-NULL
      
      for(lncv in names.cutValues){
        
      stg<-paste(c("genes_in_module","drugs_in_module","score_sad","rule_module"),lncv,sep="")
      stringtogrep<-c(stringtogrep,stg) 
      }
      
      newcol1<- colnames(filterMRAC)[which(colnames(filterMRAC)%in% stringtogrep)]
      idx_start_clinical<-tail(grep(colnames(filterMRAC),pattern="rule"),1)+1 # see the index in which start the clinical data, i known that come after the rules
      newcol2<-colnames(filterMRAC)[idx_start_clinical:ncol(filterMRAC)]
      totalCol<-c(colnames(filterMRAC)[1],newcol1,newcol2)
      
      
      ##group genes 
      subselectionGenes<-grep(totalCol,pattern="genes_in_module",value=T)
      dfgenes<-data.frame(genes_in_module=paste(as.character(filterMRAC[i,subselectionGenes]),collapse=","))
      
      ##group drugs
      subselectionDrugs<-grep(totalCol,pattern="drugs_in_module",value=T)
      dfdrugs<-data.frame(drugs_in_module=paste(as.character(filterMRAC[i,subselectionDrugs]),collapse=","))
       
      ##group score
      subselectionsad<-grep(totalCol,pattern="sad",value=T)
      dfsad<-data.frame(sad_of_module=paste(as.numeric(filterMRAC[i,subselectionsad]),collapse=","))
      ##group rule
      subselectionrule<-grep(totalCol,pattern="rule",value=T)
      dfrule<-data.frame(rule_of_module=paste(as.character(filterMRAC[i,subselectionrule]),collapse="-"))
      
      stringCV<-t(data.frame(paste(names.cutValues,collapse="-")))

      subRow<-cbind(patientID=filterMRAC[i,1],dfgenes,dfdrugs,dfrule,dfsad,n.module.cut.sad=stringCV,stringsAsFactors=F)
      #the first four columns are always the same
      simple.output.results<-rbind(simple.output.results,subRow)
      
    }
    
    write.table(simple.output.results,file=paste(paste("Analysis_GMIEC_simplified_results",as.character(unlist(list_parameter_for_custom_analysis[[4]][ts])),sep=""),".txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F) # the first row is always empty
  }

}
  
}

extractModules(MATRIX_RESULTS_ALL_CLINICAL,selection="which.max")
extractModules(MATRIX_RESULTS_ALL_CLINICAL,selection="which.min")

###
### extract data with custom user defined threshold
###

module_thr.max<-0.6
module_thr.min<- -(0.6)

selectColumns<-grep(grep(colnames(MATRIX_RESULTS_ALL_CLINICAL),pattern="#",invert=T,value=T),pattern="score_alteration",invert=T,value=T)
filterMRAC<-MATRIX_RESULTS_ALL_CLINICAL[-1,selectColumns]
idexSADcolumn<-grep(colnames(filterMRAC),pattern="sad") #index of SAD columns

idXMinValue<-apply(X=filterMRAC[,idexSADcolumn],1, FUN=function(X){
 if(length(which(as.numeric(X) < module_thr.min))!=0) 
    {which(as.numeric(X) < module_thr.min)} else {
      which.min(as.numeric(X))
      }
    }
 )

idXMaxValue<-apply(X=filterMRAC[,idexSADcolumn],1, FUN=function(X){
  if(length(which(as.numeric(X) > module_thr.max))!=0)
  {which(as.numeric(X) > module_thr.max)} else {
    which.max(as.numeric(X))
  }
}
)

