library(rlist)
library(plyr)
library(ggplot2)
library(dplyr)
library(dplR)
library("dendroTools")
library(tidyr)
library(tidyverse)
library(data.table)
library(colorRamps)
library(ggforce)
library(patchwork)
library(paletteer)
library(GGally)
library(ggstatsplot)
library(zoo)

setwd("~/Documents/Articles/Baptiste/Wood_anatomy/Elevation_dependent_response/GRA_PICE_OUTPUTS_txt")

#QWA_dataset
Site<-"Gra"
Species<-"PICE"

# initialize all the parameters
Bandwidth<- c(30)#(10,20,30,40,50,100) # Bandwidth used to divide the ring (in micron). 
Min_Mork<-1 # Minimum value of the Mork index (RTSR) to classify cells in the latewood
posEW_LW<-40 # threshold used for the detection of EW/LW boundary. posEW_LW<-40: Cells from the LW should not be found in the first 40% of the ring
qtile<-0.25 # Value of quantile (to better define)
TPer<-c(1900:2023) # period covered by the series and the chronology. Used to plot the different series over the same time period (xaxis,plt_max)
parm<-c("EWW")#("CWTRAD","CWTTAN","MRW","LA","CWTRAD",aMXD)  
parm_range<- c(1:10) # typically for CWTRAD #c(0:1500) for LA # range covered by parameter parm. Used to plot the different series over the same (yaxis, plt_max)
RTSR_th<- 75 #percentage of cells > Mork Index in a given bandwith
#storage of outputs
fileList <- list() # empty list to store all cells characteristics for each trees
fileList1 <- list() # empty list to store max of cell characteristics per year
all_cells<-list() # list containing all the cells
series_plot_list<-list() # list containing the plots of 
RRADDISTRList<- list() # empty list to store RRADDISTR outputs

# detrending and smoothing
w.length=30 # window length to compute statistics EPS, rbar...
dtd_method<-c("Spline") # Method used to detrend the TR series
nyrs_dtd=30 # Number of years used for spline detrending
f=0.5 # Frequency of spline detrending (use for dplr spline detrending)
nyrs_series<- 30 # Number of years used for spline detrending the series (for vizualisation of plt_max in in ggplot)
nyrs_chron<- 5 # Number of years used for spline smoothing the final chronology (for vizualisation in ggplot)
EPS_th<-0.85 #Thresholds for EPS and SSS

stat_matrix=matrix(NA,2,length(Bandwidth)) # matrix used to store the statistics of a chronology (rbar,eps)
colnames(stat_matrix)=Bandwidth
row.names(stat_matrix)=c("rbar","eps")

# read the output files from Roxas
path <- getwd() 
dirs<-list.dirs('.', recursive=FALSE) # list all directories in the folder
dirs<-substring(dirs,2)
treename<-substring(dirs,2) # Extract name of each tree
site<-substr(treename[1],1,8) # Extract site and species
pdf(paste(site,parm,Bandwidth,"µm_intra_annual_profile.pdf",sep="_"),onefile = TRUE, width=20, height=5)


for(bw in (1:length(Bandwidth))) { # loop to divide the rings according to Bandwith
  for(d in (1:length(dirs))) { # loop on all the trees 
  file.path<-paste(path,dirs[d],sep="")
  
print(treename[d])

#### Read cells and ring files

# cells files

masterlist_cells <- list.files(file.path, pattern = ".*\\Output_Cells.txt") # list all files containing cells for a given tree
list_cell_files <- paste(file.path, masterlist_cells, sep = "/")
con_cell_files = c(list_cell_files)
txtread_cells <- lapply(con_cell_files,function(x)read.table(x, header=T)) # list with all the cells for a given tree
txtxmerge_cells_ <-  list.rbind(txtread_cells) # matrix with all the cells for a given tree
numeric_col<-txtxmerge_cells_[,-1]   # keep only numeric columns
numeric_col <- numeric_col %>%mutate_at(vars(-YEAR),function(x){replace(x, which(x<0), NA)})
txtxmerge_cells_<-cbind(ID=txtxmerge_cells_[,1],numeric_col)
txtxmerge_cells_$aMXD<-txtxmerge_cells_$CWA/(txtxmerge_cells_$CWA+txtxmerge_cells_$LA) #  Compute amxd following Edwards et al
all_cells[[d]]<-txtxmerge_cells_

# Ring files

masterlist_rings <- list.files(file.path, pattern = ".*\\Output_Rings.txt") # list all files containing rings for a given tree
list_tree_files <- paste(file.path, masterlist_rings, sep = "/")
con_tree_files = c(list_tree_files)
txtread_rings <- lapply(con_tree_files,function(x)read.table(x, header=T))
txtxmerge_rings <-  list.rbind(txtread_rings)
txtxmerge_rings<-txtxmerge_rings[order(txtxmerge_rings$YEAR,txtxmerge_rings$CNO),]
txtxmerge_rings_<-txtxmerge_rings[!duplicated(txtxmerge_rings$YEAR , fromLast = T),] #only keep both lines for output_rings analysis
txtxmerge_cells_$TCA=rep(NA,nrow(txtxmerge_cells_))
txtxmerge_cells_$TCA<-rowSums(cbind(txtxmerge_cells_$LA,txtxmerge_cells_$CWA))

# Extract WOODID: microsection and IMAGE cointaining each year

WOODID<-unique(substring(masterlist_cells,1,12))
IMAGE<-substring(masterlist_cells,1,16)

# Extract MRW for each year and divide MRW in sectors with a given Bandwidth

MRW<-txtxmerge_rings_$MRW # MRW=mean ring width in microns for a given tree
MRW<-MRW[-1]
MRW2<-round_any(MRW,Bandwidth[bw],floor) # Round the mean ring width as an exact multiple (integer) of bandwith 
NBSEC<-(MRW2/Bandwidth[bw])+1 #Number of sectors of a given bandwith within a ring
repMRW<-as.data.frame(cbind(MRW,NBSEC))
repMRW<-rep(repMRW$MRW,repMRW$NBSEC)
RADDISTR.BAND<-do.call(c,lapply(MRW2,fun <- function(x) {
  i <- seq(0,x,Bandwidth[bw])}))

# store YEAR, WOODID, IMAGE in an output matrix with a bandwidth resolution 
# Output_matrix = empty matrix created for each tree (one matrix per tree), each line represent a portion of tree ring
# with a given bandwith

YEAR<-txtxmerge_rings_$YEAR[-1]
repYEAR<-as.data.frame(cbind(YEAR,NBSEC))
repYEAR<-rep(repYEAR$YEAR,repYEAR$NBSEC)

WOODID<-rep(WOODID,length(RADDISTR.BAND))
WOODID<-as.character(WOODID)

IMAGE<-as.data.frame(txtxmerge_rings_$ID[-1])
repIMAGE<-as.data.frame(cbind(IMAGE=IMAGE,NBSEC))
IMAGE<-as.data.frame(rep(repIMAGE[,1],repIMAGE$NBSEC))
colnames(IMAGE)<-c("IMAGE")

output_matrix<-cbind(YEAR=repYEAR,WOODID,IMAGE,RADDISTR.BAND,MRW=repMRW)
output_matrix$RRADDISTR=rep(NA,nrow(output_matrix))
output_matrix$RRADDISTR2=rep(NA,nrow(output_matrix))
output_matrix$LA=rep(NA,nrow(output_matrix))
output_matrix$TCA=rep(NA,nrow(output_matrix))
output_matrix$DRAD=rep(NA,nrow(output_matrix))
output_matrix$DTAN=rep(NA,nrow(output_matrix))
output_matrix$CWA=rep(NA,nrow(output_matrix))
output_matrix$CWTALL=rep(NA,nrow(output_matrix))
output_matrix$CWTRAD=rep(NA,nrow(output_matrix))
output_matrix$CWTTAN=rep(NA,nrow(output_matrix))
output_matrix$RTSR=rep(NA,nrow(output_matrix))
output_matrix$CTSR=rep(NA,nrow(output_matrix))
output_matrix$RRTSR=rep(NA,nrow(output_matrix))
output_matrix$EWLW.ID=rep(NA,nrow(output_matrix))
output_matrix$TB2=rep(NA,nrow(output_matrix))
output_matrix$KH=rep(NA,nrow(output_matrix))
output_matrix$N.BAND=rep(NA,nrow(output_matrix))
output_matrix$DH=rep(NA,nrow(output_matrix))
output_matrix$YR.RRADDISTR.BAND<-rep(NA,nrow(output_matrix))
output_matrix$RADDIST.CONT<-rep(NA,nrow(output_matrix))
output_matrix$IMAGE.CHANGE<-rep(NA,nrow(output_matrix))
output_matrix$IMAGE.ID<-rep(NA,nrow(output_matrix))
output_matrix$aMXD<-rep(NA,nrow(output_matrix))

# fill  output_matrix year after year

rowIndex=1


for (y in 1:length(YEAR))
  {
  YEAR_Matrix<-txtxmerge_cells_[which(txtxmerge_cells_$YEAR==YEAR[y]),] # txtxmerge_cells_ contain all the cells of one tree; YEAR_Matrix : contain all the cells of one tree for one year y
  IndSEC<-NBSEC[y] # number of band in the ring of year y
  
    for (p in 1:IndSEC)
      {

    # fill  output_matrix bandwidth after bandwith
    SubSec_Matrix<-YEAR_Matrix[which(YEAR_Matrix$RRADDISTR>=(((p-1)/IndSEC)*100) & YEAR_Matrix$RRADDISTR<=(round_any(p/IndSEC*100,0.01,ceiling))),] # isolate the cells for one band of one year (y) per tree
    output_matrix$RRADDISTR[rowIndex]<-quantile(SubSec_Matrix$RRADDISTR,qtile,na.rm=T) # relative radial position of cell center within annual ring (%)
    output_matrix$RRADDISTR2=(output_matrix$RADDISTR.BAND/output_matrix$MRW)*100 
    output_matrix$LA[rowIndex]<-quantile(SubSec_Matrix$LA,qtile,na.rm=T)
    output_matrix$TCA[rowIndex]<-quantile(SubSec_Matrix$TCA,qtile,na.rm=T)
    output_matrix$DRAD[rowIndex]<-quantile(SubSec_Matrix$DRAD,qtile,na.rm=T)
    output_matrix$DTAN[rowIndex]<-quantile(SubSec_Matrix$DTAN,qtile,na.rm=T)
    output_matrix$N.BAND[rowIndex]<-nrow(SubSec_Matrix)
    output_matrix$CWA[rowIndex]<-quantile(SubSec_Matrix$CWA,qtile,na.rm=T)
    output_matrix$CWTALL[rowIndex]<-quantile(SubSec_Matrix$CWTALL,qtile,na.rm=T)
    output_matrix$CWTTAN[rowIndex]<-quantile(SubSec_Matrix$CWTTAN,qtile,na.rm=T)
    output_matrix$CWTRAD[rowIndex]<-quantile(SubSec_Matrix$CWTRAD,qtile,na.rm=T)
    output_matrix$RTSR[rowIndex] <- SubSec_Matrix %>%
  dplyr::filter(RTSR > Min_Mork) %>%   # Filter rows where RTSR > Min_Mork
  summarise(RTSR_perc = nrow(.) / nrow(SubSec_Matrix) * 100) %>%  # Use nrow() instead of n()
  pull(RTSR_perc) # RTSR Index Mork : ratio between 4*CWTTAN and Drad
    output_matrix$CTSR[rowIndex]<-quantile(SubSec_Matrix$CTSR,qtile,na.rm=T)
    output_matrix$RRTSR[rowIndex]<-length(SubSec_Matrix$RTSR[which(SubSec_Matrix$RTSR>1)])/length(SubSec_Matrix$RTSR[which(SubSec_Matrix$RTSR>0)])*100
    output_matrix$TB2[rowIndex]<-quantile(SubSec_Matrix$TB2,qtile,na.rm=T)
    output_matrix$KH[rowIndex]<-quantile(SubSec_Matrix$KH,qtile,na.rm=T)
    output_matrix$DH[rowIndex]<-quantile(SubSec_Matrix$DH,qtile,na.rm=T)
    output_matrix$aMXD[rowIndex]<-quantile(SubSec_Matrix$aMXD,qtile,na.rm=T)


   if(!is.na(output_matrix$RTSR[rowIndex]) & output_matrix$RTSR[rowIndex]>RTSR_th & output_matrix$RRADDISTR[rowIndex]>posEW_LW)
      {
        output_matrix$EWLW.ID[rowIndex]=c("LW")
      } else 
      {output_matrix$EWLW.ID[rowIndex]=c("EW")}
      
    rowIndex=rowIndex+1
        }
     
}

output_matrix <- output_matrix %>%
  dplyr::group_by(YEAR) %>%  # Group by YEAR and band
  dplyr::mutate(EWW = (sum(EWLW.ID == "EW", na.rm = TRUE)*Bandwidth),
         LWW=MRW-EWW) %>%  # Count "EW" per group
        ungroup()%>%as.data.frame()

 # transform bandwidth in decimal year
repMRW2<-as.data.frame(cbind(MRW2,NBSEC))
repMRW2<-rep(repMRW2$MRW2,repMRW2$NBSEC)
output_matrix$YR.RRADDISTR.BAND<-(output_matrix$RADDISTR.BAND/(repMRW2+Bandwidth[bw])+output_matrix$YEAR)

output_matrix$RADDIST.CONT[1]<-Bandwidth[bw]/2
Intpos<-(which(output_matrix$YR.RRADDISTR.BAND - floor(output_matrix$YR.RRADDISTR.BAND) == 0))

for (i in 2:length(Intpos))
{
  output_matrix$RADDIST.CONT[Intpos[i]]<-cumsum(MRW)[i-1]+(Bandwidth[bw]/2)
  output_matrix$RADDIST.CONT[Intpos[i]-1]<-output_matrix$RADDIST.CONT[Intpos[i]]-Bandwidth[bw]
  
}


Napos <- which(is.na(output_matrix$RADDIST.CONT))

for (j in 1:length(Napos))
{
  output_matrix$RADDIST.CONT[Napos[j]]<-output_matrix$RADDIST.CONT[Napos[j]-1]+(Bandwidth[bw])

}
output_matrix$RADDIST.CONT[length(output_matrix$RADDIST.CONT)]<-cumsum(MRW)[length(cumsum(MRW))]-Bandwidth[bw]/2

#plot time series of cwa parameter at a bandwidth resolution

breaks.major<-round(seq(min(TPer), max(TPer), by = 5),1) #round(seq(min(output_matrix$YR.RRADDISTR.BAND), max(output_matrix$YR.RRADDISTR.BAND), by = 5),1)
breaks.minor<-round(seq(min(TPer), max(TPer), by = 1),1)

# identify changes in images used for QWA
output_matrix<-output_matrix %>% mutate(cum_MRW=cumsum(RADDISTR.BAND)) 
output_matrix$IMAGE.CHANGE[which(output_matrix$IMAGE != dplyr::lag(output_matrix$IMAGE))]<-output_matrix$cum_MRW[which(output_matrix$IMAGE != dplyr::lag(output_matrix$IMAGE))]
output_matrix$IMAGE.ID[which(output_matrix$IMAGE != dplyr::lag(output_matrix$IMAGE))]<-substr(as.character(output_matrix$IMAGE[which(output_matrix$IMAGE != dplyr::lag(output_matrix$IMAGE))]),13,15)
output_matrix$IMAGE.ID[30]<-substr(as.character(output_matrix$IMAGE[30]),13,15) # find the first image and put the label at x=15 to avoid overlap with yaxis

fileList[[d]]<-output_matrix

colid<-which(colnames(output_matrix)==parm)

# Extract the max value of the parameter in the LW for all QWA parameters except LA

output_matrix_Max_LW<-output_matrix

if (parm=="LA"){
  suppressMessages(output_matrix_Max_LW2<-output_matrix_Max_LW %>% group_by(YEAR) %>% dplyr::filter(EWLW.ID =="EW")%>% 
    dplyr::select(cum_MRW,RRADDISTR2,parm)%>%slice_max(!! rlang::sym(parm), n = 1))
}


else {
  suppressMessages(output_matrix_Max_LW2<-output_matrix_Max_LW %>% group_by(YEAR) %>% dplyr::filter(EWLW.ID =="LW")%>% 
  dplyr::select(cum_MRW,RRADDISTR2,parm)%>%slice_max(!! rlang::sym(parm), n = 1))
}

output_matrix_Max_LW2 <- as.data.frame(output_matrix_Max_LW2)%>%complete(YEAR = min(YEAR):max(YEAR))
output_matrix_Max_LW2 <- output_matrix_Max_LW2 %>% distinct(YEAR, .keep_all = TRUE)
output_matrix_Max_LW2<-as.data.frame(output_matrix_Max_LW2)
coltofill<-fill.internal.NA(as.data.frame(output_matrix_Max_LW2[,parm]))
Spline_matrix <-as.data.frame(cbind(output_matrix_Max_LW2$YEAR,output_matrix_Max_LW2$cum_MRW, ffcsaps(coltofill[,1], nyrs=nyrs_series)))
colnames(Spline_matrix)<-c("YEAR","cum_MRW","Smooth_data")
output_matrix_Max_LW2<-Spline_matrix %>% left_join(.,output_matrix_Max_LW2, by =c("YEAR","cum_MRW"))

#store output matrix for each tree in a list
fileList1[[d]]<-output_matrix_Max_LW2

########################### Plot the QWA series, the max value of QWA parameter per year and add a spline smoothing curve ###########################
breaks.major<-output_matrix_Max_LW %>%
  group_by(cumsum(YEAR != lag(YEAR, default = first(YEAR)))) %>%
  slice(1)%>%as.data.frame()%>%
  slice(which(row_number() %% 5 == 1))%>%dplyr::select(YEAR,cum_MRW)
labels<-c(breaks.major$YEAR)
breaks.major<-c(breaks.major$cum_MRW)

breaks.minor<-output_matrix_Max_LW %>%
  group_by(cumsum(YEAR != lag(YEAR, default = first(YEAR)))) %>%
  slice(1)%>%as.data.frame()%>%dplyr::select(YEAR,cum_MRW)
breaks.minor<-c(breaks.minor$cum_MRW)


plot_series<-ggplot() +
  geom_line(data=output_matrix, aes(x=cum_MRW, y=output_matrix[,colid]),size=.1,col="#B2182B")+
  scale_x_continuous(expand = expansion(),labels = labels,
                     breaks =breaks.major,minor_breaks = breaks.minor)+
  scale_y_continuous(expand = expansion(),limits = c(min(parm_range), max(parm_range)))+
  xlab(paste('Year (resol=',Bandwidth,' microns)',sep=""))+
  ylab(parm)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        plot.title = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.grid.major.x = element_line(size=.3, color="black",linetype = 2),
        panel.grid.minor.x = element_line(size=.3, color="#8ccde3",linetype = 2))+
  ggtitle(unique(output_matrix$WOODID))+
  geom_vline(xintercept = c((na.omit(output_matrix$IMAGE.CHANGE))),col="#71874d",size=1)+
  geom_label(data=output_matrix,aes(label = IMAGE.ID,angle=0,x=cum_MRW,y=max(parm_range)-.5),size=3,label.size=NA,col="#71874d",fill="white")+
  geom_point(data=fileList1[[d]], aes(x=cum_MRW, y=output_matrix_Max_LW2[,parm]),size=.4,color="#486f83")+
  geom_line(data=fileList1[[d]],aes(x=cum_MRW,y=Smooth_data),col="#486f83",alpha=.6)

print(plot_series)

}
  dev.off()
  
Max_param <- lapply(fileList1, "[", , c("YEAR",parm)) # extract a list containing the max value of param for each year and each tree
Max_param<-Reduce(function(x,y) merge(x, y, by = "YEAR", all.x = TRUE, all.y = TRUE),Max_param)
colnames(Max_param)=c("YEAR",treename)
Max_param<-Max_param[!duplicated(Max_param),]

write.table(Max_param,paste(site,"Raw_Series",parm,"q",qtile*100,Bandwidth[bw],"mu.txt",sep="_"))


# Try to remove LBM outbreaks
#tmp <- zooreg(Max_param[,2], start = min(YEAR))
#tmp2 <- apply(tmp, 2,  FUN=isat(iis=T))

#apply(tmp, 2, function(column) {
  # Fit the isat model to the column
 # isat(as.numeric(column),iis=T)
#})


#apply( tmp2,iis=TRUE,plot=TRUE) # extract a list containing the max value of param for each year and each tree

#Max_param<-isat(Max_param,iis=TRUE,plot=TRUE)
# Analysis of raw series
# create a correlation matrix for raw series
Raw_series<-Max_param[,-1]
row.names(Raw_series)<-Max_param[,1]

# Plot correlation matrix for raw series
pdf(paste(site,parm,"q",qtile*100,Bandwidth,"µm_corr_matrix.pdf",sep="_"),onefile = TRUE, width=15, height=15)


corr_matrix<-ggcorrmat(
  data     = Raw_series,
  colors   = c("#B2182B", "white", "#4D4D4D"),
  sig.level = 0.05,
  title    = paste("Correllogram for",parm,sep=" "),
  subtitle = "units: microns")
print(corr_matrix)

dev.off()

# Plot time series

####################################################################################    
############## Raw Series##########################################################    
####################################################################################    

Raw_series_long<-Raw_series%>%
  tibble::rownames_to_column(var="Year")%>%
  pivot_longer(!Year,names_to = "Tree_ID", values_to = "RW")%>%
  mutate(Year=as.numeric(Year))

sk_max<-Raw_series_long%>%
  group_by(Tree_ID)%>%top_n(5,RW)

sk_min<-Raw_series_long%>%
  group_by(Tree_ID)%>%top_n(-5,RW)

breaks.major<-round(seq(min(TPer), max(TPer), by = 5),1) #round(seq(min(output_matrix$YR.RRADDISTR.BAND), max(output_matrix$YR.RRADDISTR.BAND), by = 5),1)
breaks.minor<-round(seq(min(TPer), max(TPer), by = 1),1)

pdf(paste(site,parm,"q",qtile*100,Bandwidth,"µm_raw_series.pdf",sep="_"),onefile = TRUE, width=18, height=24)


plot_raw_series_parm<-ggplot() +
  geom_line(data=Raw_series_long,aes(x=Year,y=RW)) +
  geom_point(data=sk_max,aes(x=Year,y=RW),col="#486f83",size=1.5) +
  geom_point(data=sk_min,aes(x=Year,y=RW),col="red",size=1.5) +
  facet_wrap(~Tree_ID,nrow=d,switch = "y")+ 
  scale_x_continuous(name="Year", expand = expansion(),limits = c(min(TPer), max(TPer)),
                     breaks =breaks.major,minor_breaks = breaks.minor)+
  scale_y_continuous(name=paste(parm,"(resol.",Bandwidth,"microns)" ))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        plot.title = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.grid.major.x = element_line(size=.3, color="black",linetype = 2),
        panel.grid.minor.x = element_line(size=.3, color="lightgrey",linetype = 3),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 6),
        strip.placement = "outside")
print(plot_raw_series_parm)

####################################################################################    
############Baillie Pilcher#########################################################  
####################################################################################    

BP_dtd<-log(Raw_series/rollmean(Raw_series,5,align = "center",fill=NA,partial = T))
BP_dtd[1:2,]<-log(sweep(Raw_series%>%slice(1:2), 2, colMeans(Raw_series%>%slice(1:5)), FUN = '/')) # compute B-P for the first years

BP_series_long<-BP_dtd%>%
  tibble::rownames_to_column(var="Year")%>%
  pivot_longer(!Year,names_to = "Tree_ID", values_to = "RW")%>%
  mutate(Year=as.numeric(Year))

sk_BP_max<-BP_series_long%>%
  group_by(Tree_ID)%>%top_n(5,RW)

sk_BP_min<-BP_series_long%>%
  group_by(Tree_ID)%>%top_n(-5,RW)

plot_BP_series_parm<-ggplot() +
  geom_line(data=BP_series_long,aes(x=Year,y=RW)) +
  geom_vline(data=sk_BP_max,aes(xintercept=Year),col="#486f83",size=0.5) +
  geom_vline(data=sk_BP_min,aes(xintercept=Year),col="red",size=0.5) +
  facet_wrap(~Tree_ID,nrow=d,switch = "y")+ 
  scale_x_continuous(name="Year", expand = expansion(),limits = c(min(TPer), max(TPer)),
                     breaks =breaks.major,minor_breaks = breaks.minor)+
  scale_y_continuous(name=paste(parm,"(resol.",Bandwidth,"microns)" ))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        plot.title = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.grid.major.x = element_line(size=.3, color="black",linetype = 2),
        panel.grid.minor.x = element_line(size=.3, color="lightgrey",linetype = 3),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 6),
        strip.placement = "outside")
print(plot_BP_series_parm)

####################################################################################    
############Spline##################################################################  
####################################################################################    

NA_pos<-which(is.na(Raw_series), TRUE)
output_df_filled<-fill.internal.NA(Raw_series, fill = c("Linear"))
Spl_dtd <- detrend(output_df_filled,make.plot=F,method=c(dtd_method),nyrs=nyrs_dtd)#f=f)
Spl_dtd[NA_pos]=NA

Spl_series_long<-Spl_dtd%>%
  tibble::rownames_to_column(var="Year")%>%
  pivot_longer(!Year,names_to = "Tree_ID", values_to = "RW")%>%
  mutate(Year=as.numeric(Year))

sk_Spl_max<-Spl_series_long%>%
  group_by(Tree_ID)%>%top_n(5,RW)

sk_Spl_min<-Spl_series_long%>%
  group_by(Tree_ID)%>%top_n(-5,RW)

py_Spl<-pointer(Spl_dtd, rgv.thresh = 20, nseries.thresh = 20, round.decimals = 2)%>%
  dplyr::filter(Nature==-1)%>%dplyr::select(year=Year,Nature)

plot_Spl_series_parm<-ggplot() +
  geom_line(data=Spl_series_long,aes(x=Year,y=RW)) +
  geom_vline(data=sk_BP_max,aes(xintercept=Year),col="#486f83",size=0.5) +
  geom_vline(data=sk_BP_min,aes(xintercept=Year),col="red",size=0.5) +
  facet_wrap(~Tree_ID,nrow=d,switch = "y")+ 
  scale_x_continuous(name="Year", expand = expansion(),limits = c(min(TPer), max(TPer)),
                     breaks =breaks.major,minor_breaks = breaks.minor)+
  scale_y_continuous(name=paste(parm,"(resol.",Bandwidth,"microns)" ))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        plot.title = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.grid.major.x = element_line(size=.3, color="black",linetype = 2),
        panel.grid.minor.x = element_line(size=.3, color="lightgrey",linetype = 3),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 6),
        strip.placement = "outside")
print(plot_Spl_series_parm)


dev.off()

####################################################################################    
############Chronologies############################################################  
####################################################################################  

# Raw data

Stat_Raw_series<-rwi.stats.running(Raw_series[which(rowSums(Max_param,na.rm=T)/(rowMeans(Max_param,na.rm=TRUE))>1),],window.length = w.length,window.overlap = (w.length-1))
Raw_chron<-chron(Raw_series,biweight=T)
Raw_chron$upsd<-Raw_chron$std+apply(Raw_series,1,sd,na.rm=T)
Raw_chron$losd<-Raw_chron$std-apply(Raw_series,1,sd,na.rm=T)
Raw_chron$year=as.numeric(row.names(Raw_chron))
SSS_Raw_chron<-as.data.frame(sss(Raw_series))
SSS_Raw_chron$year<-as.numeric(row.names(SSS_Raw_chron))
colnames(SSS_Raw_chron)<-c("sss","year")
Raw_chron <- Raw_chron %>% left_join(.,SSS_Raw_chron, by ="year")
Raw_chron <- Raw_chron %>% left_join(.,Stat_Raw_series, by = join_by(year == start.year))%>% dplyr::select(year,std,samp.depth,rbar.eff,eps,sss,upsd,losd)
Raw_chron2<-fill.internal.NA(Raw_chron, fill = c("Mean")) 

Raw_chron <- Raw_chron2 %>% 
  mutate(meanSpline = ffcsaps(std, nyrs=nyrs_chron), 
         EPSmeanSpline100 = ifelse(eps<EPS_th, meanSpline, NA),
         recons_okEPS=ifelse(eps>=EPS_th, std, NA),
         recons_lowEPS=ifelse(eps<EPS_th, std, NA),
         recons_lowSSS=ifelse(sss<=EPS_th, std, NA),
         recons_okSSS=ifelse(sss>EPS_th, std, NA),
         sd_lowSSS=ifelse(sss<EPS_th, samp.depth, NA),
         sd_okSSS=ifelse(sss>=EPS_th, samp.depth, NA))

# Detrended data

Stat_Dtd_series<-rwi.stats.running(Spl_dtd[which(rowSums(Max_param,na.rm=T)/(rowMeans(Max_param,na.rm=TRUE))>1),],window.length = w.length,window.overlap = (w.length-1))
Dtd_chron<-chron(Spl_dtd,biweight=T)
Dtd_chron$upsd<-Dtd_chron$std+apply(Spl_dtd,1,sd,na.rm=T)
Dtd_chron$losd<-Dtd_chron$std-apply(Spl_dtd,1,sd,na.rm=T)
Dtd_chron$year=as.numeric(row.names(Dtd_chron))
SSS_Dtd_chron<-as.data.frame(sss(Spl_dtd))
SSS_Dtd_chron$year<-as.numeric(row.names(SSS_Dtd_chron))
colnames(SSS_Dtd_chron)<-c("sss","year")
Dtd_chron <- Dtd_chron %>% left_join(.,SSS_Dtd_chron, by ="year")
Dtd_chron <- Dtd_chron %>% left_join(.,Stat_Raw_series, by = join_by(year == start.year))%>% dplyr::select(year,std,samp.depth,rbar.eff,eps,sss,upsd,losd)
Dtd_chron2<-fill.internal.NA(Dtd_chron, fill = c("Mean")) 

Dtd_chron <- Dtd_chron2 %>% 
  mutate(meanSpline = ffcsaps(std, nyrs=nyrs_chron), 
         EPSmeanSpline100 = ifelse(eps<EPS_th, meanSpline, NA),
         recons_okEPS=ifelse(eps>=EPS_th, std, NA),
         recons_lowEPS=ifelse(eps<EPS_th, std, NA),
         recons_lowSSS=ifelse(sss<=EPS_th, std, NA),
         recons_okSSS=ifelse(sss>EPS_th, std, NA),
         sd_lowSSS=ifelse(sss<EPS_th, samp.depth, NA),
         sd_okSSS=ifelse(sss>=EPS_th, samp.depth, NA))

Dtd_chron<-Dtd_chron %>% left_join(.,py_Spl, by=("year"))%>%
  mutate(py = case_when(Nature == -1 ~ std, 
               TRUE ~ NA))


theme_set(theme_light())
theme_chron<-theme(plot.title = element_text(size = 6),
                   axis.title = element_text(size = 6),
                   axis.text = element_text(size = 6),
                   axis.title.y.left = element_text(margin = margin(0, 10, 0, 0)),
                   axis.ticks.x.bottom = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))

Plot_rbar<-Raw_chron %>% 
  ggplot() +
  geom_line(aes(year,rbar.eff),colour="#a67d1c", alpha=1,size=0.2)+
  geom_line(aes(year,eps),colour="#9F1CA6", alpha=1,size=0.2)+
  geom_hline(yintercept=mean(Raw_chron$rbar.eff,na.rm=T),size=0.4,colour="#a67d1c",linetype="dashed")+
  geom_hline(yintercept=mean(Raw_chron$eps,na.rm=T),size=0.4,colour="#9F1CA6",linetype="dashed")+
  geom_hline(yintercept=0.85,size=0.4,colour="#9F1CA6")+
  scale_x_continuous(expand = c(0,0),limits=c(min(TPer),max(TPer)),breaks = seq(min(TPer), max(TPer), by = 10),label=NULL,element_blank())+
  scale_y_continuous(expand = c(0,0),limits=c(min(c(Raw_chron$rbar.eff,Raw_chron$eps),na.rm=T),1),breaks = seq(0,1,by =.2),sec.axis = sec_axis(transform=~.*1, name="EPS"))+
  xlab("Year CE") + ylab(paste("rbar (",w.length, "years)"))+
  theme_chron

Plot_Raw_chron<-Raw_chron %>% 
  ggplot() +
  geom_ribbon(aes(ymin = losd, ymax = upsd,x=year), fill = "azure2")+
  geom_line(aes(year,std),colour="#1C8AA6", alpha=1,size=0.2)+
  geom_line(aes(year, recons_okSSS),colour="#1C8AA6", alpha=1,size=0.2)+
  geom_line(aes(year, recons_lowSSS),colour="#c60001", alpha=1,size=0.2)+
  geom_line(aes(year, meanSpline),colour="#1c5ca6", alpha=1,size=0.4)+
  scale_x_continuous(expand = c(0,0),limits=c(min(TPer),max(TPer)),breaks = seq(min(TPer), max(TPer), by = 10))+
  scale_y_continuous(expand = c(0,0),limits=c(min(Raw_chron$losd),max(Raw_chron$upsd)),
                     breaks = seq(round(min(Raw_chron$losd,na.rm=T),0), round(max(max(Raw_chron$upsd,na.rm=T)),0),
                                  by=(round(max(Raw_chron$upsd,na.rm=T),0)-round(min(Raw_chron$losd,na.rm=T),0))/2))+
  xlab("Year CE") + ylab(paste("mean",parm, "(µm)",sep=" "))+
  theme_chron

Plot_Dtd_chron<-Dtd_chron %>% 
  ggplot() +
  geom_ribbon(aes(ymin = losd, ymax = upsd,x=year), fill = "#ffdd9c")+
  geom_line(aes(year,std),colour="#AC573D", alpha=1,size=0.2)+
  geom_line(aes(year, recons_okSSS),colour="#AC573D", alpha=1,size=0.2)+
  geom_line(aes(year, recons_lowSSS),colour="#c60001", alpha=1,size=0.2)+
  geom_line(aes(year, meanSpline),colour="#5e2717", alpha=1,size=0.4)+
  geom_point(aes(x=year,y=py),col="black",size=2) +
  scale_x_continuous(expand = c(0,0),limits=c(min(TPer),max(TPer)),breaks = seq(min(TPer), max(TPer), by = 10))+
  scale_y_continuous(expand = c(0,0),limits=c(min(Dtd_chron$losd),max(Dtd_chron$upsd)),
                     breaks = seq(round(min(Dtd_chron$losd,na.rm=T),0), round(max(max(Dtd_chron$upsd,na.rm=T)),0),
                                  by=(round(max(Dtd_chron$upsd,na.rm=T),0)-round(min(Dtd_chron$losd,na.rm=T),0))/2))+
  xlab("Year CE") + ylab(paste("mean",parm, "(µm)",sep=" "))+
  theme_chron

plot_sd<-Raw_chron %>% 
  ggplot() +
  geom_bar(aes(year, sd_okSSS),stat = "identity",position="stack",color="grey",width=0.1)+
  geom_bar(aes(year, sd_lowSSS),stat = "identity",position="stack",color="#c60001",width=0.1)+
  scale_x_continuous(expand = c(0,0),limits=c(min(TPer),max(TPer)),breaks = seq(min(TPer), max(TPer), by = 10),label=NULL,element_blank())+
  scale_y_reverse(expand = c(0,0),limits=c(20,0),breaks = seq(20, 0, by = -5),position = "right")+
  ylab("sample depth")+
  theme_chron

plot_sd  / Plot_rbar/Plot_Raw_chron/Plot_Dtd_chron+ plot_layout(heights = c(1,1,3,3))
ggsave(paste(site,parm,"q",qtile*100,Bandwidth,"µm_chron.pdf",sep="_"), width = 200, height = 120, units = "mm", dpi = 320)

write.table(Raw_chron,paste(site,parm,"q",qtile*100,Bandwidth,"µm_raw_chron.txt",sep="_"),row.names=T)
write.table(Dtd_chron,paste(site,parm,"q",qtile*100,Bandwidth,"µm",dtd_method,"chron.txt",sep="_"),row.names=T)


}








