library(dplyr)       # Manipulation des données
library(tidyr)       # Restructuration des données
library(ggplot2)     # Visualisations graphiques
library(lubridate)   # Manipulation des dates
library(pheatmap)    # Création de heatmaps
library(RColorBrewer) # Palettes de couleurs
library(stringr)
library(tibble)
library(patchwork)
library(forcats)
library(psych)
library(dplR)

setwd("~/Documents/Articles/Baptiste/Wood_anatomy/Elevation_dependent_response/Monthly_correlations")

# All the shortcuts match with a real parameter DON'T CHANGE THE NAMES

Mto<-"Eobs" # Data meteo -> where does it come from
site<-c("GRA") # Data site
mto_parm<-c("tg","rr","spei") # Data meteo -> what are we using (temp, precip, drought)
Species<-c("PICE") # Species considered
Detrend<-c("Spline") # Detrending method
Raw<-c("raw") # ?
QWA_parm<-c("CWTRAD") # Parameter considered
qtile<-.75
Bandwidth<- c(30)
sig.level<-0.01
nyrs_dtd=30 # Number of years used for spline detrending
f=0.5 # Frequency of spline detrending (use for dplr spline detrending)

Beg_cor_yr<-1920
End_cor_yr<-2020

QWA_fld_name<-paste("~/Documents/Articles/Baptiste/Wood_anatomy/Elevation_dependent_response/",site,"_",Species,"_OUTPUTS_txt",sep="")
MTO_fld_name<-paste("~/Documents/Articles/Baptiste/Wood_anatomy/Elevation_dependent_response/Reconstitutions/Meteo_data/",site,"/",sep="")
RECONS_fld_name<-paste("~/Documents/Articles/Baptiste/Wood_anatomy/Elevation_dependent_response/Reconstitutions/",(paste(site,Species,"Recons/",sep="_")),sep="")

files <- list.files(QWA_fld_name,pattern = paste(Detrend,"chron.txt",sep="_"), full.names = TRUE)
QWA_parm <- str_extract(files, "(?<=_)[A-Za-z]+_[A-Za-z]+(?=_q_)")
QWA_parm[8:11] <- c("PC1","PC2","RPC1","RPC2")

Temperatures<-read.table(paste(MTO_fld_name,mto_parm[1],"_",Mto,"_",site,".txt",sep=""),header = T)%>%
  as.data.frame()%>% mutate(Year=format(as.Date(Date),"%Y"),Month=format(as.Date(Date),"%m"))%>%dplyr::filter(Year %in% c(Beg_cor_yr:End_cor_yr))

Monthly_temperatures <- Temperatures %>%
  dplyr::group_by(Year, Month) %>%
  dplyr::summarise(tg = mean(get(mto_parm[1]), na.rm = TRUE))%>%as.data.frame()

Monthly_temperatures<-Monthly_temperatures%>%
pivot_wider(names_from = Month, values_from = tg)%>%
  rename_with(~c("Year","tg_JAN","tg_FEB","tg_MAR","tg_APR","tg_MAY","tg_JUN","tg_JUL","tg_AUG","tg_SEP","tg_OCT","tg_NOV","tg_DEC"))%>%as.data.frame()

Monthly_temperatures_Spl<- detrend(Monthly_temperatures[,2:13],make.plot=F,method=c(Detrend),nyrs=nyrs_dtd,f=f,difference = T)%>%as.data.frame()
Monthly_temperatures_Spl <-data.frame(Year=Monthly_temperatures$Year,Monthly_temperatures_Spl)

Precipitation<-read.table(paste(MTO_fld_name,mto_parm[2],"_",Mto,"_",site,".txt",sep=""),header = T)%>%
  as.data.frame()%>% mutate(Year=format(as.Date(Date),"%Y"),Month=format(as.Date(Date),"%m"))%>%dplyr::filter(Year %in% c(Beg_cor_yr:End_cor_yr))

Monthly_precipitation <- Precipitation %>%
  dplyr::group_by(Year, Month) %>%
  dplyr::summarise(rr = sum(get(mto_parm[2]), na.rm = TRUE))%>%as.data.frame()

Monthly_precipitation<-Monthly_precipitation%>%
  pivot_wider(names_from = Month, values_from = rr)%>%
  rename_with(~c("Year","rr_JAN","rr_FEB","rr_MAR","rr_APR","rr_MAY","rr_JUN","rr_JUL","rr_AUG","rr_SEP","rr_OCT","rr_NOV","rr_DEC"))%>%as.data.frame()

Monthly_precipitation_Spl<- detrend(Monthly_precipitation[,2:13],make.plot=F,method=c(Detrend),nyrs=nyrs_dtd,f=f,difference = F)%>%as.data.frame() # detrend with ratio for precip, no neg values
Monthly_precipitation_Spl <-data.frame(Year=Monthly_temperatures$Year,Monthly_precipitation_Spl)

SPEI<-read.table(paste(MTO_fld_name,mto_parm[3],"_CRU_",site,".txt",sep=""),header = T)%>%
  as.data.frame()%>% mutate(Year=format(as.Date(Date),"%Y"),Month=format(as.Date(Date),"%m"))%>%dplyr::filter(Year %in% c(Beg_cor_yr:End_cor_yr))%>%
  dplyr::reframe(Year,Month,spei)

Monthly_SPEI<-SPEI%>%
  pivot_wider(names_from = Month, values_from = spei)%>%
  rename_with(~c("Year","spei_JAN","spei_FEB","spei_MAR","spei_APR","spei_MAY","spei_JUN","spei_JUL","spei_AUG","spei_SEP","spei_OCT","spei_NOV","spei_DEC"))%>%as.data.frame()

Monthly_SPEI_Spl<- detrend(Monthly_SPEI[,2:13],make.plot=F,method=c(Detrend),nyrs=nyrs_dtd,f=f,difference = T)%>%as.data.frame()
Monthly_SPEI_Spl <-data.frame(Year=Monthly_SPEI$Year,Monthly_SPEI_Spl)


Mat_cor<-matrix(NA,(length(files)),ncol(Monthly_precipitation%>% dplyr::select(!Year))*3)
colnames(Mat_cor)<-c(colnames(Monthly_temperatures[2:13]),colnames(Monthly_precipitation[2:13]),colnames(Monthly_SPEI[2:13]))
row.names(Mat_cor)<-c(QWA_parm)

Mat_sig<-matrix(NA,(length(files)),ncol(Monthly_precipitation%>% dplyr::select(!Year))*3)
colnames(Mat_sig)<-c(colnames(Monthly_temperatures[2:13]),colnames(Monthly_precipitation[2:13]),colnames(Monthly_SPEI[2:13]))
row.names(Mat_sig)<-c(QWA_parm)


for (i in 1:length(files))
  {
WA_data <- read.table(files[i],h=T)%>%dplyr::select(Year=year,std)%>%filter(Year %in% c(Beg_cor_yr:End_cor_yr))
Mat_cor[i,1:12]<-cor(WA_data$std,Monthly_temperatures_Spl%>% dplyr::select(!Year))
Mat_cor[i,13:24]<-cor(WA_data$std,Monthly_precipitation_Spl%>% dplyr::select(!Year))
Mat_cor[i,25:36]<-cor(WA_data$std,Monthly_SPEI_Spl%>% dplyr::select(!Year))
}

Mat_cor_df<-Mat_cor%>%as.data.frame()


for (i in 1:length(files))
{
  for (j in 1:(ncol(Monthly_temperatures)-1))
  {
WA_data <- read.table(files[i],h=T)%>%dplyr::select(Year=year,std)%>%filter(Year %in% c(Beg_cor_yr:End_cor_yr))
Mat_sig[i,j]<-cor.test(WA_data$std, Monthly_temperatures_Spl %>%
                         dplyr::select(-Year) %>%
                         dplyr::select(all_of(j)) %>%
                         pull(), method = "pearson")$p.value
Mat_sig[i,j+12]<-cor.test(WA_data$std, Monthly_precipitation_Spl %>%
                            dplyr::select(-Year) %>%
                            dplyr::select(all_of(j)) %>%
                            pull(), method = "pearson")$p.value
Mat_sig[i,j+24]<-cor.test(WA_data$std, Monthly_SPEI_Spl %>%
                            dplyr::select(-Year) %>%
                            dplyr::select(all_of(j)) %>%
                            pull(), method = "pearson")$p.value
  }
}

Mat_sig_df<-Mat_sig%>%as.data.frame()
Mat_sig_df_tg<-Mat_sig_df[,1:12]%>%as.data.frame()
Mat_sig_df_rr<-Mat_sig_df[,13:24]%>%as.data.frame()
Mat_sig_df_spei<-Mat_sig_df[,25:36]%>%as.data.frame()



# Define the custom order for Variable1
custom_order <- c(
                  "RPC1","PC1",
                  paste(Species,"LA",sep="_"),
                  paste(Species,"EWW",sep="_"),
                  paste(Species,"MRW",sep="_"),
                  "RPC2",
                  "PC2",
                  paste(Species,"aMXD",sep="_"),
                  paste(Species,"CWTRAD",sep="_"),
                  paste(Species,"CWTTAN",sep="_"),
                  paste(Species,"LWW",sep="_"))


# Prepare data for 'tg' variables with correct x-axis order
tg_long <- Mat_cor_df %>%
  rownames_to_column(var = "Variable1") %>%
  dplyr::select(Variable1, starts_with("tg")) %>%
  pivot_longer(-Variable1, names_to = "Variable2", values_to = "Correlation") %>%
  # Join the significance data (ensure alignment with the same Variable1 and Variable2)
  left_join(
    Mat_sig_df_tg %>%
      rownames_to_column(var = "Variable1") %>%
      dplyr::select(Variable1, starts_with("tg")) %>%
      pivot_longer(-Variable1, names_to = "Variable2", values_to = "Significant"),
    by = c("Variable1", "Variable2")
  ) %>%
  # Order Variable2 to ensure months are in chronological order
  mutate(
    Variable2 = factor(
      Variable2,
      levels = c(
        "tg_JAN", "tg_FEB", "tg_MAR", "tg_APR", "tg_MAY", "tg_JUN",
        "tg_JUL", "tg_AUG", "tg_SEP", "tg_OCT", "tg_NOV", "tg_DEC"
      )
    ),
    # Apply custom order to Variable1
    Variable1 = factor(Variable1, levels = custom_order)
  )
  


# Prepare data for 'rr' variables with correct x-axis order
rr_long <- Mat_cor_df %>%
  rownames_to_column(var = "Variable1") %>%
  dplyr::select(Variable1, starts_with("rr")) %>%
  pivot_longer(-Variable1, names_to = "Variable2", values_to = "Correlation") %>%
  # Join the significance data (ensure alignment with the same Variable1 and Variable2)
  left_join(
    Mat_sig_df_rr %>%
      rownames_to_column(var = "Variable1") %>%
      dplyr::select(Variable1, starts_with("rr")) %>%
      pivot_longer(-Variable1, names_to = "Variable2", values_to = "Significant"),
    by = c("Variable1", "Variable2")
  ) %>%
  # Order Variable2 to ensure months are in chronological order
  mutate(
    Variable2 = factor(
      Variable2,
      levels = c(
        "rr_JAN", "rr_FEB", "rr_MAR", "rr_APR", "rr_MAY", "rr_JUN",
        "rr_JUL", "rr_AUG", "rr_SEP", "rr_OCT", "rr_NOV", "rr_DEC"
      )
    ),
    # Apply custom order to Variable1
    Variable1 = factor(Variable1, levels = custom_order)
  )

# Prepare data for 'spei' variables with correct x-axis order
spei_long <- Mat_cor_df %>%
  rownames_to_column(var = "Variable1") %>%
  dplyr::select(Variable1, starts_with("spei")) %>%
  pivot_longer(-Variable1, names_to = "Variable2", values_to = "Correlation") %>%
  # Join the significance data (ensure alignment with the same Variable1 and Variable2)
  left_join(
    Mat_sig_df_spei %>%
      rownames_to_column(var = "Variable1") %>%
      dplyr::select(Variable1, starts_with("spei")) %>%
      pivot_longer(-Variable1, names_to = "Variable2", values_to = "Significant"),
    by = c("Variable1", "Variable2")
  ) %>%
  # Order Variable2 to ensure months are in chronological order
  mutate(
    Variable2 = factor(
      Variable2,
      levels = c(
        "spei_JAN", "spei_FEB", "spei_MAR", "spei_APR", "spei_MAY", "spei_JUN",
        "spei_JUL", "spei_AUG", "spei_SEP", "spei_OCT", "spei_NOV", "spei_DEC"
      )
    ),
    # Apply custom order to Variable1
    Variable1 = factor(Variable1, levels = custom_order)
  )

# Convert Significant column to logical based on p-value threshold

# Ensure Significant column is logical
tg_long <- tg_long %>%
  mutate(Significant2 = Significant < sig.level)

rr_long <- rr_long %>%
  mutate(Significant2 = Significant < sig.level)

spei_long <- spei_long %>%
  mutate(Significant2 = Significant < sig.level)



fixed_scale <- scale_fill_gradientn(
  colors = c("darkblue", "blue", "white", "orange", "darkred"),
  values = scales::rescale(seq(-0.8, 0.8, by = 0.2)),
  limits = c(-0.8, 0.8),  # Ensure consistent limits across plots
  name = "Correlation"
)

tg_heatmap <- ggplot(tg_long, aes(x = Variable2, y = Variable1, fill = Correlation)) +
  geom_tile(color = "white") +
  fixed_scale +  # Apply consistent scale
  geom_text(
    aes(label = ifelse(Significant2 == TRUE, "*", "")),
    color = "black", size = 3, fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  ) +
  labs(
    title = "Temperature",
    x = "Months",
    y = "Variables"
  ) +
  coord_fixed()

rr_heatmap <- ggplot(rr_long, aes(x = Variable2, y = Variable1, fill = Correlation)) +
  geom_tile(color = "white") +
  fixed_scale +  # Apply consistent scale
  geom_text(
    aes(label = ifelse(Significant2 == TRUE, "*", "")),
    color = "black", size = 3, fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  ) +
  labs(
    title = "Precipitation",
    x = "Months",
    y = "Variables"
  ) +
  coord_fixed()

spei_heatmap <- ggplot(spei_long, aes(x = Variable2, y = Variable1, fill = Correlation)) +
  geom_tile(color = "white") +
  fixed_scale +  # Apply consistent scale
  geom_text(
    aes(label = ifelse(Significant2 == TRUE, "*", "")),
    color = "black", size = 3, fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  ) +
  labs(
    title = "SPEI",
    x = "Months",
    y = "Variables"
  ) +
  coord_fixed()

# Combine plots
combined_plot <- tg_heatmap / rr_heatmap / spei_heatmap

# Export to PDF
# Export the combined plot to A4 portrait PDF
ggsave(path=paste(getwd(),"/",site,"_",Species,"_MO_CORREL/",sep=""),
  paste("Correl_map_",Species,"_",Beg_cor_yr,"-",End_cor_yr,".pdf",sep=""),
  plot = combined_plot,
  width = 8.27,  # A4 width in inches
  height = 11.69,  # A4 height in inches
  units = "in",
  dpi = 300
)

write.table(Mat_cor_df,paste(getwd(),"/",site,"_",Species,"_MO_CORREL/Matr_correl_",Species,Beg_cor_yr,"-",End_cor_yr,".txt",sep=""))


