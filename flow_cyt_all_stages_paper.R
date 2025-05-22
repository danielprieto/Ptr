###############################################
# Flow cytometry analysis script  
# date:2022-12-09         
# author:Daniel Prieto dprieto(at)fcien.edu.uy
###############################################
#Load libraries
library(flowCore)
library(flowDensity)
library(flowViz)
library(flowStats)
library(flowAI)
library(ggcyto)
############################
#Set visual style for plots#
############################
# Step 1: Set the colorblind-friendly theme from the effects package
trellis.par.set(effectsTheme(col = "colorblind"))

# Step 2: Apply the current lattice theme to flowViz
flowViz.par.set(theme = trellis.par.get())

#
##
# Manually set a colorblind-friendly palette
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#B5B5B5")
trellis.par.set(overlay = list(col = cb_palette))  # Apply to points

##

#Open files
setwd <- ("~/Documentos/Cris/Citometria/01062023/")
files.test1 = list.files("~/Documentos/Cris/Citometria/01062023/data/", all.files = F, full.names = TRUE)#Lee todos los archivos del directorio
fs <- read.flowSet(files.test1[1:5], transformation = F, alter.names =T, emptyValue=FALSE)#[3:25])#Elimino las dos primeras filas que son /. y /..
####################
#QC & data cleaning#
####################
flow_auto_qc(fs, html_report = T, mini_report = FALSE, fcs_QC = TRUE, fcs_highQ=TRUE, folder_results = "~/Documentos/Cris/Citometria/01062023/dataHQ")
###################################
#Load QC-filtered data to files.hq#
###################################
files.hq = list.files("~/Documentos/Cris/Citometria/01062023/dataHQ", all.files = F, full.names = TRUE)#Lee todos los archivos del directorio
fs1 <- read.flowSet(files.hq[1:7], transformation = F, alter.names =T, emptyValue=FALSE)#[3:25])#Elimino las dos primeras filas que son /. y /..
################
#Transform data#
################
bx <- logicleTransform()
bxlist <- transformList(c("FSC.A", "SSC.A", "FSC.H", "YL2.A", "YL2.H", "YL2.W", "FSC.W"), bx)
datostrans <- transform(fs1, bxlist)
#
#Check plot
plotDens (datostrans[[1]], c("FSC.A" ,"SSC.A"))
##########################################
#Limpiamos el dataset de debris (fcs/ssc)#
##########################################
clean.dt <- datostrans
for (i in 1:7) { # Loop over the length of the flowSet
  f <- datostrans[[i]]
  # First restrict the FSC - A values :
  fsc.indices <- intersect (which (exprs (f)[, "FSC.A"] < 5) , which (exprs (f) [, "FSC.A"] > 3))
  # Then restrict SSC - A values and intersect with FSC - A restriction above :
  ssc.indices <- intersect (which ( exprs ( f)[, "SSC.A"] > 2) ,
                            which ( exprs (f)[, "SSC.A"] < 5) )
  non.debris.indices <- intersect ( fsc.indices , ssc.indices )
  # Now only select the non debris cells and place the cleaned up flowFrame into the flowSet :
  f.clean <- f[non.debris.indices]
  clean.dt [[i]] <- f.clean
}
#############
#Check plots#
#############
#plotDens (clean.dt[[1]], c("FSC.A" ,"SSC.A"))
#plotDens (clean.dt[[1]], c ("YL2.H", "FSC.W"))
#plotDens (clean.dt[[1]], c ("FSC.A", "FSC.H"))
#autoplot(clean.dt[[1]]) + labs_cyto("marker")
##
##########################
#Create data-driven gates#
##########################
cells <- lymphGate(clean.dt, channels=c("SSC.A", "FSC.A"), scale=2.2, bwFac = 1, plot = F)#Create data-driven gate
xyplot(`SSC.A` ~ `FSC.A`, clean.dt, filter = cells, xlim = c(3,5), ylim = c(2, 5), smooth = T, layout = c(2,4))#Show gated data in plots
celulas <- Subset(clean.dt, cells)#Apply gate
#################################
#Singlet filtering-elipsoid gate#
#################################
singfilt1 <- lymphGate(celulas, channels = c("FSC.H", "FSC.A"), preselection = NULL, scale=4, bwFac = 1, 
                       filterId = "singGate", evaluate = T, plot =T)
#Apply gate
sing <- Subset(celulas, singfilt1)
#
#####################################
##Normalize (limit) events displayed
#####################################
limit <- function(frame, limit = 10000) {
  frame@exprs <- frame@exprs[1:min(limit, nrow(frame@exprs)), ]
  frame
}#Create limit function
filtset_sing <- fsApply(sing, limit, limit = 10000)#Filter dataset and srt limit to 8K
filtset_cells <- fsApply(celulas, limit, limit = 15000)
##
#Define names
pData(filtset_sing)$name <- c("OregonR", "St.12 control", "St.12 Ptr23c", "St.14 control", "St.14 Ptr23c",
                              "St.16 control", "St.16 Ptr23c")
pData(filtset_cells)$name <- c("OregonR", "control st.12", "Ptr23c st.12", "control st.14", "Ptr23c st.14",
                               "control st.16", "Ptr23c st.16")
#
##
#################################
#Positive region counting filter#
#################################
#mCposFilt0 <- rectangleGate("FSC.A" = c(3.5, 4.7), "YL2.A" = c(2.3, 3.5))
#mCpos <- Subset(sing, mCposFilt0)#Apply gate
mCposFilt12 <- rangeGate(filtset_sing[2], "YL2.A", plot=T, refLine=0, alpha = 0.8)#Define data-driven gate based on sample 1
mCposFilt14 <-rangeGate(filtset_sing[4], "YL2.A", plot=T, refLine=0, alpha = 0.8)#Define data-driven gate based on sample 1
mCposFilt16 <-rangeGate(filtset_sing[6], "YL2.A", plot=T, refLine=0, alpha = 0.98)#Define data-driven gate based on sample 1

#mCR <- rectangleGate(filterId ="rg" , "YL2.A" = c(2.4, 4.5), "FSC.A" = c(4.3, 4.7))#Define data-driven gate based on sample 1
#rectangulo <- filter(filtset_cells, mCR)
#mC <- filter(sing, mCposFilt)
#mCcells <- Subset(sing, mC)
#mCpos2Filt <- curv2Filter("FSC.A", "YL2.A", filterId="data-driven filter", bwFac = 3)#Create a data-driven filter object that selects high-density regions in two dimensions.
#mCpos2 <- filter(sing, mCposFilt2)
#mCpos2data <- split(sing, mCpos2)#Split area data
#require(IDPmisc)
library(viridis)#load colorblind palette
colramp <- colorRampPalette(viridis(12)) #assign a 12-color viridis-based palette for overlay
#
plot12 <- xyplot(`FSC.A` ~ `YL2.A`, data = filtset_sing[2:3], filter=mCposFilt12, smooth=F, stats = T, pos=c(0.4, 0.6), 
                 xlab = "", ylab ="",
                 par.settings = list(panel.background=list(col="transparent"), 
                           axis.text=list(col="black"), strip.border=list(col="black"), 
                           axis.line=list(col="black"), strip.background=list(col="transparent"), 
                           flow.symbol = list(cex = 0.5, alpha = 0.5), gate = list(lwd = 2)), 
       xlim=c(0, 5), ylim=c(3, 5), layout = c(1,2),
       colramp = colramp)#Plot gate
plot12
#
plot14 <- xyplot(`FSC.A` ~ `YL2.A`, data = filtset_sing[4:5], filter=mCposFilt14, smooth=F, stats = T, pos=c(0.4, 0.6),
                 xlab = "", ylab ="",
                 par.settings = list(panel.background=list(col="transparent"), 
                                     axis.text=list(col="black"), strip.border=list(col="black"), 
                                     axis.line=list(col="black"), strip.background=list(col="transparent"), 
                                     flow.symbol = list(cex = 0.5, alpha = 0.5), gate = list(lwd = 2)), 
                 xlim=c(0, 5), ylim=c(3, 5), layout = c(1,2),
                 colramp = colramp)#Plot gate
#
plot16 <- xyplot(`FSC.A` ~ `YL2.A`, data = filtset_sing[6:7], filter=mCposFilt16, smooth=F, stats = T, pos=c(0.38, 0.6),
                 xlab = "", ylab ="",
                 par.settings = list(panel.background=list(col="transparent"), 
                                     axis.text=list(col="black"), strip.border=list(col="black"), 
                                     axis.line=list(col="black"), strip.background=list(col="transparent"), 
                                     flow.symbol = list(cex = 0.5, alpha = 0.5), gate = list(lwd = 2)), 
                 xlim=c(0, 5), ylim=c(3, 5), layout = c(1,2),
                 colramp = colramp)#Plot gate
#Assemble figure
library(ggpubr) #load library
require(grid)
#
figure <- ggarrange(plot12, NULL, plot14, NULL, plot16, 
         #labels = c("A", "", "B", "", "C"),
         widths = c(1.5, -0.25, 1.5, -0.25, 1.5),
          ncol = 5)
#Annotate axes
annotate_figure(figure, left = textGrob("FSC-A", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("mCherry-A", gp = gpar(cex = 1.3)))

###
#Just more code to get stupid crap

#######################################
#Evaluate double positives by quadrants#
########################################
qg <- quadrantGate(sing, c("YL2.A", "FSC.A"), plot=T, alpha = c(2, 2), sd = c(-2.7,-0.5))#Create quadrant gate
qfs <- split(filtset_cells, qg)#split quadrant data
xyplot(`FSC.A` ~ `YL2.A`, data = filtset_cells[2:7], filter=qg, smooth=F, 
       par.settings = list(panel.background=list(col="transparent"), 
                           axis.text=list(col="black"), strip.border=list(col="black"), 
                           axis.line=list(col="black"), strip.background=list(col="transparent")), 
       xlim=c(0, 5), ylim=c(3, 5), layout = c(3,2))#Plot gate

#################################
#mCloFilt <- rectangleGate("YL2.A" = c(2.6, 3.9), "SSC.A" = c(2.8, 4.4))
#mCHFilt <- rectangleGate("FSC.A" = c(4.3, 5), "SSC.A" = c(3.9, 4.8))

#Aplico el gate de Hemocitos
#mCpos <- Subset(sing, mCposFilt)
#hemo <- Subset(sing, mCHFilt)
##


#####################################
#Trellis plots for the whole flowSet#
#####################################
xyplot(`FSC.W` ~ `YL2.A`, data = filtset_sing, filter=mCposFilt, smooth=F, 
       par.settings = list(panel.background=list(col="transparent"), axis.text=list(col="black"), 
                           strip.border=list(col="black"), axis.line=list(col="black"), 
                           strip.background=list(col="transparent"), col = "colorblind"), xlab=c("mCherry-"), ylab=c("FSC-A"), 
       stats = TRUE, xlim=c(0, 4), ylim=c(3, 5), layout = c(2,4))#scatter pĺot
xyplot(`FSC.A` ~ `YL2.A`, data = filtset_sing[2:3], smooth=F, filter = mCposFilt,  stats = T, xlim=c(0, 5), ylim=c(2, 5), layout = c(2,1), pch = 0.01)#scatter pĺot
xyplot(`FSC.A` ~ `YL2.A`, data = celulas[2:3], smooth=F, filter = mCposFilt,  stats = T, xlim=c(0, 5), ylim=c(2, 6), layout = c(2,1), pch = 0.01)#scatter pĺot
xyplot(`FSC.A` ~ `YL2.A`, data = filtset_sing[2:3], smooth=F, filter = qg,  stats = F, xlim=c(0, 5), ylim=c(2, 5), layout = c(2,1), pch = 0.01)#scatter pĺot

#
########################
#Comparative histograms#
########################
densityplot (~ `YL2.H`,  filtset_sing[6:7], layout = c(1,1), xlab=("mCherry-H"), stack = T, stats=F, xlim = c(0,5), refline = median(mCcells@frames$`5-control_st16.fcs`[,4]@exprs))#, filter = mCposFilt)#, refline = marcador)
#c1f <- curv1Filter("YL2.A", bwFac=2.5)
#c1f.results <- filter(sing, c1f)
#densityplot (~ `YL2.A`,  sing, filter = mCposFilt, layout = c(1,1), xlim = c(2,4), xlab=("mCherry-A"), stack = T)
##############################
#Pass txt info for each sample
#Get data stats
##############################
lista <- summary(qfs)
dflist <- as.data.frame(matrix(unlist(lista), nrow = length(lista), byrow = TRUE, ncol = 16))
dflist
library(data.table)
rbindlist(lista)
library(dplyr)
df<-bind_rows(lista)

library(tibble)
df2<-enframe(lista)

split(lista, rep_len(1:16, length(lista)))

lista
data_frame <- as.data.frame(do.call(rbind, lista))
data_frame#
#pData(linfossing)$name <- c("pCoPuro s/ind", "pCoPuro Cu", "pMT_YC s/ind", "pMT_YC Cu", "pCoPuro s/ind PI", "pCoPuro ind PI", "pMT_YC s/ind PI", "pMT_YC ind PI")
#ggcyto(linfossing, aes(x = `FL1.A`, fill = name)) + geom_density(alpha = 0.2)
################################################################################
#Build monoparametric density plots and send them to an manageable object
#a.k.a publishable plots                                              
################################################################################
d1 <- ggplot(sing[1:2], aes(x = `FL1.A`, fill = NULL )) + geom_density(alpha = 0.2) + xlim (0.5,0.52) + xlab("YFP") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())+theme_classic()
d2 <- ggplot(sing[5:6], aes(x = `FL2.A`, fill = NULL )) + geom_density(alpha = 0.2) + xlim(0.5,0.52) + xlab("PI") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())+theme_classic()
d3 <- ggplot(sing[3:4], aes(x = `FL1.A`, fill = NULL )) + geom_density(alpha = 0.2) + xlim(0.5,0.52) + xlab("YFP") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())+theme_classic()
d4 <- ggplot(sing[7:8], aes(x = `FL2.A`, fill = NULL )) + geom_density(alpha = 0.2) + xlim(0.5,0.52) + xlab("PI") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())+theme_classic()
#
#####################
#Multipanel figure 1#
#####################
library(ggpubr)
ggarrange(a, d2, d1, d4 + rremove("x.text"), 
         labels = c("A", "B", "C", "D"),
        ncol = 2, nrow = 2)
##
sCount <- fsApply(sing, nrow)
sCount1 <- fsApply(mCpos, nrow)
