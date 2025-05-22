###############################################
# Hemocyte distribution analysis
# date:2025-04-17         
# author:Daniel Prieto dprieto(at)fcien.edu.uy
# aim:To analyze and compare the distribution 
# of events in concentric ROIs (rings) to build
# Figure 2 of Parada et al. 2025
###############################################
#Load libraries
library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(ggsignif)
#Open files
#Read csv data
files.A <- read.csv("/home/daniel/Documentos/Cris/Paper2/Migration/resultsA.csv", header = TRUE)#Load data from group A
files.A$genotype <- "Control"#Assign genotype to control
#
files.B <- read.csv("/home/daniel/Documentos/Cris/Paper2/Migration/resultsB.csv", header = TRUE)#Load data from group B
files.B$genotype <- "Ptr23c"#Assign genotype to Ptr23c
##
migration <- rbind(files.A, files.B)#join both data groups
#write results to a table
write.csv(migration, file = "./results_migration.csv")#output directory
##
#Convert the comma-separated strings to numeric and reshape the data
migration_long <- migration %>%
  mutate(across(3:17, ~ as.numeric(gsub(",", ".", .)))) %>%
  pivot_longer(cols = 3:17, 
               names_to = "Replicate", 
               values_to = "Value")
##
#Calculate median values for each Ring-genotype combination
medians <- migration_long %>%
  group_by(Ring, genotype) %>%
  summarise(median_value = median(Value, na.rm = TRUE), .groups = 'drop')
###############################################
##Analyze the distributions
###############################################
#Mann-Whitney U test
mann_whitney_result <- wilcox.test(Value ~ genotype, data = migration_long)
print(mann_whitney_result)
#the groups are statistically different!
###############################################
###############################################
#Box plot only median lines
###############################################
plot_w_med <- ggplot(migration_long, aes(x = as.factor(Ring), y = Value, fill = genotype))+
  geom_boxplot(alpha = 0.5) +  #Violin plots with some transparency
  scale_x_discrete(breaks = seq(0,19,2))+
  ylim(c(0,6))+
  labs(x = "Ring", y = "Values", fill = "Genotype", color = "Genotype", facet = "Hemocytes per genotype by ring") +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_wrap(~ genotype, nrow = 3)  #Separate panels for each genotype
##
plot_w_med
###############################################
#PLot ring image
###############################################
library(cowplot)
library(magick)
ring_im <- ggdraw() +
  draw_image("/home/daniel/Documentos/Cris/Paper2/Migration/E24_RINGS_C1-Composite_rot.png")
##
###############################################
#Create a new df including summary stats
###############################################
summary_stats_by_ring <- migration_long %>%
  group_by(Ring, genotype) %>%
  summarise(
    count = n(),  #Number of observations
    mean = mean(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    iqr = IQR(Value, na.rm = TRUE),
    range_min = min(Value, na.rm = TRUE),
    range_max = max(Value, na.rm = TRUE),
    .groups = 'drop'
  )
#View the summary
#print(summary_stats_by_ring, n=28)
##
###############################################
#Perform a Mann-Whitney U test for each Ring to find which are statistically different
###############################################
wilcox_results <- migration_long %>%
  group_by(Ring) %>%
  summarise(
    wilcox_result = list(wilcox.test(Value ~ genotype)),
    .groups = 'drop'
  ) %>%
  mutate(wilcox_p_value = sapply(wilcox_result, function(x) x$p.value)) %>%
  select(Ring, wilcox_p_value)
#Show the Wilcoxon test results
#print(wilcox_results)
##
#Filter for significant p-values (p < 0.05)
significant_results_wilcox <- wilcox_results %>%
  filter(wilcox_p_value < 0.05)
#Show only significant results
#print(significant_results_wilcox)#show results
##
migration_long$Ring <- factor(migration_long$Ring) #Define Ring as factor
migration_long$genotype <- factor(migration_long$genotype) #Define genotype as factor
##
#Define the significant Rings
significant_rings <- c(1, 2, 17, 18, 19)#data from line 93
#Filter data to include only the significant Rings
migration_long_significant <- migration_long %>%
  filter(Ring %in% significant_rings) #create a new df with the data from rings in line 109
#Reaname rings to "Ring something"
migration_long_significant <- migration_long_significant %>% 
  mutate(Ring = paste("Ring ", Ring)) #just to improve visualization afterwards
##
###############################################
#Plot and annotate with significance layer
###############################################
#Manually set the order of Ring
migration_long_significant$Ring <- factor(migration_long_significant$Ring,
                                          levels = c("Ring  1",  "Ring  2",  "Ring  17", "Ring  18", "Ring  19"))  #Replace with your desired order
#
plot_signif <- ggplot(migration_long_significant, aes(x = genotype, y = Value, fill = genotype)) +
  geom_violin(alpha = 0.5) +
  geom_point( alpha = 0.3)+
  ylim(0,7) +
  geom_signif(
    comparisons = list(c("Control", "Ptr23c")),
    map_signif_level = TRUE, textsize = 4      #Show stars
  ) +
  labs(x = "Genotype", y = "Values") +
  theme_bw() +
  theme(legend.position = 'none')+
  facet_wrap(~ Ring, nrow = 3, drop = TRUE)  #Facet by Ring
##
###############################################
##Assign colorblind-friendly palette
###############################################
plot_w_med_npg <- plot_w_med + scale_fill_viridis_d()
plot_signif_npg <- plot_signif + scale_fill_viridis_d()
##
###############################################
#Combine the plots into a single output using patchwork
###############################################
library(patchwork)#load the library
combined_plot1 <- (ring_im | plot_w_med_npg) + 
  plot_layout(
    ncol = 2,       #columns
    nrow = 1,       #rows
    widths = c(0.4, 1),
  ) 
combined_plot1
##
###############################################
#Combine the plots into a single output using patchwork
###############################################
library(patchwork)#load the library
combined_plot <- (combined_plot1 / plot_signif_npg) + 
  plot_layout(
    ncol = 1,       #columns for the grid
    nrow = 2,       #rows
    heights = c(1, 1.1)  #Adjust widths
  ) + 
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 12, face = "bold")  #Match all tags
    )
  )
##
###############################################
#Export figure to pdf
###############################################
pdf("./Figure_2_sub.pdf")#write the figure to a pdf
combined_plot #show the plot
dev.off()#close pdf device
#EOF