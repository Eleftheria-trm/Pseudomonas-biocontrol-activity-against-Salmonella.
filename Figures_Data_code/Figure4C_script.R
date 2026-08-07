#set up for ggplot
require(ggplot2)
require(devtools)
#devtools::install_github('bbc/bbplot')
require(bbplot)
require(grid)
require(gridExtra)
require(DescTools)
require(ggpubr)
require(ggstance)
require(ggprism)
require(patchwork)
require(scales)
require(data.table)
require(pals)
require(RColorBrewer)

# CFU SEEDING NORMALISED BY WT & TIMEPOINT #####
# import data & process
setwd("~/Desktop/EleftheriaGitHub")
data <- read.csv('Figure4C_data.csv', stringsAsFactors = F)

# normalise CFU/seedling by STM50 and bio.rep
data$norm <- NA
for (i in seq(1, length(data$Strain))){
  data$norm[i] <- data$cfu.seedling[i] / mean(data$cfu.seedling[data$Strain == 'STM50' & data$bio.rep == data$bio.rep[i]])
}

# calculate mean & sd
data$mean <- NA
data$sd <- NA
for (i in seq(1, length(data$Strain))){
  data$mean[i] <- mean(data$norm[data$Strain == data$Strain[i]])
  data$sd[i] <- sqrt(var(data$norm[data$Strain == data$Strain[i]]))
}

# stats ##### 
model <- aov(norm ~ Strain, data=data)
summary(model)
TukeyHSD(model)
stats <- TukeyHSD(model)
stats <- as.data.frame(stats$Strain)
stats$strain1 <- gsub( "(.*)-(.*)", "\\1", rownames(stats))
stats$strain2 <- gsub( "(.*)-(.*)", "\\2", rownames(stats))

# rearrange labelling
stats <- subset(stats, strain1 == 'STM50' | strain2 == 'STM50')
for(i in seq(1, length(stats$diff))){
  if(stats$strain1[i] != 'STM50'){
    stats$strain2[i] <- stats$strain1[i]
    stats$strain1[i] <- 'STM50'
  }
}
# add asterisks
stats$stars <- 'ns'
for(i in seq(1, length(stats$diff))){
  if(stats$`p adj`[i] < 0.05){ stats$stars[i] <- '*' }
  if(stats$`p adj`[i] < 0.01){ stats$stars[i] <- '**' }
  if(stats$`p adj`[i] < 0.001){ stats$stars[i] <- '***' }
  if(stats$`p adj`[i] < 0.0001){ stats$stars[i] <- '****' }
}

# final plot ##### 
ggplot(data=subset(data, Timepoint == 48), aes(x=Strain, y=norm))+
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point(aes(x=Strain, y=norm, colour=Strain), size=2, position = position_jitter(width=.05)) +
  geom_boxplot(aes(x=Strain, y=norm, fill=Strain), width=.5, alpha=.25, outliers = F) +
  #geom_errorbar(aes(x=Strain, ymin=mean-sd, ymax=mean+sd, group=Strain), width=.25, alpha=.5) +
  scale_y_continuous(name=expression('Log'[10]*'-transformed CFU per seedling after 48hrs normalised to'~italic('Salmonella')), trans='log', expand=c(0,0), limits=c(10^-6,10^2),
                     breaks=c(10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10^0,10^1)) +
  scale_x_discrete(name='',limits=c("STM50","PK","trpC"),
                   labels=c(expression(italic("Salmonella")~"ctrl"),"PK307-4",expression(Delta*italic('trpC')))) +
  scale_colour_manual(values=c("#F8766D","#00BA38","#619CFF"), breaks=c("STM50","PK","trpC")) +
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"), breaks=c("STM50","PK","trpC")) +
  geom_segment(aes(x = 1, y = 15, xend = 2, yend = 15)) +
  geom_segment(aes(x = 1, y = 35, xend = 3, yend = 35)) +
  geom_text(data=subset(stats, strain2 == 'PK'), mapping = aes(x=1.5, y=21.5, label=stars), size=6) +
  geom_text(data=subset(stats, strain2 == 'trpC'), mapping = aes(x=2, y=55.5, label=stars), size=6) +
  theme_bw() +
  theme(plot.margin = unit(c(.5,0,0,.1),'cm'),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y=element_text(margin = margin(2,2,2,2), size=12),
        axis.text.x=element_text(size=12, vjust=1, hjust=1, angle=45),
        axis.line = element_line(color = "#cbcbcb"),
        axis.ticks = element_line(color = "#cbcbcb"),
        panel.grid.minor.y = element_blank(),
        legend.position='none')
