#set up for ggplot
require(ggplot2)
#install.packages('devtools')
require(devtools)
#devtools::install_github('bbc/bbplot')
require(bbplot)
require(grid)
#install.packages('gridExtra')
require(gridExtra)
#install.packages('DescTools')
require(DescTools)
#install.packages('ggpubr')
require(ggpubr)
require(scales)

setwd("~/Desktop/EleftheriaGitHub")
data <- read.csv("Figure2_data.csv", stringsAsFactors = F)

# plot 1 with ribbons = 1:1 ratio PKs mutants and STM50 #####
# import data
data1 <- subset(data, Experiment == 1)

# replace dashes in strain names with underscores
data1$Strain <- gsub('-','_',data1$Strain)

# make competition and ctrls name unique
data1$strain_competition <- paste(data1$Strain, data1$In.competition, sep='_')
# make competition and ctrls and cells name unique
data1$strain_competition_cell <- paste(data1$strain_competition, data1$Cell.counted, sep='_')
# make solo/competition column for make plots prettier...
data1$comp <- data1$Strain
data1$comp[data1$strain_competition == 'PK0307_4_STM50'] <- 'PK0307_4_STM50'

# calculate proportions of the community
for (i in seq(1, length(data1$Strain))){
  data1$rel[i] <- data1$cfu.ml[i] / sum(data1$cfu.ml[data1$Strain == data1$Strain[i] & data1$In.competition == data1$In.competition[i] &
                                                       data1$Time == data1$Time[i] & data1$tech.rep == data1$tech.rep[i]])
}

# calculate mean, sd, ci
for (i in seq(1, length(data1$Strain))){
  data1$mean[i] <- mean(data1$rel[data1$Strain == data1$Strain[i] & data1$In.competition == data1$In.competition[i] &
                                    data1$Cell.counted == data1$Cell.counted[i] & data1$Time == data1$Time[i]])
  data1$sd[i] <- sqrt(var(data1$rel[data1$Strain == data1$Strain[i] & data1$In.competition == data1$In.competition[i] &
                                      data1$Cell.counted == data1$Cell.counted[i] & data1$Time == data1$Time[i]]))
  data1$ci[i] <- qnorm(0.95)*data1$sd[i] / sqrt(length(data1$rel[data1$Strain == data1$Strain[i] & 
                                                                   data1$In.competition == data1$In.competition[i] &
                                                                   data1$Cell.counted == data1$Cell.counted[i] & data1$Time == data1$Time[i]]))
}

# plot proportions of the community over time with ribbons instead of error bars
plot1_rel2 <- ggplot(data=subset(data1, In.competition == 'STM50'), aes(x=Time,y=rel)) +
  geom_point(aes(x=Time, y=rel, colour=Cell.counted), position=position_dodge(width=.5)) +
  geom_ribbon(aes(x=Time, ymin=mean-sd, ymax=mean+sd, fill=Cell.counted), alpha=.25) +
  scale_y_continuous(name=expression("Proportion of each species in the population"),
                     limits=c(0,1), breaks=seq(0,1,.2)) +
  scale_x_continuous(name="Time (hours)", limits=c(0,72), breaks=seq(0,72,24)) +
  scale_colour_manual(name='', values=c("#F8766D","#00BA38"), limits=c('Salmonella','Pseudomonas'),
                      labels=c(expression(italic('Salmonella')),expression(italic('Pseudomonas'))))+
  scale_fill_manual(name='', values=c("#F8766D","#00BA38"), limits=c('Salmonella','Pseudomonas'),
                    labels=c(expression(italic('Salmonella')),expression(italic('Pseudomonas'))))+
  #annotate('text',x=stats$time1[stats$bug1 == 'Salmonella'],y=.45,label=stats$stars[stats$bug1 == 'Salmonella'], colour="#F8766D", size=8)+
  #annotate('text',x=stats$time1[stats$bug1 == 'Pseudomonas'],y=.55,label=stats$stars[stats$bug1 == 'Pseudomonas'], colour="#00BA38", size=8)+
  bbc_style() +
  theme(plot.margin = unit(c(.5,.5,0,0),'cm'),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14),
        axis.text.y=element_text(margin = margin(2,2,2,2), size=14),
        axis.text.x=element_text(size=14),
        axis.line = element_line(color = "#cbcbcb"),
        axis.ticks = element_line(color = "#cbcbcb"),
        panel.grid.major.y = element_blank(),
        legend.box='vertical',
        legend.title = element_text(size=14),
        legend.direction='horizontal',
        legend.position='bottom',
        legend.text=element_text(size=14, margin = margin(r = 30, unit = "pt"))) 
plot1_rel2

# stats for plot 1 = 1:1 ratio PKs mutants and STM50 #####
# anova
model <- aov(rel ~ strain_competition_cell * as.factor(Time) , data=data1)
summary(model)
stats <- TukeyHSD(model)
stats <- stats$`strain_competition_cell:as.factor(Time)`
stats <- as.data.frame(stats)

stats$strain1 <- gsub( "(.*)-(.*)", "\\1", rownames(stats))
stats$time1 <- gsub( "(.*):(.*)", "\\2", stats$strain1)
stats$strain1 <- gsub( "(.*):(.*)", "\\1", stats$strain1)
stats$bug1 <- gsub( "(.*)_(.*)", "\\2", stats$strain1)
stats$strain2 <- gsub( "(.*)-(.*)", "\\2", rownames(stats))
stats$time2 <- gsub( "(.*):(.*)", "\\2", stats$strain2)
stats$strain2 <- gsub( "(.*):(.*)", "\\1", stats$strain2)
stats$bug2 <- gsub( "(.*)_(.*)", "\\2", stats$strain2)
stats <- subset(stats, bug1 == bug2)
stats <- subset(stats, strain1 %in% c("PK0307_4_STM50_Pseudomonas","PK0307_4_STM50_Salmonella") & 
                  strain2 %in% c("PK0307_4_STM50_Pseudomonas","PK0307_4_STM50_Salmonella"))

stats$stars <- ''
stats$stars[stats$`p adj` < 0.05] <- '*'
stats$stars[stats$`p adj` < 0.01] <- '**'
stats$stars[stats$`p adj` < 0.001] <- '***'
stats$stars[stats$`p adj` < 0.0001] <- '****'

# make nicer table
stats1 <- stats[,c(7,9,6,4,11)]
rownames(stats1) <- seq(1, length(stats1$bug1))
colnames(stats1) <- c('species','time1','time2','p_value','significance')
View(stats1)

# plot 2 with ribbons = PKs applied first #####
# import data
data2 <- subset(data, Experiment == 2)

# replace dashes in strain names with underscores
data2$Strain <- gsub('-','_',data2$Strain)

# make competition and ctrls name unique
data2$strain_competition <- paste(data2$Strain, data2$In.competition, sep='_')
# make competition and ctrls and cells name unique
data2$strain_competition_cell <- paste(data2$strain_competition, data2$Cell.counted, sep='_')
# make solo/competition column for make plots prettier...
data2$comp <- data2$Strain
data2$comp[data2$strain_competition == 'PK0307_4_STM50'] <- 'PK0307_4_STM50'

# calculate proportion of each species in the community
for (i in seq(1, length(data2$Strain))){
  data2$rel[i] <- data2$cfu.ml[i] / sum(data2$cfu.ml[data2$Strain == data2$Strain[i] & data2$In.competition == data2$In.competition[i] &
                                                       data2$Time == data2$Time[i] & data2$tech.rep == data2$tech.rep[i]])
}

# calculate mean, sd, ci
for (i in seq(1, length(data2$Strain))){
  data2$mean[i] <- mean(data2$rel[data2$Strain == data2$Strain[i] & data2$In.competition == data2$In.competition[i] &
                                    data2$Cell.counted == data2$Cell.counted[i] & data2$Time == data2$Time[i]])
  data2$sd[i] <- sqrt(var(data2$rel[data2$Strain == data2$Strain[i] & data2$In.competition == data2$In.competition[i] &
                                      data2$Cell.counted == data2$Cell.counted[i] & data2$Time == data2$Time[i]]))
  data2$ci[i] <- qnorm(0.95)*data2$sd[i] / sqrt(length(data2$rel[data2$Strain == data2$Strain[i] & 
                                                                   data2$In.competition == data2$In.competition[i] &
                                                                   data2$Cell.counted == data2$Cell.counted[i] & data2$Time == data2$Time[i]]))
}

# add initial inocculum to table
x <- data2[c(12,13),]
x$Time <- 0
x$tech.rep <- NA
x$cfu <- NA
x$dil <- NA
x$cfu.ml <- NA
x$rel[x$Cell.counted == 'Pseudomonas'] <- 1
x$rel[x$Cell.counted == 'Salmonella'] <- 0
x$mean[x$Cell.counted == 'Pseudomonas'] <- 1
x$mean[x$Cell.counted == 'Salmonella'] <- 0
x$sd <- 0
x$ci <- 0

data2x <- rbind(data2, x)

# plot cfu/mL over time with ribbons instead of error bars
plot2_rel2 <- ggplot(data=subset(data2x, In.competition == 'STM50'), aes(x=Time,y=rel)) +
  geom_point(aes(x=Time, y=rel, colour=Cell.counted, group=comp), position=position_dodge(width=.5)) +
  geom_ribbon(aes(x=Time, ymin=mean-ci, ymax=mean+ci, fill=Cell.counted), alpha=.25)+
  scale_y_continuous(name=expression("Proportion of each species in the population"),
                     limits=c(0,1), breaks=seq(0,1,.2)) +
  scale_x_continuous(name="Time (hours)", limits=c(0,72), breaks=seq(0,72,24)) +
  scale_colour_manual(name='', values=c("#F8766D","#00BA38"), limits=c('Salmonella','Pseudomonas'),
                      labels=c(expression(italic('Salmonella')),expression(italic('Pseudomonas'))))+
  scale_fill_manual(name='', values=c("#F8766D","#00BA38"), limits=c('Salmonella','Pseudomonas'),
                    labels=c(expression(italic('Salmonella')),expression(italic('Pseudomonas'))))+
  #annotate('text',x=as.numeric(stats$time1[stats$bug1 == 'Salmonella']), y=4e09, label=stats$stars[stats$bug1 == 'Salmonella'], colour="#F8766D", size=8) +
  #annotate('text',x=as.numeric(stats$time1[stats$bug1 == 'Pseudomonas']), y=5e09, label=stats$stars[stats$bug1 == 'Pseudomonas'], colour="#00BA38", size=8) +
  bbc_style() +
  theme(plot.margin = unit(c(.5,.5,0,0),'cm'),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14),
        axis.text.y=element_text(margin = margin(2,2,2,2), size=14),
        axis.text.x=element_text(size=14),
        axis.line = element_line(color = "#cbcbcb"),
        axis.ticks = element_line(color = "#cbcbcb"),
        panel.grid.major.y = element_blank(),
        legend.box='vertical',
        legend.title = element_text(size=14),
        legend.direction='horizontal',
        legend.position='bottom',
        legend.text=element_text(size=14, margin = margin(r = 30, unit = "pt"))) 
plot2_rel2

# stats for plot 2 = PKs applied first #####
# anova
model <- aov(rel ~ strain_competition_cell * as.factor(Time), data=subset(data2, Time != 0))
summary(model)
stats <- TukeyHSD(model)
stats <- stats$`strain_competition_cell:as.factor(Time)`
stats <- as.data.frame(stats)

stats$strain1 <- gsub( "(.*)-(.*)", "\\1", rownames(stats))
stats$time1 <- gsub( "(.*):(.*)", "\\2", stats$strain1)
stats$strain1 <- gsub( "(.*):(.*)", "\\1", stats$strain1)
stats$bug1 <- gsub( "(.*)_(.*)", "\\2", stats$strain1)
stats$strain2 <- gsub( "(.*)-(.*)", "\\2", rownames(stats))
stats$time2 <- gsub( "(.*):(.*)", "\\2", stats$strain2)
stats$strain2 <- gsub( "(.*):(.*)", "\\1", stats$strain2)
stats$bug2 <- gsub( "(.*)_(.*)", "\\2", stats$strain2)
stats <- subset(stats, bug1 == bug2)
stats <- subset(stats, strain1 %in% c("PK0307_4_STM50_Pseudomonas","PK0307_4_STM50_Salmonella") & 
                  strain2 %in% c("PK0307_4_STM50_Pseudomonas","PK0307_4_STM50_Salmonella"))

stats$stars <- ''
stats$stars[stats$`p adj` < 0.05] <- '*'
stats$stars[stats$`p adj` < 0.01] <- '**'
stats$stars[stats$`p adj` < 0.001] <- '***'
stats$stars[stats$`p adj` < 0.0001] <- '****'

# make nicer table
stats2 <- stats[,c(7,9,6,4,11)]
rownames(stats2) <- seq(1, length(stats2$bug1))
colnames(stats2) <- c('species','time1','time2','p_value','significance')
View(stats2)

# plot 3 with ribbons = STM50 applied first #####
# import data
data3 <- subset(data, Experiment == 3)

# replace dashes in strain names with underscores
data3$Strain <- gsub('-','_',data3$Strain)

# make competition and ctrls name unique
data3$strain_competition <- paste(data3$Strain, data3$In.competition, sep='_')
# make competition and ctrls and cells name unique
data3$strain_competition_cell <- paste(data3$strain_competition, data3$Cell.counted, sep='_')
# make solo/competition column for make plots prettier...
data3$comp <- data3$Strain
data3$comp[data3$strain_competition == 'PK0307_4_STM50'] <- 'PK0307_4_STM50'

# calculate proportion of each species in the community
for (i in seq(1, length(data3$Strain))){
  data3$rel[i] <- data3$cfu.ml[i] / sum(data3$cfu.ml[data3$Strain == data3$Strain[i] & data3$In.competition == data3$In.competition[i] &
                                                       data3$Time == data3$Time[i] & data3$tech.rep == data3$tech.rep[i]])
}

# calculate mean, sd, ci
for (i in seq(1, length(data3$Strain))){
  data3$mean[i] <- mean(data3$rel[data3$Strain == data3$Strain[i] & data3$In.competition == data3$In.competition[i] &
                                    data3$Cell.counted == data3$Cell.counted[i] & data3$Time == data3$Time[i]])
  data3$sd[i] <- sqrt(var(data3$rel[data3$Strain == data3$Strain[i] & data3$In.competition == data3$In.competition[i] &
                                      data3$Cell.counted == data3$Cell.counted[i] & data3$Time == data3$Time[i]]))
  data3$ci[i] <- qnorm(0.95)*data3$sd[i] / sqrt(length(data3$rel[data3$Strain == data3$Strain[i] & 
                                                                   data3$In.competition == data3$In.competition[i] &
                                                                   data3$Cell.counted == data3$Cell.counted[i] & data3$Time == data3$Time[i]]))
}

# add initial inocculum to table
x <- data3[c(12,13),]
x$Time <- 0
x$tech.rep <- NA
x$cfu <- NA
x$dil <- NA
x$cfu.ml <- NA
x$rel[x$Cell.counted == 'Pseudomonas'] <- 0
x$rel[x$Cell.counted == 'Salmonella'] <- 1
x$mean[x$Cell.counted == 'Pseudomonas'] <- 0
x$mean[x$Cell.counted == 'Salmonella'] <- 1
x$sd <- 0
x$ci <- 0

data3x <- rbind(data3, x)

# plot proportions over time with ribbons instead of error bars
plot3_rel2 <- ggplot(data=subset(data3x, In.competition == 'STM50' | In.competition == 'None' & Time == 0), aes(x=Time,y=rel)) +
  geom_point(aes(x=Time, y=rel, colour=Cell.counted, group=comp), position=position_dodge(width=.5)) +
  geom_ribbon(aes(x=Time, ymin=mean-ci, ymax=mean+ci, fill=Cell.counted), alpha=.25) +
  scale_y_continuous(name=expression("Proportion of each species in the population"),
                     limits=c(-.04,1.04), breaks=seq(0,1,.2)) +
  scale_x_continuous(name="Time (hours)", limits=c(0,72), breaks=seq(0,72,24)) +
  scale_colour_manual(name='', values=c("#F8766D","#00BA38"), limits=c('Salmonella','Pseudomonas'),
                      labels=c(expression(italic('Salmonella')),expression(italic('Pseudomonas'))))+
  scale_fill_manual(name='', values=c("#F8766D","#00BA38"), limits=c('Salmonella','Pseudomonas'),
                    labels=c(expression(italic('Salmonella')),expression(italic('Pseudomonas'))))+
  annotate('text',x=as.numeric(stats$time1[stats$bug1 == 'Salmonella']), y=4e09, label=stats$stars[stats$bug1 == 'Salmonella'], colour="#F8766D", size=8) +
  annotate('text',x=as.numeric(stats$time1[stats$bug1 == 'Pseudomonas']), y=5e09, label=stats$stars[stats$bug1 == 'Pseudomonas'], colour="#00BA38", size=8) +
  bbc_style() +
  theme(plot.margin = unit(c(.5,.5,0,0),'cm'),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14),
        axis.text.y=element_text(margin = margin(2,2,2,2), size=14),
        axis.text.x=element_text(size=14),
        axis.line = element_line(color = "#cbcbcb"),
        axis.ticks = element_line(color = "#cbcbcb"),
        panel.grid.major.y = element_blank(),
        legend.box='vertical',
        legend.title = element_text(size=14),
        legend.direction='horizontal',
        legend.position='bottom',
        legend.text=element_text(size=14, margin = margin(r = 30, unit = "pt"))) 
plot3_rel2

# stats for plot 3 = STM50 applied first #####
# stats = anova
model <- aov(rel ~ strain_competition_cell * as.factor(Time), data=subset(data3, Time != 0))
summary(model)
stats <- TukeyHSD(model)
stats <- stats$`strain_competition_cell:as.factor(Time)`
stats <- as.data.frame(stats)

stats$strain1 <- gsub( "(.*)-(.*)", "\\1", rownames(stats))
stats$time1 <- gsub( "(.*):(.*)", "\\2", stats$strain1)
stats$strain1 <- gsub( "(.*):(.*)", "\\1", stats$strain1)
stats$bug1 <- gsub( "(.*)_(.*)", "\\2", stats$strain1)
stats$strain2 <- gsub( "(.*)-(.*)", "\\2", rownames(stats))
stats$time2 <- gsub( "(.*):(.*)", "\\2", stats$strain2)
stats$strain2 <- gsub( "(.*):(.*)", "\\1", stats$strain2)
stats$bug2 <- gsub( "(.*)_(.*)", "\\2", stats$strain2)
stats <- subset(stats, bug1 == bug2)
stats <- subset(stats, strain1 %in% c("PK0307_4_STM50_Pseudomonas","PK0307_4_STM50_Salmonella") & 
                  strain2 %in% c("PK0307_4_STM50_Pseudomonas","PK0307_4_STM50_Salmonella"))

stats$stars <- ''
stats$stars[stats$`p adj` < 0.05] <- '*'
stats$stars[stats$`p adj` < 0.01] <- '**'
stats$stars[stats$`p adj` < 0.001] <- '***'
stats$stars[stats$`p adj` < 0.0001] <- '****'

# make nicer table
stats3 <- stats[,c(7,9,6,4,11)]
rownames(stats3) <- seq(1, length(stats3$bug1))
colnames(stats3) <- c('species','time1','time2','p_value','significance')
View(stats3)

# all proportion plots together #####
ggarrange(plot1_rel2 + ggtitle(expression(bold("1:1 ratio of"~bolditalic("Pseudomonas")~"and"~bolditalic("Salmonella")))),
          plot2_rel2 + ggtitle(expression(bold(bolditalic("Salmonella")~"applied to established"~bolditalic("Pseudomonas")~"population"))),
          plot3_rel2 + ggtitle(expression(bold(bolditalic("Pseudomonas")~"applied to established"~bolditalic("Salmonella")~"population"))),
          nrow=1, ncol=3, common.legend = T, legend='bottom')
