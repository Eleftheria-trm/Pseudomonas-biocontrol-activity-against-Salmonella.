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
require(ggpubr)
#install.packages('ggpubr')
require(scales)


# data import & prep
setwd("~/Desktop/EleftheriaGitHub")
data <- read.csv('Figure1B_data.csv', stringsAsFactors = F)

# calculate proportion of cells of each colour in the biofilm
for (i in seq(1, length(data$strain))){
  data$rel[i] <- data$cfu.ml[i] / sum(data$cfu.ml[data$strain == data$strain[i] & data$time == data$time[i] &
                                                    data$bio.rep == data$bio.rep[i] & data$tech.rep == data$tech.rep[i]])
}

# calculate change in proportion from time point 0
for (i in seq(1, length(data$strain))){
  data$change[i] <- data$rel[i] - data$rel[data$strain == data$strain[i] & data$time == 0 & data$colour == data$colour[i] &
                                             data$bio.rep == data$bio.rep[i] & data$tech.rep == data$tech.rep[i]]
}

# calculate mean, sd, ci for change
for (i in seq(1, length(data$strain))){
  data$mean.chg[i] <- mean(data$change[data$strain == data$strain[i] & data$time == data$time[i] & data$colour == data$colour[i]])
  data$sd.chg[i] <- sqrt(var(data$change[data$strain == data$strain[i] & data$time == data$time[i] & data$colour == data$colour[i]]))
  data$ci.chg[i] <- qnorm(0.95) * data$sd.chg[i] / sqrt(length(data$change[data$strain == data$strain[i] & 
                                                                             data$time == data$time[i] & 
                                                                             data$colour == data$colour[i]]))
}

# replace dash with underscore in strain names
for (i in seq(1, length(data$strain))){
  data$strain[i] <- gsub('-','_',data$strain[i])
}

# run stats = anova
model <- aov(change ~ strain * as.factor(time), data=subset(data, colour=='blue'))
summary(model)
stats.chg.a <- TukeyHSD(model)

# put stats into a data frame and subset only the stuff we're interested in
stats.chg.a <- stats.chg.a$`strain:as.factor(time)`
stats.chg.a <- as.data.frame(stats.chg.a)
stats.chg.a$strain1 <- gsub( "(.*)-(.*)", "\\1", rownames(stats.chg.a))
stats.chg.a$time1 <- gsub( "(.*):(.*)", "\\2", stats.chg.a$strain1)
stats.chg.a$strain1 <- gsub( "(.*):(.*)", "\\1", stats.chg.a$strain1)
stats.chg.a$strain2 <- gsub( "(.*)-(.*)", "\\2", rownames(stats.chg.a))
stats.chg.a$time2 <- gsub( "(.*):(.*)", "\\2", stats.chg.a$strain2)
stats.chg.a$strain2 <- gsub( "(.*):(.*)", "\\1", stats.chg.a$strain2)
stats.chg.a <- subset(stats.chg.a, strain1 == strain2)
stats.chg.a <- subset(stats.chg.a, time2 == 0)

# add significance stars 
stats.chg.a$stars <- ''
stats.chg.a$stars[stats.chg.a$`p adj` < 0.05] <- '*'
stats.chg.a$stars[stats.chg.a$`p adj` < 0.01] <- '**'
stats.chg.a$stars[stats.chg.a$`p adj` < 0.001] <- '***'
stats.chg.a$stars[stats.chg.a$`p adj` < 0.0001] <- '****'

# plot 
plot.chg.a <- ggplot(data=subset(data, colour == 'blue'), aes()) +
  geom_hline(yintercept = 0, colour = "#cbcbcb")+
  geom_point(aes(x=strain, y=change, colour=as.factor(time)), position=position_dodge(width=.6)) +
  geom_errorbar(aes(x=strain, ymin=mean.chg-ci.chg, ymax=mean.chg+ci.chg, colour=as.factor(time)), position=position_dodge(width=.6), width=.5) +
  geom_errorbar(aes(x=strain, ymin=mean.chg-ci.chg, ymax=mean.chg, colour=as.factor(time)), position=position_dodge(width=.6), width=.25) +
  scale_y_continuous(name=expression(atop('Change in proportion of'~italic('Salmonella')~'CFU/mL',
                                          'in competition with'~italic('Pseudomonas')~'relative to time 0')),
                     limits=c(-1,1), breaks=seq(-1,1,.2)) +
  scale_x_discrete(name='', limits=c("none","PK0307_4","PK0311_1","SBW25"),
                   labels=c(expression(italic("Salmonella")~"ctrl"),"PK0307-4","PK0311-1","SBW25")) +
  annotate('text', x=stats.chg.a$strain1[stats.chg.a$time1 == 6], y=1, label=stats.chg.a$stars[stats.chg.a$time1 == 6], colour="#7CAE00", size=8)+
  annotate('text', x=stats.chg.a$strain1[stats.chg.a$time1 == 24], y=.95, label=stats.chg.a$stars[stats.chg.a$time1 == 24], colour="#00BFC4", size=8)+
  annotate('text', x=stats.chg.a$strain1[stats.chg.a$time1 == 48], y=.9, label=stats.chg.a$stars[stats.chg.a$time1 == 48], colour="#C77CFF", size=8)+
  scale_colour_discrete(name='Time (h)') +
  bbc_style() +
  theme(plot.margin = unit(c(.5,0,0,0),'cm'),
        axis.title=element_text(size=14),
        axis.text.y=element_text(margin = margin(2,2,2,2), size=14),
        axis.text.x=element_text(size=14, angle=45, hjust=1),
        axis.line = element_line(color = "#cbcbcb"),
        axis.ticks = element_line(color = "#cbcbcb"),
        panel.grid.major.y = element_blank(),
        legend.box='horizontal',
        legend.title = element_text(size=14),
        legend.direction='horizontal',
        legend.position='bottom',
        legend.text=element_text(size=14, margin = margin(r = 30, unit = "pt"))) 
plot.chg.a
