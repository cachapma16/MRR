setwd("D:/Google Drive/Research/Notebook/2017/06/21/")
library(reshape2)
library(corrplot)
library(Hmisc)

dat.a.ch1 <- read_csv("20170621_gaskMXP-01a_06212017/plots/MXP-01a_netShifts_ch1step_1.csv")
dat.a.ch2 <- read_csv("20170621_gaskMXP-01a_06212017/plots/MXP-01a_netShifts_ch2step_1.csv")
dat.b.ch1 <- read_csv("20170621_gaskMXP-01b_06212017/plots/MXP-01b_netShifts_ch1step_1.csv")
dat.b.ch2 <- read_csv("20170621_gaskMXP-01b_06212017/plots/MXP-01b_netShifts_ch2step_1.csv")
dat.c.ch1 <- read_csv("20170621_gaskMXP-01c_06212017/plots/MXP-01c_netShifts_ch1step_1.csv")
dat.c.ch2 <- read_csv("20170621_gaskMXP-01c_06212017/plots/MXP-01c_netShifts_ch2step_1.csv")

hold <- rbind(dat.a.ch1, dat.a.ch2, dat.c.ch1, dat.c.ch2)
ggplot(hold, aes(y = `Net Shift`, x = Target, color = interaction(Experiment, factor(Channel)))) +
        geom_boxplot() + theme_few()
head(hold)
ggplot(hold, aes(Shift.1, Shift.2, color = Target)) + geom_point() + theme_few()
hold$Id <- hold$Ring %% 4

dat.ch1 <- dcast(dat.a.ch1,
                    Id ~ Target,
                    value.var = "Net Shift",
                    fun.aggregate = mean)

dat.ch2 <- dcast(dat.a.ch2,
                 Id ~ Target,
                 value.var = "Net Shift",
                 fun.aggregate = mean)

dat.dcast <- dcast(hold,
                 Id + Experiment + Channel ~ Target,
                 value.var = "Net Shift",
                 fun.aggregate = mean)


dat.dcast$Names <- paste(dat.dcast$Experiment, dat.dcast$Id)
rownames(dat.dcast) <- dat.dcast$Names
dat.dcast$Id <- NULL
dat.dcast$Names <- NULL
dat.dcast$Experiment <- NULL
mat <- as.matrix(dat.dcast)
heatmap(mat)
mat.cor <- rcorr(mat)
corrplot(mat.cor$r)

rcorr(dat.dcast)


plot <- ggplot(dat.dcast, 
               aes(x = `pp70S6K Thr389`, 
                   y = `pS6 Ser240/4`,
                   color = Experiment))

plot + geom_point(size = 6) + theme_few()
test <- var(mat)
test
View(test)
