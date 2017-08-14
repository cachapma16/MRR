GetName <- function(){
        # get the filename from the current working directory
        directory <- basename(getwd())
        
        # directory naming from MRR: "CHIPNAME_gaskGASK_DATE"
        # extracts and returns GASK from directory name
        name <- unlist(strsplit(directory, split = "_"))[2]
        
        # define name as global variable for use in other functions
        name <<- gsub('gask','',name) # removes "gask" from name
}

AggData <- function(loc = 'plots', filename = 'groupNames_allClusters.csv') {
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # get information of chip layout from github repository
        if (!file.exists(filename)){
                git <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/"
                url <- paste0(git, filename)
                filename <- basename(url)
                download.file(url, filename)
        }
        
        # read in recipe/chip layout
        recipe <- read_csv(filename)
        colnames(recipe)[1] <- "Target" # rename column to remove byte order mark
        targets <- recipe$Target
        
        # generate list of rings to analyze (gets all *.csv files)
        rings <- list.files(directory, pattern = ".csv", recursive = FALSE)
        idfile <- grepl("group", rings)
        removeFiles <- c("comments.csv", rings[idfile])
        rings <- rings[!rings %in% removeFiles]
        
        # add data to data frame corresponding for each ring in rings
        df <- data.frame()
        for (i in rings) {
                dat <- read_csv(i, col_names = FALSE)
                ringNum <- as.numeric(strsplit(i, "\\.")[[1]][1])
                recipeCol <- which(recipe$Ring == ringNum)
                tmp <- dat[,c(1,2)] # time and shift from raw data
                what <- tmp
                tmp$ring <- ringNum
                tmp$group <- recipe$Group[recipeCol]
                tmp$groupName <- as.character(recipe$Target[[recipeCol]])
                tmp$channel <- recipe$Channel[[recipeCol]]
                tmp$run <- name
                tmp$timePoint <- seq(1:nrow(dat))
                df <- bind_rows(df, tmp)
        }
        
        # renames columns in df
        names(df) <- c("Time", "Shift", "Ring", "Group", "Target", "Channel",
                       "Experiment", "Time Point")
        
        # creates "plots" directory if one does not exist
        if (!file.exists(loc)){dir.create(loc)}
        
        # saves aggregated data with name_allRings.csv
        write_csv(df, paste(loc, '/', name, "_allRings.csv", sep=""))
        
        # returns working directory to top level
        setwd(directory)
}

SubtractControl <- function(loc = 'plots', ch, cntl){
        #load relevant libraries
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory = getwd()
        
        # get ring data and filter by channel
        dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        if (ch != "U"){dat <- filter(dat, Channel == ch)}
        dat <- filter(dat, Target != "Ignore")
        
        # get thermal control averages
        controls <- filter(dat, Target == cntl)
        ringList <- unique(controls$Ring)
        
        # gets times from first thermal control
        times <- filter(controls, Ring == ringList[1]) %>% select(Time)
        df.controls <- data.frame(times)
        
        # create dataframe with all controls
        for (i in ringList){
                ringShift <- filter(controls, Ring == i) %>% select(Shift)
                names(ringShift) <- paste('Ring', i, sep='')
                df.controls <- cbind(df.controls, ringShift)
        }
        
        # averages thermal controls
        cols <- ncol(df.controls)
        if (cols > 2) {
                df.controls$avgControls <- rowMeans(df.controls[,c(2:cols)])
        } else {
                df.controls$avgControls <- df.controls[,c(2:cols)]
        }
        avgControls <- as.vector(df.controls$avgControls)
        
        #subtracts thermal controls from each ring
        ringNames <- unique(dat$Ring)
        for(i in ringNames){
                ringDat <- filter(dat, Ring == i) %>% select(Shift)
                ringTC <- ringDat - avgControls
                dat[dat$Ring == i, "Shift"] <- ringTC
        }
        
        write_csv(dat, paste(loc,"/", name, "_", cntl, "Control", "_ch", ch, 
                             ".csv", sep = ''))   
}

PlotRingData <- function(cntl, ch, loc = 'plots', splitPlot = FALSE){
        # loads relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # use thermally controlled data if desired
        if (cntl != "raw"){
                dat <- read_csv(paste(loc, "/", name, "_", cntl, "Control", 
                                      "_ch", ch,".csv", sep=''))
        } else if (cntl == "raw") {
                dat <- read_csv(paste(loc, "/", name, "_allRings.csv", sep=''))
                if (ch != "U") {filter(dat, Channel == ch)}
        }
        
        plot_theme <- theme_few(base_size = 16)
        
        #configure plot and legend
        plots <- ggplot(dat, 
                        aes(Time, Shift, color = Target, group = Ring)) + 
                labs(x = "Time (min)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)")),
                     color = "Target") +
                plot_theme + geom_line()
        
        # alternative plots with averaged clusters
        
        dat.2 <- dat %>% group_by(Target, `Time Point`) %>% 
                summarise_each(funs(mean, sd), c(Time, Shift))
        
        plot2 <- ggplot(dat.2, aes(Time_mean, Shift_mean, color = Target)) +
                geom_line() +
                labs(x = "Time (min)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)"))) +
                plot_theme + ggtitle(name)
        
        plot3 <- plot2 + 
                geom_ribbon(aes(ymin = Shift_mean - Shift_sd,
                                ymax = Shift_mean + Shift_sd, linetype = NA),
                            fill = "slategrey", alpha = 1/8)
        
        if (splitPlot){
                plots <- plots + facet_grid(. ~ Channel)
        }
        
        #save plot, uncomment to save
        filename <- paste0(name, "_", cntl, "Control", "_ch", ch)
        filename2 <- paste0(name, "_", cntl, "Control", "_ch", ch, "_avg")
        setwd(loc)
        ggsave(plots, file = paste0(filename, ".png"), width = 10, height = 6)
        #ggsave(plots, file = paste0(filename, ".pdf"), width = 10, height = 6)
        ggsave(plot2, file = paste0(filename2, ".png"), width = 10, height = 6)
        #ggsave(plot2, file = paste0(filename2, ".pdf"), width = 10, height = 6)
        ggsave(plot3, file = paste0(filename2, "_2.png"), width = 10, height = 6)
        #ggsave(plot3, file = paste0(filename2, "_2.pdf"), width = 10, height = 6)
        setwd(directory)
}

GetNetShifts <- function(cntl, ch, loc = 'plots', time1, time2, step = 1, target){
        # load relevant libraries
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # use thermally controlled data if desired
        if (cntl != "raw"){
                dat <- read_csv(paste0(loc, "/", name, "_", cntl, "Control", 
                                       "_ch", ch, ".csv"))
        } else {
                dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        }
        
        dat <- filter(dat, Target == target)
        
        # generate list of rings and empty dataframe to store net shift data
        ringList <- unique(dat$Ring)
        dat.rings <- data.frame()
        
        # locations for each time is determined using which, min, and abs func
        for (i in ringList){
                dat.ring <- filter(dat, Ring == i)
                time1.loc <- which.min(abs(dat.ring$Time - time1))
                time1.val <- dat.ring$Shift[time1.loc]
                time2.loc <- which.min(abs(dat.ring$Time - time2))
                time2.val <- dat.ring$Shift[time2.loc]
                ring <- i
                group <- unique(dat.ring$Group)
                target <- unique(dat.ring$Target)
                experiment <- unique(dat.ring$Experiment)
                channel <- unique(dat.ring$Channel)
                tmp <- data.frame(i, group, target, time1.val, time2.val,
                                  experiment, channel, step)
                dat.rings <- rbind(dat.rings, tmp)
        }
        
        # renames dat.rings columns
        names(dat.rings) <- c("Ring", "Group", "Target", "Shift.1", "Shift.2", 
                              "Experiment", "Channel", "Step")
        
        # calculate nat shift and create new column in dataframe
        dat.rings$`Net Shift` <- dat.rings$Shift.1 - dat.rings$Shift.2
        
        is.num <- sapply(dat.rings, is.numeric)
        dat.rings[is.num] <- lapply(dat.rings[is.num], round, 3)
        
        # save net shift data
        write_csv(dat.rings, paste0(loc, "/", name, "_netShifts_ch", ch, "_step", 
                                    step, "_", target, ".csv"))
        setwd(directory)
}

PlotNetShifts <- function(cntl, ch, loc = 'plots', step = 1, target){
        # load relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # get net shift data
        if (cntl != "raw"){
                dat <- read_csv(paste0(loc, "/", name, "_netShifts_ch", ch, 
                                       "_step", step, "_", target, ".csv"))
        } else {
                dat <- read_csv(paste0(loc,"/", name, "_netShifts_chU.csv"))
        }
        
        # configure plot and legend
        dat.nothermal <- filter(dat, Target != "thermal")
        
        plots <- ggplot(dat.nothermal, 
                        aes(Target, `Net Shift`, color = Target)) +
                geom_boxplot() +
                theme_few(base_size = 16) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position="none") +
                ylab(expression(paste("Net Shift (",Delta,"pm)"))) +
                ggtitle(name)
        
        allRings <- ggplot(dat.nothermal, 
                           aes(factor(Ring), `Net Shift`, fill = Target)) +
                geom_bar(stat = "identity") +
                theme_few(base_size = 16) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                ylab(expression(paste("Net Shift (",Delta,"pm)"))) + 
                xlab("Ring")
        
        # save plot, uncomment to save
        filename <- paste0(name, "_NetShift_ch", ch, "_step", step)
        filename2 <- paste0("IndyRings_", filename)
        setwd(loc)
        ggsave(plots, file = paste0(filename, ".png"), width = 10, height = 6)
        #ggsave(plots, file = paste0(filename, ".pdf"), width = 10, height = 6)
        ggsave(allRings, file = paste0(filename2, ".png"), 
               width = 12, height = 6)
        #ggsave(allRings, file = paste0(filename2, ".pdf"), 
        #       width = 12, height = 6)
        setwd(directory)
}

CombineNetShifts <- function(){
        library(tidyverse)
        
        directory <- getwd()
        
        netList <- grep('net', 
                        list.files(pattern = '.csv', recursive = TRUE), 
                        value = TRUE)
        
        removeFiles <- grep("Combined", netList, value = TRUE)
        netList <- netList[!netList %in% removeFiles]
        
        netShifts <- lapply(netList, read_csv)
        
        netComb <- bind_rows(netShifts)
        netComb <- filter(netComb, !grepl("thermal|BL", Target))

        write_csv(netComb, path = "netShiftsCombined.csv")
        
        setwd(directory)
}

PlotCombineNetShifts <- function(){
        library(tidyverse)
        library(ggthemes)
        library(reshape2)
        
        dat <- read_csv("netShiftsCombined.csv")
        
        dat.melt <- melt(dat, id = c("Step", "Target", "Ring", "Channel"), measure = "Net Shift")
        
        ggplot(dat.melt, aes(x = Step, y = value, color = Target)) + geom_point() + facet_wrap(~Target) + theme_few(base_size = 16)
        
        plot_theme <- theme_few(base_size = 16)
        
        plot <- ggplot(dat, aes(x = Step, y = `Net Shift`, 
                                color = Target)) +
                geom_point()+
                plot_theme + ggtitle(name) + facet_wrap(~Target) +
                scale_x_log10() +
                labs(x = "Analyte Concentration (pg/mL)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)")),
                     color = "Target")
                
        
        plot

        ggsave(plot = plot, 
               filename = paste0("netShiftsCombined_ch.png"),
               width = 8, height = 6)
}

Fit <- function(loc = 'plots', ch){
        dat <- read_csv("netShiftsCombined.csv")
        
        dat$Replicate <- as.factor((dat$Ring - 1) %% 4 + 1)
        
        fit <- list()
        
        targetList <- unique(dat$Target)
        
        plot.1 <- ggplot(dat, aes(x = Step, y = `Net Shift`, color = Target)) + geom_point()
                
                
        plot.1a <- plot.1 + scale_x_log10()
        
        # plot.2 <- ggplot(dat, aes(x = Step, y = `Net Shift`, color = Target)) + 
        #         geom_point() +
        #         geom_smooth(mapping = aes(x = Step, y = `Net Shift`), method = 'nls',
        #                     formula = y ~ A.2 + (A.1-A.2)/(1 + (x/x.0)^p),
        #                     method.args = list(start = c(A.2 =max(y),
        #                                                  A.1 = min(y),
        #                                                  x.0 = mean(y),
        #                                                  p = 1)),
        #                     se = FALSE) +
        #         labs(x = "Analyte Concentration (pg/mL)", 
        #              y = expression(paste("Relative Shift (",Delta,"pm)")),
        #              color = "Target")
        # 
        # plot.2a <- plot.2 + scale_x_log10()
        # 
        # plot.2a + geom_smooth()
        
        for(i in 1:length(targetList)) {
                tar <- targetList[i]
                print(tar)
                dat.tar <- filter(dat, Target == tar)
                y <- dat.tar$`Net Shift`
                x <- dat.tar$Step
                #fit[[i]] <- nls(y ~ SSlogis(x, Asym, xmid, scal))
                fit[[i]] <- nls(formula = y ~ A.2 + (A.1-A.2)/(1 + (x/x.0)^p),
                                start = list(A.2 =max(y),
                                             A.1 = min(y),
                                             x.0 = mean(y),
                                             p = 1))
        }
        
        capture.output(fit, file = "fitInfo.txt")
        ggsave(plot.1, filename = "plot1.png", width = 8, height = 6)
        ggsave(plot.1a, filename = "plot1a.png", width = 8, height = 6)
        ggsave(plot.2, filename = "plot2.png", width = 8, height = 6)
        
}

CheckRingQuality <- function(loc = 'plots', 
                             varLevel = 500000) {
        library(tidyverse)
        
        dat <- read_csv(paste0(loc,"/", name, "_allRings.csv"))
        
        dat.avg <- dat %>% group_by(Ring) %>%
                summarise_each(funs(var), c(Shift))
        
        ringWinners <- filter(dat.avg, Shift < varLevel) %>% select(Ring)
        ringLosers <- filter(dat.avg, Shift > varLevel) %>% select(Ring)
        
        write_csv(ringWinners, paste0(loc, '/', name, "_ringWinners.csv"))
        write_csv(ringLosers, paste0(loc, '/', name, "_ringLosers.csv"))
}

AnalyzeData <- function(loc = "plots", ch, cntl = "thermal", filename = 'groupNames_allClusters.csv') {
        GetName()
        runConc <- basename(getwd())
        ch1_MCP1 <- unlist(strsplit(runConc, split = "_"))[4]
        # ch1_IL6 <- unlist(strsplit(runConc, split = "_"))[5]
        ch2_MCP1 <- unlist(strsplit(runConc, split = "_"))[6]
        # ch2_IL6 <- unlist(strsplit(runConc, split = "_"))[8]
        AggData()
        SubtractControl(ch = 1, cntl = cntl)
        SubtractControl(ch = 2, cntl = cntl)
        PlotRingData(cntl = cntl, ch = 1, splitPlot = TRUE)
        PlotRingData(cntl = cntl, ch = 2, splitPlot = TRUE)
        GetNetShifts(cntl = cntl, ch = 1, target = "MCP-1",
                     time1 = 52, time2 = 41, step = ch1_MCP1)
        # GetNetShifts(cntl = cntl, ch = 1, target = "IL-6",
        #              time1 = 52, time2 = 41, step = ch1_IL6)
        GetNetShifts(cntl = cntl, ch = 2, target = "MCP-1",
                             time1 = 52, time2 = 41, step = ch2_MCP1)
        # GetNetShifts(cntl = cntl, ch = 2, target = "IL-6",
        #                      time1 = 52, time2 = 41, step = ch2_IL6)
        PlotNetShifts(cntl = cntl, ch = 1, step = ch1_MCP1, target = "MCP-1")
        PlotNetShifts(cntl = cntl, ch = 2, step = ch2_MCP1, target = "MCP-1")
}

AnalyzeAllData <- function() {
        foldersList <- list.dirs(recursive = FALSE)
        for (i in foldersList){
                directory <- getwd()
                setwd(i)
                AnalyzeData()
                setwd(directory)
        }
        CombineNetShifts()
        PlotCombineNetShifts()
        Fit()
}
