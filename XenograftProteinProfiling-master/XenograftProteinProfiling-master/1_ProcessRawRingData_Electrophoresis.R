GetName <- function(){
        # get the filename from the current working directory
        directory <- basename(getwd())
        
        # directory naming from MRR: "CHIPNAME_gaskGASK_DATE"
        # extracts and returns GASK from directory name
        name <- unlist(strsplit(directory, split = "_"))[2]
        
        # define name as global variable for use in other functions
        name <<- gsub('gask','',name) # removes "gask" from name
}

AggData <- function(loc = 'plots', filename = 'groupNames_XPP.csv') {
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
                tmp$ring <- ringNum
                tmp$group <- recipe$Group[recipeCol]
                tmp$groupName <- as.character(recipe$Target[[recipeCol]])
                tmp$channel <- recipe$Channel[[recipeCol]]
                tmp$run <- name
                tmp$timePoint <- seq(1:nrow(dat))
                df <- rbind(df, tmp)
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
                plot_theme + geom_line() + ggtitle(name)
        
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
                            fill = "slategrey", alpha = 1/8) + ggtitle(name)
        
        if (splitPlot){
                plots <- plots + facet_grid(. ~ Channel)
        }
        
        #save plot, uncomment to save
        filename <- paste0(name, "_", cntl, "Control", "_ch", ch)
        filename2 <- paste0(name, "_", cntl, "Control", "_ch", ch, "_avg")
        setwd(loc)
        ggsave(plots, file = paste0(filename, ".png"), width = 10, height = 6)
        ggsave(plots, file = paste0(filename, ".pdf"), width = 10, height = 6)
        ggsave(plot2, file = paste0(filename2, ".png"), width = 10, height = 6)
        ggsave(plot2, file = paste0(filename2, ".pdf"), width = 10, height = 6)
        ggsave(plot3, file = paste0(filename2, "_2.png"), width = 10, height = 6)
        ggsave(plot3, file = paste0(filename2, "_2.pdf"), width = 10, height = 6)
        setwd(directory)
}

GetNetShifts <- function(cntl, ch, loc = 'plots', time1, time2, step = 1){
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
        
        # save net shift data
        write_csv(dat.rings, paste0(loc, "/", name, "_netShifts_ch", ch, "step_", 
                                    step, ".csv"))
        setwd(directory)
}

PlotNetShifts <- function(cntl, ch, loc = 'plots', step = 1){
        # load relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # get net shift data
        if (cntl != "raw"){
                dat <- read_csv(paste0(loc, "/", name, "_netShifts_ch", ch, 
                                       "step_", step, ".csv"))
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
                ylab(expression(paste("Net Shift (",Delta,"pm)")))
        
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
        ggsave(plots, file = paste0(filename, ".pdf"), width = 10, height = 6)
        ggsave(allRings, file = paste0(filename2, ".png"), 
               width = 12, height = 6)
        ggsave(allRings, file = paste0(filename2, ".pdf"), 
               width = 12, height = 6)
        setwd(directory)
}

CombineNetShifts <- function(loc = 'plots', ch){
        library(tidyverse)
        
        directory <- getwd()
        setwd(loc)
        
        netList <- grep('net', list.files(pattern = '.csv'), value = TRUE)
        
        if ( ch %in% c(1, 2) ) { 
                netList <- grep(paste0("ch", ch), netList, value = TRUE)
        }
        
        netShifts <- lapply(netList, read_csv)
        
        netComb <- bind_rows(netShifts)
        
        write_csv(x = netComb, 
                  path = paste0(name, "_netShiftsCombined_ch", ch, ".csv"))
        
        setwd(directory)
}

PlotCombineNetShifts <- function(loc = 'plots', ch){
        library(tidyverse)
        library(ggthemes)
        
        dat <- read_csv(paste0(loc, "/", name, "_netShiftsCombined_ch", 
                               ch, ".csv"))
        dat$Step <- factor(dat$Step)
        
        plot <- ggplot(dat, aes(x = Target, y = `Net Shift`, color = Step)) +
                geom_boxplot()
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

AnalyzeData <- function() {
        GetName()
        AggData(loc = "groupMeeting", filename = "groupNames_allClusters.csv")
        SubtractControl(loc = "groupMeeting", ch = 1, cntl = "thermal")
        #SubtractControl(ch = 2, cntl = "thermal")
        #SubtractControl(ch = "U", cntl = "thermal")
        #PlotRingData(cntl = "raw", ch = "U", splitPlot = FALSE)
        #PlotRingData(cntl = "thermal", ch = "U", splitPlot = TRUE)
        PlotRingData(loc = "groupMeeting", cntl = "thermal", ch = 1, splitPlot = FALSE)
        PlotRingData(loc = "groupMeeting", cntl = "raw", ch = 1, splitPlot = FALSE)
        #PlotRingData(cntl = "thermal", ch = 2, splitPlot = FALSE)
        #GetNetShifts(cntl = "thermal", ch = 1, time1 = 51, time2 = 39, step = 1)
        #GetNetShifts(cntl = "thermal", ch = 2, time1 = 51, time2 = 39, step = 1)
        #PlotNetShifts(cntl = "thermal", ch = 1, step = 1)
        #PlotNetShifts(cntl = "thermal", ch = 2, step = 1)
}

AnalyzeAllData <- function() {
        foldersList <- list.dirs(recursive = FALSE)
        for (i in foldersList){
                directory <- getwd()
                setwd(i)
                AnalyzeData()
                setwd(directory)
        }
}
