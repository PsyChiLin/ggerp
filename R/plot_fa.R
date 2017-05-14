#' Visualization of the results using adaptive factor-adjustement for multiple testing of ERP data
#'
#' the function plot fa performs factor-adjusted multiple testing of ERP curves in the linear model framework (Causeur et al, 2012) and display results the same way as that of \code{\link{plot_tete}}.
#' @export plot_fa
#' @param data a data frame with ERP data
#' @param frames The time point of the ERP data
#' @param uV The corresponding column indices within the input data frame for the ERP amplitudes
#' @param subject The corresponding column index for the subject variables (factor)
#' @param channel The corresponding column index for the channel variables (factor)
#' @param test The corresponding column index for the variables you want to compare. It could be a factor (i.e, Condition) or numeric variable (i.e, Score).
#' @param mode The options ("test","test_signal") of the functions control the kind of plots to be made.
#' @param curve.col The color for the curves. Numbers must match. (e.g. two conditions need two colors)
#' @param labs The labs of the plot
#' @param ggtheme The theme setting for ggplot2
#' @param order Reordering the plots by the similarities between the curves using TSclust package.
#' @param scalp Logical variable. Plot the curve on the scalp location (follow the coord.mat option).
#' @param coord.mat Read users' own coordinate matrix. Default is the full 10/10 system.
#' @param ylim The limits of the plot on Y axis. The setting is same for ggplot2.
#' @param significant.col The color used to indicated the significant time point.
#' @param significant.alpha The alpha of the significant.col that used to indicated the significant time point.
#' @param design Design matrix of the full model for the relationship between the ERP and the experimental variables. Typically the output of the function model.matrix (see the R package Causeur, David, Ching-Fan Sheu, and Mei-Chen Chu. "Significance analysis of ERP data." (2014).)
#' @param design0 Design matrix of the null model. Typically a submodel of the full model, obtained by removing columns from design. Default is NULL, corresponding to the model with no covariates. (see the R package Causeur, David, Ching-Fan Sheu, and Mei-Chen Chu. "Significance analysis of ERP data." (2014).)
#' @param ... Other parameters in \code{\link[ERP]{erpfatest}}.
#' @return Plot A "ggplot2" plot will automatically generate.
#' @return Test_Rst A same testing results of \code{\link[ERP]{erpfatest}}.
#' @author Chi-Lin Yu <psychilinyu@gmail.com>, Ching-Fan Sheu <csheu@mail.ncku.edu.tw>
#' @seealso \code{\link{plot_tete}} \code{\link{plot_coord}}  \code{\link[ERP]{erpfatest}}
#' @examples
#' data(DirectedForgetting)
#' time_pt <- seq(-200, 1000, 1)
#' 
#' dta_c <- DirectedForgetting %>%
#'          filter(Channel %in% c("FZ","CZ","PZ")) %>%
#'          droplevels()
#'
#' plot_fa(data = dta_c, 
#'         frames = time_pt,
#'         channel = 5, subject = 1, uV = 6:1206, test = 4,
#'         mode = "test_signal",
#'         design = (~Subject + Condition), design0 = (~Subject), nbf = 5)

plot_fa <- function(data,  # The ERP or NIRS data
                      frames = NULL, # The time point of the data
                      uV , # The column of "Data"
                      subject = NULL, # The column of "Subject"
                      channel = NULL,  # The column of "Channel"
                      test = NULL, # The column of the variable that you want to compare (e.g. "Condition")
                      mode = c("test","test_signal"), # Modes of the possible plot
                      curve.col = NULL, # the color for the curve (e.g. two conditions need two colors)
                      labs = list(x = "Time (ms)", y = "Amplitude (microvolt)"),# Labs of the plot
                      ggtheme = NULL, # set your own theme
                      order = F, # use significant time points to order the plot
                      scalp = FALSE, # Do you want to plot on the scalp ?
                      coord.mat = NULL, # read your own coordinate matrix, default will be the full 10/10 system.
                      ylim = c(-20,20), # y limits of the plot
                      significant.col = "darkgoldenrod", # when mode ==  test | mode = test_signal, the color of the significant window
                      significant.alpha = 0.1, # when mode ==  test | mode = test_signal, the alpha value of the significant window
                      design = NULL,  # the argument for erpfatest
                      design0 = NULL,...){
        # g_legend
        g_legend<-function(a.gplot){
                tmp <- ggplot_gtable(ggplot_build(a.gplot))
                leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                legend <- tmp$grobs[[leg]]
                return(legend)}
        # Some Check Functions (need More)
        options(warn=-1)
        dta <- data
        dta[,channel] <- as.factor(dta[,channel])
        dta[,subject] <- as.factor(dta[,subject])
        # Some Check Functions (need More)
        if (scalp == TRUE & order == TRUE){message("The function will not order the curve when scalp = TRUE")}
        if (is.null(frames)){frames  = 1:length(uV)}
        if (length(levels(dta[,channel])) == 1 & scalp == T){
                message("A single channel data could not be placed on scalp ! change to scalp = FALSE")
                scalp = FALSE
        }
        if (length(levels(dta[,channel])) == 1 & order == T){
                message("A single channel data could not be ordered ! change to order = FALSE ")
                order = FALSE
        }
        if (!mode %in% c("test","test_signal")){stop("mode should be 'test',or 'test_signal'")}
        if (length(mode) != 1) {stop("You should select only one mode of plot !")}
        if (!is.null(test)){
                if (!class(dta[,test]) %in% c("numeric","factor")){stop("test should be numeric or factor" )}
                if (class(dta[,test]) == "factor"){
                        if (is.null(curve.col)){
                                if (length(levels(as.factor(dta[,test]))) == 1){curve.col = c("chartreuse4")}
                                if (length(levels(as.factor(dta[,test]))) == 2){curve.col = c("chartreuse4", "firebrick")}
                                if (length(levels(as.factor(dta[,test]))) == 3){curve.col = c("chartreuse4", "firebrick","blue")}
                                # add some colors
                                message("Using Default color as the curve color")
                        }
                        if (length((levels(dta[,test]))) != length(curve.col)){stop("curve.col should be equal to levels of test !")}
                }
                if (class(dta[,test]) == "numeric"){if (is.null(curve.col)) {curve.col = c("chartreuse4")}}
        }
        if (is.null(test)){curve.col = c("chartreuse4")}
        if (is.null(test)|length(levels(dta[,test])) == 1) {testmode = "single"}
        if (length(levels(dta[,test])) > 1) {testmode = "typical"}
        if (class(dta[,test]) == "numeric") {testmode = "numer"}
        if (testmode == "numer"){
                message("The test is a numeric variable. Correlation plot will be produced.")
                mode = "Rsquared"
        }
        curve.fun <- function(x, d){return(mean(x[d]))}
        # The test function
        chan_test <- function(data,uV = NULL,channel = NULL,design_model,design0_model=NULL,...){
                dta <- data
                levelnum <- length(levels(dta[,channel]))
                if (levelnum == 1) {
                        design <- model.matrix(design_model,data=data)
                        if (is.null(design0_model)==F){
                                design0 <- model.matrix(design0_model,data=data)
                        } else {design0 = NULL}
                        test_list <- erpfatest(dta[,uV],design,design0,...)
                } else {
                        test_list=list()
                        dta_list=list()
                        for (i in 1:levelnum) {
                                dta_list[[i]] <- subset(dta,dta[,channel]==(levels(dta[,channel])[i]))
                                design <- model.matrix(design_model,data=dta_list[[i]])
                                if (is.null(design0_model)==F){
                                        design0 <- model.matrix(design0_model,data=dta_list[[i]])
                                } else {design0 = NULL}
                                test_list[[i]] <- erpfatest(dta_list[[i]][,uV],design,design0,...)
                                names(test_list)[i] <- levels(dta[,channel])[i]
                        }
                }
                return(test_list)
        }
        # Produce the test results called "tests_rst", when mode  == "test" or "test_signal"
        tests_rst <- chan_test(dta,uV,channel,design,design0,...)
        # The summarize function
        data_summarize <- function(data,uV,summary.var,fun=mean,...){
                dta <- data
                options(warn=-1)
                agglength <- length(summary.var)
                aggvar_list <- list(dta[,summary.var[1]])
                if (agglength > 1){
                        for (i in 2:agglength ){
                                aggvar_list <- append(aggvar_list,list(dta[,summary.var[2]]))
                        }
                }
                aggdata <- aggregate(dta[,uV],by=aggvar_list,
                                     fun,...)
                aggdata <- aggdata[,1:(agglength+length(uV))]
                for (i in 1: agglength){
                        colnames(aggdata)[i] <- colnames(dta)[summary.var[i]]
                }
                rownames(aggdata) <- 1:dim(aggdata)[1]
                return(aggdata)
        }
        # Draw the plot (not on 10/10 system)
        if (is.null(ggtheme)){
                if (scalp == TRUE){
                        theme_default <- function(base_size = 10, base_family = ""){
                                theme_bw(base_size = base_size, base_family = base_family) %+replace%
                                        theme(
                                                panel.border     = element_blank(),
                                                axis.line        = element_line(colour = "black"),
                                                panel.grid.major = element_line(),
                                                panel.grid.major.x = element_blank(),
                                                panel.grid.major.y = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.grid.minor.x = element_blank(),
                                                panel.grid.minor.y = element_blank(),
                                                strip.background = element_rect(color = "white", size = 0.2),
                                                legend.key       = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.title.x = element_blank(),
                                                axis.text.x = element_blank(),
                                                axis.line = element_blank(),
                                                axis.ticks =  element_blank()
                                        )
                        }
                }
                if (scalp == FALSE){
                        theme_default <- function(base_size = 10, base_family = ""){
                                theme_bw(base_size = base_size, base_family = base_family) %+replace%
                                        theme( strip.background = element_blank()
                                        )
                        }
                }
        } else {theme_default <- ggtheme}
        if (scalp == F){
                if (mode == "test") {
                        if (length(levels(dta[,channel]))==1){
                                data_test <- data.frame(signal=as.numeric(tests_rst$signal))
                                data_test$frames = frames
                                data_test$significant <-
                                        ifelse(data_test $frames %in% data_test $frames[tests_rst$significant],
                                               "sig","non-sig")
                                data_test$sign_frames <-
                                        ifelse(data_test $frames %in% data_test$frames[tests_rst$significant],
                                               data_test $frames,NA)
                                data_test$group <- rep(0,length(data_test $signal))
                                data_test$r2 <- tests_rst$r2
                                plot <- ggplot(data_test ,aes(x=frames,y=signal,group=group))+
                                                geom_vline(data=data_test ,
                                                           aes(xintercept = sign_frames),
                                                           col=significant.col,
                                                           alpha=significant.alpha)+
                                                geom_line(col = curve.col[1])+
                                                labs(labs)+
                                                theme_default()+
                                                ylim(ylim)
                        }
                        else {
                                listlen <- length(tests_rst)
                                data_list <- list()
                                for (k in 1:listlen){
                                        data_list[[k]] <- data.frame(signal=as.numeric(tests_rst[[k]]$signal))
                                        data_list[[k]]$frames = frames
                                        data_list[[k]]$significant<-ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],"sig","non-sig")
                                        data_list[[k]]$sign_frames <- ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],data_list[[k]]$frames,NA)
                                        data_list[[k]]$group <- rep(0,length(data_list[[k]]$signal))
                                        data_list[[k]]$r2 <- tests_rst[[k]]$r2
                                        data_list[[k]]$Channel <- names(tests_rst)[k]
                                }
                                data_plot <- data_list[[1]]
                                for (j in 2 : listlen){
                                        data_plot <- rbind(data_plot,data_list[[j]])
                                }
                                if (order == TRUE) {
                                        orderlist <- list()
                                        for (i in 1:length(tests_rst)){
                                                orderlist[[i]] <- length(tests_rst[[i]]$significant)
                                        }
                                        names(orderlist) <- names(tests_rst)
                                        data_plot$Channel <- factor(data_plot$Channel,
                                                                 levels = c(colnames(sort(as.data.frame(do.call(cbind,orderlist)),decreasing=T))))
                                }
                                plot <- ggplot(data_plot,aes(x=frames,y=signal,group=group))+
                                        geom_vline(data=data_plot,
                                                   aes(xintercept = sign_frames),
                                                   col=significant.col,
                                                   alpha=significant.alpha)+
                                        geom_line(col = curve.col[1])+
                                        labs(labs)+
                                        facet_wrap(~Channel)+
                                        theme_default()+
                                        theme(legend.position="none")+
                                        ylim(ylim)
                        }
                }
                if (mode == "test_signal"){
                        if (length(levels(dta[,channel]))==1){
                                if (testmode == "typical") {
                                        values <- curve.col[1]
                                        for (i in 2:length(levels(as.factor(dta[,test])))){
                                                values <- c(values,curve.col[i])
                                        }
                                        data_sign <- data.frame(frames=frames)
                                        data_sign$significant <-
                                                ifelse(data_sign$frames %in% data_sign$frames[tests_rst$significant],"sig","non-sig")
                                        data_sign$sign_frames <-
                                                ifelse(data_sign$frames %in% data_sign$frames[tests_rst$significant],data_sign$frames,NA)
                                        data_sign$group <- rep(0,length(frames))
                                        data_plot <- data_summarize(dta,uV ,summary.var = c(channel,test),fun = curve.fun)
                                        data_plot <- melt(data_plot)
                                        colnames(data_plot)[1:2] <- c("Channel","Condition")
                                        data_plot <- data_plot[order(data_plot$Channel,data_plot$Condition,data_plot$variable),]
                                        data_plot$frames <- rep(frames,dim(data_plot)[1]/length(frames))
                                        plot <- ggplot(data = data_plot,aes(x = frames, y = value))+
                                                geom_line(aes(col=Condition))+
                                                geom_vline(data = data_sign,aes(xintercept = sign_frames),
                                                           col = significant.col, alpha =  significant.alpha)+
                                                facet_wrap(~Channel)+
                                                labs(labs)+
                                                ylim(ylim)+
                                                scale_color_manual(values = values,name = colnames(dta)[test])+
                                                theme_default()
                                }
                                #if (testmode == "numer")
                        }
                        else {
                                listlen <- length(tests_rst)
                                test_list = list()
                                for (i in 1: listlen){
                                        test_list[[i]] <- data.frame(frames=frames)
                                        test_list[[i]]$significant <- ifelse(test_list[[i]]$frames %in% test_list[[i]]$frames[tests_rst[[i]]$significant],"sig","non-sig")
                                        test_list[[i]]$sign_frames <- ifelse(test_list[[i]]$frames %in% test_list[[i]]$frames[tests_rst[[i]]$significant],test_list[[i]]$frames,NA)
                                        test_list[[i]]$group <- rep(0,length(frames))
                                        test_list[[i]]$Channel <- names(tests_rst)[i]
                                }
                                test_plot <- test_list[[1]]
                                for (k in 2:listlen){
                                        test_plot <- rbind(test_plot,test_list[[k]])
                                }
                                if (testmode == "typical") {
                                        dta[,test] <- as.factor(dta[,test])
                                        values <- curve.col[1]
                                        for (i in 2:length(levels(as.factor(dta[,test])))){
                                                values <- c(values,curve.col[i])
                                        }
                                        data_fun <- data_summarize(dta,uV = uV,summary.var = c(channel,test),fun=curve.fun)
                                        data_fun_long <- melt(data_fun,id=c(colnames(dta)[channel],colnames(dta)[test]))
                                        data_fun_long <- data_fun_long[order(data_fun_long[,1],data_fun_long[,2],data_fun_long[,3]),]
                                        colnames(data_fun_long)[c(1,2,4)] <- c("Channel","Condition","FUN")
                                        data_fun_long$frames <- c(rep(frames,(dim(data_fun_long)[1]/length(frames))))
                                        data_plot <- merge(data_fun_long,test_plot,by =c("Channel","frames"))
                                        if (order == TRUE) {
                                                orderlist <- list()
                                                for (i in 1:length(tests_rst)){
                                                        orderlist[[i]] <- length(tests_rst[[i]]$significant)
                                                }
                                                names(orderlist) <- names(tests_rst)
                                                data_plot$Channel <- factor(data_plot$Channel,
                                                                            levels = c(colnames(sort(as.data.frame(do.call(cbind,orderlist)),decreasing=T))))
                                        }
                                        plot <- ggplot(data_plot,aes(x=frames,group=Condition))+
                                                geom_vline(aes(xintercept = sign_frames),
                                                           col=significant.col,alpha=significant.alpha)+
                                                geom_line(aes(y = FUN,col=Condition))+
                                                labs(labs)+# Need some changes ?
                                                facet_wrap(~Channel)+
                                                labs(labs)+
                                                ylim(ylim)+
                                                scale_color_manual(values = values,name = colnames(dta)[test])+
                                                theme_default()
                                }
                                #if (testmode == "numer"){}
                        }
                }
                if (mode == "Rsquared"){
                        if (length(levels(dta[,channel]))==1){
                                data_test <- data.frame(signal=as.numeric(tests_rst$signal))
                                data_test$frames = frames
                                data_test$significant <-
                                        ifelse(data_test $frames %in% data_test $frames[tests_rst$significant],
                                               "sig","non-sig")
                                data_test$sign_frames <-
                                        ifelse(data_test $frames %in% data_test$frames[tests_rst$significant],
                                               data_test $frames,NA)
                                data_test$group <- rep(0,length(data_test $signal))
                                data_test$r2 <- tests_rst$r2
                                data_cor <- data.frame(Channel = levels(dta[,channel]),
                                                            frames = frames,
                                                            Correlation = melt(cor(dta[,test], dta[,uV]))[,3])
                                data_test <- suppressMessages(full_join(data_cor,data_test))
                                #data_test$Correlation <- sign(data_test$signal)*sqrt(data_test$r2)
                                plot <- ggplot(data_test ,aes(x=frames,y=Correlation,group=group))+
                                        geom_vline(data=data_test ,
                                                   aes(xintercept = sign_frames),
                                                   col=significant.col,
                                                   alpha=significant.alpha)+
                                        geom_line(col = curve.col[1])+
                                        labs(labs)+
                                        theme_default()+
                                        ylim(ylim)
                        }
                        else {
                                listlen <- length(tests_rst)
                                data_list <- list()
                                for (k in 1:listlen){
                                        data_list[[k]] <- data.frame(signal=as.numeric(tests_rst[[k]]$signal))
                                        data_list[[k]]$frames = frames
                                        data_list[[k]]$significant<-ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],"sig","non-sig")
                                        data_list[[k]]$sign_frames <- ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],data_list[[k]]$frames,NA)
                                        data_list[[k]]$group <- rep(0,length(data_list[[k]]$signal))
                                        data_list[[k]]$r2 <- tests_rst[[k]]$r2
                                        data_list[[k]]$Channel <- names(tests_rst)[k]
                                }
                                data_plot <- data_list[[1]]
                                for (j in 2 : listlen){
                                        data_plot <- rbind(data_plot,data_list[[j]])
                                }
                                if (order == TRUE) {
                                        orderlist <- list()
                                        for (i in 1:length(tests_rst)){
                                                orderlist[[i]] <- length(tests_rst[[i]]$significant)
                                        }
                                        names(orderlist) <- names(tests_rst)
                                        data_plot$Channel <- factor(data_plot$Channel,
                                                                    levels = c(colnames(sort(as.data.frame(do.call(cbind,orderlist)),decreasing=T))))
                                }
                                #data_plot$Correlation <- sign(data_plot$signal)*sqrt(data_plot$r2)
                                numerlist <- list()
                                datalist <- list()
                                for (i in 1:length(levels(dta[,channel]))){
                                        numerlist[[i]] <- filter(dta,dta[,channel] == levels(dta[,channel])[i])
                                        datalist[[i]] <- data.frame(Channel = levels(dta[,channel])[i],
                                                                    frames = frames,
                                                                    Correlation = melt(cor(numerlist[[i]][,test], numerlist[[i]][,uV]))[,3])
                                }
                                data_cor <- do.call("rbind",datalist)
                                data_plot <- suppressMessages(full_join(data_cor,data_plot))
                                plot <- ggplot(data_plot,aes(x=frames,y=Correlation,group=group))+
                                        geom_vline(data=data_plot,
                                                   aes(xintercept = sign_frames),
                                                   col=significant.col,
                                                   alpha=significant.alpha)+
                                        geom_line(col = curve.col[1])+
                                        labs(labs)+
                                        facet_wrap(~Channel)+
                                        theme_default()+
                                        theme(legend.position="none")+
                                        ylim(ylim)
                        }
                }
                print(plot)
        }
        # Draw the plot ( on 10/10 system)
        if (scalp == T){ # do.call  grid.arrange problem
                message("Label of X-axis will be suppressed.")
                get_legend<-function(myggplot){
                        tmp <- ggplot_gtable(ggplot_build(myggplot))
                        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                        legend <- tmp$grobs[[leg]]
                        return(legend)
                }
                if (is.null(coord.mat)){
                        message("Use default coord.mat")
                        erplay <- rbind(c(NA,NA,NA,NA,"Fp1","Fpz","Fp2",NA,NA,NA,NA),
                                        c(NA,NA,NA,"AF7","AF3","AFz","AF4","AF8",NA,NA,NA),
                                        c(NA,"F7","F5","F3","F1","Fz","F2","F4","F6","F8",NA),
                                        c("FT9","FT7","FC5","FC3","FC1","FCz","FC2","FC4","FC6","FT8","FT10"),
                                        c("T9","T7","C5","C3","C1","CZ","C2","C4","C6","T8","T10"),
                                        c("TP9","TP7","CP5","CP3","CP1","CPz","CP2","CP4","CP6","TP8","TP10"),
                                        c("P9","P7","P5","P3","P1","PZ","P2","P4","P6","P8","P10"),
                                        c(NA,NA,"PO9","PO7","PO3","POz","PO4","PO8","PO10",NA,NA),
                                        c(NA,NA,NA,NA,"O1","Oz","O2",NA,NA,NA,NA),
                                        c(NA,NA,NA,NA,"I1","Iz","I2",NA,NA,NA,NA))
                } else {erplay <- coord.mat}
                plotlist <- list()
                if (mode == "test") {
                        listlen <- length(tests_rst)
                        data_list <- list()
                        for (k in 1:listlen){
                                data_list[[k]] <- data.frame(signal=as.numeric(tests_rst[[k]]$signal))
                                data_list[[k]]$frames = frames
                                data_list[[k]]$significant<-ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],"sig","non-sig")
                                data_list[[k]]$sign_frames <- ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],data_list[[k]]$frames,NA)
                                data_list[[k]]$group <- rep(0,length(data_list[[k]]$signal))
                                data_list[[k]]$r2 <- tests_rst[[k]]$r2
                                data_list[[k]]$Channel <- names(tests_rst)[k]
                        }
                        data_plot <- data_list[[1]]
                        for (j in 2 : listlen){
                                data_plot <- rbind(data_plot,data_list[[j]])
                        }
                        data_plot$Channel <- as.factor(data_plot$Channel)
                                for (i in 1:length(levels(data_plot$Channel))){
                                        dta_sub <- subset(data_plot,data_plot$Channel == levels(data_plot$Channel)[i])
                                        plotlist[[i]] <- ggplot(dta_sub,aes(x=frames,y=signal,group=group))+
                                                geom_vline(data=dta_sub,
                                                           aes(xintercept = sign_frames),
                                                           col=significant.col,
                                                           alpha=significant.alpha)+
                                                geom_line(col = curve.col[1])+
                                                theme_default()+
                                                ylim(ylim)+
                                                labs(labs(list(title=levels(data_plot$Channel)[i])))+
                                                theme(legend.position="none")
                                        names(plotlist)[i] <-levels(data_plot$Channel)[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                                string <-as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        print(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        }
                if (mode == "test_signal"){
                        listlen <- length(tests_rst)
                        test_list = list()
                        for (i in 1: listlen){
                                test_list[[i]] <- data.frame(frames=frames)
                                test_list[[i]]$significant <- ifelse(test_list[[i]]$frames %in% test_list[[i]]$frames[tests_rst[[i]]$significant],"sig","non-sig")
                                test_list[[i]]$sign_frames <- ifelse(test_list[[i]]$frames %in% test_list[[i]]$frames[tests_rst[[i]]$significant],test_list[[i]]$frames,NA)
                                test_list[[i]]$group <- rep(0,length(frames))
                                test_list[[i]]$Channel <- names(tests_rst)[i]
                        }
                        test_plot <- test_list[[1]]
                        for (k in 2:listlen){
                                test_plot <- rbind(test_plot,test_list[[k]])
                        }
                        if (testmode == "typical")  {
                                values <- curve.col[1]
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                data_fun <- data_summarize(dta,uV = uV,summary.var = c(channel,test),fun=curve.fun)
                                data_fun_long <- melt(data_fun,id=c(colnames(dta)[channel],colnames(dta)[test]))
                                data_fun_long <- data_fun_long[order(data_fun_long[,1],data_fun_long[,2],data_fun_long[,3]),]
                                colnames(data_fun_long)[c(1,2,4)] <- c("Channel","Condition","FUN")
                                data_fun_long$frames <- c(rep(frames,(dim(data_fun_long)[1]/length(frames))))
                                data_plot <- merge(data_fun_long,test_plot,by =c("Channel","frames"))
                                dta_sub <- subset(data_plot,data_plot[,1] == levels(data_plot[,1])[1])
                                flegend <- ggplot(data = dta_sub,aes(x = frames, y = FUN))+
                                        geom_line(aes(col = Condition))+
                                        geom_vline(aes(xintercept = sign_frames),
                                                   col=significant.col,alpha=significant.alpha)+
                                        labs(labs(list(title=levels(data_plot[,1])[1])))+
                                        ylim(ylim)+
                                        scale_color_manual(values = values, name = colnames(dta)[test])+
                                        theme_default()
                                legend <- get_legend(flegend)
                                for (i in 1:length(levels(data_plot[,1]))){
                                        dta_sub <- subset(data_plot,data_plot[,1] == levels(data_plot[,1])[i])
                                        plotlist[[i]] <- ggplot(data = dta_sub,aes(x = frames, y = FUN))+
                                                geom_vline(aes(xintercept = sign_frames),
                                                           col=significant.col,alpha=significant.alpha)+
                                                geom_line(aes(col = Condition))+
                                                labs(labs(list(title=levels(data_plot[,1])[i])))+
                                                ylim(ylim)+
                                                scale_color_manual(values = values, name = colnames(dta)[test])+
                                                theme_default()+
                                                theme(legend.position = "none")
                                        names(plotlist)[i] <-levels(data_plot[,1])[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                                string <- as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                printlist[[length(printlist)+1]] <- legend
                                erplaynum[1,ncol(erplaynum)] <- length(printlist)
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        }
                        #if (testmode == "numer")
                }
                if (mode == "Rsquared"){
                        listlen <- length(tests_rst)
                        data_list <- list()
                        for (k in 1:listlen){
                                data_list[[k]] <- data.frame(signal=as.numeric(tests_rst[[k]]$signal))
                                data_list[[k]]$frames = frames
                                data_list[[k]]$significant<-ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],"sig","non-sig")
                                data_list[[k]]$sign_frames <- ifelse(data_list[[k]]$frames %in% data_list[[k]]$frames[tests_rst[[k]]$significant],data_list[[k]]$frames,NA)
                                data_list[[k]]$group <- rep(0,length(data_list[[k]]$signal))
                                data_list[[k]]$r2 <- tests_rst[[k]]$r2
                                data_list[[k]]$Channel <- names(tests_rst)[k]
                        }
                        data_plot <- data_list[[1]]
                        for (j in 2 : listlen){
                                data_plot <- rbind(data_plot,data_list[[j]])
                        }
                        data_plot$Channel <- as.factor(data_plot$Channel)
                        numerlist <- list()
                        datalist <- list()
                        for (i in 1:length(levels(dta[,channel]))){
                                numerlist[[i]] <- filter(dta,dta[,channel] == levels(dta[,channel])[i])
                                datalist[[i]] <- data.frame(Channel = levels(dta[,channel])[i],
                                                            frames = frames,
                                                            Correlation = melt(cor(numerlist[[i]][,test], numerlist[[i]][,uV]))[,3])
                        }
                        data_cor <- do.call("rbind",datalist)
                        data_plot <- suppressMessages(full_join(data_cor,data_plot))
                        #data_plot$Correlation <- sign(data_plot$signal)*sqrt(data_plot$r2)
                        for (i in 1:length(levels(data_plot$Channel))){
                                dta_sub <- subset(data_plot,data_plot$Channel == levels(data_plot$Channel)[i])
                                plotlist[[i]] <- ggplot(dta_sub,aes(x=frames,y=Correlation,group=group))+
                                        geom_vline(data=dta_sub,
                                                   aes(xintercept = sign_frames),
                                                   col=significant.col,
                                                   alpha=significant.alpha)+
                                        geom_line(col = curve.col[1])+
                                        theme_default()+
                                        ylim(ylim)+
                                        labs(labs(list(title=levels(data_plot$Channel)[i])))+
                                        theme(legend.position="none")
                                names(plotlist)[i] <-levels(data_plot$Channel)[i]
                        }
                        erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                        string <-as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                        if (!identical(string, character(0))){
                                message("The following electrodes will not be plotted :")
                                print(string)
                        }
                        printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                        erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                        for (i in 1:nrow(erplay)){
                                for (j in 1:ncol(erplay)){
                                        if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                        }
                                }
                        }
                        plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                }
        }
        re_list <- list()
        re_list[[1]] <- plot
        re_list[[2]] <- tests_rst
        names(re_list)[1:2] <- c("Plot","Test_Rst")
        return(re_list)
}


