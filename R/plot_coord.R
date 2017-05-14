#' Visualization of significant differences in ERPs using gganimation
#'
#' The plot coord takes the results of \code{\link{plot_fa}} to create a animated GIF plot for visualizing significant differences in ERPs across locations on the scalp over time.
#' @param tests_rst the testing result from \code{\link{plot_fa}}.
#' @param frames The time point of the ERP data
#' @param show Which frames to show in the animation.
#' @param point_size The size of the points (channels) on the scalp.
#' @param high_col The color for "highly" significant time point (smaller p value).
#' @param low_col The color for "low or no significant" time point (bigger p value).
#' @param na_col The color for no specified channels on full 10-10 system.
#' @param text Logical variable. Whether to show the channel label or not.
#' @param text_size The text size of the channel label (when text = TRUE).
#' @param text_col The text color of the channel label (when text = TRUE).
#' @param circle Logical variable. Whether to show the scalp circle or not.
#' @param nose Logical variable.  Whether to show the nose or not.
#' @param cir_nose_col The color of the nose (when nose = TRUE).
#' @param loop The parameter for adjusting animations in \code{\link[animation]{ani.options}}. Whether to iterate or not (default TRUE to iterate for infinite times)
#' @param ... Further parameters in \code{\link[animation]{ani.options}} to adjust the animation (i.e.  interval = 0.5), and  \code{\link[gganimate]{gg_animate}} to save the animation.
#' @author Chi-Lin Yu <psychilinyu@gmail.com>, Ching-Fan Sheu <csheu@mail.ncku.edu.tw>
#' @seealso \code{\link{plot_tete}} \code{\link{plot_fa}}
#' @export plot_coord
#' @examples
#' data(DirectedForgetting)
#' time_pt <- seq(-200, 1000, 1)
#'
#' test_res <- plot_fa(data = DirectedForgetting, frames = time_pt,
#'                     channel = 5, subject = 1, uV = 6:1206, test = 4,
#'                     mode = "test_signal",
#'                     design = (~Subject + Condition), design0 = (~Subject), nbf = 5)
#'
#' plot_coord(tests_rst = test_res$Test_Rst,
#'            frames = time_pt,
#'            show = seq(200, 500, by = 1),
#'            loop = 1,
#'            interval = 0.1,
#'            filename = "GIFfile.gif")


plot_coord <- function(tests_rst,
                       frames,
                       show,
                       point_size = 14,
                       high_col= "firebrick",
                       low_col = "forestgreen",
                       na_col = "grey70",
                       text = T,
                       text_size = 3,
                       text_col = "black",
                       circle = T,
                       nose = T,
                       cir_nose_col="black",
                       loop = 1,
                       ... # option in "ani.options"
                       ){
        # getlegend
        g_legend<-function(a.gplot){
                tmp <- ggplot_gtable(ggplot_build(a.gplot))
                leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                legend <- tmp$grobs[[leg]]
                return(legend)}
        # some check function
        ## Log scale : only for pval and corectedpval
        #elect_coord=readRDS("Elect_Location.Rdata")
        data(Elect_Location)
        num <- length(tests_rst)
        rstlist <- list()
        for (i in 1:num){
                rstlist[[i]] <- data.frame(pval = tests_rst[[i]]$pval)
                rstlist[[i]]$correctedpval <- tests_rst[[i]]$correctedpval
                rstlist[[i]]$frame <- frames
        }
        for (i in 1:num) {
                names(rstlist)[i] <- names(tests_rst)[i]
                rstlist[[i]]$Electrode <- names(rstlist)[i]
        }
        sum_rstlist <- list()
        for (i in 1:num){
                sum_rstlist[[i]] <- rstlist[[i]][rstlist[[i]]$frame %in% show,]
        }
        sum_rstlist_fun <- list()
        sum_rstlist_fun <- sum_rstlist
        plotdta <- data.frame(sum_rstlist_fun[[1]])
        for (i in 2:num){
                plotdta <- rbind(plotdta,as.data.frame(sum_rstlist_fun[[i]]))
        }
        print("The following electrodes will not be plot :")
        print(as.character(plotdta$Electrode[!toupper(plotdta$Electrode) %in% toupper(elect_coord$Electrode)])[!duplicated(as.character(plotdta$Electrode[!toupper(plotdta$Electrode) %in% toupper(elect_coord$Electrode)]))])
        currdta <- subset(plotdta,toupper(plotdta$Electrode) %in% toupper(elect_coord$Electrode))
        currdta$Electrode <- toupper(currdta$Electrode)
        currele <- elect_coord
        currele$Electrode <- toupper(currele$Electrode)
        currdtaele <- merge(currele,currdta,all=T)
        circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
                r = diameter / 2
                tt <- seq(0,2*pi,length.out = npoints)
                xx <- center[1] + r * cos(tt)
                yy <- center[2] + r * sin(tt)
                return(data.frame(x = xx, y = yy))
        }
        diameter <- max(dist(currdtaele[,2:3]))+point_size*0.03
        datcir <- circleFun(c(mean(as.numeric(currdtaele$x)),
                              mean(as.numeric(currdtaele$y))),
                            diameter,npoints = 100)
        colnames(datcir)[1:2] <- c("cirx","ciry")
        datnose <- data.frame(nosex=c(datcir[23,1],0,-datcir[23,1]),
                              nosey=c(datcir[23,2],datcir[23,2]+0.3,datcir[23,2]))
        currdtaele <- as.data.frame(currdtaele)
        currdtaele$x <- as.numeric(currdtaele$x)
        currdtaele$y <- as.numeric(currdtaele$y)
        currdtaele$logp <- -log10(currdtaele$correctedpval)
        currdtaele = na.omit(currdtaele)
        maxlim <- max(currdtaele$logp)
        eledta <- currdtaele[,1:3]
        eledta <- eledta[!duplicated(eledta),]
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
                                axis.text.y = element_blank(),
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                #axis.line = element_blank(),
                                axis.ticks =  element_blank()
                        )
        }
        rename <- function(x){
                if (x < 10) {
                        return(name <- paste('000',i,'plot.png',sep=''))
                }
                if (x < 100 && i >= 10) {
                        return(name <- paste('00',i,'plot.png', sep=''))
                }
                if (x >= 100) {
                        return(name <- paste('0', i,'plot.png', sep=''))
                }
        }
        p <- ggplot(data = currdtaele,aes(x=x,y=y))+
                        geom_point(aes(fill=logp,frame = frame),col=cir_nose_col,size=point_size,pch=21)+
                        geom_text(data = eledta, aes(label=Electrode),
                                  size=text_size,
                                  alpha = ifelse(text==T,1,0),
                                  col=text_col)+
                        geom_path(data = datcir, aes(x = cirx,y = ciry),
                                  alpha = ifelse(circle==T,1,0),
                                  col = cir_nose_col)+
                        geom_path(data = datnose, aes(nosex,nosey),
                                  alpha = ifelse(nose==T,1,0),
                                  col = cir_nose_col)+
                        scale_fill_gradient2(high = high_col,
                                             low = low_col,
                                             midpoint = -log(0.05),
                                             na.value = na_col,
                                             limit=c(0,maxlim),
                                             guide_colorbar(title = "p-value"))+
                        theme_default()
        print(plot <- gganimate(p,...))
        return(plot)
}

