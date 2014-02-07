bcv = function(y) {
    # a prettier version of edgeRs BCV plot
#       Check y
        if(!is(y,"DGEList")) stop("y must be a DGEList.")

#       Compute AveLogCPM if not found in y
        A <- aveLogCPM(y$counts,offset=getOffset(y))
        #logCPM <- object$logCPM

#       Points to determine y axis limits
        disp <- sqrt(getDispersion(y))
        if(is.null(disp)) stop("No dispersions to plot")
        if(attr(disp,"type")=="common") disp <- rep(disp,length=length(A))
        df = data.frame(A, disp)

#       Make plot
        p = ggplot(df, aes(A, disp)) + geom_point(alpha=1/10)
        labels <- cols <- NULL
         if(!is.null(y$tagwise.dispersion)) {
           df$tagwise.dispersion = sqrt(y$tagwise.dispersion)
           p = p + geom_point(data=df, aes(A, tagwise.dispersion), alpha=1/10)
         }
        if(!is.null(y$common.dispersion)) {
           df$common.dispersion = sqrt(y$common.dispersion)
           p = p + geom_line(data=df, aes(A, common.dispersion), color="red")
        }
        if(!is.null(y$trended.dispersion)) {
           o <- order(A)
           df$trended.dispersion = y$trended.dispersion
#           x = df[o]
           p = p + geom_line(data=df, aes(A, sqrt(trended.dispersion)), color="blue")
}
        p + xlab("average log CPM") +
            ylab("biological coefficient of variation")
}
