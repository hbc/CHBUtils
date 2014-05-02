#' Plot MDS
#'
#' @param counts  matrix
#' @param condition  vector indicating the group of each sample
#' @param k  number of dimension
#' @param d distance method (euclidian, cor)
#' @param xi PC for x-axis (1 to k)
#' @param yi PC for y-axis (1 to k)
#' @export
mds = function(counts, condition=NA,k=6,d="euclidian",xi=1,yi=2) {
    require(ggplot2)
    nprobes = nrow(counts)
    nsamples = ncol(counts)
    if (d=="euclidian"){
     distances = dist(t(counts))
    }else if (d=="cor"){
      distances = as.dist(1-cor(counts))
    }
    fit = cmdscale(distances, eig=TRUE, k=k)
    #names = c("one", "two", "three", "four", "five", "six")
    #colnames(fit$points) = names[1:k]
    eigs = data.frame(variance_explained=fit$eig / sum(fit$eig))
    xnames = paste0("PC",1:k," ",round(eigs[1:k,1],digits=2)*100,"%")
    df = as.data.frame(fit$points[,c(xi,yi)])
    names(df) = c("one", "two")
    df$label = rownames(df)
    if(!is.na(condition)) { 
        df$condition = condition
        p = ggplot(df, aes(one, two, label=label, color=condition)) +
            geom_text(aes(one, two, label=label), size=3) +
            labs(list(x=xnames[xi],y=xnames[yi])) +
            scale_x_continuous(expand=c(0.3, 0.3))
    }
    else {
        p = ggplot(df, aes(one, two)) +
            geom_text(aes(one, two, label=label), size=3) +
            labs(list(x=xnames[xi],y=xnames[yi])) + scale_x_continuous(expand=c(0.3, 0.3))
            
    }
        
    return(p)
}

#' Plot PC importance
#'
#' @param counts  matrix
#' @param k  number of dimension
#' @param d distance method (euclidian, cor)
#' @export
variance_by_component = function(counts,k=6,d="euclidian") {
    nsamples = ncol(counts)
    if (d=="euclidian"){
      distances = dist(t(counts))
    }else if (d=="cor"){
      distances = as.dist(1-cor(counts))
    }
    fit = cmdscale(distances, eig=TRUE, k=k)
    eigs = data.frame(variance_explained=fit$eig / sum(fit$eig))
    eigs$component = factor(rownames(eigs), levels=rownames(eigs))
    p = ggplot(eigs, aes(component, variance_explained)) + geom_point() +
            ylab("variance explained") +
            xlab("principal component")
    return(p)
}
