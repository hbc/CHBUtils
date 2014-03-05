mds = function(counts, condition,k=6,d="euclidian",xi=1,yi=2) {
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
    df$condition = condition
    p = ggplot(df, aes(one, two, color=condition)) +
        geom_point()+
      labs(list(x=xnames[xi],y=xnames[yi]))
    return(p)
}

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
