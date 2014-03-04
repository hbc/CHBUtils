mds = function(counts, condition,k=6) {
    require(ggplot2)
    nprobes = nrow(counts)
    nsamples = ncol(counts)
    distances = dist(t(counts))
    fit = cmdscale(distances, eig=TRUE, k=k)
    names = c("one", "two", "three", "four", "five", "six")
    colnames(fit$points) = names[1:k]
    df = as.data.frame(fit$points)
    df$label = rownames(df)
    df$condition = condition
    p = ggplot(df, aes(one, two, color=condition)) +
        geom_point()
    return(p)
}

variance_by_component = function(counts,k=6) {
    nsamples = ncol(counts)
    distances = dist(t(counts))
    fit = cmdscale(distances, eig=TRUE, k=k)
    eigs = data.frame(variance_explained=fit$eig / sum(fit$eig))
    eigs$component = factor(rownames(eigs), levels=rownames(eigs))
    p = ggplot(eigs, aes(component, variance_explained)) + geom_point() +
            ylab("variance explained") +
            xlab("principal component")
    return(p)
}
