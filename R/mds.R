mds = function(counts, condition) {
    require(ggplot2)
    require(gridExtra)
    nprobes = nrow(counts)
    nsamples = ncol(counts)
    distances = dist(t(counts))
    fit = cmdscale(distances, eig=TRUE, k=6)
    colnames(fit$points) = c("one", "two", "three", "four", "five", "six")
    df = as.data.frame(fit$points)
    df$label = rownames(df)
    df$condition = condition
    p1 = ggplot(df, aes(one, two, color=condition)) +
        geom_point() +
        theme(text=element_text(family="Gill Sans"))

    eigs = data.frame(variance_explained=fit$eig / sum(fit$eig))
    eigs$component = factor(rownames(eigs), levels=rownames(eigs))
    quartz()
    p2 = ggplot(eigs, aes(component, variance_explained)) + geom_point() +
        theme(text=element_text(family="Gill Sans")) +
            ylab("Proportion of variance explained") +
            xlab("Principal component")     
    grid.arrange(p1, p2)
}
