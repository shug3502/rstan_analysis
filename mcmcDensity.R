mcmcDensity <- function(fit, pars = names(fit), byChain = FALSE, nParPerPage = 16, prior = NULL){
    require(dplyr)
    require(tidyr)
    require(ggplot2)
    simsTable <- getSimsTable(fit, pars = pars)  
    posterior <- simsTable %>%
      gather(key = parameter, value = value, -chain, -iteration)
    simsTable <- getSimsTable(fit, pars = pars)

    parameters <- sort(unique(posterior$parameter))
    nParameters <- length(parameters)
    nPages <- ceiling(nParameters / nParPerPage)
    parameters <- data.frame(parameter = parameters,
                             page = sort(rep(1:nPages, length = nParameters)),
                             stringsAsFactors = FALSE)
    posterior <- posterior %>% left_join(parameters)
    posterior$chain <- as.factor(posterior$chain)

    if(!is.null(prior)) prior <- prior %>% left_join(parameters)
    
    for(i in 1:nPages){
        xplot <- subset(posterior, page == i)
        p1 <- ggplot(xplot, aes(x = value))
        if(byChain) p1 <- p1 + aes(color = chain)
        p1 <- p1 + geom_density() + 
                  labs(x = "value", y = "density") +
                      theme(text = element_text(size = 12), axis.text = element_text(size = 8),
                            legend.position = "none", strip.text = element_text(size = 8)) +
                                facet_wrap(~ parameter, ncol = 4, nrow = 4, scales = "free")
        if(!is.null(prior))
            p1 <- p1 + geom_line(data = subset(prior, page == i), aes(x = value, y = density),
                                 color = "red")
        print(p1)
    }
    NULL
}

getSimsTable <- function(x, ...){
    require(dplyr)
    nChains <- dim(x)[2]
    nPost <- dim(x)[1]
    x %>%
        as.data.frame(...) %>%
            mutate(chain = rep(1:nChains, ea = nPost),
                   iteration = rep(1:nPost, nChains))
}

