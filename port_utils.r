
library(rpart)
library(rpart.plot)

#' Violation definition
#'
#' Produce a row of results from a problematic subgroup for the user output (table of results).
#'
#' @param subgrp the subgroup produced in port()
#' @param var Variables names defining the subgroup
#' @param type_var Type of the variable (continuous or categorical)
#' @param data Dataset
#' @param pruning Boolean, should we remove subgroup like >30 & <40?
#' @param type_A Type of the exposure (binary, continuous, or categorical)
#' @param mediation Boolean, are we in mediation analysis? Only for table labeling
#' @param beta Beta hyperparameter
#' @param A Column name of the exposure
#' @param check.side Which side of the positivity assumption should be checked? P(A|C), beta is 'zero', P(A|C)>1-beta is 'one', both is 'both' (default)
#'
#' @return Table with the violations (subgroup name | A=1 probability | A level | subgroup size (N) | subgroup size (%))
#'
#' @examples
define_cutoff <- function(subgrp, var, type_var, data, pruning, type_A="b", mediation=FALSE, beta, A, check.side="both"){
  n <- subgrp$n[nrow(subgrp)]
  pourcent <- round(n/subgrp$n[1], digits = 3) * 100
  if(type_A!="b"){
    all_prob <- subgrp$proba[nrow(subgrp),]
    if(check.side=="both"){
      pb_prob <- which(all_prob<=beta | all_prob>=(1-beta))
    }else{
      if(check.side=="zero"){
        pb_prob <- which(all_prob<=beta)
      }else{
        pb_prob <- which(all_prob>beta)
      }
    }
    proba <- round(all_prob[pb_prob], digits = 3) |> setNames(nm=NULL)
    treat <- levels(data[,A])[pb_prob]
  }else{
    proba <- round(subgrp$proba[nrow(subgrp)], digits = 3)
    treat <- "exposed"
  }
  newcut <- subgrp
  
  for (v in var) {
    if (type_var[var == v] == "continuous") {
      imcs <- newcut$var2[grepl(v, newcut$var2, fixed = TRUE)]
      loc_imcs <- which(newcut$var2 %in% imcs)
      imcs <- sapply(imcs, function(i) gsub(v, "", i))
      imcs_sup <- imcs[grepl("<", imcs)]
      loc_sup <- loc_imcs[grepl("<", imcs)]
      imcs_inf <- imcs[grepl(">", imcs)]
      loc_inf <- loc_imcs[grepl(">", imcs)]
      if (pruning) {
        if (length(imcs_sup) > 0 && length(imcs_inf) >
            0) {
          newcut <- NULL
          return(list(subgrp = newcut, table=data.frame()))
        }
      }
      imcs_sup <- as.numeric(vapply(imcs_sup, function(i) substr(i, start = 3, stop = nchar(i)), FUN.VALUE = character(1)))
      imcs_inf <- as.numeric(vapply(imcs_inf, function(i) substr(i, start = 3, stop = nchar(i)), FUN.VALUE = character(1)))
      if (length(imcs_sup) > 1 || length(imcs_inf) > 1) {
        imc_sup <- ifelse(length(imcs_sup) > 0, min(imcs_sup), imcs_sup)
        imc_inf <- ifelse(length(imcs_inf) > 0, max(imcs_inf), imcs_inf)
        keep <- c(loc_sup[which(imcs_sup == imc_sup)],
                  loc_inf[which(imcs_inf == imc_inf)])
        newcut <- newcut[-setdiff(loc_imcs, keep),]
      }
      if ("root" %in% newcut$var2 && length(newcut$var2) > 1) {
        newcut <- newcut[-1, ]
      }
    }
    if (type_var[var == v] == "categorical") {
      index <- grep(v, newcut$var2, fixed = TRUE)
      l <- index[-length(index)]
      if (length(l) > 0) {
        newcut <- newcut[-l, ]
      }
      if ("root" %in% newcut$var2) {
        newcut <- newcut[-1, ]
      }
    }
  }
  
  if(mediation){
    tab <- data.frame(subgroup=paste(unique(newcut$var2), collapse=" & ") |> rep(length(treat)),
                      proba.mediator=unlist(proba),
                      mediator=treat,
                      subgroup.size=n |> rep(length(treat)),
                      subgroup.rel.size=pourcent |> rep(length(treat)),
                      row.names = NULL)
  }else{
    tab <- data.frame(subgroup=paste(unique(newcut$var2), collapse=" & ") |> rep(length(treat)),
                      proba.exposure=unlist(proba),
                      exposure=treat,
                      subgroup.size=n  |> rep(length(treat)),
                      subgroup.rel.size=pourcent  |> rep(length(treat)),
                      row.names = NULL)
  }
  
  
  return(list(data = newcut, table=tab))
}


#' Gruber bound
#'
#' @param data Dataset used for computation
#'
#' @return The bound used to define the non-positivity
gruber <- function(data){
  5/(sqrt(nrow(data))*log(nrow(data)))
}


#' To obtain the cutoff / class
#'
#' @param object rpart object
#' @param data data set
#' @param digits
#' @param minlength
#' @param pretty
#' @param collapse
#' @param ...
#'
#' @return return the rpart object with the cutoffs (< or >) included
#'
#' @examples
labels.rpart <- function(object, data, digits = 4, minlength = 1L,
                         pretty, collapse = TRUE, ...) {
  if (missing(minlength) && !missing(pretty)) {
    minlength <- if (is.null(pretty))
      1L
    else if (is.logical(pretty)) {
      if (pretty)
        4L
      else 0L
    }
    else 0L
  }
  ff <- object$frame
  n <- nrow(ff)
  if (n == 1L)
    return("root")
  whichrow <- !(ff$var == "<leaf>")
  vnames <- ff$var[whichrow]
  index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + whichrow))
  irow <- index[c(whichrow, FALSE)]
  ncat <- object$splits[irow, 2L] #1 if >=, -1 if <
  lsplit <- rsplit <- character(length(irow))
  if (any(ncat < 2L)) {
    jrow <- irow[ncat < 2L]
    splits <- object$splits[jrow, 4L]
    
    if(length(splits) == 1) names(splits) <- row.names(object$splits)[jrow]#vnames
    
    cutpoint <- sapply(1:length(splits), function(s) findInterval(splits[s], sort(unique(data[,unique(names(splits[s]))]))))
    
    temp1 <- (ifelse(ncat < 0, "< ", ">="))[ncat < 2L]
    temp2 <- (ifelse(ncat < 0, ">=", "< "))[ncat < 2L]
    
    formatg <- function(x, digits = getOption("digits"), 
                        format = paste0("%.", digits, "g")) {
      if (!is.numeric(x)) 
        stop("'x' must be a numeric vector")
      temp <- sprintf(format, x)
      if (is.matrix(x)) 
        matrix(temp, nrow = nrow(x))
      else temp
    }
    
    lsplit[ncat < 2L] <- sapply(1:length(splits), function(s) 
      paste0(temp1[s], formatg(sort(unique(data[,unique(names(splits[s]))]))[cutpoint[s]+1], digits=digits))
    )
    rsplit[ncat < 2L] <- sapply(1:length(splits), function(s) 
      paste0(temp2[s], formatg(sort(unique(data[,unique(names(splits[s]))]))[cutpoint[s]+1], digits = digits))
    )
  }
  if (any(ncat > 1L)) {
    xlevels <- attr(object, "xlevels")
    jrow <- seq_along(ncat)[ncat > 1L]
    crow <- object$splits[irow[ncat > 1L], 4L]
    cindex <- (match(vnames, names(xlevels)))[ncat > 1L]
    if (minlength == 1L) {
      if (any(ncat > 52L))
        warning("more than 52 levels in a predicting factor, truncated for printout", domain = NA)
    }
    else if (minlength > 1L)
      xlevels <- lapply(xlevels, abbreviate, minlength, ...)
    for (i in seq_along(jrow)) {
      j <- jrow[i]
      splits <- object$csplit[crow[i], ]
      cl <- if (minlength == 1L)
        ""
      else ","
      lsplit[j] <- paste((xlevels[[cindex[i]]])[splits == 1L], collapse = ",")
      rsplit[j] <- paste((xlevels[[cindex[i]]])[splits ==  3L], collapse = ",")
    }
  }
  if (!collapse) {
    ltemp <- rtemp <- rep("<leaf>", n)
    ltemp[whichrow] <- lsplit
    rtemp[whichrow] <- rsplit
    return(cbind(ltemp, rtemp))
  }
  lsplit <- paste0(ifelse(ncat < 2L, "", "="), lsplit)
  rsplit <- paste0(ifelse(ncat < 2L, "", "="), rsplit)
  varname <- (as.character(vnames))
  node <- as.numeric(row.names(ff))
  parent <- match(node%/%2L, node[whichrow])
  odd <- (as.logical(node%%2L))
  labels <- character(n)
  labels[odd] <- paste0(varname[parent[odd]], rsplit[parent[odd]])
  labels[!odd] <- paste0(varname[parent[!odd]], lsplit[parent[!odd]])
  labels[1L] <- "root"
  labels
}

#' Found parent node
#'
#' @param x Frame in a rpart object
#'
#' @return return the number of the parent node
#'
#' @examples
parent <- function(x) {
  if (x[1] != 1)
    c(Recall(if (x%%2 == 0L) x/2 else (x - 1)/2),
      x)
  else x
}



#' Positivity Regression trees
#'
#' Check the positivity assumption and identify problematic subgroups/covariables.
#'
#' @param A Label of the column related to the exposure
#' @param type_A Type of the exposure ('b' for binary, 'c' for continuous, or 'n' for nominal)
#' @param cov.quanti Columns' labels of the quantitative variables in the adjustment set
#' @param cov.quali Columns' labels of the qualitative variables in the adjustment set
#' @param data Dataset, must be a data.frame object
#' @param alpha Subgroup minimal size (as a proportion of the whole sample)
#' @param beta Threshold for non-positivity (i.e., extreme exposure's probability).
#' @param gamma Maximal number of variables to define a subgroup
#' @param check.side Which side of the positivity assumption should be checked? P(A|C)<beta is 'zero', P(A|C)>1-beta is 'one', both is 'both' (default)
#' @param mediation Boolean, is the exposure a mediator?
#' @param graph Provide the trees as additional output?
#' @param pruning Boolean, should we keep the 'internal' violation (e.g., >30 & <40)
#' @param minbucket Rpart hyperparameter, minimum number of individual in the leaves
#' @param minsplit Rpart hyperparameter, minimum number of individual in a node to be split
#' @param maxdepth Rpart hyperparameter, maximum number of successive nodes
#' @param tweak Text size for the graphs
#'
#' @return Data.frame with all positivity violations (subgroup name | A=1 probability | A level | subgroup size (N) | subgroup relative size (%))
#' @export
#'
#' @examples
port <- function (A, type_A="b", cov.quanti, cov.quali, data, alpha = 0.05, beta = 'gruber', gamma = 2, check.side="both", mediation=FALSE, graph="none", pruning = FALSE, minbucket = 6, minsplit = 20, maxdepth = 30, tweak=1){
  
  if (!(is.data.frame(data) | is.matrix(data))){
    stop("The argument \'data\' need to be a data.frame or a matrix")
  }
  if(nrow(data)==0) return("Empty dataset.") #use return instead of stop for sport(). 
  if (alpha > 0.5 | alpha < 0)
    stop("The argument \'alpha\' must be a proportion (e.g., 0.05 for 5%).")
  if(beta=='gruber') beta <- gruber(data)
  if (beta > 1 | beta <= 0)
    stop("The argument \'beta\' must be a non-null proportion (e.g., 0.05 for 5%).")
  if(type_A=="b"){
    if (!all(names(table(data[, A])) == c("0","1"))) {
      stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the argument \'A\' when the argument \'type_A\' is \'b\' (i.e., binary).")
    }
  }
  if (type_A=="c") { # categorisation according to the quartiles
    data[, A] <- as.factor(cut(data[, A], breaks = quantile(data[, A], seq(0, 1, by = 0.25), na.rm = FALSE)))
  }
  if(type_A=="n"){ # nominal
    data[, A] <- as.factor(data[,A])
    if(length(levels(data[,A])) < 3){
      stop("At least three modalities are required in the argument \'A\' when the argument \'type_A\' is \'n\' (i.e., nominal).")
    }
  }
  check_prop_a <- prop.table(table(data[,A]))
  if(any(check_prop_a<=beta) & check.side!="one") return("The whole sample presents at least one exposure modality's prevalence lower than beta.")
  if(any(check_prop_a >= 1-beta) & check.side!="zero") return("The whole sample presents at least one exposure modality's prevalence higher than 1-beta.")
  
  if(length(cov.quali)>1){
    data[, cov.quali] <- apply(data[, cov.quali], 2, as.factor)
  }else{
    if(length(cov.quali)==1) data[, cov.quali] <- as.factor(data[, cov.quali])
  }
  
  covariates <- c(cov.quanti, cov.quali)
  m <- length(covariates)
  up <- sapply(1:gamma, function(x) ncol(combn(m,x)))
  savegraph <- problem_covariates <- problem_cutoffs <- list()
  
  combi <- sapply(1:gamma, function(x) combn(m, x) |> split(f=col(combn(m, x))) |> unname(), simplify=F ) |> unlist(recursive = FALSE)
  
  for (q in 1:length(combi)) {
    if (q %in% (up+1)) { # remove problematic covariates already identified when gamma is updated
      if (length(problem_cutoffs) > 0) {
        var_prob <- problem_covariates |> unlist() |> unique()
        bad_cov <- which(covariates %in% var_prob)
      }
    }
    if (exists("bad_cov") && any(bad_cov %in% combi[[q]])) { #pass combination if violation found in a predictor involved inside
      next
    }
    covariables <- covariates[combi[[q]]]
    cart_max <- rpart::rpart(reformulate(covariables, A), data = data, cp = 0, minbucket = minbucket, method = ifelse(type_A=="b", "anova", "class"), minsplit = minsplit, maxdepth = maxdepth)
    frame <- cart_max$frame
    frame$var2 <- labels.rpart(cart_max, data = data)
    if (nrow(frame) > 1) {
      for (j in 2:nrow(frame)) {
        for (cov in covariables) {
          if (grepl(cov, frame$var2[j])) {
            frame$var[j] <- cov
          }
        }
      }
    }
    problematic_nodes <- numeric(0)
    if(type_A=="b"){
      for (i in 1:nrow(frame)) { #read tree with alpha & beta -> save problematic nodes 
        if(check.side=="both"){
          if ((frame$yval[i] >= 1 - beta || frame$yval[i] <=
               beta) && frame$n[i] >= nrow(data) * alpha) { 
            problematic_nodes <- c(problematic_nodes, as.numeric(rownames(frame)[i]))
          }
        }else{
          if(check.side=="zero"){
            if (frame$yval[i] <= beta && frame$n[i] >= nrow(data) * alpha) { 
              problematic_nodes <- c(problematic_nodes, as.numeric(rownames(frame)[i]))
            }
          }else{
            if (frame$yval[i] >= 1 - beta && frame$n[i] >= nrow(data) * alpha) { 
              problematic_nodes <- c(problematic_nodes, as.numeric(rownames(frame)[i]))
            }
          }
        }
      }
    }else{
      p_a <- data.frame(matrix(frame$yval2[,(length(table(data[,A]))+2):(dim(frame$yval2)[2]-1)], ncol = length(levels(data[,A])))) |> setNames(nm = levels(data[,A])) #keep proba of each exposure modality
      for (i in 1:nrow(frame)) { #read tree with alpha & beta -> save problematic nodes
        if(check.side=="both"){
          if (any( (p_a[i,] >= 1 - beta | p_a[i,] <=
                    beta) & frame$n[i] >= nrow(data) * alpha )) {
            problematic_nodes <- c(problematic_nodes, as.numeric(rownames(frame)[i]))
          }
        }else{
          if(check.side=="zero"){
            if (any( p_a[i,] <= beta & frame$n[i] >= nrow(data) * alpha )) {
              problematic_nodes <- c(problematic_nodes, as.numeric(rownames(frame)[i]))
            }
          }else{
            if (any( p_a[i,] >= 1 - beta & frame$n[i] >= nrow(data) * alpha )) {
              problematic_nodes <- c(problematic_nodes, as.numeric(rownames(frame)[i]))
            }
          }
        }
        
      }
    }
    for (n in problematic_nodes) { #check if a node is an ancestor of another, if yes keep the ancestor.
      for (p in parent(n)[-length(parent(n))]) {
        if (p %in% problematic_nodes) {
          problematic_nodes <- problematic_nodes[!problematic_nodes == n]
        }
      }
    }
    
    problematic_path <- list()
    for (i in problematic_nodes) { # save the path from root to the pb node
      problematic_path[[as.character(i)]] <- frame[which(rownames(frame) %in% parent(i)), c("var", "var2", "n")]
      if(type_A=="b"){
        problematic_path[[as.character(i)]]$proba <- frame$yval[which(rownames(frame) %in%  parent(i))]
      }else{
        problematic_path[[as.character(i)]]$proba <- p_a[which(rownames(frame) %in%  parent(i)),]
      }
      
    }
    problem_names <- 1
    for (i in names(problematic_path)) {
      problem_names <- c(rownames(problematic_path[[i]]), problem_names)
    }
    problem_names <- unique(problem_names)
    
    ####
    if(graph!="none"){
      cols<-ifelse(as.numeric(row.names(cart_max$frame)) %in% problem_names,"red","blue")
      rpart.plot::rpart.plot(cart_max, branch.col=cols, tweak=tweak) |> suppressWarnings()
      savegraph[[q]] <- list(tree=recordPlot(load=c("rpart", "rpart.plot")), problem=ifelse(length(problematic_path) != 0, TRUE, FALSE))
    }
    ###
    
    if (length(problematic_path) > 0) { # save all new prb path in problem_cutoff (can thereby reuse problematic_path for the next iteration)
      for (k in 1:length(problematic_path)){
        if (length(problem_cutoffs) == 0 || !any(sapply(problem_cutoffs, problematic_path[k][[names(problematic_path[k])]], FUN = identical)) ) {
          problem_cutoffs[[as.character(length(problem_cutoffs) + 1)]] <- problematic_path[k][[names(problematic_path[k])]]
          problem_covariates[[as.character(length(problem_covariates) + 1)]] <- covariables
        }
      }
    }
  } # end of the "q" for-loop
  
  res <- data.frame()
  if (length(problem_cutoffs) > 0) {
    for (i in 1:length(problem_cutoffs)) {
      vars <- unique(problem_cutoffs[[i]]$var)
      type <- character(length(vars))
      for (k in 1:length(vars)) {
        type[k] <- ifelse(vars[k] %in% cov.quanti, "continuous", "categorical")
      }
      problem_covariates[[i]] <- paste(problem_covariates[[i]],  collapse = ";")
      cut <- define_cutoff(problem_cutoffs[[i]], vars, type_var=type, data, pruning = pruning, type_A, mediation, beta, A, check.side)
      
      if (is.null(cut$data)) {
        problem_cutoffs[[i]] <- NULL
        problem_covariates[[i]] <- NULL
      }else {
        problem_cutoffs[[i]] <- cut$data
      }
      
      
      res <- rbind(res, cut$table)
      
    }
    
    #enregistrer les variables plutot que les graphs et ajouter une fonction pour retrouver les arbres entiers a partir de l'objet port.
    if(graph!="none") savegraph <- savegraph[which(sapply(lapply(savegraph, '[[', 'problem'), function(x) !is.null(x)))] #cut null result from tree not built because an one-variable violation was found
    
    output <- switch(graph,
                     all=list(
                       table = res,
                       problematic.graphs = lapply(savegraph[sapply(savegraph, '[[', 'problem')], '[[', 'tree'),
                       correct.graphs = lapply(savegraph[!sapply(savegraph, '[[', 'problem')], '[[', 'tree')
                     ),
                     none=res,
                     prob=list(
                       table = res,
                       problematic.graphs = lapply(savegraph[sapply(savegraph, '[[', 'problem')], '[[', 'tree')
                     ),
                     correct=list(
                       table = res,
                       correct.graphs = lapply(savegraph[!sapply(savegraph, '[[', 'problem')], '[[', 'tree')
                     )
    )
    
  }else{
    
    output <- "No problematic subgroup was identified."
    
  }
  
  return(output)
}


#' Sequential positivity Regression Trees
#'
#' Check the sequential positivity assumption for longitudinal/panel data to identify problematic subgroups/covariables over time.
#' This function is a wrapper around port() and deals with wide or long data.frames.
#'
#' @param A Intervention indicators. Must be a list of each intervention columns' labels.
#' @param D.bar Treatement strategy history's column names. Must be a matrix (n*times).
#' @param time Column name for the time when pooled=T, vector of time when pooling=F (optional).
#' @param static Boolean, is the rule static (TRUE, default) or dynamic (FALSE).
#' @param monotony Boolean, is the treatment pattern is monotone in data? FALSE by default. See Details.
#' @param add.subset Column name related to an indicator of additionnal subseting (1 to keep, 0 to remove), See Details.
#' @param lag Lagged values to consider for covariates (default is 1). See Details
#' @param type_A Type of the intervention ('b' for binary, 'c' for continuous, or 'n' for nominal)
#' @param cov.quanti Columns' labels of the quantitative variables in the adjustment set. Must be a list with the same length than A if pooled=FALSE. See Details
#' @param cov.quali Columns' labels of the qualitative variables in the adjustment set. Must be a list with the same length than A if pooled=FALSE. See Details
#' @param data Dataset, either in the wide format (if pooled=FALSE) or long format (if pooled=TRUE)
#' @param pooling Boolean, should the positivity be checked at each time (FALSE, default) or globally (TRUE)?
#' @param alpha Subgroup minimal size (as a proportion of the whole sample)
#' @param beta Threshold for non-positivity (i.e., extreme exposure's probability).
#' @param gamma Maximal number of variables to define a subgroup
#' @param minbucket Rpart hyperparameter, minimum number of individual in the leaves
#' @param minsplit Rpart hyperparameter, minimum number of individual in a node to be split
#' @param maxdepth Rpart hyperparameter, maximum number of successive nodes
#' @param pruning Boolean, should we keep the 'internal' violation (e.g., 20<cov<30)
#'
#' @return List of regimen-specific (and time-specific if pooled=FALSE) data.frames with the potential positivity violations identified.
#' @export
#' 
sport <- function(A, D.bar, static=TRUE, monotony=FALSE, add.subset=NULL, time=NULL, lag=0, type_A="b", cov.quanti, cov.quali, data, pooling=FALSE, alpha = 0.05, beta = 'gruber', gamma = 2, minbucket = 6, minsplit = 20, maxdepth = 30, pruning = FALSE){
  
  ###
  # add checks here
  ###
 
  tps <- length(A)
  
  if(!pooling){ #need data in wide format
    
    if(!pooling & monotony){
      if( data[data[,D.bar[length(D.bar)]]==0,D.bar] |> rowSums(, na.rm=T) |> sum() != 0){ stop("In case of treatment monotony, the rule values of the variables passed in the argument \'D.bar\' must be monotone too (i.e., cannot be 1 then 0).") }
    } 
    
    if(is.null(time)) time <- paste0("T", 1:tps)
    
    
    res <- sapply(1:tps,
                  function(t){
                    
                    if(monotony & t!=1){
                        subdata <- subset(data, get(A[t-1])==0)
                    }else{
                        subdata <- data
                    }
                      
                    
                    
                    if(!is.null(add.subset)) subdata <- subset(subdata, get(add.subset[t])==1)
                      
                    subdata <- list(Aeq1 = subset(subdata, get(D.bar[t])==1),
                                    Aeq0 = subset(subdata, get(D.bar[t])==0))
                    
                    
                    qlcov <- unique(unlist(cov.quali[max(0,(t-lag)):t])) 
                    qtcov <-  unique(unlist(cov.quanti[max(0,(t-lag)):t]))
                    
                    print(paste("time", t)) 

                    
                    res.t <- list(Aeq1 = port(A[t], type_A=type_A, cov.quanti=qtcov, cov.quali=qlcov, data=subdata[[1]], alpha=alpha, beta=beta, gamma=gamma, check.side="zero", mediation=FALSE, pruning=pruning, minbucket=minbucket, minsplit=minsplit, maxdepth=maxdepth),
                                  Aeq0 = port(A[t], type_A=type_A, cov.quanti=qtcov, cov.quali=qlcov, data=subdata[[2]], alpha=alpha, beta=beta, gamma=gamma, check.side="one", mediation=FALSE, pruning=pruning, minbucket=minbucket, minsplit=minsplit, maxdepth=maxdepth)
                      )
                    
                    
                    return(res.t)
                    
                  },
                  simplify=FALSE
    )
    
    names(res) <- time #voir ici
    
  }else{ #if pool=T, need data in long_format
    
    
    #for monotony, should pass by add.subset (or, little trick, pass the name of the previous A in monotony)
    
    
    if(!is.null(add.subset)){
      subdata <- subset(data, get(add.subset)==1)
    }else{
      subdata <- data
    }
    
    if(is.character(monotony)) subdata <- subset(data, get(monotony)==0)
    
    subdata <- list(Aeq1 = subset(subdata, get(D.bar)==1),
                    Aeq0 = subset(subdata, get(D.bar)==0))
      
      res <- list(Aeq1 = port(A=A, type_A=type_A, cov.quanti=cov.quanti, cov.quali=cov.quali, data=subdata[[1]], alpha = alpha, beta = beta, gamma = gamma, check.side="zero", mediation=FALSE, pruning = pruning, minbucket = minbucket, minsplit = minsplit, maxdepth = maxdepth),
                  Aeq0 = port(A=A, type_A=type_A, cov.quanti=cov.quanti, cov.quali=cov.quali, data=subdata[[2]], alpha = alpha, beta = beta, gamma = gamma, check.side="one", mediation=FALSE, pruning = pruning, minbucket = minbucket, minsplit = minsplit, maxdepth = maxdepth))
    
    
  }
  
  return(res)
  
  
}
