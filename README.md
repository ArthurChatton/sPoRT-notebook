Sequential Positivity-Regression Trees (sPoRT) notebook
================
Arthur Chatton
2024-12-02

## Goal

This is the notebook related to the paper “Is checking for sequential
positivity violations getting you down? Try sPoRT!” by Chatton,
Schomaker, Luque-Fernandez, Platt and Schnitzer (Submitted).

This notebook aims to illustrate the sPoRT algorithm’s use with
illustrations on a simulated dataset. The algorithm is implemented in R
but we assume only basic knowledge of R.

## Set-up: Import data and functions

    library(rpart)
    library(rpart.plot)

    simdata <- read.csv('sim_data_sport.csv', stringsAsFactors = TRUE)
    source('port_utils.r')

The simulated dataset comes from Schomaker et al. (2019).

    summary(simdata[,1:13])

    ##         W1            W2             W3             L1_0       
    ##  southern:3768   female:2426   Min.   :12.00   Min.   :   0.0  
    ##  west    :1232   male  :2574   1st Qu.:18.00   1st Qu.: 400.0  
    ##                                Median :30.00   Median : 640.0  
    ##                                Mean   :32.65   Mean   : 650.3  
    ##                                3rd Qu.:42.00   3rd Qu.: 880.0  
    ##                                Max.   :54.00   Max.   :1920.0  
    ##                                                                
    ##       L2_0             L3_0             Y_0               L1_1       
    ##  Min.   :0.0000   Min.   :-5.000   Min.   :-10.000   Min.   :   0.0  
    ##  1st Qu.:0.1000   1st Qu.:-2.500   1st Qu.: -4.000   1st Qu.: 440.0  
    ##  Median :0.1500   Median :-2.000   Median : -3.000   Median : 680.0  
    ##  Mean   :0.1384   Mean   :-1.762   Mean   : -3.061   Mean   : 697.2  
    ##  3rd Qu.:0.2000   3rd Qu.:-1.000   3rd Qu.: -2.000   3rd Qu.: 940.0  
    ##  Max.   :0.4000   Max.   : 1.500   Max.   :  3.000   Max.   :2020.0  
    ##                                                                      
    ##       L2_1             L3_1              A_1              C_1       
    ##  Min.   :0.0000   Min.   :-10.000   Min.   :0.0000   Min.   :0.000  
    ##  1st Qu.:0.1000   1st Qu.: -3.000   1st Qu.:0.0000   1st Qu.:0.000  
    ##  Median :0.1500   Median : -2.000   Median :1.0000   Median :0.000  
    ##  Mean   :0.1529   Mean   : -1.935   Mean   :0.5284   Mean   :0.085  
    ##  3rd Qu.:0.2000   3rd Qu.: -1.000   3rd Qu.:1.0000   3rd Qu.:0.000  
    ##  Max.   :0.4500   Max.   :  2.000   Max.   :1.0000   Max.   :1.000  
    ##                                                                     
    ##       Y_1         
    ##  Min.   :-10.000  
    ##  1st Qu.: -4.000  
    ##  Median : -3.000  
    ##  Mean   : -2.992  
    ##  3rd Qu.: -2.000  
    ##  Max.   :  3.000  
    ##  NA's   :425

The dataset `simdata` is composed of 5000 individuals and 43 variables.

1.  `Wx` are baseline confounders
2.  `Lx_t` are time-varying confouders at time t (t=0 for baseline
    value)
3.  `A_t` are treatment value at time t (1 if treated)
4.  `C_t` are censoring indicator at time t (1 if censored)
5.  `Y_t` are outcome value at time t

Note that the values are slightly discretized to save computational
time. Use clinically meaningful thresholds in practice.

## sPoRT with stratification on time

### Construction of treatment strategies

Let four treatment strategies:

1.  Initiate treatment immediately (`d1`)
2.  Never initiate treatment (`d2`)
3.  Initiate treatment when `L1` drops below 750 or when `L2` drops
    below 0.25 (`d3`)
4.  Initiate treatment when `L1` drops below 350 or when `L2` drops
    below 0.15 (`d4`)

<!-- -->

    #define number of time-points
    tps <- 1:6

    for(t in tps){
      
      #Creation of the strategies' dummy variables
      if(t==tps[1]){ #first time-point, no censoring yet
        
        simdata[,paste0('d1_',t)] <- 1
        simdata[,paste0('d2_',t)] <- 0
        simdata[,paste0('d3_',t)] <- 1*(simdata[,paste0('L1_',t)] < 750 | simdata[,paste0('L2_',t)] < 0.25)
        simdata[,paste0('d4_',t)] <- 1*(simdata[,paste0('L1_',t)] < 350 | simdata[,paste0('L2_',t)] < 0.15)
        
      }else{ #other time-points, censoring may occurs and should be considered
        
        simdata[,paste0('d1_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1)
        simdata[,paste0('d2_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 0)
        simdata[,paste0('d3_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1*(simdata[,paste0('L1_',t)] < 750 | simdata[,paste0('L2_',t)] < 0.25 | simdata[,paste0('d3_',tps[which(t==tps)-1])] == 1))
        simdata[,paste0('d4_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1*(simdata[,paste0('L1_',t)] < 350 | simdata[,paste0('L2_',t)] < 0.15 | simdata[,paste0('d4_',tps[which(t==tps)-1])] == 1))
        
      }
      
    }

### sPoRT’s main inputs

    ## Input 1: the treatment strategy columns' name. We can put several strategies through a list
    D.bar <- list(d1=grep("d1", names(simdata), value=T),
                 d2=grep("d2", names(simdata), value=T),
                 d3=grep("d3", names(simdata), value=T),
                 d4=grep("d4", names(simdata), value=T)
                 )

    ## Input 2a: qualitative confounders columns' name, one vector per time-point
    cov.quali <- rep(list(c("W1", "W2")), 6)

    ## Input 2b: quantative confounders columns' name, one vector per time-point
    cov.quanti <- list(
      t1 = c("W3", grep("_0", names(simdata), value=T), "L1_1", 'L2_1', 'L3_1'),
      t2 = c("W3", grep("_0", names(simdata), value=T), "L1_2", 'L2_2', 'L3_2'),
      t3 = c("W3", grep("_0", names(simdata), value=T), "L1_3", 'L2_3', 'L3_3'),
      t4 = c("W3", grep("_0", names(simdata), value=T), "L1_4", 'L2_4', 'L3_4'),
      t5 = c("W3", grep("_0", names(simdata), value=T), "L1_5", 'L2_5', 'L3_5'),
      t6 = c("W3", grep("_0", names(simdata), value=T), "L1_6", 'L2_6', 'L3_6')
    )

    # Note that we put baseline confounders in each list because we want them each time when stratifying on time

    ## Input 3: Actual treatment columns' name and variable type
    treat <- grep("A_",names(simdata),value=T)
    type <- "b" #binary treatment

    ## Input 4: sPoRT hyperparameters
    alpha <- 0.05 #Personal preference, other values are meaningful as well.
    beta <- "gruber" #to use Gruber's bound defined in Gruber et al. (2022) Adapt to the sample size at time t.
    gamma <- 2 #maximum complexity level of the subgroups (here combination of two confounders max.)

    ## Input 5: sPoRT hyperparameters
        monotony <- TRUE #there is a monotone treatment pattern in the data. The function will automatically check positivity on untreated individuals.
        static <- c(TRUE, TRUE, FALSE, FALSE) #d1 and d2 are static, while d3 and d4 are dynamic. When dynamic, sPoRT will check the probability of being treated among those following the rule at time t AND the probability of being to be untreated among those not following the rule.

Now, we can run the `sport` function.

### sPoRT running and output

    #let's begin with stratification on time (argument pooling=FALSE)
    #we did not consider rule 1 since the immediate initiation make only sense at time 1.

        d2_strat_results <- sport(D.bar=D.bar[[2]], A=treat, lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooling=F, beta=beta, alpha=alpha, gamma=gamma, monotony=monotony, static=static[2])

    ## [1] "time 1"
    ## [1] "time 2"
    ## [1] "time 3"
    ## [1] "time 4"
    ## [1] "time 5"
    ## [1] "time 6"

        d3_strat_results <- sport(D.bar=D.bar[[3]], A=treat, lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooling=F, beta=beta, alpha=alpha, gamma=gamma, monotony=monotony, static=static[3])

    ## [1] "time 1"
    ## [1] "time 2"
    ## [1] "time 3"
    ## [1] "time 4"
    ## [1] "time 5"
    ## [1] "time 6"

        d4_strat_results <- sport(D.bar=D.bar[[4]], A=treat, lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooling=F, beta=beta, alpha=alpha, gamma=gamma, monotony=monotony, static=static[4])

    ## [1] "time 1"
    ## [1] "time 2"
    ## [1] "time 3"
    ## [1] "time 4"
    ## [1] "time 5"
    ## [1] "time 6"

The output of sPoRT is a list of table. Each table represents the
violations identified at a certain time for a certain strategy. For
instance:

    d2_strat_results

    ## $T1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=980          0.004  exposed           915              18.3
    ## 2 L1_1>=1040          0.003  exposed           896              17.9
    ## 
    ## $T2
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=940          0.006  exposed          1063              45.3
    ## 2  L1_2>=960          0.010  exposed          1419              60.4
    ## 3 L2_2>=0.35          0.008  exposed           130               5.5
    ## 
    ## $T3
    ##                              subgroup proba.exposure exposure subgroup.size
    ## 1                           L1_0>=940          0.007  exposed          1057
    ## 2                           L1_3>=960          0.013  exposed          1621
    ## 3                          L2_3>=0.35          0.009  exposed           221
    ## 4 Y_0< -2 & W3>=36 & Y_0>=-8 & W3< 42          0.007  exposed           149
    ## 5   L2_0>=0.25 & Y_0< -2 & L2_0< 0.35          0.008  exposed           244
    ## 6  L2_0>=0.25 & L2_0< 0.3 & W2=female          0.008  exposed           126
    ##   subgroup.rel.size
    ## 1              51.7
    ## 2              79.3
    ## 3              10.8
    ## 4               7.3
    ## 5              11.9
    ## 6               6.2
    ## 
    ## $T4
    ##                               subgroup proba.exposure exposure subgroup.size
    ## 1                            L1_0>=900          0.014  exposed          1185
    ## 2                            L2_0>=0.3          0.015  exposed           136
    ## 3                  L3_0>=0 & L3_0< 0.5          0.009  exposed           107
    ## 4                           L1_4>=1040          0.012  exposed          1534
    ## 5                           L2_4>=0.35          0.010  exposed           315
    ## 6                     Y_0>=-2 & W3< 36          0.014  exposed           346
    ## 7  Y_0>=-6 & W3>=42 & W3< 54 & Y_0< -3          0.015  exposed           136
    ## 8          L3_4>=-1 & W3< 36 & L3_4< 0          0.011  exposed           284
    ## 9           L3_4>=-1 & W3< 48 & W3>=42          0.008  exposed           123
    ## 10                   W1=west & Y_0>=-2          0.006  exposed           179
    ##    subgroup.rel.size
    ## 1               62.3
    ## 2                7.2
    ## 3                5.6
    ## 4               80.7
    ## 5               16.6
    ## 6               18.2
    ## 7                7.2
    ## 8               14.9
    ## 9                6.5
    ## 10               9.4
    ## 
    ## $T5
    ##                                    subgroup proba.exposure exposure
    ## 1                                L1_0>=1120          0.000  exposed
    ## 2                                L1_5>=1060          0.011  exposed
    ## 3                    L2_5>=0.35 & L2_5< 0.5          0.014  exposed
    ## 4                        W3< 48 & L2_0>=0.3          0.000  exposed
    ## 5       W3>=24 & Y_0< -3 & Y_0>=-4 & W3< 42          0.008  exposed
    ## 6                          L3_5>=0 & W3< 42          0.013  exposed
    ## 7              L2_0>=0.2 & Y_0>=-2 & Y_0< 0          0.015  exposed
    ## 8                       L3_5>=0 & L2_0>=0.2          0.005  exposed
    ## 9                      W1=west & L2_0>=0.25          0.008  exposed
    ## 10 L3_5< 2 & L3_0>=-1 & L3_5>=1 & L3_0< 0.5          0.009  exposed
    ## 11                        Y_0>=-3 & L3_5>=1          0.015  exposed
    ## 12              Y_0>=-6 & W1=west & Y_0< -3          0.005  exposed
    ##    subgroup.size subgroup.rel.size
    ## 1            511              28.0
    ## 2           1487              81.6
    ## 3            351              19.3
    ## 4             99               5.4
    ## 5            131               7.2
    ## 6            314              17.2
    ## 7            261              14.3
    ## 8            220              12.1
    ## 9            121               6.6
    ## 10           106               5.8
    ## 11           201              11.0
    ## 12           182              10.0
    ## 
    ## $T6
    ##                                     subgroup proba.exposure exposure
    ## 1                                 L1_0>=1080          0.005  exposed
    ## 2                                 L1_6>=1140          0.009  exposed
    ## 3                                  L2_6>=0.4          0.011  exposed
    ## 4  L2_0>=0.15 & W3< 24 & L2_0< 0.25 & W3>=18          0.015  exposed
    ## 5                           L3_6>=1 & W3>=36          0.014  exposed
    ## 6 L3_6>=0 & L3_6< 2 & L2_0>=0.2 & L2_0< 0.25          0.009  exposed
    ## 7                         L3_6>=0 & L3_0< -1          0.010  exposed
    ## 8                 L3_6>=1 & Y_0< 0 & Y_0>=-2          0.007  exposed
    ##   subgroup.size subgroup.rel.size
    ## 1           607              35.1
    ## 2          1301              75.2
    ## 3           182              10.5
    ## 4           137               7.9
    ## 5           146               8.4
    ## 6           115               6.7
    ## 7           199              11.5
    ## 8           136               7.9

are the violations identified at each time for the strategy d2.

These violations are defined in terms of subgroups, such as the
individuals with L1\_0&gt;=980 have 0.004% chance of being `A_1=1`,
which represents 915 individuals (18.3% of the whole sample size). The
`exposure` column seems unnecessary here, but it is useful for
categorical treatment (not implemented yet for sPoRT).

However, the current use of sPoRT assume smoothing over treatment
history. One can stratify on both time and treatment history through the
argument `add.subset` as follow:

    # the `add.subset` argument takes a vector of column's names of `data` (one per time-point). Rows should be coded 1 if kept, and 0 if dropped.

    for(i in tps){
      simdata[,paste0("addsub_", tps)] <- 1*(simdata[,paste0("A_", tps)] == simdata[,paste0("d3_", tps)])
    }

    add.subset <- grep("addsub_", names(simdata), value=T)

    sport(D.bar=D.bar[[3]], A=treat, lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooling=F, beta=beta, alpha=alpha, gamma=gamma, monotony=monotony, static=static[3], add.subset=add.subset)

    ## [1] "time 1"
    ## [1] "time 2"

    ## Error in port(A[t], type_A = type_A, cov.quanti = qtcov, cov.quali = qlcov, : Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the argument 'A' when the argument 'type_A' is 'b' (i.e., binary).

Here, positivity is checked at each time for the individuals that should
initiate at this time but having not initiated nor should have initiated
before. The error is present because at one time there is only treated
or untreated individuals in this subset, which is a common occurence
when stratifying on treatment history.

## sPoRT with pooling over time-points

### Pivoting to long format

When we want an analysis at the person-time level, one must check
positivity at the same scale. The sPoRT algorithm can handle long-format
datasets as well. But first, we need to define the lagged values of the
treatment (for subseting on A\_t-1=0).

    simdata$Alag_1 <- 1
        simdata$Alag_2 <-simdata$A_1
        simdata$Alag_3 <-simdata$A_2
        simdata$Alag_4 <-simdata$A_3
        simdata$Alag_5 <-simdata$A_4
        simdata$Alag_6 <-simdata$A_5
        
    pool <- reshape(simdata, 
                    varying = list(
                      grep("L1", names(simdata), value=T)[-1],
                      grep("L2", names(simdata), value=T)[-1],
                      grep("L3", names(simdata), value=T)[-1],
                      grep("A_", names(simdata), value=T),
                      grep("Alag_", names(simdata), value=T),
                      grep("Y_", names(simdata), value=T)[-1],
                      grep("C_", names(simdata), value=T),
                      grep("d2", names(simdata), value=T),
                      grep("d3", names(simdata), value=T),
                      grep("d4", names(simdata), value=T)
                    ),
                    v.names = c("L1", "L2", "L3", "A", "Alag", "Y", "C", "d2", "d3", "d4"),
                    timevar = "T", 
                    times = tps, 
                    direction = "long")

### sPoRT’s main inputs

        ## Input 1: the treatment strategy columns' name. We can put several strategies through a vector
        D.bar <- c("d1", "d2", "d3", "d4")

        ## Input 2: qualitative and quantitative confounders
        cov.quanti <- c("Y_0", "L1_0",  "L2_0", "L3_0", "W3", "L1", "L2", "L3") 
        cov.quali <- c("W1", "W2") 

        ## Input 3: Actual treatment column's name and variable type
        treat <- "A"
        type <- "b" #binary treatment

        ## Input 4: sPoRT hyperparameters
        alpha <- 0.05 #Personal preference, other values are meaningful as well.
        beta <- "gruber" #to use Gruber's bound defined in Gruber et al. (2022) Adapt to sample size.
        gamma <- 2 #maximum complexity level of the subgroups (here, a combination of two confounders max.)

        ## Input 5: sPoRT hyperparameters, monotony must be the column name of the lagged A when pooling=TRUE due to implementation issue (only if there is a monotone treatment pattern in the data). 
        monotony <- "Alag" #The function will automatically check positivity on individuals with Alag=0
        static <- c(TRUE, TRUE, FALSE, FALSE) #d1 is static, d3 and d4 are dynamic. 

Now, we can run the `sport` function, with the argument `pooling=TRUE`.

### sPoRT running and output

        d2_pool_results <- sport(D.bar=D.bar[[2]], A="A", time='T', lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=pool, pooling=TRUE, alpha=alpha, beta=beta, gamma=gamma, static=static[2], monotony=monotony)

        d3_pool_results <- sport(D.bar=D.bar[[3]], A="A", time='T', lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=pool, pooling=TRUE, alpha=alpha, beta=beta, gamma=gamma, static=static[3], monotony=monotony)

        d4_pool_results <- sport(D.bar=D.bar[[4]], A="A", time='T', lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=pool, pooling=TRUE, alpha=alpha, beta=beta, gamma=gamma, static=static[4], monotony=monotony)

Results’ presentation is similar; but one sees only one table per
strategy due to pooling over time.

        d2_pool_results 

    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.002  exposed          2775              28.2
    ## 2   L1>=1160          0.004  exposed          5149              52.3

## Other inputs and miscenallous

The `sport`function can takes several other arguments. First, `lag` is
used to add lagged effects (e.g., `lag=1` means that we add
automatically the counfounders from the one time to the next, such as
`L1_1` in the second time). Second, `pruning` can be set to `TRUE` when
one focus only on structural violations, because this remove violations
bounded between two values from the output for clarity (e.g., this may
remove violations such as `L1_0>=960 & L1_0< 860`). Last, the `rpart`
hyperparameters can also be set (Therneau, 2019). By default, the
function use `minbucket = 6`, `minsplit = 20`, and `maxdepth = 30` to
allows small subgroups. Other values can be set to split only
smaller/bigger nodes.

### References

Chatton A, Schomaker M, Luque-Fernandez M-A, Platt RW, Schnitzer ME. Is
checking for sequential positivity violations getting you down? Try
sPoRT! *Submitted*. 2024.

Gruber S, Phillips RV, Lee H, van der Laan MJ. Data-Adaptive Selection
of the Propensity Score Truncation Level for
Inverse-Probability-Weighted and Targeted Maximum Likelihood Estimators
of Marginal Point Treatment Effects. *Am J Epidemiol.*
2022;191(9):1640–51.

Schomaker M, Luque-Fernandez MA, Leroy V, Davies MA. Using longitudinal
targeted maximum likelihood estimation in complex settings with dynamic
interventions. *Statistics in Medicine.* 2019;38(24):4888–911.

Therneau TM, Atkinson B. rpart: Recursive Partitioning and Regression
Trees. 2019. <https://CRAN.R-project.org/package=rpart> (20 October
2021, date last accessed).

### Updates history

First update planned for the release of the `port` package.

