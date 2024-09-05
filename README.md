Sequential Positivity-Regression Trees (sPoRT) notebook
================
Arthur Chatton
2024-07-18

## Goal

This is the notebook related to the paper “Is checking for sequential
positivity violations getting you down? Try sPoRT!” by Chatton,
Schomaker, Luque-Fernandez, Platt, and Schnitzer (Submitted).

This notebook aims to illustrate the sPoRT algorithm’s use with
illustrations on a simulated dataset. The algorithm is implemented in R,
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
2.  `Lx_t` are time-varying confounders at time t (t=0 for baseline
    value)
3.  `A_t` are treatment value at time t (1 if treated)
4.  `C_t` are censoring indicator at time t (1 if censored)
5.  `Y_t` are outcome value at time t

Note that the values are slightly discretized to save computational
time. Use clinically meaningful thresholds in practice.

## sPoRT with stratification on time

### Construction of treatment strategies / intervention rules

Let four intervention rules:

1.  Initiate treatment immediately (rule 1)
2.  Never initiate treatment (S2)
3.  Initiate treatment when `L1` drops below 750 or when `L2` drops
    below 0.25 (S3)
4.  Initiate treatment when `L1` drops below 350 or when `L2` drops
    below 0.15 (S4)

<!-- -->

    #define number of time-points
    tps <- 1:6

    for(t in tps){
      
      #Creation of the rules' dummy variables
      if(t==tps[1]){ #first time-point, no censoring yet
        
        simdata[,paste0('d1_',t)] <- 1
        simdata[,paste0('d3_',t)] <- 1*(simdata[,paste0('L1_',t)] < 750 | simdata[,paste0('L2_',t)] < 0.25)
        simdata[,paste0('d4_',t)] <- 1*(simdata[,paste0('L1_',t)] < 350 | simdata[,paste0('L2_',t)] < 0.15)
        
      }else{ #other time-points, censoring may occur and should be considered
        
        simdata[,paste0('d1_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1)
        simdata[,paste0('d3_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1*(simdata[,paste0('L1_',t)] < 750 | simdata[,paste0('L2_',t)] < 0.25 | simdata[,paste0('d3_',tps[which(t==tps)-1])] == 1))
        simdata[,paste0('d4_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1*(simdata[,paste0('L1_',t)] < 350 | simdata[,paste0('L2_',t)] < 0.15 | simdata[,paste0('d4_',tps[which(t==tps)-1])] == 1))
        
      }
      
    }

Note that we didn’t code the rule d2 (never initiate) because this rule
will be checked using rule 1 (but with a probability close to one).

### sPoRT’s main inputs

    ## Input 1: the treatment strategy columns' name. 
    D.bar <- list(d1=grep("d1", names(simdata), value=T),
                 d3=grep("d3", names(simdata), value=T),
                 d4=grep("d4", names(simdata), value=T)
                 )

    ## Input 2a: qualitative confounders columns' name, one vector per time-point
    cov.quali <- rep(list(c("W1", "W2")), 6)

    ## Input 2b: quantitative confounders columns' name, one vector per time-point
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

    ## Input 4: PoRT hyperparameters
    alpha <- 0.05 #Personal preference, other values are meaningful as well.
    beta <- "gruber" #to use Gruber's bound defined in Gruber et al. (2022) Adapt to sample size.
    gamma <- 2 #maximum complexity level of the subgroups (here, a combination of two confounders max.)

    ## Input 5: sPoRT hyperparameters
    monotony <- TRUE #there is a monotone treatment pattern in the data. The function will automatically check positivity on untreated individuals.
    static <- c(TRUE, FALSE, FALSE) #d1 is static, d3 and d4 are dynamic. When dynamic, sPoRT will check the probability of being treated among those following the rule at time t AND the probability of being to be untreated among those not following the rule.

Now, we can run the `sport` function.

### sPoRT running and output

    # let's begin with stratification on time (argument pooling=FALSE)

    d1d2_strat_results <- sport(D.bar=D.bar[[1]], A=treat, lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooling=F, beta=beta, alpha=alpha, gamma=gamma, monotony=monotony, static=static[1])

    ## [1] "time 1"
    ## [1] "time 2"
    ## [1] "time 3"
    ## [1] "time 4"
    ## [1] "time 5"
    ## [1] "time 6"

    d3_strat_results <- sport(D.bar=D.bar[[2]], A=treat, lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooling=F, beta=beta, alpha=alpha, gamma=gamma, monotony=monotony, static=static[2])

    ## [1] "time 1"
    ## [1] "time 2"
    ## [1] "time 3"
    ## [1] "time 4"
    ## [1] "time 5"
    ## [1] "time 6"

    d4_strat_results <- sport(D.bar=D.bar[[3]], A=treat, lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooling=F, beta=beta, alpha=alpha, gamma=gamma, monotony=monotony, static=static[3])

    ## [1] "time 1"
    ## [1] "time 2"
    ## [1] "time 3"
    ## [1] "time 4"
    ## [1] "time 5"
    ## [1] "time 6"

The output of sPoRT is a list of tables. Each table represents the
violations identified at a certain time for a certain strategy. For
instance:

    d1d2_strat_results

    ## $T1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=980          0.004  exposed           915              18.3
    ## 2  L1_0< 400          0.993  exposed          1223              24.5
    ## 3 L1_1>=1040          0.003  exposed           896              17.9
    ## 4  L1_1< 460          0.995  exposed          1288              25.8
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
    ##                              subgroup proba.exposure exposure subgroup.size
    ## 1                           L1_0>=900          0.014  exposed          1185
    ## 2                           L2_0>=0.3          0.015  exposed           136
    ## 3                 L3_0>=0 & L3_0< 0.5          0.009  exposed           107
    ## 4                          L1_4>=1040          0.012  exposed          1534
    ## 5                          L2_4>=0.35          0.010  exposed           315
    ## 6                    Y_0>=-2 & W3< 36          0.014  exposed           346
    ## 7  Y_0>=-6 & W3>=18 & W3>=42 & W3< 54          0.015  exposed           136
    ## 8         L3_4>=-1 & W3< 48 & L3_4< 0          0.011  exposed           284
    ## 9          L3_4>=-1 & W3< 48 & W3>=36          0.008  exposed           123
    ## 10                  W1=west & Y_0>=-2          0.006  exposed           179
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
    ## 5      Y_0>=-7 & Y_0>=-6 & Y_0>=-4 & W3< 42          0.008  exposed
    ## 6                          L3_5>=0 & W3< 42          0.013  exposed
    ## 7              Y_0>=-7 & L2_0>=0.2 & Y_0< 0          0.015  exposed
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
    ##                                        subgroup proba.exposure exposure
    ## 1                                    L1_0>=1080          0.005  exposed
    ## 2                                    L1_6>=1140          0.009  exposed
    ## 3                                     L2_6>=0.4          0.011  exposed
    ## 4     L2_0>=0.15 & W3< 24 & L2_0< 0.25 & W3>=18          0.015  exposed
    ## 5                     L3_6>=1 & W3>=24 & W3>=36          0.014  exposed
    ## 6 L3_6>=0 & L2_0>=0.15 & L2_0>=0.2 & L2_0< 0.25          0.009  exposed
    ## 7                            L3_6>=0 & L3_0< -1          0.010  exposed
    ## 8                    L3_6>=1 & Y_0< 0 & Y_0>=-2          0.007  exposed
    ##   subgroup.size subgroup.rel.size
    ## 1           607              35.1
    ## 2          1301              75.2
    ## 3           182              10.5
    ## 4           137               7.9
    ## 5           146               8.4
    ## 6           115               6.7
    ## 7           199              11.5
    ## 8           136               7.9

are the violations identified at time 1 for the rule d1 (other results
at the end of this document).

These violations are defined in terms of subgroups, such as the
individuals with L1\_0&gt;=980 have 0.004% chance of being `A_1=1`,
which represents 915 individuals (18.3% of the sample size at time t).
The `exposure` column seems unnecessary here, but it is useful for
categorical treatment.

## sPoRT with pooling over time-points

### Pivoting to long format

When we want an analysis at the person-time level, one must check
positivity at the same scale. The sPoRT algorithm can handle long format
datasets as well. But first, we need to define the lagged values of the
treatment (for subseting on A<sub>t-1</sub>=0).

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
                      grep("Y_", names(simdata), value=T)[-1],
                      grep("C_", names(simdata), value=T),
                      grep("d1", names(simdata), value=T),
                      grep("d3", names(simdata), value=T),
                      grep("d4", names(simdata), value=T),
                      grep("Alag_", names(simdata), value=T)
                    ),
                    v.names = c("L1", "L2", "L3", "A", "Y", "C", "d1", "d3", "d4", "Alag"),
                    timevar = "T", 
                    times = 1:6, 
                    direction = "long")

### sPoRT’s main inputs

    ## Input 1: the treatment strategy columns' name. We can put several strategies through a vector
    D.bar <- c("d1", "d3", "d4")

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
    static <- c(TRUE, FALSE, FALSE) #d1 is static, d3 and d4 are dynamic. 

Now, we can run the `sport` function, with the argument `pooling=TRUE`.

### sPoRT running and output

    d1d2_pool_results <- sport(D.bar=D.bar[[1]], A="A", time='T', lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=pool, pooling=TRUE, alpha=alpha, beta=beta, gamma=gamma, static=static[1], monotony=monotony)

    d3_pool_results <- sport(D.bar=D.bar[[2]], A="A", time='T', lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=pool, pooling=TRUE, alpha=alpha, beta=beta, gamma=gamma, static=static[2], monotony=monotony)

    d4_pool_results <- sport(D.bar=D.bar[[3]], A="A", time='T', lag=0, type_A=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=pool, pooling=TRUE, alpha=alpha, beta=beta, gamma=gamma, static=static[3], monotony=monotony)

Results’ presentation is similar, but one sees only one table per rule
due to pooling over time.

    d1d2_pool_results

    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.002  exposed          2775              28.2
    ## 2   L1>=1160          0.004  exposed          5149              52.3

## Other inputs and miscenallous

The `sport`function can takes several other arguments. First, `lag` is
used to add lagged effects (e.g., `lag=1` means that we automatically
add the confounders from one time to the next, such as `L1_1` in the
second time). Second, `add.subset` can be specified when one wants to
subset more (e.g., A<sub>t-1</sub>=d<sub>t-1</sub> with the longitudinal
targeted maximum likelihood estimator). Third, `pruning` can be set to
`TRUE` when one focuses only on structural violations because this
removes violations bounded between two values from the output for
clarity (e.g., this may remove violations such as
`L1_0>=960 & L1_0< 860`). Last, the `rpart` hyperparameters can also be
set (Therneau, 2019). By default, the function use `minbucket = 6`,
`minsplit = 20`, and `maxdepth = 30` to allow small subgroups. Other
values can be set to split only smaller/bigger nodes.

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

### Supplementary results

Note that sPoRT for dynamic rules checks positivity in two subsets: the
individuals following the rule (`Aeq1`) and those not following the rule
(`Aeq0`). For tyhe first subset, a lack of positivity is present when
some subgroups have a tiny probability to be treated. In contrast, the
probability to be treated must be high in the second subset for
diagnosing positivity issues.

Stratified on time, rule 3:

    d3_strat_results

    ## $T1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=980          0.004  exposed           915              18.3
    ## 2  L1_0< 400          0.993  exposed          1223              24.5
    ## 3 L1_1>=1040          0.003  exposed           896              17.9
    ## 4  L1_1< 460          0.995  exposed          1288              25.8
    ## 
    ## $T2
    ## $T2$Aeq1
    ##    subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=940          0.008  exposed           735              40.4
    ## 2 L1_2>=960          0.012  exposed          1010              55.5
    ## 
    ## $T2$Aeq0
    ## [1] "The whole sample presents an exposure prevalence higher than 'beta'."
    ## 
    ## 
    ## $T3
    ## $T3$Aeq1
    ##                              subgroup proba.exposure exposure subgroup.size
    ## 1                           L1_0>=940          0.008  exposed           736
    ## 2                           L1_3>=960          0.016  exposed          1161
    ## 3 W3>=36 & Y_0< -2 & Y_0>=-8 & W3< 42          0.009  exposed           112
    ##   subgroup.rel.size
    ## 1              47.9
    ## 2              75.6
    ## 3               7.3
    ## 
    ## $T3$Aeq0
    ## [1] "The whole sample presents an exposure prevalence higher than 'beta'."
    ## 
    ## 
    ## $T4
    ## $T4$Aeq1
    ##                                        subgroup proba.exposure exposure
    ## 1                                     L1_0>=900          0.014  exposed
    ## 2                           L3_0>=0 & L3_0< 0.5          0.012  exposed
    ## 3                                    L1_4>=1000          0.018  exposed
    ## 4                                    L2_4>=0.25          0.014  exposed
    ## 5                   L2_0>=0.2 & W3>=36 & W3>=42          0.016  exposed
    ## 6                              W3< 54 & Y_0>=-2          0.016  exposed
    ## 7                   L3_4>=-1 & W3< 54 & L3_4< 0          0.015  exposed
    ## 8           L3_4>=-1 & W3< 54 & W3>=42 & W3< 48          0.013  exposed
    ## 9                              W1=west & W3>=30          0.018  exposed
    ## 10                L2_0>=0.15 & Y_0< 0 & Y_0>=-2          0.017  exposed
    ## 11          L2_0< 0.15 & L2_0>=0.05 & L2_0< 0.1          0.013  exposed
    ## 12                         L2_0>=0.15 & L3_4>=1          0.013  exposed
    ## 13 L3_4>=-1 & L2_0>=0.15 & L2_0< 0.25 & L3_4< 0          0.011  exposed
    ## 14                         W1=west & L2_0>=0.15          0.013  exposed
    ## 15                            W1=west & Y_0>=-2          0.008  exposed
    ##    subgroup.size subgroup.rel.size
    ## 1            839              59.5
    ## 2             83               5.9
    ## 3           1195              84.8
    ## 4            553              39.2
    ## 5            122               8.7
    ## 6            256              18.2
    ## 7            201              14.3
    ## 8             79               5.6
    ## 9            222              15.8
    ## 10           117               8.3
    ## 11            79               5.6
    ## 12            75               5.3
    ## 13            89               6.3
    ## 14           235              16.7
    ## 15           131               9.3
    ## 
    ## $T4$Aeq0
    ## [1] "The whole sample presents an exposure prevalence higher than 'beta'."
    ## 
    ## 
    ## $T5
    ## $T5$Aeq1
    ##                                            subgroup proba.exposure exposure
    ## 1                                        L1_0>=1120          0.000  exposed
    ## 2                                        L1_5>=1040          0.017  exposed
    ## 3                                           L3_5>=1          0.018  exposed
    ## 4                      L2_0>=0.15 & W3< 48 & W3< 30          0.011  exposed
    ## 5       W3< 48 & L3_0< -0.5 & L3_0< -1 & L3_0>=-1.5          0.014  exposed
    ## 6               W3< 48 & Y_0< -1 & W3>=24 & Y_0>=-4          0.015  exposed
    ## 7               L2_0>=0.15 & L3_0>=-0.5 & L2_0< 0.2          0.011  exposed
    ## 8                     L2_0>=0.15 & Y_0>=-1 & Y_0< 0          0.011  exposed
    ## 9  L2_5>=0.25 & L2_0< 0.25 & L2_0>=0.15 & L2_5< 0.3          0.016  exposed
    ## 10               L2_5< 0.2 & L2_0< 0.1 & L2_5>=0.15          0.013  exposed
    ## 11                L3_0>=0 & W1=southern & L3_0< 0.5          0.014  exposed
    ## 12                W2=female & L3_0>=-1 & L3_0< -0.5          0.018  exposed
    ## 13                             L2_5>=0.25 & Y_0>=-2          0.018  exposed
    ## 14                                Y_0< -3 & W1=west          0.015  exposed
    ##    subgroup.size subgroup.rel.size
    ## 1            347              25.8
    ## 2           1116              83.0
    ## 3            169              12.6
    ## 4             93               6.9
    ## 5             72               5.4
    ## 6            131               9.7
    ## 7             87               6.5
    ## 8             94               7.0
    ## 9            247              18.4
    ## 10            79               5.9
    ## 11            72               5.4
    ## 12           112               8.3
    ## 13           219              16.3
    ## 14           130               9.7
    ## 
    ## $T5$Aeq0
    ## [1] "The whole sample presents an exposure prevalence higher than 'beta'."
    ## 
    ## 
    ## $T6
    ## $T6$Aeq1
    ##                                          subgroup proba.exposure exposure
    ## 1                                      L1_0>=1080          0.007  exposed
    ## 2                                      L1_6>=1120          0.016  exposed
    ## 3        L2_0>=0.15 & L2_0< 0.2 & W3< 54 & W3>=18          0.000  exposed
    ## 4                    L3_0< -1.5 & W3>=48 & W3< 54          0.015  exposed
    ## 5        L2_6>=0.2 & W3< 24 & W3>=18 & L2_6< 0.35          0.017  exposed
    ## 6                     L2_6>=0.2 & W3>=24 & W3>=48          0.014  exposed
    ## 7                       L3_6>=1 & W3>=24 & W3>=36          0.018  exposed
    ## 8           L3_6>=-4 & W3< 54 & W3>=48 & L3_6< -1          0.015  exposed
    ## 9   L2_0>=0.15 & L2_0< 0.2 & L3_0< 0.5 & L3_0>=-2          0.000  exposed
    ## 10 L2_0< 0.2 & L2_0>=0.1 & L2_6>=0.3 & L2_0>=0.15          0.014  exposed
    ## 11  L3_6>=1 & L2_0< 0.25 & L2_0>=0.1 & L2_0>=0.15          0.018  exposed
    ## 12               L3_0>=-3 & L2_6>=0.3 & L3_0< 0.5          0.019  exposed
    ## 13                           L3_6>=1 & L3_0< -0.5          0.000  exposed
    ## 14                L3_6< 1 & L3_0< -1.5 & L3_6>=-1          0.019  exposed
    ## 15     L2_6>=0.2 & Y_0>=-8 & Y_0>=-2 & L2_6< 0.25          0.015  exposed
    ## 16                     L3_6>=1 & Y_0< 0 & Y_0>=-2          0.010  exposed
    ##    subgroup.size subgroup.rel.size
    ## 1            416              32.9
    ## 2            967              76.6
    ## 3             70               5.5
    ## 4             68               5.4
    ## 5            119               9.4
    ## 6             72               5.7
    ## 7            111               8.8
    ## 8             67               5.3
    ## 9             73               5.8
    ## 10           139              11.0
    ## 11           109               8.6
    ## 12           103               8.2
    ## 13            89               7.0
    ## 14           157              12.4
    ## 15            67               5.3
    ## 16           101               8.0
    ## 
    ## $T6$Aeq0
    ## [1] "The whole sample presents an exposure prevalence higher than 'beta'."

Stratified on time, rule 4:

    d4_strat_results

    ## $T1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=980          0.004  exposed           915              18.3
    ## 2  L1_0< 400          0.993  exposed          1223              24.5
    ## 3 L1_1>=1040          0.003  exposed           896              17.9
    ## 4  L1_1< 460          0.995  exposed          1288              25.8
    ## 
    ## $T2
    ## $T2$Aeq1
    ##                            subgroup proba.exposure exposure subgroup.size
    ## 1                         L1_0>=960          0.000  exposed           210
    ## 2                         L1_2>=940          0.023  exposed           351
    ## 3 L3_0< -2 & L3_0>=-2.5 & L2_0< 0.1          0.029  exposed            34
    ## 4                 W2=male & L3_0>=0          0.026  exposed            38
    ##   subgroup.rel.size
    ## 1              31.7
    ## 2              53.0
    ## 3               5.1
    ## 4               5.7
    ## 
    ## $T2$Aeq0
    ##                           subgroup proba.exposure exposure subgroup.size
    ## 1                        L1_0>=920          0.006  exposed           881
    ## 2                        L1_2>=960          0.007  exposed          1096
    ## 3                       L2_2>=0.35          0.008  exposed           130
    ## 4             L2_0>=0.3 & L3_0>=-2          0.010  exposed            99
    ## 5 L2_0>=0.25 & Y_0< -3 & L2_0< 0.3          0.009  exposed           108
    ## 6              Y_0>=-3 & L2_0>=0.3          0.011  exposed            94
    ##   subgroup.rel.size
    ## 1              52.3
    ## 2              65.0
    ## 3               7.7
    ## 4               5.9
    ## 5               6.4
    ## 6               5.6
    ## 
    ## 
    ## $T3
    ## $T3$Aeq1
    ##                                                       subgroup proba.exposure
    ## 1                                                    L1_0>=920          0.015
    ## 2                                                    L1_3>=920          0.018
    ## 3                                                    L2_3>=0.2          0.022
    ## 4                                 W3< 24 & L2_0>=0.05 & W3>=18          0.024
    ## 5                                 W3>=24 & L2_0>=0.05 & W3< 30          0.028
    ## 6                                   L3_0< -1 & W3< 48 & W3>=36          0.027
    ## 7                               L3_0< -1 & W3< 48 & L3_0>=-1.5          0.026
    ## 8                                             W3< 24 & Y_0>=-2          0.021
    ## 9                                   W3< 24 & Y_0< -2 & Y_0< -3          0.032
    ## 10                                   Y_0< -2 & W3>=36 & W3< 42          0.000
    ## 11                                   Y_0< -2 & W3>=48 & W3< 54          0.024
    ## 12                                           W3< 24 & L3_3>=-2          0.031
    ## 13                                  L3_3< -2 & W3< 54 & W3>=36          0.028
    ## 14                                            W3< 24 & W1=west          0.000
    ## 15                                          W3< 24 & W2=female          0.015
    ## 16                                   W2=male & W3>=36 & W3< 42          0.026
    ## 17 L2_0>=0.05 & L3_0< -1.5 & L3_0>=-2 & L2_0< 0.15 & L2_0< 0.1          0.033
    ## 18             L2_0>=0.05 & L3_0>=-3 & L3_0< -2.5 & L2_0< 0.15          0.033
    ## 19                               L3_0< -1 & Y_0< -3 & L3_0>=-2          0.000
    ## 20                              Y_0>=-3 & L3_0>=-2.5 & Y_0>=-2          0.031
    ## 21                 L3_0< -1.5 & L3_3>=-2 & L3_0>=-2 & L3_3< -1          0.029
    ## 22                         L3_0< -1.5 & L3_0>=-2 & W1=southern          0.030
    ## 23                             L3_0>=-3.5 & L3_0< -2 & W1=west          0.022
    ## 24                           L3_0< -1.5 & W2=female & L3_0>=-2          0.023
    ## 25                         L3_0>=-3.5 & W2=female & L3_0< -2.5          0.028
    ## 26                                          L3_3>=-2 & Y_0< -4          0.000
    ## 27                                 Y_0< -1 & L3_3>=0 & Y_0>=-2          0.024
    ## 28                                 Y_0>=-3 & Y_0< -2 & W1=west          0.033
    ## 29                               Y_0>=-5 & W2=female & Y_0< -4          0.031
    ##    exposure subgroup.size subgroup.rel.size
    ## 1   exposed           267              48.6
    ## 2   exposed           438              79.8
    ## 3   exposed            46               8.4
    ## 4   exposed            42               7.7
    ## 5   exposed            36               6.6
    ## 6   exposed            74              13.5
    ## 7   exposed            39               7.1
    ## 8   exposed            47               8.6
    ## 9   exposed            31               5.6
    ## 10  exposed            40               7.3
    ## 11  exposed            42               7.7
    ## 12  exposed            98              17.9
    ## 13  exposed            36               6.6
    ## 14  exposed            35               6.4
    ## 15  exposed            65              11.8
    ## 16  exposed            38               6.9
    ## 17  exposed            30               5.5
    ## 18  exposed            30               5.5
    ## 19  exposed            64              11.7
    ## 20  exposed            32               5.8
    ## 21  exposed            34               6.2
    ## 22  exposed            67              12.2
    ## 23  exposed            45               8.2
    ## 24  exposed            44               8.0
    ## 25  exposed            36               6.6
    ## 26  exposed            39               7.1
    ## 27  exposed            42               7.7
    ## 28  exposed            30               5.5
    ## 29  exposed            32               5.8
    ## 
    ## $T3$Aeq0
    ##                                   subgroup proba.exposure exposure
    ## 1                                L1_0>=940          0.005  exposed
    ## 2                                L1_3>=960          0.011  exposed
    ## 3                               L2_3>=0.35          0.009  exposed
    ## 4                                  L3_3>=1          0.009  exposed
    ## 5 L2_0>=0.25 & L2_0< 0.3 & W3< 54 & W3< 24          0.013  exposed
    ## 6        L2_0>=0.25 & L2_0< 0.3 & L3_0>=-1          0.000  exposed
    ## 7                     L2_0>=0.25 & Y_0< -2          0.015  exposed
    ## 8         L2_0>=0.25 & L2_0< 0.3 & W1=west          0.013  exposed
    ## 9       L2_0>=0.25 & L2_0< 0.3 & W2=female          0.008  exposed
    ##   subgroup.size subgroup.rel.size
    ## 1           811              54.3
    ## 2          1220              81.7
    ## 3           221              14.8
    ## 4           109               7.3
    ## 5            80               5.4
    ## 6            90               6.0
    ## 7           269              18.0
    ## 8            80               5.4
    ## 9           126               8.4
    ## 
    ## 
    ## $T4
    ## $T4$Aeq1
    ##                      subgroup proba.exposure exposure subgroup.size
    ## 1                  L1_0>=1000          0.000  exposed           189
    ## 2       L1_0< 800 & L1_0>=760          0.000  exposed            31
    ## 3                  L2_0< 0.05          0.000  exposed            28
    ## 4                     L3_0>=0          0.019  exposed            53
    ## 5           Y_0>=-2 & Y_0< -1          0.033  exposed           121
    ## 6                      Y_0>=0          0.034  exposed            29
    ## 7                  L1_4>=1000          0.021  exposed           422
    ## 8                   L2_4>=0.2          0.024  exposed           125
    ## 9  L3_4>=-2 & W3< 48 & W3>=36          0.021  exposed            97
    ## 10          L3_4>=-2 & W3< 48          0.029  exposed           103
    ## 11 W3< 48 & W3< 36 & L3_4>=-1          0.033  exposed            60
    ## 12 W3>=48 & L3_4>=0 & L3_4< 1          0.036  exposed            28
    ## 13       W3< 18 & W1=southern          0.026  exposed            38
    ## 14           W3>=30 & W1=west          0.033  exposed            61
    ## 15           W3< 24 & W2=male          0.031  exposed            65
    ## 16         W2=female & W3< 18          0.032  exposed            31
    ##    subgroup.rel.size
    ## 1               37.4
    ## 2                6.1
    ## 3                5.5
    ## 4               10.5
    ## 5               24.0
    ## 6                5.7
    ## 7               83.6
    ## 8               24.8
    ## 9               19.2
    ## 10              20.4
    ## 11              11.9
    ## 12               5.5
    ## 13               7.5
    ## 14              12.1
    ## 15              12.9
    ## 16               6.1
    ## 
    ## $T4$Aeq0
    ##                 subgroup proba.exposure exposure subgroup.size
    ## 1        W3< 30 & W3>=24          0.018  exposed           170
    ## 2              L1_0>=840          0.014  exposed          1040
    ## 3              L2_0>=0.3          0.015  exposed           136
    ## 4  L3_0>=-2.5 & L3_0< -2          0.005  exposed           200
    ## 5    L3_0< 0.5 & L3_0>=0          0.014  exposed            70
    ## 6                Y_0>=-1          0.013  exposed           230
    ## 7             L1_4>=1060          0.010  exposed          1110
    ## 8             L2_4>=0.35          0.010  exposed           315
    ## 9                L3_4>=1          0.007  exposed           134
    ## 10               W1=west          0.018  exposed           393
    ##    subgroup.rel.size
    ## 1               12.2
    ## 2               74.5
    ## 3                9.7
    ## 4               14.3
    ## 5                5.0
    ## 6               16.5
    ## 7               79.5
    ## 8               22.6
    ## 9                9.6
    ## 10              28.2
    ## 
    ## 
    ## $T5
    ## $T5$Aeq1
    ##                                          subgroup proba.exposure exposure
    ## 1                                          W3< 18          0.019  exposed
    ## 2                                       L1_0>=880          0.029  exposed
    ## 3                                      L1_5>=1040          0.021  exposed
    ## 4                                      L2_5>=0.25          0.032  exposed
    ## 5                                         L3_5>=1          0.015  exposed
    ## 6                                        L3_5< -3          0.033  exposed
    ## 7  L2_0< 0.1 & L2_0>=0.05 & L3_0>=-2.5 & L3_0< -2          0.000  exposed
    ## 8              L3_0< -0.5 & L2_0>=0.15 & L3_0>=-2          0.000  exposed
    ## 9                             L2_0< 0.1 & Y_0< -3          0.031  exposed
    ## 10               L2_0< 0.1 & L2_0>=0.05 & W1=west          0.025  exposed
    ## 11                   Y_0< -3 & L3_0< -1 & Y_0>=-4          0.000  exposed
    ## 12      Y_0>=-5 & L3_0< -1 & Y_0< -4 & L3_0>=-2.5          0.028  exposed
    ## 13                             W1=west & L3_0< -2          0.022  exposed
    ## 14              W1=southern & L3_0< 0.5 & L3_0>=0          0.030  exposed
    ## 15               W2=female & L3_0>=-1.5 & L3_0< 0          0.031  exposed
    ## 16                  W2=male & L3_0< -1 & L3_0>=-3          0.034  exposed
    ## 17                              W1=west & Y_0< -3          0.000  exposed
    ##    subgroup.size subgroup.rel.size
    ## 1             53              11.1
    ## 2            308              64.3
    ## 3            389              81.2
    ## 4             31               6.5
    ## 5             65              13.6
    ## 6             30               6.3
    ## 7             24               5.0
    ## 8             26               5.4
    ## 9             65              13.6
    ## 10            40               8.4
    ## 11            46               9.6
    ## 12            36               7.5
    ## 13            45               9.4
    ## 14            33               6.9
    ## 15            98              20.5
    ## 16            29               6.1
    ## 17            42               8.8
    ## 
    ## $T5$Aeq0
    ##                                    subgroup proba.exposure exposure
    ## 1                                 L1_0>=860          0.015  exposed
    ## 2                                   L3_0>=0          0.000  exposed
    ## 3                                L1_5>=1000          0.018  exposed
    ## 4                                L2_5>=0.35          0.017  exposed
    ## 5                                   L3_5>=0          0.015  exposed
    ## 6                        W3< 48 & L2_0>=0.2          0.000  exposed
    ## 7   W3< 48 & L2_0>=0.2 & L2_0< 0.3 & W3>=42          0.013  exposed
    ## 8   L2_0>=0.2 & L2_0< 0.3 & W3>=24 & W3< 30          0.014  exposed
    ## 9  L2_0< 0.2 & L2_0>=0.15 & W3< 36 & W3>=24          0.010  exposed
    ## 10      Y_0>=-6 & W3< 48 & W3>=24 & Y_0< -3          0.019  exposed
    ## 11       W3< 48 & Y_0< 1 & W3>=24 & Y_0>=-1          0.019  exposed
    ## 12            W3>=24 & W3< 36 & W1=southern          0.017  exposed
    ## 13              W3< 48 & W2=female & W3>=30          0.009  exposed
    ## 14                W2=male & W3>=24 & W3< 30          0.012  exposed
    ## 15                      L2_0>=0.2 & Y_0>=-2          0.017  exposed
    ## 16         L2_0< 0.2 & L2_0>=0.15 & Y_0>=-1          0.013  exposed
    ## 17                      L2_0>=0.2 & W1=west          0.016  exposed
    ## 18           Y_0< 1 & Y_0>=-1 & W1=southern          0.019  exposed
    ## 19              Y_0>=-4 & Y_0< -3 & W1=west          0.011  exposed
    ## 20             Y_0< 1 & Y_0>=-1 & W2=female          0.010  exposed
    ##    subgroup.size subgroup.rel.size
    ## 1            972              72.3
    ## 2            107               8.0
    ## 3           1228              91.4
    ## 4            362              26.9
    ## 5            388              28.9
    ## 6             99               7.4
    ## 7             78               5.8
    ## 8             74               5.5
    ## 9             96               7.1
    ## 10           212              15.8
    ## 11           108               8.0
    ## 12           237              17.6
    ## 13           223              16.6
    ## 14            82               6.1
    ## 15           297              22.1
    ## 16            77               5.7
    ## 17           244              18.2
    ## 18           160              11.9
    ## 19            88               6.5
    ## 20           100               7.4
    ## 
    ## 
    ## $T6
    ## $T6$Aeq1
    ##                                subgroup proba.exposure exposure subgroup.size
    ## 1                       W3< 54 & W3>=42          0.037  exposed           107
    ## 2                             L1_0>=940          0.033  exposed           246
    ## 3                 L3_0>=-2.5 & L3_0< -2          0.033  exposed            61
    ## 4                            L1_6>=1120          0.019  exposed           324
    ## 5                L2_6>=0.2 & L2_6< 0.25          0.037  exposed           108
    ## 6                               L3_6>=1          0.028  exposed            72
    ## 7      L2_0>=0.15 & Y_0>=-4 & L2_0< 0.2          0.026  exposed            38
    ## 8            Y_0>=-5 & Y_0< 1 & Y_0< -1          0.034  exposed            58
    ## 9  L2_0>=0.15 & L2_0< 0.2 & W1=southern          0.038  exposed            53
    ## 10               L2_0>=0.15 & W2=female          0.024  exposed            42
    ## 11          Y_0< -2 & W1=west & Y_0>=-3          0.000  exposed            23
    ##    subgroup.rel.size
    ## 1               23.9
    ## 2               54.9
    ## 3               13.6
    ## 4               72.3
    ## 5               24.1
    ## 6               16.1
    ## 7                8.5
    ## 8               12.9
    ## 9               11.8
    ## 10               9.4
    ## 11               5.1
    ## 
    ## $T6$Aeq0
    ##                                               subgroup proba.exposure exposure
    ## 1                                            L1_0>=940          0.014  exposed
    ## 2                                           L1_6>=1160          0.005  exposed
    ## 3                                            L2_6>=0.4          0.011  exposed
    ## 4             W3< 24 & L2_0< 0.3 & L2_0>=0.15 & W3>=18          0.016  exposed
    ## 5                                  W3< 24 & L3_0< -1.5          0.012  exposed
    ## 6                  W3< 24 & Y_0>=-5 & Y_0< -2 & W3>=18          0.017  exposed
    ## 7                                     L3_6>=0 & W3>=36          0.018  exposed
    ## 8                            L3_6< 0 & W3< 24 & W3>=18          0.016  exposed
    ## 9                                     W3< 24 & W1=west          0.010  exposed
    ## 10                                L2_0>=0.3 & L3_0>=-2          0.011  exposed
    ## 11                                 L3_6>=0 & L2_0>=0.2          0.017  exposed
    ## 12                             L2_0>=0.3 & W1=southern          0.012  exposed
    ## 13 Y_0>=-5 & L3_0< 0 & L3_0>=-4 & L3_0< -1.5 & Y_0< -3          0.012  exposed
    ## 14                                  L3_6>=0 & L3_0< -1          0.007  exposed
    ## 15          L3_6< 0 & L3_0< -1.5 & L3_6>=-2 & L3_0>=-2          0.008  exposed
    ## 16                     L3_0>=-2 & L3_0< -1.5 & W1=west          0.015  exposed
    ## 17                     L3_0>=-2 & L3_0< -1.5 & W2=male          0.018  exposed
    ## 18                                   Y_0>=-2 & L3_6>=1          0.008  exposed
    ##    subgroup.size subgroup.rel.size
    ## 1            776              60.6
    ## 2            950              74.2
    ## 3            182              14.2
    ## 4            127               9.9
    ## 5            165              12.9
    ## 6            120               9.4
    ## 7            227              17.7
    ## 8            125               9.8
    ## 9             96               7.5
    ## 10            90               7.0
    ## 11           242              18.9
    ## 12            85               6.6
    ## 13            86               6.7
    ## 14           148              11.6
    ## 15           131              10.2
    ## 16            68               5.3
    ## 17           113               8.8
    ## 18           129              10.1

Pooled over time, rule 3:

    d3_pool_results

    ## $Aeq1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1080          0.003  exposed          2066              28.0
    ## 2   L1>=1140          0.006  exposed          3820              51.8
    ## 
    ## $Aeq0
    ##                 subgroup proba.exposure exposure subgroup.size
    ## 1                Y_0>=-2          0.009  exposed           866
    ## 2              L1_0>=980          0.004  exposed          1361
    ## 3 L2_0>=0.3 & L2_0< 0.35          0.012  exposed           501
    ## 4             L3_0>=-0.5          0.004  exposed           455
    ## 5                 W3< 24          0.010  exposed           684
    ## 6               L1>=1020          0.008  exposed          2082
    ## 7                L2>=0.3          0.013  exposed          1988
    ## 8                  L3>=0          0.002  exposed           622
    ##   subgroup.rel.size
    ## 1              35.0
    ## 2              55.0
    ## 3              20.3
    ## 4              18.4
    ## 5              27.7
    ## 6              84.2
    ## 7              80.4
    ## 8              25.2

Pooled over time, rule 4:

    d4_pool_results

    ## $Aeq1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1060          0.006  exposed           680              25.7
    ## 2   L1>=1120          0.010  exposed          1307              49.5
    ## 
    ## $Aeq0
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.001  exposed          2214              30.7
    ## 2   L1>=1140          0.005  exposed          4243              58.9

### Updates history

First update planned for the release of the `port` package.
