Sequential Positivity-Regression Trees (sPoRT) notebook
================
Arthur Chatton
2024-07-18

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

1.  Initiate treatment immediately (S1)
2.  Never initiate treatment (S2)
3.  Initiate treatment when `L1` drops below 750 or when `L2` drops
    below 0.25 (S3)
4.  Initiate treatment when `L1` drops below 350 or when `L2` drops
    below 0.15 (S4)

<!-- -->

    #define number of time-points
    tps <- 1:6

    for(t in tps){
      
      #Creation of the strategies' dummy variables
      if(t==tps[1]){ #first time-point, no censoring yet
        
        simdata[,paste0('S1_',t)] <- 1
        simdata[,paste0('S3_',t)] <- 1*(simdata[,paste0('L1_',t)] < 750 | simdata[,paste0('L2_',t)] < 0.25)
        simdata[,paste0('S4_',t)] <- 1*(simdata[,paste0('L1_',t)] < 350 | simdata[,paste0('L2_',t)] < 0.15)
        
      }else{ #other time-points, censoring may occurs and should be considered
        
        simdata[,paste0('S1_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1)
        simdata[,paste0('S3_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1*(simdata[,paste0('L1_',t)] < 750 | simdata[,paste0('L2_',t)] < 0.25 | simdata[,paste0('S3_',tps[which(t==tps)-1])] == 1))
        simdata[,paste0('S4_',t)] <- ifelse(simdata[,paste0('C_',tps[which(t==tps)-1])]==1, NA, 1*(simdata[,paste0('L1_',t)] < 350 | simdata[,paste0('L2_',t)] < 0.15 | simdata[,paste0('S4_',tps[which(t==tps)-1])] == 1))
        
      }
      
    }

Note that we didn’t code the strategy S2 because this strategy will
presents the same positivity issues than S1. Here, we pool over
treatment strategies, but we can add more stratification by including
the previous treatment in `cov.quali`. For sustained treatment
strategies as in our example, it doesn’t make sense but this can be
sometimes welcomed.

### sPoRT’s main inputs

    ## Input 1: the treatment strategy columns' name. We can put several strategies through a list
    regi <- list(S1=grep("S1", names(simdata), value=T),
                 S3=grep("S3", names(simdata), value=T),
                 S4=grep("S4", names(simdata), value=T)
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
    alpha <- 0.0233 #to have at least 100 individuals by subgroups at the last time for the static strategy S1. Personal preference, other values are meaningful as well.
    beta <- "gruber" #to use Gruber's bound defined in Gruber et al. (2022) Adapt to sample size.
    gamma <- 2 #maximum complexity level of the subgroups (here combination of two confounders max.)

Now, we can run the `sport` function.

### sPoRT running and output

    results <- sport(regimen=regi, exposure=treat, lag=0, type_expo=type, cov.quanti=cov.quanti, cov.quali=cov.quali, data=simdata, pooled=F, beta=beta, alpha=alpha, gamma=gamma)

    ## [1] "time 1 and groupe 1"
    ## [1] "time 2 and groupe 1"
    ## [1] "time 3 and groupe 1"
    ## [1] "time 4 and groupe 1"
    ## [1] "time 5 and groupe 1"
    ## [1] "time 6 and groupe 1"
    ## [1] "time 1 and groupe 2"
    ## [1] "time 2 and groupe 2"
    ## [1] "time 3 and groupe 2"
    ## [1] "time 4 and groupe 2"
    ## [1] "time 5 and groupe 2"
    ## [1] "time 6 and groupe 2"
    ## [1] "time 1 and groupe 3"
    ## [1] "time 2 and groupe 3"
    ## [1] "time 3 and groupe 3"
    ## [1] "time 4 and groupe 3"
    ## [1] "time 5 and groupe 3"
    ## [1] "time 6 and groupe 3"

    results

    ## $S1
    ## $S1$T1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=960          0.004  exposed           915              18.3
    ## 2  L1_0< 400          0.993  exposed          1223              24.5
    ## 3 L1_1>=1020          0.003  exposed           896              17.9
    ## 4  L1_1< 460          0.995  exposed          1288              25.8
    ## 
    ## $S1$T2
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=960          0.007  exposed           915              20.0
    ## 2  L1_0< 520          0.992  exposed          1377              30.1
    ## 3 L1_2>=1080          0.007  exposed           867              19.0
    ## 4  L1_2< 580          0.999  exposed          1185              25.9
    ## 
    ## $S1$T3
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1000          0.006  exposed           787              17.8
    ## 2  L1_0< 520          0.996  exposed          1219              27.6
    ## 3 L1_3>=1160          0.006  exposed           830              18.8
    ## 4  L1_3< 700          0.998  exposed          1299              29.4
    ## 
    ## $S1$T4
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1160          0.000  exposed           397               9.1
    ## 2  L1_0< 540          0.994  exposed          1288              29.6
    ## 3 L1_4>=1280          0.003  exposed           683              15.7
    ## 4  L1_4< 800          0.997  exposed          1405              32.3
    ## 
    ## $S1$T5
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.004  exposed           513              11.9
    ## 2  L1_0< 520          0.997  exposed          1145              26.5
    ## 3 L1_5>=1300          0.006  exposed           682              15.8
    ## 4  L1_5< 840          0.998  exposed          1480              34.2
    ## 
    ## $S1$T6
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.008  exposed           513              11.9
    ## 2  L1_0< 560          0.995  exposed          1335              31.1
    ## 3 L1_6>=1360          0.003  exposed           595              13.9
    ## 4  L1_6< 920          0.994  exposed          1742              40.5
    ## 
    ## 
    ## $S3
    ## $S3$T1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=960          0.004  exposed           915              18.3
    ## 2  L1_0< 400          0.993  exposed          1223              24.5
    ## 3 L1_1>=1020          0.003  exposed           896              17.9
    ## 4  L1_1< 460          0.995  exposed          1288              25.8
    ## 
    ## $S3$T2
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.000  exposed           333               8.3
    ## 2  L1_0< 520          0.992  exposed          1377              34.4
    ## 3 L1_2>=1120          0.004  exposed           488              12.2
    ## 4  L1_2< 580          0.999  exposed          1185              29.6
    ## 
    ## $S3$T3
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1000          0.009  exposed           542              14.0
    ## 2  L1_0< 520          0.996  exposed          1219              31.6
    ## 3 L1_3>=1160          0.005  exposed           551              14.3
    ## 4  L1_3< 700          0.998  exposed          1299              33.6
    ## 
    ## $S3$T4
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.000  exposed           344               9.0
    ## 2  L1_0< 520          0.997  exposed          1172              30.7
    ## 3 L1_4>=1280          0.002  exposed           447              11.7
    ## 4  L1_4< 800          0.997  exposed          1405              36.8
    ## 
    ## $S3$T5
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.000  exposed           346               9.1
    ## 2  L1_0< 580          0.990  exposed          1457              38.4
    ## 3 L1_5>=1300          0.007  exposed           451              11.9
    ## 4  L1_5< 880          0.992  exposed          1680              44.2
    ## 
    ## $S3$T6
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.006  exposed           347               9.2
    ## 2  L1_0< 560          0.995  exposed          1335              35.4
    ## 3 L1_6>=1360          0.005  exposed           392              10.4
    ## 4  L1_6< 920          0.995  exposed          1731              45.9
    ## 
    ## 
    ## $S4
    ## $S4$T1
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=960          0.004  exposed           915              18.3
    ## 2  L1_0< 400          0.993  exposed          1223              24.5
    ## 3 L1_1>=1020          0.003  exposed           896              17.9
    ## 4  L1_1< 460          0.995  exposed          1288              25.8
    ## 
    ## $S4$T2
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=940          0.010  exposed           196              10.2
    ## 2  L1_0< 540          0.996  exposed           979              51.0
    ## 3 L2_0>=0.15          0.987  exposed            79               4.1
    ## 4 L1_2>=1100          0.007  exposed           141               7.3
    ## 5  L1_2< 660          0.996  exposed          1017              53.0
    ## 6  L2_2>=0.2          1.000  exposed            60               3.1
    ## 
    ## $S4$T3
    ##                subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1            L1_0>=1040          0.000  exposed           127               6.8
    ## 2             L1_0< 660          0.987  exposed          1134              61.0
    ## 3            L1_3>=1180          0.007  exposed           145               7.8
    ## 4             L1_3< 780          0.998  exposed          1041              56.0
    ## 5    L2_3>=0.2 & W3>=18          1.000  exposed            61               3.3
    ## 6 L2_3>=0.2 & L2_0>=0.1          0.986  exposed            74               4.0
    ## 
    ## $S4$T4
    ##                 subgroup proba.exposure exposure subgroup.size
    ## 1              L1_0>=980          0.011  exposed           188
    ## 2              L1_0< 660          0.988  exposed          1094
    ## 3             L1_4>=1280          0.008  exposed           128
    ## 4              L1_4< 900          0.987  exposed          1126
    ## 5 L2_4>=0.2 & L2_0>=0.15          1.000  exposed            54
    ##   subgroup.rel.size
    ## 1              10.2
    ## 2              59.4
    ## 3               6.9
    ## 4              61.1
    ## 5               2.9
    ## 
    ## $S4$T5
    ##                 subgroup proba.exposure exposure subgroup.size
    ## 1             L1_0>=1040          0.007  exposed           137
    ## 2              L1_0< 660          0.993  exposed          1074
    ## 3             L1_5>=1280          0.007  exposed           140
    ## 4              L1_5< 960          0.985  exposed          1187
    ## 5 L2_0>=0.1 & L2_5>=0.25          1.000  exposed            54
    ##   subgroup.rel.size
    ## 1               7.5
    ## 2              58.7
    ## 3               7.7
    ## 4              64.9
    ## 5               3.0
    ## 
    ## $S4$T6
    ##                 subgroup proba.exposure exposure subgroup.size
    ## 1              L1_0< 660          0.992  exposed          1059
    ## 2             L1_6>=1340          0.008  exposed           118
    ## 3              L1_6< 960          0.992  exposed          1121
    ## 4 L2_6>=0.2 & L2_0>=0.15          1.000  exposed            48
    ##   subgroup.rel.size
    ## 1              58.0
    ## 2               6.5
    ## 3              61.4
    ## 4               2.6

The output of sPoRT is a list of table. Each table represents the
violations identified at a certain time for a certain strategy. For
instance:

    results$S1$T1

    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1  L1_0>=960          0.004  exposed           915              18.3
    ## 2  L1_0< 400          0.993  exposed          1223              24.5
    ## 3 L1_1>=1020          0.003  exposed           896              17.9
    ## 4  L1_1< 460          0.995  exposed          1288              25.8

are the four violations identified at time 1 for the strategy S1.

These violations are defined in terms of subgroups, such as the
individuals with L1\_0&gt;=960 have 0.4% chance of being `A_1=1`, which
represents 915 individuals (18.3% of the whole sample size). The
`exposure` column seems unnecessary here, but it is useful for
categorical treatment.

## sPoRT with pooling over time-points

### Pivoting to long format

When we want an analysis at the person-time level, one must checks
positivity at the same scale. The sPoRT algorithm can handle long format
dataset as well.

    pool <- reshape(simdata, 
                    varying = list(
                      grep("L1", names(simdata), value=T)[-1],
                      grep("L2", names(simdata), value=T)[-1],
                      grep("L3", names(simdata), value=T)[-1],
                      grep("A_", names(simdata), value=T),
                      grep("Y_", names(simdata), value=T)[-1],
                      grep("C_", names(simdata), value=T),
                      grep("S1", names(simdata), value=T),
                      grep("S3", names(simdata), value=T),
                      grep("S4", names(simdata), value=T)
                    ),
                    v.names = c("L1", "L2", "L3", "A", "Y", "C", "S1", "S3", "S4"),
                    timevar = "T", 
                    times = 1:6, 
                    direction = "long")

### sPoRT’s main inputs

    ## Input 1: the treatment strategy columns' name. We can put several strategies through a vector
    regi <- c("S1", "S3", "S4")

    ## Input 2: qualitative and quantitative confounders
    cov.quanti <- c("Y_0", "L1_0",  "L2_0", "L3_0", "W3", "L1", "L2", "L3") 
    cov.quali <- c("W1", "W2") 

    ## Input 3: Actual treatment column's name and variable type
    treat <- "A"
    type_expo = "b" #binary treatment

    ## Input 4: sPoRT hyperparameters
    alpha <- 0.01 #to have at least 300 person-time by subgroups for strategy S1. Personal preference, other values are meaningful as well.
    beta <- "gruber" #to use Gruber's bound defined in Gruber et al. (2022) Adapt to sample size.
    gamma <- 2 #maximum complexity level of the subgroups (here combination of two confounders max.)

Now, we can run the `sport` function, with the argument `pooled=TRUE`.

### sPoRT running and output

    pool_results <- sport(regimen=c("S1", "S3", "S4"), exposure="A", time='T', id="id", lag=0, type_expo="b", cov.quanti=cov.quanti, cov.quali=cov.quali, data=pool, pooled=TRUE, alpha=alpha, beta=beta, gamma=gamma)
    pool_results

    ## [[1]]
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.003  exposed          3078              11.4
    ## 2  L1_0< 400          0.998  exposed          4667              17.3
    ## 3   L1>=1360          0.001  exposed          2367               8.8
    ## 4    L1< 460          0.998  exposed          2780              10.3
    ## 
    ## [[2]]
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.001  exposed          2224               9.2
    ## 2  L1_0< 400          0.998  exposed          4667              19.2
    ## 3   L1>=1360          0.002  exposed          1583               6.5
    ## 4    L1< 460          0.998  exposed          2780              11.5
    ## 
    ## [[3]]
    ##     subgroup proba.exposure exposure subgroup.size subgroup.rel.size
    ## 1 L1_0>=1100          0.002  exposed          1016               7.1
    ## 2  L1_0< 400          0.998  exposed          3915              27.4
    ## 3   L1>=1380          0.002  exposed           464               3.2
    ## 4    L1< 460          0.998  exposed          2665              18.7

Results’ presentation is similar; but one sees only one table per
strategy due to polling over time. The results are very close between
the strategies in this example, but we can see that the subgroups varies
across strategies due to the various number of individuals following the
strategy at the previous time.

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
