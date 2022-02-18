Group Project1\_Monte Carlo Simulation Design.
================
hw2849
P8160 Advanced Statistical Computing

## Project 2: Design a simulation study to assess the three hypothsis testings

Design a set of distributions/models under both proportional-hazard and
non-proportional-hazard assumptions, and carry out a simulation study to
compare performance of those hypothesis tests in those models. Based on
your numerical investigations, write a practical recommendation for
general users to choose a suitable testing tool.

### Non-proportional

Under non-proportional test assumption, we generate data set of a sample
size n = 1000, from two different distributions, exponential
distribution (p = 500) and Weibull distribution (p = 500). We have a
time variable, treatment and control group, and an event variable
indicating censoring.

``` r
# sampling survival times from exponential distribution & Weibull distribution
set.seed(8160)
p = 500
treatment = rep(1, 500)
control = rep(0, 500)

test_data = data.frame(
  time = c(rexp(p), rweibull(p, 0.3, 0.1)), 
  group = c(treatment, control)
) %>% 
  mutate(
    event = case_when(time >= 4 ~ 0, 
                      time < 4 ~ 1)
  )

plot(survfit(Surv(time, event) ~ group, data = test_data), 
     xlab = "Time", 
     ylab = "Overall survival probability")
```

![](Group_Project_1_Simulation_Studies_1_files/figure-gfm/non-proportional%20data%20simulation-1.png)<!-- -->

``` r
# Exponential
model_exp = survreg(Surv(time, event) ~ group, test_data, dist = "exponential")
summary(model_exp)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(time, event) ~ group, data = test_data, 
    ##     dist = "exponential")
    ##               Value Std. Error     z       p
    ## (Intercept) -0.3098     0.0457 -6.77 1.3e-11
    ## group        0.3915     0.0644  6.08 1.2e-09
    ## 
    ## Scale fixed at 1 
    ## 
    ## Exponential distribution
    ## Loglik(model)= -857.8   Loglik(intercept only)= -876.2
    ##  Chisq= 36.73 on 1 degrees of freedom, p= 1.4e-09 
    ## Number of Newton-Raphson Iterations: 8 
    ## n= 1000

``` r
## Weibull
model_weibull = survreg(Surv(time, event) ~ group, test_data, dist = "weibull", scale = 1)
summary(model_weibull)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(time, event) ~ group, data = test_data, 
    ##     dist = "weibull", scale = 1)
    ##               Value Std. Error     z       p
    ## (Intercept) -0.3098     0.0457 -6.77 1.3e-11
    ## group        0.3915     0.0644  6.08 1.2e-09
    ## 
    ## Scale fixed at 1 
    ## 
    ## Weibull distribution
    ## Loglik(model)= -857.8   Loglik(intercept only)= -876.2
    ##  Chisq= 36.73 on 1 degrees of freedom, p= 1.4e-09 
    ## Number of Newton-Raphson Iterations: 8 
    ## n= 1000

#### Non-proportional: log-rank and weighted log-rank tests

``` r
## log-rank test
logrank.test(test_data$time, test_data$event, test_data$group, rho = 0, gamma = 0)
```

    ## Call:
    ## logrank.test(time = test_data$time, event = test_data$event, 
    ##     group = test_data$group, rho = 0, gamma = 0)
    ## 
    ##     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## 1 500      478      323      74.0       115
    ## 2 500      488      643      37.2       115
    ## 
    ##  Chisq= 115  on 1 degrees of freedom, p= <2e-16
    ##  rho   =  0 gamma =  0

``` r
## weighted log-rank test for an early effect 
logrank.test(test_data$time, test_data$event, test_data$group, rho = 1, gamma = 0)
```

    ## Call:
    ## logrank.test(time = test_data$time, event = test_data$event, 
    ##     group = test_data$group, rho = 1, gamma = 0)
    ## 
    ##     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## 1 500      478      323      74.0       318
    ## 2 500      488      643      37.2       318
    ## 
    ##  Chisq= 305  on 1 degrees of freedom, p= <2e-16
    ##  rho   =  1 gamma =  0

``` r
## weighted log-rank test for a late effect 
logrank.test(test_data$time, test_data$event, test_data$group, rho = 0, gamma = 1)
```

    ## Call:
    ## logrank.test(time = test_data$time, event = test_data$event, 
    ##     group = test_data$group, rho = 0, gamma = 1)
    ## 
    ##     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## 1 500      478      323      74.0       374
    ## 2 500      488      643      37.2       374
    ## 
    ##  Chisq= 0.2  on 1 degrees of freedom, p= 0.7
    ##  rho   =  0 gamma =  1

``` r
# take maximum over the previous three log-rank with multiple comparison control
logrank.maxtest(test_data$time, test_data$event, test_data$group)
```

    ## Call:
    ## logrank.maxtest(time = test_data$time, event = test_data$event, 
    ##     group = test_data$group)
    ## 
    ##  Two sided p-value = 0 (Bonferroni corrected: 0)
    ## 
    ##  Individual weighted log-rank tests:
    ##   Test      z     p
    ## 1    1 10.726 0.000
    ## 2    2  0.406 0.685
    ## 3    3 17.453 0.000

**Decisions**

*Log-rank test:* The Chi-Squared test statistic is 115 with 1 degree of
freedom and the corresponding p-value is 0. Since this p-value is less
than .05, we reject the null hypothesis.

*Weighted log-rank test for an early effect:* The Chi-Squared test
statistic is 305 with 1 degree of freedom and the corresponding p-value
is 0. Since this p-value is less than 0.05, we reject the null
hypothesis.

*Weighted log-rank test for a late effect:* The Chi-Squared test
statistic is 0.2 with 1 degree of freedom and the corresponding p-value
is 0.7. Since this p-value is greater than 0.05, we fail to reject the
null hypothesis.

### Proportional

Under proportional assumption, we generate survival data with
`sim.survdata` with variables `group`, `time`, and `fail`.

-   group: 1 = treatment, 0 = control

``` r
set.seed(8160)

## simulating data
sim_data = sim.survdata(N = 500, T = 100, xvars = 1, censor = .4, num.data.frames = 2) 

## data frame
test_df = data.frame( #trt = 1, ctl = 0
  time_1 = sim_data[[1]]$data$y,
  time_0 = sim_data[[2]]$data$y,
  fail_1 = sim_data[[1]]$data$failed,
  fail_0 = sim_data[[2]]$data$failed
) %>% 
  pivot_longer( # make df readable
    time_1:fail_0,
    names_to = c(".value", "group"),
    names_sep = "_"
  ) %>% 
  mutate(group = as.numeric(group))

## plot survival curves for each group
plot(survfit(Surv(time, fail) ~ group, data = test_df), 
     xlab = "Time", 
     ylab = "Overall survival probability")
```

![](Group_Project_1_Simulation_Studies_1_files/figure-gfm/proportional%20simulation-1.png)<!-- -->

#### Proportional: log-rank tests

``` r
## log-rank test
logrank.test(test_df$time, test_df$fail, test_df$group, rho = 0, gamma = 0)
```

    ## Call:
    ## logrank.test(time = test_df$time, event = test_df$fail, group = test_df$group, 
    ##     rho = 0, gamma = 0)
    ## 
    ##     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## 1 500      277      194      35.0      61.7
    ## 2 500      301      384      17.8      61.7
    ## 
    ##  Chisq= 61.7  on 1 degrees of freedom, p= 4e-15
    ##  rho   =  0 gamma =  0

``` r
## weighted log-rank test for an early effect 
logrank.test(test_df$time, test_df$fail, test_df$group, rho = 1, gamma = 0)
```

    ## Call:
    ## logrank.test(time = test_df$time, event = test_df$fail, group = test_df$group, 
    ##     rho = 1, gamma = 0)
    ## 
    ##     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## 1 500      277      194      35.0       119
    ## 2 500      301      384      17.8       119
    ## 
    ##  Chisq= 62.1  on 1 degrees of freedom, p= 3e-15
    ##  rho   =  1 gamma =  0

``` r
## weighted log-rank test for a late effect 
logrank.test(test_df$time, test_df$fail, test_df$group, rho = 0, gamma = 1)
```

    ## Call:
    ## logrank.test(time = test_df$time, event = test_df$fail, group = test_df$group, 
    ##     rho = 0, gamma = 1)
    ## 
    ##     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## 1 500      277      194      35.0       374
    ## 2 500      301      384      17.8       374
    ## 
    ##  Chisq= 28.8  on 1 degrees of freedom, p= 8e-08
    ##  rho   =  0 gamma =  1

``` r
# take maximum over the previous three log-rank with multiple comparison control
logrank.maxtest(test_df$time, test_df$fail, test_df$group)
```

    ## Call:
    ## logrank.maxtest(time = test_df$time, event = test_df$fail, group = test_df$group)
    ## 
    ##  Two sided p-value = 4.44e-15 (Bonferroni corrected: 9.99e-15)
    ## 
    ##  Individual weighted log-rank tests:
    ##   Test    z        p
    ## 1    1 7.85 4.00e-15
    ## 2    2 5.37 7.96e-08
    ## 3    3 7.88 3.33e-15

**Decisions**

*Log-rank test:* The Chi-Squared test statistic is 61.7 with 1 degree of
freedom and the corresponding p-value is 0. Since this p-value is less
than .05, we reject the null hypothesis.

*Weighted log-rank test for an early effect:* The Chi-Squared test
statistic is 62.1 with 1 degree of freedom and the corresponding p-value
is 0. Since this p-value is less than 0.05, we reject the null
hypothesis.

*Weighted log-rank test for a late effect:* The Chi-Squared test
statistic is 28.8 with 1 degree of freedom and the corresponding p-value
is 0. Since this p-value is greater than 0.05, we reject the null
hypothesis.

type I: false positive type II: false negative\*

repeat tests under the same setting
