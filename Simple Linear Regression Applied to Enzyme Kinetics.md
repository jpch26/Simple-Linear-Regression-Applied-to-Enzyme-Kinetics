R in the Lab: Simple Linear Regression Applied to Enzyme Kinetics
================

## Problem

You have performed a
<a href="https://en.wikipedia.org/wiki/Enzyme_kinetics" target="_blank">enzyme
kinetics</a> experiment. You are supposing a
<a href="https://www.ncbi.nlm.nih.gov/books/NBK22430/" target="_blank">Michaelis-Menten</a>
kinetics and your aim is to determine the constants Vmax, Km and some
other related.

![](images/enzyme_experiment.jpg) The first ten observations in your
data would look like this:

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>

<thead>

<tr>

<th style="text-align:right;">

S (M)

</th>

<th style="text-align:right;">

v (M/min)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

8.0e-06

</td>

<td style="text-align:right;">

1.2e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

8.0e-06

</td>

<td style="text-align:right;">

1.5e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

8.0e-06

</td>

<td style="text-align:right;">

1.1e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

1.8e-05

</td>

<td style="text-align:right;">

2.5e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

1.8e-05

</td>

<td style="text-align:right;">

2.7e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

1.8e-05

</td>

<td style="text-align:right;">

2.4e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

2.8e-05

</td>

<td style="text-align:right;">

3.2e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

2.8e-05

</td>

<td style="text-align:right;">

3.2e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

2.8e-05

</td>

<td style="text-align:right;">

3.3e-09

</td>

</tr>

<tr>

<td style="text-align:right;">

3.8e-05

</td>

<td style="text-align:right;">

3.9e-09

</td>

</tr>

</tbody>

</table>

The substrate units are
<a href="https://en.wikipedia.org/wiki/Molar_concentration" target="_blank">molarity</a>
expressed simply as M (mol/L), and the velocity units are molarity over
minute since you are measuring the change in the substrate concentration
per time unit.

By the way, you can see that the numbers are expressed as , e.g.,
8.0e-6. This is the way that R uses to represent scientific notation, so
the previous number can be read as 8.0x10<sup>-6</sup>.

The Michelis-Menten equation describes the reaction velocity as a
function of the substrate concentration:

![](images/equation%201.png)

## Solution

A regular approach to obtain Vmax and Km is to convert the
Michaelis-Menten equation into a “linear form” like:

![](images/equation%202.png)

Once you obtained your experimental data the next step is to fit a
linear model using data transformations (the velocity and the substrate
concentration inverse), then take the model coefficients and perform
simple calculations to obtain the constants.

Let’s do this and a nice plot. I also want to show you an alternative
approach that not require a linear transformation and, hence, not any
data transformations.

As always, in the first place I need to show you how I have organized my
analysis project:

![](images/files_organization.png)

  - The analysis folder contains all my R scripts.  
  - The data folder contains the experimental data and all CSV and TXT
    files products of my analysis scripts.  
  - The graph folder contains all the plots results of my analysis
    scripts.

These three folders are contained in a main folder that can be named as,
e.g., “enzyme\_experiment”, or you can specified the date of your
analysis like “enzime\_experiment\_4FEB2021”.

### Linear regression

The R analysis script for this analysis:

``` r
# Linear regression analysis to enzyme kinetics ------------------------------

# Packages 
if (!"ggplot2" %in% .packages()) library(ggplot2)
if (!"dplyr" %in% .packages()) library(dplyr)

# 1 Import data -----------------------------------------------------------

kinetic_data <- read.csv("data/kinetic_data.csv")

# 2 Data transformations ---------------------------------------------------

kinetic_data <- kinetic_data %>% 
  mutate(
    S_i = 1/S,
    v_i = 1/v
    )

# 2.1 Save data transformations
write.csv(kinetic_data, "data/kinetic_data_trans.csv", row.names = FALSE)

# 3 Fitting linear model -------------------------------------------------

# 3.1 Fit linear model 
model_lm <- lm(v_i ~ S_i, data = kinetic_data)

# 3.2 Model summary
model_sum <- summary(model_lm)

# 3.3 Model summary
capture.output(file = "data/model_summary.txt", model_sum)

# 4 Results report -------------------------------------------------------

# 4.1 Obtain the model coefficients 
model_coeff <- coef(model_lm)

# 4.2 Calculate Vmax
v_max <- unname(model_coeff[1]^-1)

# 4.3 Calculate km 
km <- unname(model_coeff[2] * v_max)

# 4.4 Enzyme concentration
enz = 4.5e-8

# 4.5 Calculate kcat and kcat/km 
kcat <- v_max /enz
kcat_km <- kcat / km

# 4.6 Make a data frame with the results
kin_report <- data.frame(
  Vmax = v_max,
  Km   = signif(km, 2),
  Enz  = enz,
  kcat = signif(kcat, 2),
  kcat_Km = signif(kcat_km, 2) 
)

# 4.7 Kinetic report
write.csv(kin_report, "data/kinetic_report.csv", row.names = FALSE)

# 5 Plot for transformed data and linear model -------------------------------

# 5.1 Make the plot usin ggplot2
model_plot_lm <- ggplot(kinetic_data, aes(S_i, v_i)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  xlab("1/[S] (1/M)") +
  ylab("1/v (min/M)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15)
  )

# 5.2 Save Plot
ggsave("graphs/kinetic_plot.jpeg", model_plot_lm)
```

The next code checks if the necessary packages are already loaded and,
if they don’t, then loaded them.

``` r
# Packages 
if (!"ggplot2" %in% .packages()) library(ggplot2)
if (!"dplyr" %in% .packages()) library(dplyr)
```

You can omit the if statement while working on your own.

Then, I imported the experimental data and performed a data
transformation using the `dplyr` function `mutate()`:

``` r
# 1 Import data -----------------------------------------------------------

kinetic_data <- read.csv("data/kinetic_data.csv")

# 2 Data transformations ---------------------------------------------------

kinetic_data <- kinetic_data %>% 
  mutate(
    S_i = 1/S,
    v_i = 1/v
    )
```

The transformation was simply the inverse of both the concentration and
the velocity reaction. The following table shows you the first ten rows
of the previous data.

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>

<thead>

<tr>

<th style="text-align:right;">

S

</th>

<th style="text-align:right;">

v

</th>

<th style="text-align:right;">

S\_i

</th>

<th style="text-align:right;">

v\_i

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

8.0e-06

</td>

<td style="text-align:right;">

1.2e-09

</td>

<td style="text-align:right;">

125000.00

</td>

<td style="text-align:right;">

828359716

</td>

</tr>

<tr>

<td style="text-align:right;">

8.0e-06

</td>

<td style="text-align:right;">

1.5e-09

</td>

<td style="text-align:right;">

125000.00

</td>

<td style="text-align:right;">

648934642

</td>

</tr>

<tr>

<td style="text-align:right;">

8.0e-06

</td>

<td style="text-align:right;">

1.1e-09

</td>

<td style="text-align:right;">

125000.00

</td>

<td style="text-align:right;">

873355158

</td>

</tr>

<tr>

<td style="text-align:right;">

1.8e-05

</td>

<td style="text-align:right;">

2.5e-09

</td>

<td style="text-align:right;">

55555.56

</td>

<td style="text-align:right;">

401078092

</td>

</tr>

<tr>

<td style="text-align:right;">

1.8e-05

</td>

<td style="text-align:right;">

2.7e-09

</td>

<td style="text-align:right;">

55555.56

</td>

<td style="text-align:right;">

365033479

</td>

</tr>

<tr>

<td style="text-align:right;">

1.8e-05

</td>

<td style="text-align:right;">

2.4e-09

</td>

<td style="text-align:right;">

55555.56

</td>

<td style="text-align:right;">

418003857

</td>

</tr>

<tr>

<td style="text-align:right;">

2.8e-05

</td>

<td style="text-align:right;">

3.2e-09

</td>

<td style="text-align:right;">

35714.29

</td>

<td style="text-align:right;">

310241761

</td>

</tr>

<tr>

<td style="text-align:right;">

2.8e-05

</td>

<td style="text-align:right;">

3.2e-09

</td>

<td style="text-align:right;">

35714.29

</td>

<td style="text-align:right;">

312616064

</td>

</tr>

<tr>

<td style="text-align:right;">

2.8e-05

</td>

<td style="text-align:right;">

3.3e-09

</td>

<td style="text-align:right;">

35714.29

</td>

<td style="text-align:right;">

307573858

</td>

</tr>

<tr>

<td style="text-align:right;">

3.8e-05

</td>

<td style="text-align:right;">

3.9e-09

</td>

<td style="text-align:right;">

26315.79

</td>

<td style="text-align:right;">

255222351

</td>

</tr>

</tbody>

</table>

Next, I fit the linear model on the transformed data following the
linear relation which I showed you previously (equation 2):

``` r
# 3 Fitting linear model -------------------------------------------------

# 3.1 Fit linear model 
model_lm <- lm(v_i ~ S_i, data = kinetic_data)

# 3.2 Model summary
model_sum <- summary(model_lm)

# 3.3 Model summary
capture.output(file = "data/model_summary.txt", model_sum)
```

I made a nice summary using the function `summary()`. You can type
`model_sum`to see these results (or look for the “model\_summary.txt”
file in the data folder):

``` r
model_sum
```

    ## 
    ## Call:
    ## lm(formula = v_i ~ S_i, data = kinetic_data)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -127152864   -3710613    1211299    4143150   97267652 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 1.219e+08  3.572e+06   34.12   <2e-16 ***
    ## S_i         5.234e+03  1.091e+02   47.98   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 22990000 on 61 degrees of freedom
    ## Multiple R-squared:  0.9742, Adjusted R-squared:  0.9738 
    ## F-statistic:  2302 on 1 and 61 DF,  p-value: < 2.2e-16

In the fourth step, I calculated the Vmax and Km constants, and using
the enzyme concentration in the experiment I also calculated the
constants kcat and kcat/Km that can be useful in you results discussion.
I put all these constants into a data frame and saved it in the data
folder as a CSV file (“kinetic\_report.csv”).

``` r
# 4 Results report -------------------------------------------------------

# 4.1 Obtain the model coefficients 
model_coeff <- coef(model_lm)

# 4.2 Calculate Vmax
v_max <- unname(model_coeff[1]^-1)

# 4.3 Calculate km 
km <- unname(model_coeff[2] * v_max)

# 4.4 Enzyme concentration
enz = 4.5e-8

# 4.5 Calculate kcat and kcat/km 
kcat <- v_max /enz
kcat_km <- kcat / km

# 4.6 Make a data frame with the results
kin_report <- data.frame(
  Vmax = v_max,
  Km   = signif(km, 2),
  Enz  = enz,
  kcat = signif(kcat, 2),
  kcat_Km = signif(kcat_km, 2) 
)

# 4.7 Kinetic report
write.csv(kin_report, "data/kinetic_report.csv", row.names = FALSE)
```

The results are these:

``` r
kin_report_lm
```

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>

<thead>

<tr>

<th style="text-align:right;">

Vmax

</th>

<th style="text-align:right;">

Km

</th>

<th style="text-align:right;">

Enz

</th>

<th style="text-align:right;">

kcat

</th>

<th style="text-align:right;">

kcat\_Km

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

8.2e-09

</td>

<td style="text-align:right;">

4.3e-05

</td>

<td style="text-align:right;">

4.5e-08

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:right;">

4200

</td>

</tr>

</tbody>

</table>

To calculate Vmax and Km firstly I took the model coefficients:

``` r
coef(model_lm)
```

    ##  (Intercept)          S_i 
    ## 1.219002e+08 5.233499e+03

From equation 2 the coefficients are \((Intercept) = 1/Vmax\) and
\(S_i = Km/Vmax\) (the slope of the line).

Finally I also plotted the transformed data and represented the fitted
model as a blue line. For this I used the `ggplot2` package.

``` r
# 5 Plot for transformed data and linear model -------------------------------

# 5.1 Make the plot usin ggplot2
model_plot_lm <- ggplot(kinetic_data, aes(S_i, v_i)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  xlab("1/[S] (1/M)") +
  ylab("1/v (min/M)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15)
  )

# 5.2 Save Plot
ggsave("graphs/kinetic_plot.jpeg", model_plot)
```

I saved this plot as “kinetic\_plot.jpeg” in my graph folder.

``` r
model_plot_lm
```

![](Simple%20Linear%20Regression%20Applied%20to%20Enzyme%20Kinetics_files/figure-gfm/show%20linear%20model%20plot-1.png)<!-- -->

### Nonlinear Regression Analysis

Instead of transforming the equation 1 into a linear form, we can obtain
the constants in a straight way using the function `nls()`.

``` r
# Linear regression analysis to enzyme kinetics ------------------------------

# Packages 
if (!"ggplot2" %in% .packages()) library(ggplot2)
if (!"dplyr" %in% .packages()) library(dplyr)

# 1 Import data -----------------------------------------------------------
kinetic_data <- read.csv("data/kinetic_data.csv")

# 2 Fitting non linear model -------------------------------------------------

# 2.1 Fit non linear model
model_nls <- nls(
  v ~ (Vmax * S)/(Km + S), 
  data = kinetic_data, 
  start = c(Vmax = 0.01, Km = 0)
)

# 2.2 Model summary
model_sum_nls <- summary(model_nls)

# 2.3 Save model summary
capture.output(file = "data/model_summary_nls.txt", model_sum_nls)

# 3 Kinetic report ---------------------------------------------------------

# 3.1 Enzyme concentration (M)
enz <- 4.5e-8

# 3.2 Model coefficients
coeffs_nls <- unname(coef(model_nls))

# 3.3 Calculate kcat and kcat/km 
kcat <- coeffs_nls[1]/enz
kcat_km <- kcat / coeffs_nls[2]

# 3.4 Make a data frame with the results
kin_report_nls <- data.frame(
  Vmax = signif(coeffs_nls[1], 2), 
  Km = signif(coeffs_nls[2], 2),
  Enz = enz, 
  kcat = signif(kcat, 2), 
  kcat_km = signif(kcat_km, 2)
  )

# 3.5 Save kinetic report
write.csv(kin_report_nls, "data/kinetic_report_nls.csv", row.names = FALSE)

# 4 Plot for transformed data and linear model -------------------------------

# 4.1 Make plot with ggplot2
model_plot_nls <- ggplot(kinetic_data, aes(S, v)) +
  geom_point() +
  geom_smooth(method = "nls", se = FALSE, formula = y ~ a*x/(b+x),
              method.args = list(start = c(a = 0.01, b = 0))) +
  xlab("[S] (M)") +
  ylab("v (M/min)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15)
  )

# 4.2 Save plot
ggsave("graphs/kinetic_plot_nls.jpeg", model_plot_nls)
```

The steps I followed were similar than those in the linear model script.
First, I checked and loaded the necessary packages, I imported the data
and in the second step I fitted the non linear model:

``` r
# 2 Fitting non linear model -------------------------------------------------

# 2.1 Fit non linear model
model_nls <- nls(
  v ~ (Vmax * S)/(Km + S), 
  data = kinetic_data, 
  start = c(Vmax = 0.01, Km = 0)
)

# 2.2 Model summary
model_sum_nls <- summary(model_nls)

# 2.3 Save model summary
capture.output(file = "data/model_summary_nls.txt", model_sum_nls)
```

Note that I specified the model as the equation 1 (`v ~ (Vmax * S)/(Km +
S)`). The function `nls()` also needs start values for both
coefficients. You can see more information about this functions by
typing `?nls`.

I saved the results in my data folder as “model\_summary\_nls.txt”. You
can simply type `model_sum_nls` to see this summary:

``` r
model_sum_nls
```

    ## 
    ## Formula: v ~ (Vmax * S)/(Km + S)
    ## 
    ## Parameters:
    ##       Estimate Std. Error t value Pr(>|t|)    
    ## Vmax 8.038e-09  7.324e-11  109.75   <2e-16 ***
    ## Km   4.037e-05  1.291e-06   31.26   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.525e-10 on 61 degrees of freedom
    ## 
    ## Number of iterations to convergence: 7 
    ## Achieved convergence tolerance: 1.262e-06

As you can see, this is straight way to obtain the Vmax and Km
constants.

On the third step I made a data frame with all the results and saved it
as “kinetic\_report\_nls.csv” in my data folder:

``` r
# 3 Kinetic report ---------------------------------------------------------

# 3.1 Enzyme concentration (M)
enz <- 4.5e-8

# 3.2 Model coefficients
coeffs_nls <- unname(coef(model_nls))

# 3.3 Calculate kcat and kcat/km 
kcat <- coeffs_nls[1]/enz
kcat_km <- kcat / coeffs_nls[2]

# 3.4 Make a data frame with the results
kin_report_nls <- data.frame(
  Vmax = signif(coeffs_nls[1], 2), 
  Km = signif(coeffs_nls[2], 2),
  Enz = enz, 
  kcat = signif(kcat, 2), 
  kcat_km = signif(kcat_km, 2)
  )

# 3.5 Save kinetic report
write.csv(kin_report_nls, "data/kinetic_report_nls.csv", row.names = FALSE)
```

This is the table with all the calculated constants:

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>

<thead>

<tr>

<th style="text-align:right;">

Vmax

</th>

<th style="text-align:right;">

Km

</th>

<th style="text-align:right;">

Enz

</th>

<th style="text-align:right;">

kcat

</th>

<th style="text-align:right;">

kcat\_km

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

8e-09

</td>

<td style="text-align:right;">

4e-05

</td>

<td style="text-align:right;">

4.5e-08

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:right;">

4400

</td>

</tr>

</tbody>

</table>

If you compare this table with the obtained in the linear regression
section you will notice that there are not big differences. Later I will
compare both methods taking into consideration their prediction
performance.

Finally I made a plot using the data and represented the non linear
model as a blue line, note that this way allow us to represent a model
completely related to the original data.

``` r
# 4 Plot for transformed data and linear model -------------------------------

# 4.1 Make plot with ggplot2
model_plot_nls <- ggplot(kinetic_data, aes(S, v)) +
  geom_point() +
  geom_smooth(method = "nls", se = FALSE, formula = y ~ a*x/(b+x),
              method.args = list(start = c(a = 0.01, b = 0))) +
  xlab("[S] (M)") +
  ylab("v (M/min)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15)
  )

# 4.2 Save plot
ggsave("graphs/kinetic_plot_nls.jpeg", model_plot_nls)
```

In the function `geom_smooth()` I designated the `method`, the `formula`
and the `method.args` arguments. This was similar to the code where I
fitted the non linear model using the function `nls()`.

``` r
model_plot_nls
```

![](Simple%20Linear%20Regression%20Applied%20to%20Enzyme%20Kinetics_files/figure-gfm/model%20plot%20nls-1.png)<!-- -->

I saved this plot in my graph folder as “kinetic\_plot\_nls.jpeg”.

### Comparing both Methods

First, we can simply compare the tables with the constant results.

Results from the linear regression approach:

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>

<thead>

<tr>

<th style="text-align:right;">

Vmax

</th>

<th style="text-align:right;">

Km

</th>

<th style="text-align:right;">

Enz

</th>

<th style="text-align:right;">

kcat

</th>

<th style="text-align:right;">

kcat\_Km

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

8.2e-09

</td>

<td style="text-align:right;">

4.3e-05

</td>

<td style="text-align:right;">

4.5e-08

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:right;">

4200

</td>

</tr>

</tbody>

</table>

Results from the nonlinear regression approach:

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>

<thead>

<tr>

<th style="text-align:right;">

Vmax

</th>

<th style="text-align:right;">

Km

</th>

<th style="text-align:right;">

Enz

</th>

<th style="text-align:right;">

kcat

</th>

<th style="text-align:right;">

kcat\_km

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

8e-09

</td>

<td style="text-align:right;">

4e-05

</td>

<td style="text-align:right;">

4.5e-08

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:right;">

4400

</td>

</tr>

</tbody>

</table>

Yeah, there is no a big difference, but what would you think if I tell
you that the real values for Vmax and Km are exactly those obtained in
the nonlinear regression? I mean, I stated these values when I simulated
the experimental data.

Other way to compare the models is by determining how good they are to
make predictions, namely, how far the reaction velocity predictions are
from the experimental values?

I made this evaluation with the next script:

``` r
# Compare linear and non linear models performance ------------------------

# 1 Data ------------------------------------------------------------------
data <- read.csv("data/kinetic_data_trans.csv")

# 2 Linear model ----------------------------------------------------------

# 2.1 Fit model
model_lm <- lm(v_i ~ S_i, data)

# 2.2 Predicted values for linear model
predicts_lm <- predict(model_lm) 
predicts_lm <- predicts_lm^-1

# 2.3 Average sum of squared residuals
squared_res_lm <- (data$v - predicts_lm) ^ 2
sum_res_prom_lm <- sum(squared_res_lm) / (nrow(data)-length(coef(model_lm)))

# 2.4 RSME for linear model
rsme_lm <- sqrt(sum_res_prom_lm)

# 3 Non linear model --------------------------------------------------------

# 3.1 Fit non linear model
model_nlm <- model_nls <- nls(
  v ~ (Vmax * S)/(Km + S), 
  data = kinetic_data, 
  start = c(Vmax = 0.01, Km = 0)
)

# 3.2 Predicted values for non linear model 
predicts_nlm <- predict(model_nlm)

# 3.3 Average sum of squared residuals
squared_res_nlm <- (data$v - predicts_nlm) ^ 2
sum_res_prom_nlm <- sum(squared_res_nlm) / (nrow(data)-length(coef(model_nlm)))

# 3.4 RSME for non linear model
rsme_nlm <- sqrt(sum_res_prom_nlm)
```

In both cases I fitted the model. In the linear model I had to transform
the predictions to the original units with the code:

``` r
# 2.2 Predicted values for linear model
predicts_lm <- predict(model_lm) 
predicts_lm <- predicts_lm^-1
```

Remember, to fit the linear model I’ve used the equation 2, so firstly I
transformed the substrate concentrations and the reaction velocity data.

Next, I calculated the residuals (real value - predicted value), I
squared these values, added them up, and divided by the degrees of
freedom to obtain an average of these deviations.

``` r
# 2.3 Average sum of squared residuals
squared_res_lm <- (data$v - predicts_lm) ^ 2
sum_res_prom_lm <- sum(squared_res_lm) / (nrow(data)-length(coef(model_lm)))
```

The residual standard error or RMSE is the square root of the previous
value, so we have a measurement of the average deviations in the
original response units.

``` r
# 2.4 RSME for linear model
rsme_lm <- sqrt(sum_res_prom_lm)
```

I followed the same pattern for the nonlinear regression:

``` r
# 3.1 Fit non linear model
model_nlm <- model_nls <- nls(
  v ~ (Vmax * S)/(Km + S), 
  data = kinetic_data, 
  start = c(Vmax = 0.01, Km = 0)
)

# 3.2 Predicted values for non linear model 
predicts_nlm <- predict(model_nlm)

# 3.3 Average sum of squared residuals
squared_res_nlm <- (data$v - predicts_nlm) ^ 2
sum_res_prom_nlm <- sum(squared_res_nlm) / (nrow(data)-length(coef(model_nlm)))

# 3.4 RSME for non linear model
rsme_nlm <- sqrt(sum_res_prom_nlm)
```

I must point out that you can see the residual standard error by looking
at each model summary, but for more clarity I decided to show you the
way you can calculate it, besides the residuals for the lineal fitting
were not in the original velocity data units.

To compare both RSME you can simply divide them:

``` r
rsme_lm / rsme_nlm
```

    ## [1] 1.040217

So we can conclude that both models have a similar performance
estimating reaction velocity values. Anyway, I think that the nonlinear
regression is a better approach since offer a straight way to estimate
Vmax and Km without any data transformation.

That’s it\! You can clone the repository with the code and results of
this R tutorial, including the experimental data simulation:

<a href="https://github.com/jpch26/Simple-Linear-Regression-Applied-to-Enzyme-Kinetics" target="_blank">Simple
Linear Regression Applied to Enzyme Kinetics</a>

Try to reproduce the analysis step by step, modify or improve the code.
It’s all yours\!
