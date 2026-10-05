# SCIM311 Chapter 1: R code of the chapter
# Run from top to bottom to reproduce every figure in the chapter.
# Part of the SCIM311 lecture notes.

## ----students-----------------------------------------------------------------
students <- data.frame(
  student = paste0("S", 1:10),
  height  = c(165, 170, 160, 175, 168, 172, 158, 180, 169, 162),  # cm
  weight  = c(58, 65, 55, 72, 60, 68, 54, 75, 62, 57),            # kg
  gpa     = c(2.8, 3.5, 2.9, 3.8, 3.1, 3.6, 2.7, 3.9, 3.3, 3.0)
)




## ----cross-sectional, fig.cap = "Cross-sectional data: ten students measured once."----
plot(students$weight, students$gpa, pch = 19, col = "darkgreen",
     main = "Cross-sectional data: weight vs GPA", ylim = c(2.6, 4.0),
     xlab = "Weight (kg)", ylab = "GPA")
text(students$weight, students$gpa, students$student, pos = 3, cex = 0.7)


## ----s1-gpa-------------------------------------------------------------------
s1 <- data.frame(semester = 1:8,
                 gpa = c(2.8, 3.0, 3.2, 3.4, 3.5, 3.6, 3.7, 3.8))




## ----time-series, fig.cap = "Time series data: one student (S1) measured repeatedly."----
plot(s1$semester, s1$gpa, type = "o", pch = 19, col = "blue",
     main = "Time series data: student S1 GPA",
     xlab = "Semester", ylab = "GPA")


## ----regression, fig.cap = "Regression method: sales of product B against purchases of product A (simulated)."----
set.seed(11)
purchases_A <- seq(0, 10, length.out = 30)
sales_B <- 1 + 2.3 * purchases_A + rnorm(30, sd = 4)
fit_AB <- lm(sales_B ~ purchases_A)

plot(purchases_A, sales_B, pch = 19, col = "blue",
     main = "Regression method: two product sales",
     xlab = "Purchases of product A", ylab = "Sales of product B")
abline(fit_AB, col = "red", lwd = 2)
legend("topleft", bty = "n", pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2),
       col = c("blue", "red"),
       legend = c("Observed data",
                  sprintf("Fitted line: y = %.2f + %.2f x",
                          coef(fit_AB)[1], coef(fit_AB)[2])))


## ----components, fig.height = 6, fig.cap = "Four series with different components. Data sets built into R: austres, co2, sunspot.year, lh."----
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Trend: quarterly Australian residents (thousands), 1971-1993
plot(austres, type = "o", pch = 20, main = "Trend",
     xlab = "Year", ylab = "Residents (1000s)")

# Trend + seasonality: monthly CO2 at Mauna Loa (ppm), 1959-1997
plot(co2, main = "Trend + seasonality",
     xlab = "Year", ylab = expression(CO[2] ~ "(ppm)"))

# Cycles: yearly mean sunspot numbers, 1700-1988
plot(sunspot.year, main = "Cycles",
     xlab = "Year", ylab = "Sunspot number")

# Irregular only: luteinizing hormone, 48 samples at 10-minute intervals
plot(lh, type = "o", pch = 20, main = "Irregular fluctuation",
     xlab = "Sample (10-minute intervals)", ylab = "Hormone level")

par(mfrow = c(1, 1))


## ----additive-multiplicative, fig.height = 5.5, fig.cap = "The same trend and seasonal pattern combined additively and multiplicatively."----
t <- 1:60
trend <- 10 + 2 * t
season <- sin(2 * pi * t / 12)

Y_add <- trend + 5 * season          # additive model
Y_mul <- trend * (1 + 0.3 * season)  # multiplicative model

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(t, Y_add, type = "l", col = "blue", lwd = 2,
     main = "Additive model", xlab = "Time", ylab = "Y (additive)")
plot(t, Y_mul, type = "l", col = "darkorange", lwd = 2,
     main = "Multiplicative model", xlab = "Time", ylab = "Y (multiplicative)")
par(mfrow = c(1, 1))


## ----iid-noise, fig.cap = "Simulated iid N(0, 1) noise."----------------------
set.seed(1)
n <- 200
x <- rnorm(n, mean = 0, sd = 1)  # iid N(0, 1)
plot(ts(x), main = "IID noise: N(0, 1)", xlab = "t", ylab = expression(x[t]))
abline(h = 0, col = "grey50", lty = 2)


## ----random-walk, fig.cap = "A simulated random walk started at zero."--------
set.seed(2)
n <- 200
eps <- rnorm(n, mean = 0, sd = 1)
S <- cumsum(eps)  # random walk: S_t = S_{t-1} + eps_t
plot(ts(S), main = "Random walk", xlab = "t", ylab = expression(S[t]))
abline(h = 0, col = "grey50", lty = 2)




## ----lake-huron, fig.cap = "Level of Lake Huron (feet above 570 ft) with the least-squares trend line."----
level <- as.numeric(LakeHuron) - 570  # feet above 570 ft
t <- seq_along(level)                 # t = 1, ..., 98
fit_lake <- lm(level ~ t)
coef(fit_lake)

plot(1875:1972, level, type = "o", pch = 0,
     main = "Level of Lake Huron, 1875-1972",
     xlab = "Year", ylab = "Level (feet above 570 ft)")
lines(1875:1972, fitted(fit_lake), lwd = 2)


## ----rw-three-trends, fig.cap = "Three simulated realisations of the same random walk; each seems to have its own trend."----
set.seed(31)
walks <- replicate(3, cumsum(rnorm(60)))   # three random walks, t = 1..60
matplot(walks, type = "l", lty = 1, lwd = 2,
        col = c("black", "steelblue", "firebrick"),
        main = "Three realisations of one random walk",
        xlab = "t", ylab = expression(S[t]))
abline(h = 0, col = "grey50", lty = 2)


## ----nottem-trend, fig.cap = "Monthly air temperature at Nottingham, 1920-1939 (first eight years shown), with the fitted cosine trend (red) and the seasonal means (blue)."----
t <- as.numeric(time(nottem))   # 1920, 1920.083, ... (years)
fit_cos <- lm(nottem ~ cos(2 * pi * t) + sin(2 * pi * t))  # 3 parameters
month <- factor(cycle(nottem))  # month 1, ..., 12
fit_means <- lm(nottem ~ month - 1)                         # 12 parameters

round(coef(fit_cos), 2)
round(c(cosine = summary(fit_cos)$sigma,
        means = summary(fit_means)$sigma), 2)   # residual SD

show <- t < 1928
plot(t[show], nottem[show], type = "p", pch = 20,
     main = "Nottingham temperature with two seasonal trends",
     xlab = "Year", ylab = "Temperature (F)")
lines(t[show], fitted(fit_cos)[show], col = "red", lwd = 2)
lines(t[show], fitted(fit_means)[show], col = "blue", lty = 2)


## ----plot-template, eval = FALSE----------------------------------------------
## # my_ts_data is any time series (ts) object
## plot(my_ts_data, main = "Time plot of my data", xlab = "Time", ylab = "Value")


## ----air-base, fig.cap = "Time plot of AirPassengers drawn with base R."------
my_ts_data <- AirPassengers

plot(my_ts_data,
     main = "Time plot of AirPassengers",
     xlab = "Year", ylab = "Passengers (1000s)")


## ----sales-ts, fig.cap = "Monthly sales, January 2020 to December 2022."------
sales_data <- c(10, 12, 15, 13, 17, 20, 18, 22, 25, 23, 28, 30,
                12, 14, 17, 15, 19, 22, 20, 24, 27, 25, 30, 32,
                15, 16, 19, 18, 21, 24, 23, 26, 29, 28, 33, 35)

# Monthly data starting in January 2020, three years
sales_ts <- ts(sales_data, start = c(2020, 1), frequency = 12)
print(sales_ts)
plot(sales_ts, main = "Monthly sales", xlab = "Year", ylab = "Sales")


## ----decompose, fig.height = 6, fig.cap = "Classical additive decomposition of the monthly sales."----
decomposed_sales_add <- decompose(sales_ts, type = "additive")
plot(decomposed_sales_add)

# If a multiplicative model seems more appropriate:
# decomposed_sales_mult <- decompose(sales_ts, type = "multiplicative")
# plot(decomposed_sales_mult)


## ----stl, fig.height = 6, fig.cap = "STL decomposition of the monthly sales."----
stl_sales <- stl(sales_ts, s.window = "periodic")
plot(stl_sales)

# stl() is part of base R (package stats); no extra package is needed.




## ----sensor-stream, fig.cap = "A simulated sensor stream with a missing block (minutes 45-49) and a spike (minute 85)."----
set.seed(311)
n <- 120
t <- 1:n
sensor <- 20 + 0.01 * t + sin(2 * pi * t / 24) + rnorm(n, 0, 0.15)

sensor_with_events <- sensor
sensor_with_events[45:49] <- NA                # communication failure
sensor_with_events[85] <- sensor[85] + 4       # unusual physical event

plot(sensor_with_events, type = "o", pch = 16,
     main = "One stream, two different problems",
     xlab = "Minute", ylab = "Sensor value")
abline(v = c(45, 49, 85), col = c("steelblue", "steelblue", "firebrick"),
       lty = 2)

