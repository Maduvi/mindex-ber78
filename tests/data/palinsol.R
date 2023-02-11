##
## use R library "palinsol" to produce test data
## Maduvi
##

library(palinsol)
library(ggplot2)
library(cowplot)

## time vector
time <- seq(-500e3, 0, 1e3)

## astronomical insolation during summer/winter
N65sum <- sapply(time, function(x) c(tim = x, Insol_l1l2(ber78(x), l1=0, l2=pi, lat=65 * pi /180, S0=1361.371, avg=TRUE)))
N65win <- sapply(time, function(x) c(tim = x, Insol_l1l2(ber78(x), l1=pi, l2=2*pi, lat=65 * pi /180, S0=1361.371, avg=TRUE)))
S65sum <- sapply(time, function(x) c(tim = x, Insol_l1l2(ber78(x), l1=0, l2=pi, lat=-65 * pi /180, S0=1361.371, avg=TRUE)))
S65win <- sapply(time, function(x) c(tim = x, Insol_l1l2(ber78(x), l1=pi, l2=2*pi, lat=-65 * pi /180, S0=1361.371, avg=TRUE)))

## soslstices
N65sso <- sapply(time, function(x) c(tim = x, Insol(ber78(x), long=pi/2, lat=65 * pi /180, S0=1361.371)))
N65wso <- sapply(time, function(x) c(tim = x, Insol(ber78(x), long=3*pi/2, lat=65 * pi /180, S0=1361.371)))
S65sso <- sapply(time, function(x) c(tim = x, Insol(ber78(x), long=pi/2, lat=-65 * pi /180, S0=1361.371)))
S65wso <- sapply(time, function(x) c(tim = x, Insol(ber78(x), long=3*pi/2, lat=-65 * pi /180, S0=1361.371)))

## astrosummer/winter length
Ts <- sapply(time, function(x) c(tim = x, l2day(ber78(x), pi) - l2day(ber78(x), 0)))
Tw <- sapply(time, function(x) c(tim = x, 360 - l2day(ber78(x), pi) + l2day(ber78(x), 0)))

## gather dataframe
df <- data.frame(time=N65sum["tim", 0:-1],
                 N65sum=N65sum["ecc", 0:-1],
                 N65win=N65win["ecc", 0:-1],
                 S65sum=S65sum["ecc", 0:-1],
                 S65win=S65win["ecc", 0:-1],
                 N65sso=N65sso[2, 0:-1],
                 N65wso=N65wso[2, 0:-1],
                 S65sso=S65sso[2, 0:-1],
                 S65wso=S65wso[2, 0:-1],
                 Ts=Ts["varpi", 0:-1],
                 Tw=Tw["varpi", 0:-1])

write.csv(df, "./palinsol.csv", row.names=FALSE)
