library(stringr)
library(dplyr)
library("ggplot2")

times = read.csv("times.txt")
times
avg_times = aggregate( Time ~ Nodes:Alg, data=times, FUN=mean)
avg_times
avg_times$sd = aggregate( Time ~ Nodes:Alg, data=times, FUN=sd)$Time
avg_times


avg_lagrange = avg_times %>%
select(Nodes, Alg, Time, sd) %>%
filter(str_detect(Alg, "Lagrange"))

avg_newton = avg_times %>%
select(Nodes, Alg, Time, sd) %>%
filter(str_detect(Alg, "Newton"))

avg_gsl = avg_times %>%
select(Nodes, Alg, Time, sd) %>%
filter(str_detect(Alg, "GSL"))

avg_akima = avg_times %>%
select(Nodes, Alg, Time, sd) %>%
filter(str_detect(Alg, "Akima"))

avg_spline = avg_times %>%
select(Nodes, Alg, Time, sd) %>%
filter(str_detect(Alg, "Spline"))

#wykres + SD - Lagrange
ggplot() + 
geom_errorbar(data=avg_lagrange, mapping=aes(x = Nodes, ymin=Time-sd, ymax=Time+sd), color="black") +
geom_line(data = avg_lagrange, aes( y = Time, x = Nodes), size = 1, color = "red" ) +
geom_point(data=avg_lagrange, mapping=aes(x = Nodes, y = Time), color="black") + 

labs(title = "Lagrange interpolation efficiency\n", x = "Number of nodes", y = "Execution time [s]", color = NULL) +
theme_bw() +
geom_point(size=5)

#wykres + SD - Newton
ggplot() + 
geom_errorbar(data=avg_newton, mapping=aes(x = Nodes, ymin=Time-sd, ymax=Time+sd), color="black") +
geom_line(data = avg_newton, aes( y = Time, x = Nodes), size = 1, color ="red" ) +
geom_point(data=avg_newton, mapping=aes(x = Nodes, y = Time), color="black") + 

labs(title = "Newton interpolation efficiency\n", x = "Number of nodes", y = "Execution time [s]", color = NULL) +
theme_bw() +
geom_point(size=5)

#wykres + SD - Spline
ggplot() + 
geom_errorbar(data=avg_spline, mapping=aes(x = Nodes, ymin=Time-sd, ymax=Time+sd), color="black") +
geom_line(data = avg_spline, aes( y = Time, x = Nodes), size = 1, color = "red" ) +
geom_point(data=avg_spline, mapping=aes(x = Nodes, y = Time), color="black") + 

labs(title = "Spline interpolation efficiency\n", x = "Number of nodes", y = "Execution time [s]", color = NULL) +
theme_bw() +
geom_point(size=5)

#wykres + SD - Akima
ggplot() + 
geom_errorbar(data=avg_akima, mapping=aes(x = Nodes, ymin=Time-sd, ymax=Time+sd), color="black") +
geom_line(data = avg_akima, aes( y = Time, x = Nodes), size = 1, color = "red" ) +
geom_point(data=avg_akima, mapping=aes(x = Nodes, y = Time), color="black") + 

labs(title = "Akima interpolation efficiency\n", x = "Number of nodes", y = "Execution time [s]", color = NULL) +
theme_bw() +
geom_point(size=5)

#wykres + SD - GSL
ggplot() + 
geom_errorbar(data=avg_gsl, mapping=aes(x = Nodes, ymin=Time-sd, ymax=Time+sd), color="black") +
geom_line(data = avg_gsl, aes( y = Time, x = Nodes), size = 1, color = "red" ) +
geom_point(data=avg_gsl, mapping=aes(x = Nodes, y = Time), color="black") + 

labs(title = "GSL interpolation efficiency\n", x = "Number of nodes", y = "Execution time [s]", color = NULL) +
theme_bw() +
geom_point(size=5)

#czasy na 1 wykresie - bez SD
ggplot() + 
geom_line(data = avg_lagrange, aes( y = Time, x = Nodes, color="blue"), size = 1 ) +
geom_line(data = avg_newton, aes( y = Time, x = Nodes, color="red"), size = 1 ) +
geom_line(data = avg_spline, aes( y = Time, x = Nodes, color="green"), size = 1 ) +
geom_line(data = avg_akima, aes( y = Time, x = Nodes, color="black"), size = 1 ) +
geom_line(data = avg_gsl, aes( y = Time, x = Nodes, color="grey"), size = 1 ) +

scale_color_discrete(name = "Algorithms", labels = c("Akima", "Lagrange", "Spline", "GSL", "Newton")) +
labs(title = "Interpolation efficiency\n", x = "Number of nodes", y = "Execution time [s]", color = "Legend\n") +
theme_bw() +      #obramowka

geom_point(data=avg_lagrange, mapping=aes(x = Nodes, y = Time), color="black") + 
geom_point(data=avg_newton, mapping=aes(x = Nodes, y = Time), color="black") + 
geom_point(data=avg_spline, mapping=aes(x = Nodes, y = Time), color="black") + 
geom_point(data=avg_akima, mapping=aes(x = Nodes, y = Time), color="black") + 
geom_point(data=avg_gsl, mapping=aes(x = Nodes, y = Time), color="black") + 
geom_point(size=5) 