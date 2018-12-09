library(stringr)
library(dplyr)
library("ggplot2")

result = read.csv("result.txt")
result

result_lagrange = result %>%
select(Type, X, Y) %>%
filter(str_detect(Type, "Lagrange"))

result_newton = result %>%
select(Type, X, Y) %>%
filter(str_detect(Type, "Newton"))

result_gsl = result %>%
select(Type, X, Y) %>%
filter(str_detect(Type, "GSL"))

result_akima = result %>%
select(Type, X, Y) %>%
filter(str_detect(Type, "Akima"))

result_spline = result %>%
select(Type, X, Y) %>%
filter(str_detect(Type, "Spline"))

nodes = result %>%
select(Type, X, Y) %>%
filter(str_detect(Type, "Generated_nodes"))

#Lagrange
ggplot() + 
geom_line(data = result_lagrange, aes( y = Y, x = X), size = 1, color = "red" ) +
geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 

labs(   title = "Lagrange interpolation polynomial\n", 
        subtitle = "Black dots are generated nodes\n", 
        x = "x", 
        y = "y", 
        color = NULL) +
theme_bw() +
geom_point(size=5)

#Newton
ggplot() + 
geom_line(data = result_newton, aes( y = Y, x = X), size = 1, color = "red" ) +
geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 

labs(   title = "Newton interpolation polynomial\n", 
        subtitle = "Black dots are generated nodes\n", 
        x = "x", 
        y = "y", 
        color = NULL) +
theme_bw() +
geom_point(size=5)

#Spline
ggplot() + 
geom_line(data = result_spline, aes( y = Y, x = X), size = 1, color = "red" ) +
geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 

labs(   title = "Spline interpolation polynomial\n", 
        subtitle = "Black dots are generated nodes\n", 
        x = "x", 
        y = "y", 
        color = NULL) +
theme_bw() +
geom_point(size=5)

#Akima
ggplot() + 
geom_line(data = result_akima, aes( y = Y, x = X), size = 1, color = "red" ) +
geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 

labs(   title = "Akima interpolation polynomial\n", 
        subtitle = "Black dots are generated nodes\n", 
        x = "x", 
        y = "y", 
        color = NULL) +
theme_bw() +
geom_point(size=5)

#GSL
ggplot() + 
geom_line(data = result_gsl, aes( y = Y, x = X), size = 1, color = "red" ) +
geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 

labs(   title = "GSL interpolation polynomial\n", 
        subtitle = "Black dots are generated nodes\n", 
        x = "x", 
        y = "y", 
        color = NULL) +
theme_bw() +
geom_point(size=5)

#Lagrange + Newton + GSL na jednym wykresie
ggplot() + 

geom_line(data = result_gsl, aes( y = Y, x = X, color = "red"), size = 1 ) +
geom_line(data = result_newton, aes( y = Y, x = X, color = "green"), size = 1 ) +
geom_line(data = result_lagrange, aes( y = Y, x = X, color = "blue"), size = 1 ) +

scale_color_discrete(name = "Algorithms", labels = c("Lagrange", "GSL", "Newton")) +
labs(   title = "Algorithm interpolation result comparison\n", 
        subtitle = "Black dots are generated nodes\n", 
        x = "x", 
        y = "y", 
        color = NULL) +
theme_bw() +

geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 
geom_point(size=5) +

coord_cartesian(xlim = c(40,50), ylim = c(-100, 100) ) #ograniczenie dziedziny i przeciwdziedziny

#Nodes scatter
ggplot() +
geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 
geom_point(size=5)

#Spline vs. Akima vs. Lagrange
ggplot() + 

geom_line(data = result_spline, aes( y = Y, x = X, color = "red"), size = 1 ) +
geom_line(data = result_akima, aes( y = Y, x = X, color = "green"), size = 1 ) +
geom_line(data = result_lagrange, aes( y = Y, x = X, color = "blue"), size = 1 ) +

scale_color_discrete(name = "Algorithms", labels = c("Lagrange", "Akima", "Spline")) +
labs(   title = "Spline vs. Akima vs. Lagrange\n", 
        subtitle = "Black dots are generated nodes\n", 
        x = "x", 
        y = "y", 
        color = NULL) +
theme_bw() +

geom_point(data = nodes, mapping=aes(x = X, y = Y), color="black") + 
geom_point(size=5) +

coord_cartesian(xlim = c(15,80), ylim = c(-1200, 1200) )