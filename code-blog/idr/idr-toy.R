library(idr)
library(tidyverse)
library(ggplot2)

data("simu.idr")
x <- cbind(-simu.idr$x, -simu.idr$y)
x_r <- cbind(-simu.idr$x, sample(-simu.idr$y, nrow(x)))

expand.grid(seq(-1,10,3), seq(0.4,1,0.3), seq(0.1,1,0.3), seq(0.1,1,0.3))

comb = expand.grid(seq(-1,10,3), seq(0.4,1.5,0.3))

good = apply(comb, 1, function(c){
    #print(class(c[1]))
    idr_x = est.IDR(x, mu = c[1], sigma = c[2], rho = 0.8, p = 0.7)
    bind_rows(
        select.IDR(x, idr_x$idr, 0.01)[["x"]] %>% 
            as.data.frame() %>% 
            mutate(mu = c[1], sd = c[2], keep = 1),
        as.data.frame(x) %>% 
            mutate(mu = c[1], sd = c[2], keep = 0)
    )
}) %>% bind_rows()

ggplot(filter(good, keep == 0), aes(V1, V2)) +
    geom_point(color = "grey") +
    geom_point(data = filter(good, keep == 1), color = "orange") +
    facet_grid(mu ~ sd) +
    xlab("sim rep1") +
    ylab("sim rep2") +
    ggtitle("mu values versus sd values. rho = 0.8, proportion = 0.7") +
    ggsave("good_param_testing.png")

bad = apply(comb, 1, function(c){
    idr_x = est.IDR(x_r, mu = c[1], sigma = c[2], rho = 0.8, p = 0.7)
    bind_rows(
        select.IDR(x_r, idr_x$idr, 0.01)[["x"]] %>% 
            as.data.frame() %>% 
            mutate(mu = c[1], sd = c[2], keep = 1),
        as.data.frame(x) %>% 
            mutate(mu = c[1], sd = c[2], keep = 0)
    )
}) %>% bind_rows()

ggplot(filter(bad, keep == 0), aes(V1, V2)) +
    geom_point(color = "grey") +
    geom_point(data = filter(bad, keep == 1), color = "orange") +
    facet_grid(mu ~ sd) +
    ggsave("good_param_testing.png")

idr_xr = est.IDR(x_r, mu = 10, sigma = 1, rho = 0.4, p <- 0.1)
plot(x)
points(select.IDR(x, idr_x$idr, 0.01)[["x"]], co = "red")
plot(x_r)
points(select.IDR(x_r, idr_xr$idr, 0.01)[["x"]], co = "red")


summary(-log10(idr_x$idr))
summary(-log10(idr_xr$idr))
summary(-log10(idr_xr3$idr))
select.IDR(x, idr_x$IDR, 0.01)[["n"]]
select.IDR(x_r, idr_xr$IDR, 0.01)[["n"]]
select.IDR(x_r3, idr_xr3$IDR, 0.01)[["n"]]
