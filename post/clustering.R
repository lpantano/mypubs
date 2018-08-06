library(dplyr)
library(ggplot2)

x <- 1:10

# going UP
y <- x - x^2 + x^3 + 10
plot(x,y) 
noise <- rnorm(length(x), mean=5, sd=20)
lm <- lm(y + noise~ poly(x,3))
sim1 <- simulate(lm, 5)
df = sim1 %>% as.data.frame() %>% mutate(point=1:10) %>%
    reshape::melt.data.frame(id.vars = "point")
ggplot(df, aes(x=point, y=value, group=variable)) + geom_point()

# going DOWN
y <- 200 - x - x^2 
plot(x,y) 
noise <- rnorm(length(x), mean=5, sd=20)
lm <- lm(y + noise~ poly(x,3))
sim2 <- simulate(lm, 5)
df = sim2 %>% as.data.frame() %>% mutate(point=1:10) %>%
    reshape::melt.data.frame(id.vars = "point")
ggplot(df, aes(x=point, y=value, group=variable)) + geom_point()


# merge matrix
exp = rbind(t(sim1),t(sim2))
colnames(exp) <- paste0("s", colnames(exp))
row.names(exp) <- paste0("r", 1:10)

df = t(exp) %>% as.data.frame() %>% mutate(point=1:10) %>%
    reshape::melt.data.frame(id.vars = "point")
ggplot(df, aes(x=point, y=value, group=variable)) + geom_point()

# sim colData
de = data.frame(row.names=colnames(exp), 
                group = paste0("g", c(1,1,2,2,3,3,4,4,5,5)) )

res = degPatterns(exp, de, minc=1, time="group", col=NULL)

# bysteps

counts_group = t(sapply(rownames(exp), function(g){
    sapply(levels(de[,"group"]), function(i){
        idx = which(de[,"group"] == i)
        mean(exp[g, idx], na.rm=TRUE)
    })
}))
colnames(counts_group) = unique(de[,"group"])

library(cluster)
m = (1-cor(t(counts_group), method = "kendall"))
d = as.dist(m^2)
c = diana(d, diss = TRUE, stand = FALSE)

plot(as.dendrogram(c))
c$dc
table(cutree(as.hclust(c), h = c$dc))

# real data
library(rio)
p <- "~/repos/mypubs/post"
keep = import(file.path(p,"clusters_genes.tsv"))
keep_genes = keep[keep$cluster %in% c(1,2,3), "genes"]

ma = read.csv(file.path(p,"rlog_counts.tsv"),row.names=1)
ma_keep = as.matrix(ma[keep_genes, 1:15])

ma_de = data.frame(row.names=colnames(ma_keep), 
                   group=factor(gsub("_[0-9]$", "", colnames(ma_keep))))

res = degPatterns(ma_keep, ma_de, minc=1, time="group", col=NULL,
                  reduce = TRUE)


counts_group = t(sapply(rownames(ma_keep), function(g){
    sapply(levels(ma_de[,"group"]), function(i){
        idx = which(ma_de[,"group"] == i)
        mean(ma_keep[g, idx], na.rm=TRUE)
    })
}))
colnames(counts_group) = unique(ma_de[,"group"])

m = (1-cor(t(counts_group), method = "kendall"))
d = as.dist(m^2)
c = diana(d, diss = TRUE, stand = FALSE)

plot(as.dendrogram(c),leaflab = "none")
c$dc
group <- cutree(as.hclust(c), h = c$dc)
table(group)

# reduce function
ngroup <- unique(group)[1:4]
cor <- lapply(ngroup, function(nc1){
    sapply(ngroup, function(nc2){
        g1 = colMeans(counts_group[names(group[group==nc1]),])
        g2 = colMeans(counts_group[names(group[group==nc2]),])
        (1-cor.test(g1, g2)$estimate)^2
    })
})
cor <- do.call(rbind, cor)
colnames(cor) <- ngroup
rownames(cor) <- ngroup
h <- hclust(as.dist(cor), method = "ward.D2")
plot(h)
c <- cutree(h, h = (1-0.7)^2)
new <- c[as.character(group)]
names(new) <- names(group)
new
