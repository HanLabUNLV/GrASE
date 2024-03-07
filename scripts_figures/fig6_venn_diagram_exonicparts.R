library(eulerr)
#DEXSeq exonic unit
fit <- euler(c("rMATS detected" = 138997-95503-166,
               "DEXSeq sig" = 3647-817-166-289,
               "DEXSeq sig&rMATS detected" = 166,
               "rMATS detected&rMATS tested" = 138997-95503-7432,
               "rMATS detected&rMATS tested&rMATS sig" = 7432-289,
               "rMATS detected&rMATS tested&DEXSeq sig" = 917-289,
               "rMATS detected&rMATS tested&rMATS sig&DEXSeq sig" = 289))

plot(fit,
     fills = list(fill = c("#77bca2", "#e1926b","#a09cc8", "#bae1ff"), alpha = 0.7),
     labels = list(col = "black", font = 4, cex=1.5),quantities = F)
