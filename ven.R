library(VennDiagram)
setwd('Raf2236Venn.csv')

csv = read.csv("venn2.csv", header = T, na.strings = c(''))

ABA_up= csv$ABA_up
ABA_down = csv$ABA_down
DKO_up = csv$DKO_up
DKO_down = csv$DKO_down

ABA_up= ABA_up[!is.na(ABA_up)]    
ABA_down = ABA_down[!is.na(ABA_down)]
DKO_up = DKO_up[!is.na(DKO_up)]
DKO_down = DKO_down[!is.na(DKO_down)]

data = list(ABA_Up = ABA_up, DKO_Up = DKO_up, ABA_Down = ABA_down, DKO_Down = DKO_down)
venn.diagram(data,filename = "Venn2.png")



