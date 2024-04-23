# tetst
setwd("/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/analysis/testsamples/DGE")
sayHello <- function(){
   print('hello')
}

sayHello()

extraction = rep(c("test", "tempus"), 5)
print(extraction)
write(extraction, file = "extraction2.txt")