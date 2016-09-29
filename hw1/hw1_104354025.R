# read PAM1 from data
pam1 <-  read.table(file.choose(),h=T)

# check PAM1 data
dim(pam1)
str(pam1)

ppam1 <- as.matrix(pam1/colSums(pam1))

# construct PAM250 from PAM1
pamn <- function(x) {
  ppamn <- ppam1
    for (i in 1:(x-1)){
      ppamn <- ppamn %*% ppam1
    }
  ppamn
}

pam250 <- round(pamn(250)*100)

# output PAM250 as a file
write.table(pam250,"E:/git/nccubio/hw1/pam250.txt")