###########################################
#Choose file "test.txt"
Seqs <- scan(file.choose(),"")

laboA <- unlist(strsplit(Seqs[2],""))
lycsB <- unlist(strsplit(Seqs[4],""))

#Choose file "PAM100.txt"
PAM100 <- read.table(file.choose(),T)
#Gap represents "-"
colnames(PAM100)[24] <- "-"
rownames(PAM100)[24] <- "-"

#Choose file "PAM250.txt"
PAM250 <- read.table(file.choose(),T)
colnames(PAM250)[24] <- "-"
rownames(PAM250)[24] <- "-"

###############################################
#Sequence alignment function
#Input variable
#sequence A is vector 
#sequence B is vector
#PAM is 100 or 250 ,they represent PAM100 amd PAM250
#Type is alignment type is "Global" or "Local" 
#OpenExtension is a vector which is (OE,GG,GP)
#OE is TRUE or FALSE, OE = TRUE is that scores of open-gap and extension-gap are different.
#GG is the score of gap match gap (But this will never occur in the case), GP is the score of gap match protein.

SeqAln <- function(laboA, lycsB, PAMs=250, Type="Global",
          OpenExtension = TRUE, OG = -10, EG = -2){

#Choose score scheme
  if(PAMs == 100){
    PAM <- PAM100
  }else{PAM <- PAM250}
  Index <- colnames(PAM)
  dPAM <- dim(PAM)
#OpenGap matrix
  PAM[,dPAM[2]] <- rep(OG,dPAM[1])
  PAM[dPAM[1],] <- rep(OG,dPAM[2])
#ExtensionGap matrix
  PAME <- PAM
#Let extension gap match non-gap default is -2 
#gap matcg gap default is 0
  if(OpenExtension == TRUE){
    dPAM <- dim(PAM)
    PAME[,dPAM[2]] <- rep(EG,dPAM[1])
    PAME[dPAM[1],] <- rep(EG,dPAM[2])
  }


#Initial local score matrix
  ILPS <- matrix(numeric(),ncol = length(laboA)+1,nrow = length(lycsB)+1)
   colnames(ILPS) <- c("-",laboA)
    rownames(ILPS) <- c("-",lycsB)
  r <- dim(ILPS)
    for(i in 1:r[2]){
      for(j in 1:r[1]){
        ILPS[j,i] <- PAM[which(rownames(ILPS)[j] == Index),which(colnames(ILPS)[i] == Index)]
      }
    }
#Decide alignment type   
  if(Type == "Global"){
    From <- c(2,2)
    Bound <- -Inf
  }else{
    From <- c(r[1],r[2])
    Bound <- -1
  }

  LocalScore <-matrix(numeric(),ncol = length(laboA)+1,nrow = length(lycsB)+1)
   colnames(LocalScore) <- c("-",laboA)
    rownames(LocalScore) <- c("-",lycsB)



#Sequence alignment
  maxSP <- -Inf
  LocalAlnPath <- c()
  for(i in 1:From[2]){
    for(j in 1:From[1]){
      if(i == 1 & j == 1){
      }else{
        x0 <- x <- i
        y0 <- y <- j
        aln <- matrix(c(colnames(ILPS)[x0],rownames(ILPS)[y0]),2,1)
        if( sum(aln[,dim(aln)[2]] == "-") == 0 ){
          ScoreMatrix <- PAM
        }else{
          ScoreMatrix <- PAME
        }
        score <- ILPS[y0,x0]
        while(is.na(LocalScore[y0,x0])){

#To decide the moving direction,W=1 is right , W=2 is right and down, W=3 is down.
          if(y < r[1] & x < r[2]){
            D <- c(ScoreMatrix[24,which(Index == colnames(ILPS)[x+1])],ScoreMatrix[which(Index == rownames(ILPS)[y+1]),which(Index == colnames(ILPS)[x+1])],ScoreMatrix[which(Index == rownames(ILPS)[y+1]),24])
          }else{
            if(y == r[1] & x < r[2]){
              D <- c(ScoreMatrix[24,which(Index == colnames(ILPS)[x+1])],NA,NA)
            }else{
              if(y < r[1] & x == r[2]){
                D <- c(NA,NA,ScoreMatrix[which(Index == rownames(ILPS)[y+1]),24])
              }else{D <- c(NA,NA,NA)
              }              
            }
          }
          W <- which( D == max(D,na.rm=T))

#Different directions are the same scores. Randomly choice direction
          if(length(W) == 0 ){
            LocalScore[y0,x0] <- sum(score)
            LocalAlnPath <- aln
          }else{
            if(length(W) != 1){W <- sample(W,1)}
#Move to next point and record the match sequence
            if(D[W] > Bound){
              score <- c(score,D[W])
              if(W == 1){
                y <- y ; x <- x+1 ; aln <- cbind(aln,c(colnames(ILPS)[x],"-"))
              }else{
                if(W == 2){
                  y <- y+1 ; x <- x+1 ; aln <- cbind(aln,c(colnames(ILPS)[x],rownames(ILPS)[y]))
                }else{
                  y <- y+1 ; x <- x ; aln <- cbind(aln,c("-",rownames(ILPS)[y]))
                }
              }
            }else{
            LocalScore[y0,x0] <- sum(score)
            LocalAlnPath <- aln
            }    
          }
        }
#Output 
      if(Type == "Global"){
        if(maxSP[1] < LocalScore[y0,x0]){
          maxSP <- LocalScore[y0,x0]
           maxAlnPath <- list(LocalAlnPath)
        }else{
          if(maxSP[1] == LocalScore[y0,x0]){
           maxSP <- c(maxSP,LocalScore[y0,x0])
            maxAlnPath <-list(maxAlnPath ,LocalAlnPath)
          }
        }
      }else{
        if(maxSP[1] < LocalScore[y0,x0]){
          maxSP <- LocalScore[y0,x0]
           maxAlnPath <- list(LocalAlnPath)
        }else{
          if(maxSP[1] == LocalScore[y0,x0]){
            maxSP <- c(maxSP,LocalScore[y0,x0])
             maxAlnPath <-list(maxAlnPath ,LocalAlnPath)
          }
        }
      }
    }
  }
}
#Output sequence alignment list with Sequence_alignment, SP-score ,alignment type
SA <- list("Sequence_alignment" = maxAlnPath,"SP-score" = maxSP ,"Type" = Type,
"OpenAndExtend"=OpenExtension, "CostOpenGap"=PAM[dPAM[1],1],"CostExtentGap"=PAME[dPAM[1],1])

#Combine the sequence order
for(L in 1:length(SA[[2]])){
  ss <- matrix(unlist(SA[[1]][L]),nrow=2)
  sss <- c()
  for(k in 1:dim(ss)[2]){
    sss <- paste(sss,ss[,k],sep="")
  }
  SA$Sequence_alignment[[L]]<- matrix(sss,2,1)
}

print(SA)

}


##################################################
######################Exsample####################
#########PAM250
#Local Opengap and Extensiongap Alignment
saLOP_250 <- SeqAln(laboA,lycsB,250,"Local",TRUE,-10,-2)
#Global Opengap and Extensiongap Alignment
saGOP_250 <- SeqAln(laboA,lycsB,250,"Global",TRUE,-10,-2)
#Local Opengap Alignment
saLO_250 <- SeqAln(laboA,lycsB,250,"Local",FALSE,-10,-2)
#Global Opengap Alignment
saGO_250 <- SeqAln(laboA,lycsB,250,"Global",FALSE,-10,-2)


########PAM100
#Local Opengap and Extensiongap Alignment
saLOP_100 <- SeqAln(laboA,lycsB,100,"Local",TRUE,-10,-2)
#Global Opengap and Extensiongap Alignment
saGOP_100 <- SeqAln(laboA,lycsB,100,"Global",TRUE,-10,-2)
#Local Opengap Alignment
saLO_100 <- SeqAln(laboA,lycsB,100,"Local",FALSE,-10,-2)
#Global Opengap Alignment
saGO_100 <- SeqAln(laboA,lycsB,100,"Global",FALSE,-10,-2)


#Output FASTA file
install.packages("seqinr")
library("seqinr")  
SeqName<- gsub(">","",Seqs)[c(1,3)]
#Choose the sequence alignment
#Example choose "saLOP_250" 
write.fasta(as.list(unlist(saLOP_250[[1]])),
rep(SeqName,times=length(saLOP_250[[2]])),"result.fasta")


