# sequence, solvent accessibility, and secondary structure element for 1311512 and 138284
SEQ.131612  <- "QSFSEEDSFKKCSSEVEAKNKIEELLASLLNRVCQDGRKPHTVRLIIRRYSSEKHYGRESRQCPIPSHVIQVMTPMVDILMKLFRNMTLLSVCFCNLK"
SA.131612   <- "AAAAAAAAAAAAAAAAAAAAAaAAAaAAaAAAAAAAAAAaAAaAaAAAAAAAAAAAAAAAAAAAaAAAAAAAAAAAAAAaaAAAAAaAAaAaAaAAaA"
SSE.131612  <- "CEEEEEECPCPCCCPPCCCHHHHHHHHHHHHHHHHHPCEEEEEEEEECCCCCCCCCCCEEEEEECCHHHHCCHHHHHHHHHHHHHCCCCEEEEEEPEE"

SEQ.138284  <- "VRKSIGRIVTMKRNSRNLEEIKPYLFRAIEESYYKLDKRIPKAIHVVAVTEDLDIVSRGRTFPHGISKETAYSESVKLLQKILEEDERKIRRIGVRFSKFI"
SA.138284   <- "AAAAAaAAAAaAAAAAAAAAaAAAaAAaaAAaAAAaAAAAaAaaAaAaaaAAAAAAAAaAAAAAAaAAAAaAAAaAAaaAAaAAAAAAAaAAaaaAaAAAA"
SSE.138284  <- "CCCEEEEEEEEEEEECCHHHHHHHHHHHHHHHHHHHPCCCEEEEEEEEEECCPCEEEEEEECCCCCCHHHHHHHHHHHHHHHHHHCCCCEEEEEEEEEPEC"

# structurally aligned sequences of 1311512 and 138284 using TM-align (length = 89, RMSD = 2.04, TM-score = 0.73731, ID = 0.157)
SAlignedSEQ.131612  <- "Q-----SFSEEDSFKKCS-S-EVEAKNKIEELLASLLNRVCQDGRKPHTVRLI-IRRYSSEKHYGRESRQCPIPSHVIQ----VMTPMVDILMKLFR---N-----MTLLSVCFCNLK-"
SAlignedSEQ.138284  <- "-VR-K-SIGRIVTMKRNSRNLE-EIKPYLFRAIEESYYKLDK--RIPKAIHVVAVTE--DLD---IVSRGRTFP-----HGIS-KETAYSESVKLLQKILEEDERKIRRIGVRFSKF-I"

RefAlign <- c(SAlignedSEQ.131612,SAlignedSEQ.138284);

# a function for global alignment
align <- function(seq1,seq2,score_matrix=NA,gap_open=-1,gap_ext=-1,match_score=1,mismatch_score=-1) {
    aas1 <- strsplit(seq1,"")[[1]];
    aas2 <- strsplit(seq2,"")[[1]];
    len1 <- length(aas1);
    len2 <- length(aas2);

    score <- mat.or.vec(len2+1,len1+1);
    point <- mat.or.vec(len2+1,len1+1);
    direc <- factor(1:4,labels=c("U", "L", "D", "N"));

    point[1,1] <- "N";
    point[1,2:(len1+1)] <- "L";
    point[2:(len2+1),1] <- "U";

    score[1, 2] = gap_open;
    score[1, 3:(len1+1)] = gap_ext*1:(len1-1);
    score[2, 1] = gap_open;
    score[3:(len2+1), 1] = gap_ext*1:(len2-1);

    for (i in 1:(len2)) {
        aa2 <- aas2[i];
        for (j in 1:(len1)) {
            aa1 <- aas1[j];

            diag_score  <- score[i, j] + ifelse(aa1 == aa2, match_score, mismatch_score);
            up_score    <- score[i, j + 1] + gap;
            left_score  <- score[i + 1, j] + gap;

            if (diag_score >= up_score) {
                if (diag_score >= left_score) {
                    score[i+1, j+1] = diag_score;
                    point[i+1, j+1] = "D";
                } else {
                    score[i+1, j+1] = left_score;
                    point[i+1, j+1] = "L";
                }
            } else {
                if (up_score > left_score) {
                    score[i+1, j+1] = up_score;
                    point[i+1, j+1] = "U";
                } else {
                    score[i+1, j+1] = left_score;
                    point[i+1, j+1] = "L";
                }
            }
        }
    }

    ali1    <- "";
    ali2    <- "";
    ii      <- dim(point)[1];
    jj      <- dim(point)[2];
    keepgo  <- TRUE;

    while(keepgo) {
        p = point[ii,jj]
        if (p == "N") { break; }
        s = score[ii,jj];
        if (p == "D") {
            ali1 <- paste(aas1[jj-1],ali1,sep="")
            ali2 <- paste(aas2[ii-1],ali2,sep="")
            jj <- jj - 1;
            ii <- ii - 1;
        } else if (p == "L") {
            ali1 <- paste(aas1[jj-1],ali1,sep="");
            ali2 <- paste(ali2,"-",sep="");
            jj <- jj - 1;
        } else if (p == "U") {
            ali1 <- paste(ali1,"-",sep="");
            ali2 <- paste(ali2,aas2[ii-1],sep="");
            ii <- ii - 1;
        } else {
            stop("Something wrong!\n");
        }
    }
    return(c(ali1, ali2));
}

# a function for PID calculation
#calcPID <- function(s1, s2) {
    #aas1    <- strsplit(s1, "")[[1]];
    #aas2    <- strsplit(s2, "")[[1]];
    #len     <- length(aas1);
    #ali     <- 0;
    #gap     <- 0;
    #idt     <- 0;
    #for (i in 1:len) {
        #if (aas1[i] != "-" & aas2[i] != "-") {
            #ali <- ali + 1;
            #if (aas1[i] == aas2[i]) {
                #idt <- idt + 1;
            #}
        #} else if ((aas1[i] == "-" &  aas2[i] != "-") | (aas1[i] != "-" &  aas2[i] == "-")) {
            #gap <- gap + 1;
        #}
    #}
    ##return(idt / (ali + gap));
    #return(idt / (ali));
#}

#calcPID(SAlignedSEQ.131612, SAlignedSEQ.138284)
