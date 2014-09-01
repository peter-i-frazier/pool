## Dependencies: NB_utility.R

# Generate feature vectors of all peptides, and verify their nterm and cterm add up to the peptide sequence. 

# :param data.org:
#		A data.frame. Each row corresponds to a peptide. The matrix has n+3 columns. The first n columns contain labels. The last # 3 columns nterm, cterm, sequence
#		E.g. we have 2 enzymes, and the first peptide 'AGDVSPTLG' works with Enzyme_1('sfp') but not with Enzyme_2('acph'), then the 1st row of data.org could be:
#			sfp	PfAcpH   nterm	cterm   sequence
#		1	1	0	     'AGDV'	'PTLG'  'AGDVSPTLG'
# :param classlist:
#  A vector specifying the class each amino-acids belongs to. 
#		A	R	N	D	C	Q	...
#		5	4	2	1	8	2	...
# :param nL:
#  Maximum number of amino acids at the left of the serene in the feature vector.
# :param nR:
#  Maximum number of amino acids at the right of the serene in the feature vector.
# :return feature:
#  A matrix. Each row corresponds to the feature vector of a peptide. The matrix has nL+nR+n columns, where nL is the number of 
#features to the left of the serene and nR is the number of features to the right of the serene. The first (nL + nR) columns 
#are the feature vector of a peptide. If a peptide has more amino acids to the left of the serene than nL, only the nL amino 
#acids that are closest to the serene will be transformed to features according to the class mapping list. 
#  If a peptide has less than nL amino-acids to the left of the serene then the leftmost positions of its feature vector will be 
#filled with '-1'. The same goes with the right side. The last n columns indicates whether a peptide works with the n enzymes 
#respectively. E.g. we have 2 enzymes(n=2), nL = 5, nR =5, the first peptide is 'AGDVSPTLG', and it works with Enzyme_1('sfp') 
#but not Enzyme_2('acph'), then the 1st row of
#       	feature could be: 
#			1	2	3	...	9	10	sfp	PfAcpH
#		1	-1	5	2	...	1	-1	1	0
#		'-1' appears in position 1 and 10 because there're only 4 amino acids at the left or right of the serene 'S' respectively. If whether a peptide works with an 
#		enzyme is not specified, the corresponding column will be filled with 'NA'.

getFeatures <- function(data.org, classlist, nL, nR)
{
    error <- 0
    #nVal: number of values each feature can take
    nVal <- length(unique(as.numeric(classlist)))
    #nOUTCOME: number of outcome values
    nOUTCOME <- dim(data.org)[2]-3
    if (nOUTCOME <= 0) {
        print ("no labels")
        error <- error + 1
    }
    feature <- c()
    for (r in 1:dim(data.org)[1]) {
        if (data.org[r,'nterm'] == '') {
            next
        }
        # verify nterm and cterm add up to the whole sequence
        if (paste(data.org[r,'nterm'],'S',data.org[r,'cterm'],sep='') != data.org[r,'sequence']) {
            print ("add up verification failed")
            print (r)
            error <- error + 1
        }
        one_feature <- rep(-1, nL+nR+nOUTCOME)
        sequence <- unlist(strsplit(data.org[r, 'nterm'],''))
        l.seq  <- length(sequence)
        l <- min(nL,l.seq)
        for (i in 1:l) {
            one_feature[nL+1-i] <- classlist[1, sequence[l.seq+1-i]]
	    }
        sequence <- unlist(strsplit(data.org[r, 'cterm'],''))
        l.seq  <- length(sequence)
        c <- min(nR,l.seq)
        for (i in 1:c) {
            one_feature[nL+i] <- classlist[1, sequence[i]]
	    }
        for (i in 1:nOUTCOME) {
            if (is.na(data.org[r,i])) {
                one_feature[nL+nR+i] <- -1
            } else {
                one_feature[nL+nR+i] <- data.org[r,i]
            }
        }
        feature <- rbind(feature, one_feature)
    }    
    #outcome.names : name of outcome values
    outcome.names <- colnames(data.org)[1:(nOUTCOME)]
    feature <- as.data.frame(feature)
    colnames(feature) <- c(paste('L',nL:1,sep=""), paste('R',1:nR,sep=""), outcome.names)
    if (error > 0) {
        print ("there is error when getting features!")
    }
    return( feature )
}

# Calculate area under ROC curve given tunable parameters: gamma_0, gamma_1, prior
# 
# :param X: each row is feature vector for a peptide in the training dataset
# :type X: matrix
# :param Y: label of peptides in the training dataset
# :type Y: vector
# :param AAclass: table for reduced AA alphabet
# :type AAclass: vector
# :param S.Pos: position of Serine
# :type S.Pos: int
# :param nL: number of AAs to the left of Serine, usually equals to S.Pos
# :type nL: int
# :param nR: number of AAs to the right of Serine
# :type nR: int
# :param gamma_0: scaling model parameter
# :type gamma_0: float
# :param gamma_1: scaling model parameter
# :type gamma_1: float
# :param prior: prior probability of positive label
# :type prior: float
# :param itr: iteration for MC version of estimating probability of hit
# :type itr: int
# :return auc: area under ROC curve
# :rtyupe: float

loocv <- function(X, Y, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr) {
    num_rows <- dim(X)[1]
    prob <- rep(0, num_rows)
    for (n in 1:num_rows) {
        print (n)
        train_X <- X[-n,]
        train_Y <- Y[-n]
        prob[n] <- Naive_Bayes(train_X, train_Y, X[n,], AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr)
    }
    xy <- ROC_xy_output(prob, Y)
    auc <- AUC(rev(xy$x), rev(xy$y))
    return (auc)
}

fold_cv <- function(X, Y, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr, num_fold) {
    num_rows <- dim(X)[1]
    group_size <- floor(num_rows / (num_fold - 1))
    prob <- c()
    for (n in 1:(num_fold-1)) {
        print (n)
        train_X <- X[-c(((n-1)*group_size+1):(n*group_size)),]
        train_Y <- Y[-c(((n-1)*group_size+1):(n*group_size))]
        test_X <- X[c(((n-1)*group_size+1):(n*group_size)),]
        part_prob <- Naive_Bayes(train_X, train_Y, test_X, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr)
        prob <- c(prob, part_prob)
    }
    train_X <- X[-c(((num_fold-1)*group_size+1):num_rows),]
    train_Y <- Y[-c(((num_fold-1)*group_size+1):num_rows)]
    test_X <- X[c(((num_fold-1)*group_size+1):num_rows),]
    part_prob <- Naive_Bayes(train_X, train_Y, test_X, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr)
    prob <- c(prob, part_prob)

    xy <- ROC_xy_output(prob, Y)
    auc <- AUC(rev(xy$x), rev(xy$y))
    return (auc)
}

cv <- function(X, Y, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr) {
    num_rows <- dim(X)[1]
    prob <- Naive_Bayes(X, Y, X, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr)
    xy <- ROC_xy_output(prob, Y)
    auc <- AUC(rev(xy$x), rev(xy$y))
    return (auc)
}

generate_recommendation_MAP_old <- function(X, Y, AAclass, S.Pos, max_L, max_R, min_L, min_R, gamma_0, gamma_1, add_ins, num_recom) {
	num_features <- dim(X)[2]
	rec <- c()	
	train_x <- X
	train_y <- Y
	num_peptides <- 0
    count_repeated_recom <- 0

	while (num_peptides < num_recom) {		
        ratio <- get_ratio_old_method(train_x, train_y, AAclass, gamma_0, gamma_1)
		best_class <- as.numeric(apply(ratio, 2, which.max))
		length_left <- ceiling(runif(1, min = min_L-1, max = max_L))
		length_right  <- ceiling(runif(1, min = min_R-1, max = max_R))
		best_peptide <- rep(-1, num_features)
		best_peptide[(S.Pos-length_left+1):(S.Pos+length_right)] <- best_class[(S.Pos-length_left+1):(S.Pos+length_right)]
		if(isNew(train_x, best_peptide)){
			rec <- rbind(rec, best_peptide)
			num_peptides = num_peptides + 1
            print ("nth peptide")
            print (num_peptides)
		} else {
            count_repeated_recom <- count_repeated_recom + 1
            print ("repeated count")
            print (count_repeated_recom)
        }
        for (i in add_ins) {
            train_x <- rbind(train_x, best_peptide)
            train_y <- c(train_y, 0)
        }
	}
	colnames(rec) <- colnames(X)
	rownames(rec) <- c(1:dim(rec)[1])
	return (rec)
}

generate_recommendation_MAP_new <- function(matrix_x_label_prefer, vector_y_label_prefer, matrix_x_label_unprefer, vector_y_label_unprefer, matrix_x_unlabel, vector_y_unlabel, classlist, serine_position, num_recom, max_left, max_right, min_left, min_right, gamma_0_label_prefer, gamma_1_label_prefer, gamma_0_label_unprefer, gamma_1_label_unprefer, gamma_0_unlabel, gamma_1_unlabel, add_ins) {
	num_features <- dim(matrix_x_label_prefer)[2]
	rec <- c()	
	train_x_label_prefer <- matrix_x_label_prefer 
	train_y_label_prefer <- vector_y_label_prefer
	train_x_label_unprefer <- matrix_x_label_unprefer 
	train_y_label_unprefer <- vector_y_label_unprefer
	train_x_unlabel <- matrix_x_unlabel
	train_y_unlabel <- vector_y_unlabel
	num_peptides <- 0
    count_repeated_recom <- 0

	while (num_peptides < num_recom) {		
        ratio <- GetRatioNewMethod(classlist, train_x_label_prefer, train_y_label_prefer, gamma_0_label_prefer, gamma_1_label_prefer, train_x_label_unprefer, train_y_label_unprefer, gamma_0_label_unprefer, gamma_1_label_unprefer, train_x_unlabel, train_y_unlabel, gamma_0_unlabel, gamma_1_unlabel)
		best_class <- as.numeric(apply(ratio, 2, which.max))
		length_left <- ceiling(runif(1, min = min_left-1, max = max_left))
		length_right  <- ceiling(runif(1, min = min_right-1, max = max_right))
		best_peptide <- rep(-1, num_features)
		best_peptide[(serine_position-length_left+1):(serine_position+length_right)] <- best_class[(serine_position-length_left+1):(serine_position+length_right)]
		if(isNew(rbind(train_x_label_prefer, train_x_label_unprefer, train_x_unlabel), best_peptide)){
			rec <- rbind(rec, best_peptide)
			num_peptides = num_peptides + 1
            print ("nth peptide")
            print (num_peptides)
		} else {
            count_repeated_recom <- count_repeated_recom + 1
            print ("repeated count")
            print (count_repeated_recom)
        }
        for (i in add_ins) {
            train_x_label_prefer <- rbind(train_x_label_prefer, best_peptide)
            train_y_label_prefer <- c(train_y_label_prefer, 0)
            train_x_label_unprefer <- rbind(train_x_label_unprefer, best_peptide)
            train_y_label_unprefer <- c(train_y_label_unprefer, 1)
            train_x_unlabel <- rbind(train_x_unlabel, best_peptide)
            train_y_unlabel <- c(train_y_unlabel, 0)
        }
	}
	colnames(rec) <- colnames(matrix_x_label_prefer)
	rownames(rec) <- c(1:dim(rec)[1])
	return (rec)
}
