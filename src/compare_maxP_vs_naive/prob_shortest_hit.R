prob_shortest_hit <- function(peptides, prob, b) {
#================================================================================
# Calculate P(shortest hit <= b)
	idx <- c()
	no <- 0
	if (is.vector(peptides)) {
		if (get_length(peptides)<=b){
			result <- prob
		} else {
			result <- 0
		}
	} else {
		for (n in 1:dim(peptides)[1]) {
			if (get_length(peptides[n,])<=b) {
				idx <- c(idx, n)
				no <- no +1
			}
		}
		if (no == 0) {
			result <- 0
		} else {
			prod <- 1
			for (i in 1:length(idx)) {
				prod <- prod * (1.0-prob[idx[i]])
			}
			result <- 1.0-prod
		}
	}
	return (result)
}

get_length <- function(peptide) {
	len <- 1
	for (i in 1:length(peptide)) {
		if (peptide[i] != -1){
			len <- len + 1
		}
	}
	return (len)
}

