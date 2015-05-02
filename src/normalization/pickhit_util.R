threshold <- function(org.data, treat.1, tr.1, treat.2, tr.2, treat.3=NA, tr.3=NA) {
    sub.data <- org.data[org.data[, treat.1] > tr.1 & org.data[, treat.2] < tr.2,]
    if (!is.na(treat.3))
        sub.data <- sub.data[sub.data[, treat.3] < tr.3,]
    print (nrow(sub.data))
    return (sub.data)
}

to.log <- function(v)
    return (log(0.1 * pmax((v - median(v)) / (max(v) - median(v)), 0) + 1e-8)) 

log.rank <- function(org.data, m1, coef1, m2, coef2, m3=NA, coef3=NA) {
    d1 <- to.log(org.data[, m1])
    d2 <- to.log(org.data[, m1] - org.data[, m2])
    measure <- coef1 * d1 + coef2 * d2
    if (!is.na(m3)) {
        d3 <- to.log(org.data[, m1] - org.data[, m3])
        measure <- measure + coef3 * d3
    }
    org.data <- cbind(org.data, measure)
    return (org.data[order(measure, decreasing=T),])
}

slice.table <- function(org.table, unit.length, dest.folder, filename) {
    org.table <- org.table[org.table[, 'TS'] != 2,]
    if (floor(nrow(org.table) / unit.length) > 0) {
        for (i in 0:(floor(nrow(org.table) / unit.length) - 1))
            write.csv(org.table[(i * unit.length + 1):((i + 1) * unit.length),], sprintf("%s/%s_%d_to_%d.csv", dest.folder, filename, i * unit.length + 1, (i+1) * unit.length), row.names=F)
        write.csv(org.table[(floor(nrow(org.table) / unit.length) * unit.length + 1):nrow(org.table),], sprintf("%s/%s_%d_to_%d.csv", dest.folder, filename, floor(nrow(org.table) / unit.length) * unit.length + 1, nrow(org.table)), row.names=F)
    }
    write.csv(org.table, sprintf("%s/%s.csv", dest.folder, filename), row.names=F)
}
