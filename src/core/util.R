library(RMySQL)

# constants
USER <- "jialeiwang"
PASSWORD <- "wangjialei123"
HOST <- "work.cxcjqzn7ydtp.us-east-1.rds.amazonaws.com"

QuerySql <- function(dbname, query) {
  # Submit query to db, and return retrived data
  #
  # Args:
  #   dbname: string that specifies name of db
  #   query: string
  # 
  # Returns:
  #   data.frame
  con <- dbConnect(MySQL(), user=USER, password=PASSWORD, host=HOST, dbname=dbname)
  data <- dbGetQuery(conn=con, statement=query)
  dbDisconnect(con)
  return (data)
}

get_con <- function(dbname)
  con <- dbConnect(MySQL(), user=USER, password=PASSWORD, host=HOST, dbname=dbname)
