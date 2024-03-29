testthat::context("Testing simpers")

## Folder path
extdata <- system.file("extdata", package = "Simpers")

## Testing elementary mode w/ files
x <- file.path(extdata, "Elementary_Mode", "iDC.txt")
sm <- file.path(extdata, "Elementary_Mode", "operations.txt")
dgm <- Simpers::simpers(x = x, sm = sm, elementary_mode = TRUE)

## Testing elementary mode w/ simplextree + input vectors
st <- simplextree::simplex_tree()
st$insert(as.list(as.integer(readLines(x)[-1])))
sm_vec <- readLines(sm)
dgm2 <- Simpers::simpers(x = st, sm = sm_vec, elementary_mode = TRUE)

## Test they are equal
testthat::expect_equal(dgm, dgm2)


x <- simplextree::simplex_tree()
x$insert(as.list(seq(0, 24)))
sm <- c("i 23 24", "# 1.000000", "i 0 1", "# 2.000000", "i 22 24", "# 3.000000", "i 0 2", "# 4.000000", "i 102",
  "# 5.000000", "i 2 102", "# 6.000000", "i 103", "# 7.000000", "i 103 22", "# 8.000000", "i 1 17",
  "# 9.000000", "i 8 23", "# 10.000000", "i 104", "# 11.000000", "i 104 23", "# 12.000000",
  "i 105", "# 13.000000", "i 1 105", "# 14.000000", "i 7 22", "# 15.000000", "i 2 18",
  "# 16.000000", "i 106", "# 17.000000", "i 2 106", "# 18.000000", "i 107", "# 19.000000", "i 107 22",
  "# 20.000000", "i 1 16", "# 21.000000", "i 9 23", "# 22.000000", "i 108", "# 23.000000", "i 108 23",
  "# 24.000000", "i 109", "# 25.000000", "i 1 109", "# 26.000000", "i 6 22", "# 27.000000", "i 2 19",
  "# 28.000000", "i 110", "# 29.000000", "i 2 110", "# 30.000000", "i 111", "# 31.000000", "i 111 22",
  "# 32.000000", "i 1 15", "# 33.000000", "i 10 23", "# 34.000000", "i 112", "# 35.000000", "i 112 23",
  "# 36.000000", "i 113", "# 37.000000", "i 1 113", "# 38.000000", "i 5 22", "# 39.000000", "i 2 20",
  "# 40.000000", "i 114", "# 41.000000", "i 2 114", "# 42.000000", "i 115", "# 43.000000", "i 115 22",
  "# 44.000000", "i 116", "# 45.000000", "i 116 1", "# 46.000000", "i 1 14", "# 47.000000", "i 11 23",
  "# 48.000000", "i 117", "# 49.000000", "i 117 23", "# 50.000000", "i 118", "# 51.000000", "i 1 118",
  "# 52.000000", "c 23 t 22", "# 53.000000", "i 4 22", "# 54.000000", "c 1 t 2", "# 55.000000", "i 2 21",
  "# 56.000000", "i 119", "# 57.000000", "i 2 119", "# 58.000000", "i 120", "# 59.000000", "i 120 22",
  "# 60.000000", "i 2 13", "# 61.000000", "i 12 22", "# 62.000000", "i 121", "# 63.000000", "i 121 22",
  "# 64.000000", "i 122", "# 65.000000", "i 2 122", "# 66.000000", "i 3 22", "# 67.000000")
dgm <- simpers(x, sm)

x <- simplextree::simplex_tree()
x$insert(as.list(seq(0, 24)))
sm <- c("i 23 24", "# 0.038914", "i 0 1", "# 0.038915", "i 22 24", "# 0.042033", "i 0 2", "# 0.042034", "i 123", "# 0.044425", "i 2 123", "# 0.044425", "i 124", "# 0.053772", "i 124 22", "# 0.053773", "i 1 17", "# 0.109929", "i 8 23", "# 0.109930", "i 125", "# 0.114699", "i 125 23", "# 0.114700", "i 126", "# 0.129938", "i 1 126", "# 0.129938", "i 7 22", "# 0.138352", "i 2 18", "# 0.138352", "i 127", "# 0.170217", "i 2 127", "# 0.170217", "i 128", "# 0.177982", "i 128 22", "# 0.177982", "i 1 16", "# 0.201078", "i 9 23", "# 0.201078", "i 129", "# 0.209284", "i 129 23", "# 0.209285", "i 130", "# 0.221994", "i 1 130", "# 0.221994", "i 6 22", "# 0.233548", "i 2 19", "# 0.233548", "i 131", "# 0.261159", "i 2 131", "# 0.261160", "i 132", "# 0.266905", "i 132 22", "# 0.266905", "i 1 15", "# 0.269601", "i 10 23", "# 0.269602", "i 133", "# 0.282695", "i 133 23", "# 0.282695", "i 134", "# 0.306195", "i 1 134", "# 0.306195", "i 5 22", "# 0.312552", "i 2 20", "# 0.312553", "i 135", "# 0.321497", "i 2 135", "# 0.321498", "i 136", "# 0.326445", "i 136 22", "# 0.326445", "i 137", "# 0.340823", "i 137 1", "# 0.340823", "i 1 14", "# 0.356993", "i 11 23", "# 0.356994", "i 138", "# 0.365615", "i 138 23", "# 0.365615", "i 139", "# 0.366602", "i 1 139", "# 0.366602", "c 23 t 22", "# 0.371543", "i 4 22", "# 0.371543", "c 1 t 2", "# 0.395785", "i 2 21", "# 0.395785", "i 140", "# 0.398950", "i 2 140", "# 0.398950", "i 141", "# 0.400866", "i 141 22", "# 0.400867", "i 2 13", "# 0.406149", "i 12 22", "# 0.406150", "i 142", "# 0.423526", "i 142 22", "# 0.423526", "i 143", "# 0.429924", "i 2 143", "# 0.429924", "i 3 22", "# 0.432740")
dgm2 <- simpers(x, sm)




## Reduce birth/death codes to a list
bd_idx <- which(sapply(si_ops, is.numeric))
si_filt <- mapply(function(i, j){ unlist(si_ops[i:(j-1)]) }, c(1, bd_idx[-length(bd_idx)]+1L), bd_idx, SIMPLIFY = FALSE)
si_filt <- lapply(si_filt, function(si_ops){ unlist(lapply(si_ops, function(op){ parse(text=op) })) })

x <- simplextree::simplex_tree()
x$insert(list(1:3, 4:5, 6))
sm <- c("i 4 5 6", "# 1", "i 7", "i 8", "# 2", "i 8 9 10 11", "# 3")
validate_sm(x, sm)


