cs <- 1:10
ys <- c(12, 90, 80, 5, 63, 15, 67, 22, 56, 33)
ns <- c(50, 150, 63, 10, 63, 8, 56, 19, 63, 19)

dat <- data.frame(cs, ys, ns)
colnames(dat) <- c("county", "approve", "disapprove")

dat$prop_approve <- dat$approve / (dat$approve + dat$disapprove)
