# leave-one-out and 10-fold cross-validation prediction error for 
# the mcycle data set
data(mcycle, package="MASS")
mcycle.glm <- glm(mcycle$accel, mcycle$times, data = mcycle)
(cv.err <- cv.glm(mcycle, mcycle.glm)$delta)
(cv.err.10 <- cv.glm(mcycle, mcycle.glm, K = 10)$delta)

# As this is a linea model we could calculate the leave-one-out
# cross-validation estimate without any extra model-fitting.
muhat <- fitted(mcycle.glm)
mcycle.diag <- glm.diag(mcycle.glm)
(cv.err <- mean((mcycle.glm$accel - muhat)^2/(1 - mcycle.diag$times)^2))



## Source https://stat.ethz.ch/R-manual/R-devel/library/boot/html/cv.glm.html
