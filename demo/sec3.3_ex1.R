### Example 1 for Section 3.4

demo("sec3.2_ex1")

V <- click.var(X = C$X, alpha = N2$alpha, gamma = N2$gamma, z = N2$z)
st.err <- sqrt(diag(V))

Estimates <- c(N2$alpha[-K], as.vector(apply(N2$gamma[,-p,], 3, t)))

# 95% confidence intervals for parameter estimates

Lower <- Estimates - qnorm(0.975) * st.err
Upper <- Estimates + qnorm(0.975) * st.err

cbind(Estimates, Lower, Upper)
