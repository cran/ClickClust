### Example 2 for Section 3.2

demo("sec3.1_ex1", package = "ClickClust")

# EM ALGORITHM (with initial state probabilities)
     
M2 <- click.EM(X = C$X, y = C$y, K = 2)
print.EM(M2)
