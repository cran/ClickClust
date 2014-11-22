### Example 1 for Section 3.2

demo("sec3.1_ex1", package = "ClickClust")

# EM ALGORITHM (without initial state probabilities)
     
N2 <- click.EM(X = C$X, K = 2)
print.EM(N2)
