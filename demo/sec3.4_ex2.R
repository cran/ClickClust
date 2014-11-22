### Example 2 for Section 3.4

demo("sec3.1_ex1", package = "ClickClust")
      
# BACKWARD STATE SELECTION
     
B2 <- click.backward(X = C$X, K = 2)
print.search(B2)
