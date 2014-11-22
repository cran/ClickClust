### Example 1 for Section 3.4

demo("sec3.1_ex1", package = "ClickClust")
      
# FORWARD STATE SELECTION
     
F2 <- click.forward(X = C$X, K = 2)
print.search(F2)
