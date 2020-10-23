# Linear Algebra Applications
The goal of this mini project is to demonstrate some of the applications 
of computational linear algebra in image processing. Some of the topics
covered were topics used in CS 475/675 as motivation for the concepts being taught.
When covering these applications for the course we used MATLAB for most of the implementation.
So this is my way to learn the LinearAlgebra, Image and SparseArrays packages and get
further practice with the Julia language.

## Low Rank Approximations
Low rank approximations are used to be able to approximate a matrix of values.
So what makes the low rank approximation worthwhile?
Well, if we can determine the SVD of a matrix A then we have a way to represent
A as a multiplication of the left singular values, U, singular values, S and 
the right singular values, V.
Now by theorem we can represent this multiplication of U*S*transpose(T) as a sum
of r rank-one matrices. So,
A = U*S*transpose(T) 
= summation of j = 1:r, S[j]*U[j,:]*transpose(V[j,:]) 

So now we can form an approximation of A by taking only the first k <= r of the summation.
So to demonstrate this we first take an image and represent it in 3 matrices that represent
the R,G,B values respectively. We then take the SVD of each matrix and choose a k-th rank 
approximation to take for each matrix. We then recreate a coloured image from the 
approximate R,G,B matrices. In the samples we have a k = 5,10 and 50 approximations for the 
three images lake.tiff, mandrill.tiff and waterfall.tiff.

So we see that as k is around 50 our approximations are quite accurate.
So what makes it worthwhile? Well in the case of the lake.tiff image it was a 900x997 image
so we require a matrix with 879300 elements that has a R,G,B value stored in each
to represent this image. Yet if we consider the k=50 approximation of the image we only need
to store 50*1 + 50*900 + 50*977 elements each with a R,G,B value stored in it and this totals only
45000 + 48850 + 50 = 93900 elements!

So we can then do the same process on other data sets that are able to be represented as a matrix
and we can look at the "most important" data points.





