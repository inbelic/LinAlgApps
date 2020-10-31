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
Now by theorem we can represent this multiplication of U * S * transpose(T) as a sum
of r rank-one matrices. So,
A = U * S * transpose(T) 
= summation of j = 1:r, S[j] * U[j,:] * transpose(V[j,:]) 

So now we can form an approximation of A by taking only the first k <= r of the summation.
So to demonstrate this we first take an image and represent it in 3 matrices that represent
the R,G,B values respectively. We then take the SVD of each matrix and choose a k-th rank 
approximation to take for each matrix. We then recreate a coloured image from the 
approximate R,G,B matrices. In the samples we have a k = 5,10 and 50 approximations for the 
three images lake.tiff, mandrill.tiff and waterfall.tiff.

So we see that as k is around 50 our approximations are quite accurate.
So what makes it worthwhile? Well in the case of the waterfall.tiff image it was a 337x337 image
so we require a matrix with 113569 elements that has a R,G,B value stored in each
to represent this image. Yet if we consider the k=50 approximation of the image we only need
to store 50 * 1 + 50 * 337 + 50 * 337 elements each with a R,G,B value stored in it and this 
totals only 50 + 16850 + 16850 = 33750 elements!

So we can then do the same process on other data sets that are able to be represented as a matrix
and we can look at the "most important" data points.


## Image Denoising
Image denoising is and will only become more and more prevalent in the future. Currently,
one of the biggest uses of image denoising is in medical imaging and going into the future
we see that synthetic images generated from raytracing will produce noise if we don't run
the raytracing for a long period of time. So as these simulations get increasingly complex
we can instead raytrace for a small time and then clean up the synthetic image using image denoising.

So I am presenting a denoising function that uses a Lapacian Regularization in the denoising 
function. However, we could easily modify the code to take into account Tikhonov Regularization 
or Total Variation Regularization in place of the create_lapacian function to alter our results. 

Our create_lapacian function creates a Lapacian graph of the matrix representation U with the
equation 4U(i,i) - U(i,i+1) - U(i+1,i) - U(i,i-1) - U(i-1,i). 

Then the "core" of the denoising algorithm is to repeatedly solve the system U^(k+1)*x = U^k
where x iteratively becomes closer to the denoised image. 




