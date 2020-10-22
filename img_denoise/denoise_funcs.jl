using LinearAlgebra, SparseArrays, Images, ImageView

function load_to_float(img_name)
	img = load(img_name)
	img = Gray.(img)
	img = convert(Array{Float64,2}, img)
	return img
end

# Adds random noise to an image for demonstrative purposes
function make_noisy(img, strength)
	m, n = size(img)
	noise = rand(m, n) - rand(m, n)
	noisy_img = img + (noise * strength)
	return map(clamp01nan, noisy_img)
end


# Approximates the solution to the linear solution Ax=b using
# the SOR (Succesive Over Relaxation) algorithm.
# omega, A,b and an inital guess for x are required
# a maximum number of iteration is defaulted to 20000
# and a tolerance of an a solution is defualted to 0.01
#
# returns the approximate solution x, and the iterations taken
function SOR(omega, A, b, x, max_iter = 20000, tol = 0.01)
	m, n = size(A)
	D = sparse(Diagonal(A))
	L = sparse(LowerTriangular(A)) - D
	M = (1/omega)*D + L
	
	global iters_taken = 0
	while iters_taken < max_iter
		global iters_taken = iters_taken + 1
		y = M*x + (b-A*x)
		x = M \ y

		r = b - A*x
		
		if norm(r,2) < tol
			break
		end
	end
	return x, iters_taken 
end


# ######### CG IS CURRENTLY NOT WORKING ################
# Approximates the solution to the linear solution Ax=b using
# the Conjugate Gradiant algorithm.
# A,b and an inital guess for x are required
# a maximum number of iteration is defaulted to 20000
# and a tolerance of an a solution is defualted to 0.01
#
# returns the approximate solution x, and the iterations taken
function CG(A, b, x, max_iter = 20000, tol = 0.01)
	m, n = size(A)
	global r = b - A*x     # r is the residual of our linear system
	
	iter = 0
	global r_rT
	global r_last
	global p
	global alpha
	global A_p

	for dim in 0:(n-1)     # looping through the n dimensional space of A
		iter = iter + 1
	
		if dim == 0 
			beta = 0
			p = r 
		else
			beta = r_rT / (transpose(r_last)*r_last)[1]
		end
		
		# compute values that will be used more than once
		r_rT = (transpose(r)*r)[1]

		p = r + beta*p	
		A_p = A*p
	
		alpha = r_rT / (transpose(p)*A_p)[1]

		x = x + alpha*p

		r_last = r
		r = r - alpha*A_p

		if norm(r,2) < tol*norm(b,2)
			break
		end 
		if iter == max_iter
			break
		end
	end

	return x, iter
end



# Used for scaling the values of A
alpha_h(val, alpha, h) = (val = alpha*val/(h*h)) 

# Creates a finite difference approximation of the Laplacian 
# from a given array U
# 
# returns L a non-normalized lapacian graph
function create_lapacian(U, alpha)
	m, n = size(U)
	len = m*n
	h = m	
	L = spzeros(len, len)
		
	# Ugly for, if loop to place all the pixels weighting in W
	for (U_col, U_row) in Iterators.product(1:n,1:m)
		row = (U_row - 1)*n + U_col			
		
		x_pos = [ U_row U_col ]
		x_int = U[U_row,U_col]
		
		L[row, row] = alpha_h(4*x_int, alpha, h)

		if U_col > 1
			y_pos = [U_row U_col-1]
			added_weight = -1*U[y_pos[1], y_pos[2]]
			L[row,row-1] = alpha_h(added_weight, alpha, h)
		end 
		if U_row > 1
			y_pos = [U_row-1 U_col]
			added_weight = -1*U[y_pos[1], y_pos[2]]
			L[row,row-n] = alpha_h(added_weight, alpha, h)
		end
		if U_col < n
			y_pos = [U_row U_col+1]
			added_weight = -1*U[y_pos[1], y_pos[2]]
			L[row,row+1] = alpha_h(added_weight, alpha, h)
		end 
		if U_row < m
			y_pos = [U_row+1 U_col]
			added_weight = -1*U[y_pos[1], y_pos[2]] 
			L[row,row+n] = alpha_h(added_weight, alpha, h) 
		end
	end

	return L
end


# Flattens our image into a 1-dimensional version of the image
# where the indices
function flatten_img(U)
	m, n = size(U)
	len = m*n
	flat_U = zeros(len,1)
	for (U_col, U_row) in Iterators.product(1:n,1:m)
		row = (U_row - 1)*n + U_col
		flat_U[row] = U[U_row, U_col]
	end	
	return flat_U
end




