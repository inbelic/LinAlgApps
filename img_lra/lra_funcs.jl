using LinearAlgebra, Images

# Returns three arrays, r_sep, g_sep, b_sep that contain the 
# R, G, B values of the original img respectively
function colour_seperation(img)
	m, n = size(img)
	r_sep = Float64[]
	g_sep = Float64[]
	b_sep = Float64[]

	for pixel in img
		push!(r_sep, pixel.r)
		push!(g_sep, pixel.g)
		push!(b_sep, pixel.b)
	end

	r_sep = reshape(r_sep, m, n)	
	g_sep = reshape(g_sep, m, n)
	b_sep = reshape(b_sep, m, n)

	return r_sep, g_sep, b_sep
end

# Combines three arrays representing the r,g,b values
# of an image and mixes the arrays to create the 
# coloured image of the three arrays combined
function mix_colours(r_sep, g_sep, b_sep)
	m, n = size(r_sep) # arbitrary choice of r,g,b_sep
	mixed_img = RGBA{N0f8}[];
	for (col, row) in Iterators.product(1:n, 1:m)
		r = clamp01nan(r_sep[row,col])
		g = clamp01nan(g_sep[row,col])
		b = clamp01nan(b_sep[row,col])
		pixel = RGBA{N0f8}(r,g,b,1.0)
		push!(mixed_img, pixel)
	end
	return transpose(reshape(mixed_img, n, m))
end


# Determine the k-th rank approximation of a matrix
# using the SVD of the matrxi
function low_rank_approx(matrix, k)
	r = rank(matrix)
	U, S, V = svd(matrix)
	S = Diagonal(S)

	return U[:,1:k]*S[1:k,1:k]*transpose(V[:,1:k])
end


