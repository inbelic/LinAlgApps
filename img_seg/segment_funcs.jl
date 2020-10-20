using SparseArrays, LinearAlgebra, Arpack, Clustering, Images

# Merges the segmentation image array and the original
# grayscale array to create
function comp_merge(U, img_seg)
	comp = Gray.(U) + img_seg	
	return map(clamp01nan, comp)		
end


# Turns an image into an array of Float64's that represent
# its pixel intensity
function img_to_array(img_name)
	img = load(img_name)
	img = Gray.(img)
	U = convert(Array{Float64,2}, img)
	return U
end

# Returns a unique color given the
# current index, k
function get_coloring(k)
	r = (1/(k+1))*mod(k,3)
	b = (1/(k+1))*mod(k+1,3)
	g = (1/(k+1))*mod(k+2,3)
	return r, b, g
end

# Turns the index array of Int64, index,
# given from the kmeans function into an
# RGBA array, img_array, used to create a 
# new image
# s denotes the size of the original image
# K determines the number of clusters
function array_to_img(index, s, K)
	colours = Array{RGBA{N0f8},1}() #initialize colouring array
	for k in 1:K
		r, b, g = get_coloring(k)
		push!(colours, RGBA{N0f8}(r,b,g,0.85))
	end
	
	img_array = Array{RGBA{N0f8},1}()
	for i in 1:(s[1]*s[2])
		push!(img_array, colours[index[i]])
	end
	
	return transpose(reshape(img_array, s[2], s[1]))
end

# Returns the real eigenvalues of a lapacian matrix
# (returns the k smallest magnitude eigenvalues
# and their corresponding eigenvectors)
function get_real_eigs(lapacian, K)
	complex_D, complex_V = eigs(lapacian, nev=K, which=:SM)
	return real(complex_D), real(complex_V)
end

# Returns an index array of the pixels after clustering
function get_index(V, K)
	m, n = size(V)
	P = zeros(m, n)
	for row in 1:m
		P[row,:] = V[row,:]/norm(V[row,:])
	end
	R = kmeans(transpose(P), K; maxiter=2000, display=:none)
	return assignments(R)
end


# Provides a sample weighting for segmentation based
# on pixel intensity
# x_pos, y_pos represent the x,y position in the image
# respectively
# x_int, y_int represent the x,y pixel intensity
# respectively, both within [0,1] for grayscale images
#
# returns weight used for assessment in spectral clustering
function sample_weight(x_pos, y_pos, x_int, y_int)
	dist_weighting = 150
	int_weighting  = 1/500
	
	dist_val = norm(x_pos - y_pos)^2
	dist_val = exp(-dist_val/dist_weighting)

	int_val = norm(x_int - y_int)^2
	int_val = exp(-int_val/int_weighting)

	return dist_val * int_val
end

# Creates a lapacian graph from the original image 
# that is given as an array, U, and a weight function
# weight, where the type of array can be any datatype
# compatible with the given function weight.
#
# returns a normalized lapacian graph of size m*n x m*n
function create_lapacian(U, weight)
	m, n = size(U)
	len = m*n
	
	W = spzeros(len,len)
	D = spzeros(len,len)
	
	# Ugly for, if loop to place all the pixels weighting in W
	for (U_col, U_row) in Iterators.product(1:n,1:m)
		row = (U_row - 1)*n + U_col			
		d_sum = 0
		
		x_pos = [ U_row U_col ]
		x_int = U[U_row,U_col]
	
		if U_col > 1
			y_pos = [U_row U_col-1]
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]])
			W[row,row-1] = added_weight
			d_sum = d_sum + added_weight
		end 
		if U_row > 1 && U_col < n
			y_pos = [U_row-1 U_col+1]
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]])
			W[row,row-n+1] = added_weight
			d_sum = d_sum + added_weight
		end
		if U_row > 1
			y_pos = [U_row-1 U_col]
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]])
			W[row,row-n] = added_weight
			d_sum = d_sum + added_weight
		end
		if U_row > 1 && U_col > 1 
			y_pos = [U_row-1 U_col-1]
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]])
			W[row,row-n-1] = added_weight
			d_sum = d_sum + added_weight
		end 
		if U_col < n
			y_pos = [U_row U_col+1]
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]])
			W[row,row+1] = added_weight
			d_sum = d_sum + added_weight
		end 
		if U_row < m && U_col > 1
			y_pos = [U_row+1 U_col-1]
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]])
			W[row,row+n-1] = added_weight
			d_sum = d_sum + added_weight
		end
		if U_row < m
			y_pos = [U_row+1 U_col];
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]]); 
			W[row,row+n] = added_weight;
			d_sum = d_sum + added_weight;
		end 
		if U_row < m && U_col < n
			y_pos = [U_row+1 U_col+1];
			added_weight = weight(x_pos, y_pos, x_int, U[y_pos[1], y_pos[2]]);
			W[row,row+n+1] = added_weight;
			d_sum = d_sum + added_weight;
		end
		
		D[row,row] = d_sum;
	end

	# Calculate inverse of D manually to avoid using built in inverse
	for row in 1:len
		D[row,row] = D[row,row]^(1/2)
		D[row,row] = 1/D[row,row]
	end
	
	
	# Compute the normalized lapacian from the equation
	# NL = Iden - DWD
	# We note that Iden, D, W are all of the same dimensions
	
	Iden = spzeros(len, len) + I #Initalize the identity matrix
	WD = sparse(W*D)
	DWD = sparse(D*WD)
	NL = Iden - DWD
	NL = sparse(NL)
	
	return NL;
end


