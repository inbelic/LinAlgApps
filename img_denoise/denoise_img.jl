include("denoise_funcs.jl")

img_name = ARGS[1]
number_iters = ARGS[2]
alpha = ARGS[3]

img = load_to_float(img_name)

noisy_img = make_noisy(img, 0.25)

global U = noisy_img
global x = flatten_img(U)

m, n = size(U)
len = m * n

for iter in 1:number_iters
	U = transpose(reshape(x, m, n))

	cur_diff = norm(U-img,2)
	println(cur_diff)
	println()
	
	A = create_lapacian(U,alpha) + sparse(I, len, len) 
	b = x
	global x, iters_taken = SOR(omega,A,b,x)
	
	println(iter)
	println(iters_taken)
end

imshow(img)
imshow(U)
imshow(transpose(reshape(x,m,n)))
