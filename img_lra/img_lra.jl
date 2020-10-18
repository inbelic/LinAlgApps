include("lra_funcs.jl")

img_name = ARGS[1]
k = ARGS[2]

U = load(img_name)

r_sep, g_sep, b_sep = colour_seperation(U)

r_lra = low_rank_approx(r_sep, k)
g_lra = low_rank_approx(g_sep, k)
b_lra = low_rank_approx(b_sep, k)

lra_img = mix_colours(r_lra, g_lra, b_lra)

output = string(SubString(img_name, 1, length(img_name)-5),"-lra.tiff")

save(output, lra_img)


