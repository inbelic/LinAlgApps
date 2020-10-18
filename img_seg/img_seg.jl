include("segment_funcs.jl")

img_name = ARGS[1]
K = ARGS[2]

U = img_to_array(img_name)
s = size(U)


lapacian = create_lapacian(U, sample_weight)

D, V = get_real_eigs(lapacian, K)

index = get_index(V, K)

img_seg = array_to_img(index, s, K) 

save("seg.tiff", img_seg)

comp = comp_merge(U, img_seg)
save("comp.tiff", comp)

