using PIC

using PyPlot
M = PyPlot

function imshaw(image)
    M.imshow(transpose(image) ;
             origin="lower", extent=[0.5, size(image,1)+0.5, 0.5, size(image,2)+0.5])
end

# you need a FitResult for YJ and a wavelamp file and a specpos file
# for example from test/antoine-test.jl

i = 123

d = readfits("mean-wavelamp.fits")
imshaw(d)
colorbar()

box = result.lenses_boxes[i]
box_xs = [ box.xmin-0.5, box.xmin-0.5, box.xmax+0.5, box.xmax+0.5, box.xmin-0.5]
box_ys = [ box.ymin-0.5, box.ymax+0.5, box.ymax+0.5, box.ymin-0.5, box.ymin-0.5]
plot(box_xs, box_ys; color="red")

poly_xs = PIC.polynome_with_reference(result.λ0, result.wavelamps_fits_cx[:,i]).(PIC.IFS_λRANGE)
poly_ys = PIC.polynome_with_reference(result.λ0, result.wavelamps_fits_cy[:,i]).(PIC.IFS_λRANGE)
plot(poly_xs, poly_ys; color="yellow")

ctrs_x = PIC.polynome_with_reference(result.λ0, result.wavelamps_fits_cx[:,i]).(WAVELAMPS_λLASERS_YJ)
ctrs_y = PIC.polynome_with_reference(result.λ0, result.wavelamps_fits_cy[:,i]).(WAVELAMPS_λLASERS_YJ)
scatter(ctrs_x, ctrs_y; s=8, color="red")

for (j, (ctr_x, ctr_y)) in enumerate(zip(ctrs_x, ctrs_y))

    fwhm = result.wavelamps_fits_fwhm[j,i]
    circle_xs = ctr_x .+ [cos(ang) * fwhm for ang in range(0,2π, 1000)]
    circle_ys = ctr_y .+ [sin(ang) * fwhm for ang in range(0,2π, 1000)]
    plot(circle_xs, circle_ys; color="red")
end

figure()
d = readfits("mean-specpos.fits")
imshaw(d)
colorbar()

box = result.lenses_boxes[i]
box_xs = [ box.xmin-0.5, box.xmin-0.5, box.xmax+0.5, box.xmax+0.5, box.xmin-0.5]
box_ys = [ box.ymin-0.5, box.ymax+0.5, box.ymax+0.5, box.ymin-0.5, box.ymin-0.5]
plot(box_xs, box_ys; color="red")

xs    = PIC.polynome_with_reference(result.λ0, result.specpos_fits_cx[:,i]).(result.wavelamps_λvals[box][3,:])
fwhms = PIC.polynome_with_reference(result.λ0, result.specpos_fits_cλ[:,i]).(result.wavelamps_λvals[box][3,:])
scatter(xs, axes(box,2); s=8, color="red")

for (i, (x, fwhm)) in enumerate(zip(xs, fwhms))
    plot([x - fwhm, x + fwhm], [box.ymin + i - 1, box.ymin + i - 1]; color="orange")
end
