function exporte(filepath, A)
    
    (; nlens, nλ, lasers_λs, λref, lens_dx_lower, lens_dx_upper, lens_dy_lower, lens_dy_upper, bbox_width, bbox_height, lasers_order, lamp_order, nrows_lamp_amplitudes, assigned_lenslets, lenslets_models, lamp_amplitudes, lasers_cxs, lasers_cys, lasers_fwhms, lasers_amplitudes, lasers_pixels_dists, lasers_pixels_λs) = A
    
    FitsFile(filepath, "w!") do fits
    
        write(fits, FitsHeader("COMMENT" => "DATA is in a Table HDU"), [0;;])
    
        hdu = FitsTableHDU(fits,
            "ASSIGNED" => Bool,
            "BBOX_XMIN" => Int,
            "BBOX_XMAX" => Int,
            "BBOX_YMIN" => Int,
            "BBOX_YMAX" => Int,
            "LASERS_CXS" => (Float64, lasers_order+1),
            "LASERS_CYS" => (Float64, lasers_order+1),
            "LASERS_FWHMS" => (Float64, nλ),
            "LASERS_AMPLITUDES" => (Float64, nλ),
            "LASERS_PIXELS_DISTS" => (Float64, (bbox_width, bbox_height)),
            "LASERS_PIXELS_LAMBDAS" => (Float64, (bbox_width, bbox_height)),
            "PROFILE_CLAMBDAS" => (Float64, lamp_order+1),
            "PROFILE_CXS" => (Float64, lamp_order+1),
            "LAMP_AMPLITUDE" => (Float64, nrows_lamp_amplitudes))
            
        hdu["EXTNAME"] = "PIC_DATA"
        hdu["PIC_PACKAGE_VERSION"] = string(pkgversion(PIC))
        hdu["NLENS"] = nlens
        hdu["NLAMBDA"] = nλ
        for (i,λ) in enumerate(lasers_λs)
            hdu["LASER_LAMBDA_$i"] = λ
        end
        hdu["LAMBDAREF"] = λref
        hdu["LENS_DX_LOWER"] = lens_dx_lower
        hdu["LENS_DX_UPPER"] = lens_dx_upper
        hdu["LENS_DY_LOWER"] = lens_dy_lower
        hdu["LENS_DY_UPPER"] = lens_dy_upper
        hdu["BBOX_WIDTH"] = bbox_width
        hdu["BBOX_HEIGHT"] = bbox_height
        hdu["LASERS_ORDER"] = lasers_order
        hdu["PROFILE_ORDER"] = lamp_order
        hdu["NROWS_LAMP_AMPLITUDE"] = nrows_lamp_amplitudes
        
        mapdef(f,V,default) = map(eachindex(V)) do i; isassigned(V,i) ? f(V[i]) : default end
        vects_to_mat(vs) = reduce(hcat, vs)
        mats_to_cub(ms) = reduce(ms) do m1, m2; cat(m1, m2; dims=3) end
        
        write(hdu, "ASSIGNED" => Vector{Bool}(assigned_lenslets))
        write(hdu, "BBOX_XMIN" => mapdef(lenslets_models, -1) do lm; lm.bbox.xmin end)
        write(hdu, "BBOX_XMAX" => mapdef(lenslets_models, -1) do lm; lm.bbox.xmax end)
        write(hdu, "BBOX_YMIN" => mapdef(lenslets_models, -1) do lm; lm.bbox.ymin end)
        write(hdu, "BBOX_YMAX" => mapdef(lenslets_models, -1) do lm; lm.bbox.ymax end)
        write(hdu, "LASERS_CXS" => lasers_cxs)
        write(hdu, "LASERS_CYS" => lasers_cys)
        write(hdu, "LASERS_FWHMS" => lasers_fwhms)
        write(hdu, "LASERS_AMPLITUDES" => lasers_amplitudes)
        write(hdu, "LASERS_PIXELS_DISTS" => lasers_pixels_dists)
        write(hdu, "LASERS_PIXELS_LAMBDAS" => lasers_pixels_λs)
        def = fill(NaN64, lamp_order+1)
        write(hdu, "PROFILE_CLAMBDAS" => vects_to_mat(mapdef(lenslets_models, def) do lm
            lm.lamp_model.cλs
        end))
        write(hdu, "PROFILE_CXS" => vects_to_mat(mapdef(lenslets_models, def) do lm
            lm.lamp_model.cxs
        end))
        write(hdu, "LAMP_AMPLITUDE" => lamp_amplitudes)
    end
    
    nothing
end

function importe(filepath)
    FitsFile(filepath) do fits

        hdu = fits["PIC_DATA"]

        nlens = hdu["NLENS"].integer
        nλ = hdu["NLAMBDA"].integer
        lasers_λs = [ hdu["LASER_LAMBDA_$i"].float for i in 1:nλ ]
        λref = hdu["LAMBDAREF"].float
        lens_dx_lower = hdu["LENS_DX_LOWER"].integer
        lens_dx_upper = hdu["LENS_DX_UPPER"].integer
        lens_dy_lower = hdu["LENS_DY_LOWER"].integer
        lens_dy_upper = hdu["LENS_DY_UPPER"].integer
        bbox_width = hdu["BBOX_WIDTH"].integer
        bbox_height = hdu["BBOX_HEIGHT"].integer
        lasers_order = hdu["LASERS_ORDER"].integer
        lamp_order = hdu["PROFILE_ORDER"].integer
        nrows_lamp_amplitudes = hdu["NROWS_LAMP_AMPLITUDE"].integer

        D = read(hdu)

        assigned_lenslets = BitVector(D["ASSIGNED"])
        lenslets_models = Vector{LensletModel}(undef, nlens)

        lasers_cxs = D["LASERS_CXS"]
        lasers_cys = D["LASERS_CYS"]
        lasers_fwhms = D["LASERS_FWHMS"]
        lasers_amplitudes = D["LASERS_AMPLITUDES"]
        lasers_pixels_dists = D["LASERS_PIXELS_DISTS"]
        lasers_pixels_λs = D["LASERS_PIXELS_LAMBDAS"]

        lasers_dists      = fill(NaN64, 2048, 2048)
        λmap              = fill(NaN64, 2048, 2048)
        lamp_amplitudes   = fill(NaN64, nrows_lamp_amplitudes, nlens)
        
        for i in 1:nlens
            assigned_lenslets[i] || continue
            bbox = BoundingBox{Int}(; xmin=D["BBOX_XMIN"][i], xmax=D["BBOX_XMAX"][i],
                                      ymin=D["BBOX_YMIN"][i], ymax=D["BBOX_YMAX"][i])
            lamp_model = ProfileModel(
                λref, lamp_order, D["PROFILE_CLAMBDAS"][:,i], D["PROFILE_CXS"][:,i])
            lenslets_models[i] = LensletModel(bbox, lamp_model)
            lamp_amplitudes[:,i] .= D["LAMP_AMPLITUDE"][:,i]
        end

        (; nlens, nλ, lasers_λs, λref, lens_dx_lower, lens_dx_upper, lens_dy_lower, lens_dy_upper,
           bbox_width, bbox_height, lasers_order, lamp_order, nrows_lamp_amplitudes,
           assigned_lenslets, lenslets_models, lamp_amplitudes,
           lasers_cxs, lasers_cys, lasers_fwhms, lasers_amplitudes, lasers_pixels_dists,
           lasers_pixels_λs)
    end
end

function compar(A, B)
    eq = true

    if A.nlens != B.nlens
        @warn "different number of lenses"
        return false
    end
    if A.nλ != B.nλ
        @warn "different number of lasers"
        return false
    end
    if A.lasers_order != B.lasers_order
        @warn "different lasers order"
        return false
    end
    if A.lamp_order != B.lamp_order
        @warn "different profile order"
        return false
    end
    if (A.bbox_width, A.bbox_height) != (B.bbox_width, B.bbox_height)
        @warn "different bbox size"
        return false
    end

    nlens = A.nlens
    nλ = A.nλ
    nrows_lampAmplitude = size(A.lamp_amplitudes,1)
    
    errprint = 0
    for i in 1:nlens
        if A.assigned_lenslets[i] != B.assigned_lenslets[i]
            if A.assigned_lenslets[i]
                @warn "lens $i assigned in left but not in right"
            else
                @warn "lens $i assigned in right but not in left"
            end
            eq = false
            errprint += 1
        end
        if errprint >= 10
            @warn "too many errors for assigned_lenslets, stopping"
            break
        end
    end
    
    bothassigned = A.assigned_lenslets .& B.assigned_lenslets
    
    errprint = 0
    for i in 1:nlens
        bothassigned[i] || continue
        if isassigned(A.lenslets_models, i) & isassigned(B.lenslets_models, i)
            bboxA = A.lenslets_models[i].bbox
            bboxB = B.lenslets_models[i].bbox
            if bboxA != bboxB
                @warn "lens $i different bboxs: $bboxA $bboxB"
               eq = false
               errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for bboxs, stopping"
                break
            end
        end
    end
    
    errprint = 0
    for i in 1:nlens
        bothassigned[i] || continue
        
        for j in 1:(A.lasers_order+1)
            if !isapprox(A.lasers_cxs[j,i], B.lasers_cxs[j,i]; rtol=0.05, atol=2)
                @warn "lens $i different lasers cxs[$j]: ($(A.lasers_cxs[j,i]) != $(B.lasers_cxs[j,i]))"
                eq = false
                errprint += 1
            end
            if !isapprox(A.lasers_cys[j,i], B.lasers_cys[j,i]; rtol=0.05, atol=2)
                @warn "lens $i different lasers cys[$j]: ($(A.lasers_cys[j,i]) != $(B.lasers_cys[j,i]))"
                eq = false
                errprint += 1
            end
        end
        for l in 1:nλ
            if !isapprox(A.lasers_fwhms[l,i], B.lasers_fwhms[l,i]; atol=0.05)
                @warn "lens $i different lasers fwhms[$l]: ($(A.lasers_fwhms[l,i]) != $(B.lasers_fwhms[l,i]))"
                eq = false 
                errprint += 1
            end
        end
        for l in 1:nλ
            if !isapprox(A.lasers_amplitudes[l,i], B.lasers_amplitudes[l,i]; atol=2)
                @warn "lens $i different lasers amplitudes[$l]: ($(A.lasers_amplitudes[l,i]) != $(B.lasers_amplitudes[l,i]))"
                eq = false
                errprint += 1
            end
            if errprint >= 10
                break
            end
        end
        for y in 1:A.bbox_height, x in 1:A.bbox_width
            if !isapprox(A.lasers_pixels_dists[x,y,i], B.lasers_pixels_dists[x,y,i]; atol=0.05, nans=true)
                @warn "lens $i lasers_pixels_dists x $x y $y ($(A.lasers_pixels_dists[x,y,i]) != $(B.lasers_pixels_dists[x,y,i]))"
                eq = false 
                errprint +=1
            end
            if errprint >= 10
                break
            end
        end
        for y in 1:A.bbox_height, x in 1:A.bbox_width
            if !isapprox(A.lasers_pixels_λs[x,y,i], B.lasers_pixels_λs[x,y,i]; atol=0.0001, nans=true)
                @warn "lens $i lasers_pixels_λs x $x y $y ($(A.lasers_pixels_λs[x,y]) != $(B.lasers_pixels_λs[x,y]))"
                eq = false 
                errprint += 1
            end
            if errprint >= 10
                break
            end
        end
    
        if errprint >= 30
            @warn "too many errors for lasers, stopping"
            break
        end
    end
    
    errprint = 0
    for i in 1:nlens
        if errprint >= 10
            @warn "too many errors for profile, stopping searching them"
            break
        end
        if isassigned(A.lenslets_models, i) & isassigned(B.lenslets_models, i)
            profileA = A.lenslets_models[i].lamp_model
            profileB = B.lenslets_models[i].lamp_model
            if !isapprox(profileA.λref, profileB.λref)
                @warn "different profile λref lens $i"
                eq = false
                errprint += 1
                continue
            end
            if !(profileA.order == profileB.order)
                @warn "different profile order lens $i"
                eq = false
                errprint += 1
                continue
            end
            for j in 1:(profileA.order+1)
                if !isapprox(profileA.cλs[j], profileB.cλs[j]; rtol=0.05, atol=2)
                    @warn "different profile cλs lens $i coeff $j ($(profileA.cλs[j]) != $(profileB.cλs[j]))"
                    eq = false
                    errprint += 1
                end
                if !isapprox(profileA.cxs[j], profileB.cxs[j]; rtol=0.05, atol=2)
                    @warn "different profile cxs lens $i coeff $j ($(profileA.cxs[j]) != $(profileB.cxs[j]))"
                    eq = false
                    errprint += 1
                end
            end
        elseif !isassigned(A.lenslets_models, i) & !isassigned(B.lenslets_models, i)
            # nothing to do
        else
            continue # already warned in previous tests
        end
    end
    
    if (nrows_lampAmplitude,nlens) == size(A.lamp_amplitudes) == size(B.lamp_amplitudes)
        errprint = 0
        for i in 1:nlens
            bothassigned[i] || continue
            for r in 1:nrows_lampAmplitude
                if !isapprox(A.lamp_amplitudes[r,i], B.lamp_amplitudes[r,i]; atol=1, rtol=0.01, nans=true)
                    @warn "lamp_amplitudes lens $i row $r ($(A.lamp_amplitudes[r,i]) != $(B.lamp_amplitudes[r,i]))"
                    eq = false
                    errprint += 1
                end
                if errprint >= 10
                    break
                end
            end
            if errprint >= 10
                @warn "too many errors for lamp_amplitudes, stopping searching them"
                break
            end
        end
    else
        @warn "different lamp_amplitudes sizes"
        eq = false
    end
    eq
end
