function exporte(filepath::String, A::NamedTuple) ::Nothing
    
    (; nλ, lasers_λs, λref, lasers_order, lamp_order, assigned_lenslets, lasers_cxs, lasers_cys,
      lasers_fwhms, lasers_amplitudes, lasers_pixels_dists, lasers_pixels_λs,
      lamp_backs, lamp_amplitudes, bboxes, lamp_cfwhms, lamp_cxs) = A
    
    FitsFile(filepath, "w!") do fits
    
        write(fits, FitsHeader("COMMENT" => "DATA is in a Table HDU"), [0;;])
    
        hdu = FitsTableHDU(fits,
            "ASSIGNED" => Bool,
            "BBOXES" => (Int, 4),
            "LASERS_CXS" => (Float64, lasers_order+1),
            "LASERS_CYS" => (Float64, lasers_order+1),
            "LASERS_FWHMS" => (Float64, nλ),
            "LASERS_AMPLITUDES" => (Float64, nλ),
            "LASERS_PIXELS_DISTS" => (Float64, (BBOX_WIDTH, BBOX_HEIGHT)),
            "LASERS_PIXELS_LAMBDAS" => (Float64, (BBOX_WIDTH, BBOX_HEIGHT)),
            "LAMP_CFWHMS" => (Float64, lamp_order+1),
            "LAMP_CXS" => (Float64, lamp_order+1),
            "LAMP_BACKS" => Float64,
            "LAMP_AMPLITUDES" => (Float64, BBOX_HEIGHT))
            
        hdu["EXTNAME"] = "PIC_DATA"
        hdu["PIC_PACKAGE_VERSION"] = string(pkgversion(PIC))
        hdu["NLENS"] = NLENS
        hdu["NLAMBDA"] = nλ
        for (i,λ) in enumerate(lasers_λs)
            hdu["LASER_LAMBDA_$i"] = λ
        end
        hdu["LAMBDAREF"] = λref
        hdu["LASERS_ORDER"] = lasers_order
        hdu["PROFILE_ORDER"] = lamp_order
        
        write(hdu, "ASSIGNED" => Vector{Bool}(assigned_lenslets))
        write(hdu, "BBOXES" => [ bboxes[i][c] for c in 1:4, i in 1:NLENS ])
        write(hdu, "LASERS_CXS" => lasers_cxs)
        write(hdu, "LASERS_CYS" => lasers_cys)
        write(hdu, "LASERS_FWHMS" => lasers_fwhms)
        write(hdu, "LASERS_AMPLITUDES" => lasers_amplitudes)
        write(hdu, "LASERS_PIXELS_DISTS" => lasers_pixels_dists)
        write(hdu, "LASERS_PIXELS_LAMBDAS" => lasers_pixels_λs)
        write(hdu, "LAMP_CFWHMS" => lamp_cfwhms),
        write(hdu, "LAMP_CXS" => lamp_cxs),
        write(hdu, "LAMP_BACKS" => lamp_backs)
        write(hdu, "LAMP_AMPLITUDES" => lamp_amplitudes)
    end
    
    nothing
end

function importe(filepath)
    FitsFile(filepath) do fits

        hdu = fits["PIC_DATA"]

        nλ = hdu["NLAMBDA"].integer
        lasers_λs = [ hdu["LASER_LAMBDA_$i"].float for i in 1:nλ ]
        λref = hdu["LAMBDAREF"].float
        lasers_order = hdu["LASERS_ORDER"].integer
        lamp_order = hdu["PROFILE_ORDER"].integer

        D = read(hdu)

        assigned_lenslets = BitVector(D["ASSIGNED"])
        bboxes = [ BoundingBox{Int}(D["BBOXES"][:,i]...) for i in 1:NLENS ]
        lasers_cxs = D["LASERS_CXS"]
        lasers_cys = D["LASERS_CYS"]
        lasers_fwhms = D["LASERS_FWHMS"]
        lasers_amplitudes = D["LASERS_AMPLITUDES"]
        lasers_pixels_dists = D["LASERS_PIXELS_DISTS"]
        lasers_pixels_λs = D["LASERS_PIXELS_LAMBDAS"]
        lamp_backs = D["LAMP_BACKS"]
        lamp_amplitudes = D["LAMP_AMPLITUDES"]

        (; nλ, lasers_λs, λref, lasers_order, lamp_order,
           assigned_lenslets, bboxes,
           lasers_cxs, lasers_cys, lasers_fwhms, lasers_amplitudes, lasers_pixels_dists,
           lasers_pixels_λs, lamp_backs, lamp_amplitudes)
    end
end

function compar(A, B)
    eq = true

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
    
    bbox_width  = size(first(A.bboxes),1)
    bbox_height = size(first(A.bboxes),2)

    nλ = A.nλ
    
    errprint = 0
    for i in 1:NLENS
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
    for i in 1:NLENS
        bothassigned[i] || continue
        
        if A.bboxes[i] != B.bboxes[i]
            @warn "lens $i different bboxes: ($(A.bboxes[i]) != $(B.bboxes[i]))"
            eq = false
            errprint += 1
        end
        
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
        
        for y in 1:bbox_height, x in 1:bbox_width
            if !isapprox(A.lasers_pixels_dists[x,y,i], B.lasers_pixels_dists[x,y,i]; atol=0.05, nans=true)
                @warn "lens $i lasers_pixels_dists x $x y $y ($(A.lasers_pixels_dists[x,y,i]) != $(B.lasers_pixels_dists[x,y,i]))"
                eq = false 
                errprint +=1
            end
            if errprint >= 10
                break
            end
        end
        
        for y in 1:bbox_height, x in 1:bbox_width
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
    
    if size(A.lamp_backs) == size(B.lamp_backs)
        errprint = 0
        for i in 1:NLENS
            bothassigned[i] || continue
            if !isapprox(A.lamp_backs[i], B.lamp_backs[i]; atol=1, rtol=0.01, nans=true)
                @warn "lamp_backs lens $i row $r ($(A.lamp_backs[i]) != $(B.lamp_backs[i]))"
                eq = false
                errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for lamp_backs, stopping searching them"
                break
            end
        end
    else
        @warn "different lamp_backs sizes"
        eq = false
    end
    
    if size(A.lamp_amplitudes) == size(B.lamp_amplitudes)
        errprint = 0
        for i in 1:NLENS
            bothassigned[i] || continue
            for r in 1:bbox_height
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
