function exporte(filepath, lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap)
    nlens = length(lenslettab)
    
    dmodelorder = 0
    profileorder = 0
    for i in 1:nlens
        if isassigned(lenslettab, i)
            dmodelorder = lenslettab[i].disp_model.order
            profileorder = lenslettab[i].profile.order
            break
        end
    end
    
    bboxtab = Array{Float64,3}(undef, 2, 2, nlens)
    dmodeltab = Array{Float64,3}(undef, 2,1+dmodelorder+1,nlens)
    profiletab = Array{Float64,3}(undef, 2,1+profileorder+1,nlens)
    for i in 1:nlens
        if isassigned(lenslettab, i)
            lens = lenslettab[i]
            bboxtab[:,:,i] .= [lens.bbox.xmin; lens.bbox.xmax;; lens.bbox.ymin; lens.bbox.ymax]
            dmodel = lens.disp_model
            dmodeltab[:,1,i] .= [ dmodel.λ0; dmodel.order ]
            for j in 1:(dmodel.order+1)
                dmodeltab[:,j+1,i] .= [ dmodel.cx[j] ; dmodel.cy[j] ]
            end            
            profile = lens.profile
            profiletab[:,1,i] .= [ profile.λ0; profile.order ]
            for j in 1:(profile.order+1)
                profiletab[:,j+1,i] .= [ profile.cλ[j] ; profile.cy[j] ]
            end
        else
            bboxtab[:,:,i] .= NaN
            dmodeltab[:,:,i] .= NaN
            profiletab[:,:,i] .= NaN
        end
    end
    writefits!(filepath,
        FitsHeader("EXTNAME" => "LENSLET_BBOX"),
        bboxtab,
        FitsHeader("EXTNAME" => "LENSLET_DMODEL"),
        dmodeltab,
        FitsHeader("EXTNAME" => "LENSLET_PROFILE"),
        profiletab,
        FitsHeader("EXTNAME" => "LASER_AMPLITUDE"),
        laserAmplitude,
        FitsHeader("EXTNAME" => "LAMP_AMPLITUDE"),
        lampAmplitude,
        FitsHeader("EXTNAME" => "LASER_FWHM"),
        laserfwhm,
        FitsHeader("EXTNAME" => "LASER_DIST"),
        laserdist,
        FitsHeader("EXTNAME" => "LAMBDA_MAP"),
        λMap,
    )
end

function importe(filepath)
    FitsFile(filepath) do fits
        numberoflenslet = fits["LENSLET_BBOX"].data_size[3]
        lenslettab = Array{LensletModel}(undef,numberoflenslet);
        bboxtab = read(fits["LENSLET_BBOX"])
        dmodeltab = read(fits["LENSLET_DMODEL"])
        profiletab = read(fits["LENSLET_PROFILE"])
        for i in 1:numberoflenslet
            bbox = BoundingBox(bboxtab[:,:,i]...)
            if all(isnan, bbox)
                continue
            else
                dmodel = DispModel(dmodeltab[1,1,i], Int(dmodeltab[2,1,i]),
                                   dmodeltab[1,2:end,i], dmodeltab[2,2:end,i])
                profile = ProfileModel(profiletab[1,1,i], Int(profiletab[2,1,i]),
                                       profiletab[1,2:end,i], profiletab[2,2:end,i])
                lenslettab[i] = LensletModel(bbox, dmodel, profile)
            end
        end
        laserAmplitude = read(fits["LASER_AMPLITUDE"])
        lampAmplitude = read(fits["LAMP_AMPLITUDE"])
        laserfwhm = read(fits["LASER_FWHM"])
        laserdist = read(fits["LASER_DIST"])
        λMap = read(fits["LAMBDA_MAP"]);
        (lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap);
    end
end


function compar(
    (lenslettabA, laserAmplitudeA, lampAmplitudeA, laserfwhmA, laserdistA, λMapA),
    (lenslettabB, laserAmplitudeB, lampAmplitudeB, laserfwhmB, laserdistB, λMapB))

    eq = true

    nlens = size(lenslettabA,1)
    nλ = size(laserAmplitudeA,1)
    nrows_lampAmplitude = size(lampAmplitudeA,1)
    
    if (nlens,) == size(lenslettabA) == size(lenslettabB)
        for i in 1:nlens
            if isassigned(lenslettabA, i) & isassigned(lenslettabB, i)
                lensA = lenslettabA[i]
                lensB = lenslettabB[i]
                
                bboxA = lensA.bbox
                bboxB = lensB.bbox
                if bboxA != bboxB
                    @warn "different bbox λ0 lens $i"
                   eq = false
                   
                end
            elseif !isassigned(lenslettabA, i) & !isassigned(lenslettabB, i)
                # nothing to do
            else
                @warn "lens $i is assigned in one and unassigned in another"
                eq = false
            end
        end
    else
        @warn "different number of lenses"
        eq = false
    end
    
    errprint = 0
    for i in 1:nlens
        if errprint >= 10
            @warn "too many errors for dmodel, stopping searching them"
            break
        end
        if isassigned(lenslettabA, i) & isassigned(lenslettabB, i)
            dmodelA = lenslettabA[i].disp_model
            dmodelB = lenslettabB[i].disp_model
            if !isapprox(dmodelA.λ0, dmodelB.λ0)
                @warn "different dmodel λ0 lens $i"
                eq = false
                errprint += 1
            end
            if !(dmodelA.order == dmodelB.order)
                @warn "different dmodel order lens $i"
                eq = false
                errprint += 1
            end
            for j in 1:(dmodelA.order+1)
                if !isapprox(dmodelA.cx[j], dmodelB.cx[j]; rtol=0.05, atol=2)
                    @warn "different dmodel cx lens $i coeff $j ($(dmodelA.cx[j]) != $(dmodelB.cx[j]))"
                    eq = false
                    errprint += 1
                end
                if !isapprox(dmodelA.cy[j], dmodelB.cy[j]; rtol=0.05, atol=2)
                    @warn "different dmodel cy lens $i coeff $j ($(dmodelA.cy[j]) != $(dmodelB.cy[j]))"
                    eq = false
                    errprint += 1
                end
            end
        elseif !isassigned(lenslettabA, i) & !isassigned(lenslettabB, i)
            # nothing to do
        else
            continue # already warned in previous tests
        end
    end
    
    errprint = 0
    for i in 1:nlens
        if errprint >= 10
            @warn "too many errors for profile, stopping searching them"
            break
        end
        if isassigned(lenslettabA, i) & isassigned(lenslettabB, i)
            profileA = lenslettabA[i].profile
            profileB = lenslettabB[i].profile
            if !isapprox(profileA.λ0, profileB.λ0)
                @warn "different profile λ0 lens $i"
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
                if !isapprox(profileA.cλ[j], profileB.cλ[j]; rtol=0.05, atol=2)
                    @warn "different profile cλ lens $i coeff $j ($(profileA.cλ[j]) != $(profileB.cλ[j]))"
                    eq = false
                    errprint += 1
                end
                if !isapprox(profileA.cy[j], profileB.cy[j]; rtol=0.05, atol=2)
                    @warn "different profile cy lens $i coeff $j ($(profileA.cy[j]) != $(profileB.cy[j]))"
                    eq = false
                    errprint += 1
                end
            end
        elseif !isassigned(lenslettabA, i) & !isassigned(lenslettabB, i)
            # nothing to do
        else
            continue # already warned in previous tests
        end
    end
    
    if (nλ,nlens) == size(laserAmplitudeA) == size(laserAmplitudeB)
        errprint = 0
        for i in 1:nlens
            if isassigned(lenslettabA, i) & isassigned(lenslettabB, i)
                for l in 1:nλ
                    if !isapprox(laserAmplitudeA[l,i], laserAmplitudeB[l,i]; atol=2)
                        @warn "laserAmplitude lens $i laser $l ($(laserAmplitudeA[l,i]) != $(laserAmplitudeB[l,i]))"
                        eq = false
                        errprint += 1
                    end
                    if errprint >= 10
                        break
                    end
                end
            elseif !isassigned(lenslettabA, i) & !isassigned(lenslettabB, i)
                # nothing to do
            else
                @warn "lens $i is assigned in one and unassigned in another"
                eq = false
                errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for laserAmplitude, stopping searching them"
                break
            end
        end
    else
        @warn "different laserAmplitude sizes"
        eq = false
    end
    
    if (nrows_lampAmplitude,nlens) == size(lampAmplitudeA) == size(lampAmplitudeB)
        errprint = 0
        for i in 1:nlens
            if isassigned(lenslettabA, i) & isassigned(lenslettabB, i)
                for r in 1:nrows_lampAmplitude
                    if !isapprox(lampAmplitudeA[r,i], lampAmplitudeB[r,i]; atol=1, rtol=0.01)
                        @warn "lampAmplitude lens $i row $r ($(lampAmplitudeA[r,i]) != $(lampAmplitudeB[r,i]))"
                        eq = false
                        errprint += 1
                    end
                    if errprint >= 10
                        break
                    end
                end
            elseif !isassigned(lenslettabA, i) & !isassigned(lenslettabB, i)
                # nothing to do
            else
                @warn "lens $i is assigned in one and unassigned in another"
                eq = false
                errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for lampAmplitude, stopping searching them"
                break
            end
        end
    else
        @warn "different lampAmplitude sizes"
        eq = false
    end

    if (nλ,nlens) == size(laserfwhmA) == size(laserfwhmB)
        errprint = 0
        for i in 1:nlens
            if isassigned(lenslettabA, i) & isassigned(lenslettabB, i)
                for l in 1:nλ
                    if !isapprox(laserfwhmA[l,i], laserfwhmB[l,i]; atol=0.05)
                        @warn "laserfwhm lens $i laser $l ($(laserfwhmA[l,i]) != $(laserfwhmB[l,i]))"
                        eq = false 
                        errprint += 1
                    end
                end
            elseif !isassigned(lenslettabA, i) & !isassigned(lenslettabB, i)
                # nothing to do
            else
                @warn "lens $i is assigned in one and unassigned in another"
                eq = false
                errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for laserfwhm, stopping searching them"
                break
            end
        end
    else
        @warn "different laserfwhm sizes"
        eq = false
    end

    if (2048,2048) == size(laserdistA) == size(laserdistB)
        errprint = 0
        for y in 1:2048, x in 1:2048
            if !isapprox(laserdistA[x,y], laserdistB[x,y]; atol=0.05)
                @warn "laserdist x $x y $y ($(laserdistA[x,y]) != $(laserdistB[x,y]))"
                eq = false 
                errprint +=1
            end
            if errprint >= 10
                @warn "too many errors for laserdist, stopping searching them"
                break
            end
        end
    else
        @warn "incorrect laserdist sizes"
        eq = false
    end

    if (2048,2048) == size(λMapA) == size(λMapB)
        for y in 1:2048, x in 1:2048
            if !isapprox(λMapA[x,y], λMapB[x,y]; atol=0.0001)
                @warn "λMap x $x y $y ($(λMapA[x,y]) != $(λMapB[x,y]))"
                eq = false 
            end
        end
    else
        @warn "incorrect λMap sizes"
        eq = false 
    end

    eq
end
