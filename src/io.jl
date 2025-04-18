function exporte(filepath, t)

    (lenslets_models, lasers_fwhms, lasers_amplitudes, lasers_dists, λmap, lamp_amplitudes) = t
    
    nlens = length(lenslets_models)
    
    dmodelorder = 0
    profileorder = 0
    for i in 1:nlens
        if isassigned(lenslets_models, i)
            dmodelorder = lenslets_models[i].disp_model.order
            profileorder = lenslets_models[i].profile_model.order
            break
        end
    end
    
    bboxtab = Array{Float64,3}(undef, 2, 2, nlens)
    dmodeltab = Array{Float64,3}(undef, 2,1+dmodelorder+1,nlens)
    profiletab = Array{Float64,3}(undef, 2,1+profileorder+1,nlens)
    for i in 1:nlens
        if isassigned(lenslets_models, i)
            lens = lenslets_models[i]
            bboxtab[:,:,i] .= [lens.bbox.xmin; lens.bbox.xmax;; lens.bbox.ymin; lens.bbox.ymax]
            dmodel = lens.disp_model
            dmodeltab[:,1,i] .= [ dmodel.λref; dmodel.order ]
            for j in 1:(dmodel.order+1)
                dmodeltab[:,j+1,i] .= [ dmodel.cxs[j] ; dmodel.cys[j] ]
            end            
            profile = lens.profile_model
            profiletab[:,1,i] .= [ profile.λref; profile.order ]
            for j in 1:(profile.order+1)
                profiletab[:,j+1,i] .= [ profile.cλs[j] ; profile.cxs[j] ]
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
        FitsHeader("EXTNAME" => "LASER_FWHM"),
        laserfwhm,
        FitsHeader("EXTNAME" => "LASER_AMPLITUDE"),
        laserAmplitude,
        FitsHeader("EXTNAME" => "LASER_DIST"),
        lasers_dists,
        FitsHeader("EXTNAME" => "LAMBDA_MAP"),
        λmap,
        FitsHeader("EXTNAME" => "LAMP_AMPLITUDE"),
        lampAmplitude
    )
end

function importe(filepath)
    FitsFile(filepath) do fits
        nlens = fits["LENSLET_BBOX"].data_size[3]
        lenslets_models = Array{LensletModel}(undef,nlens);
        bboxtab = read(fits["LENSLET_BBOX"])
        dmodeltab = read(fits["LENSLET_DMODEL"])
        profiletab = read(fits["LENSLET_PROFILE"])
        for i in 1:nlens
            bbox = BoundingBox(bboxtab[:,:,i]...)
            if all(isnan, bbox)
                continue
            else
                dmodel = DispersionModel(dmodeltab[1,1,i], Int(dmodeltab[2,1,i]),
                                   dmodeltab[1,2:end,i], dmodeltab[2,2:end,i])
                profile = ProfileModel(profiletab[1,1,i], Int(profiletab[2,1,i]),
                                       profiletab[1,2:end,i], profiletab[2,2:end,i])
                lenslets_models[i] = LensletModel(bbox, dmodel, profile)
            end
        end
        lasers_fwhms = read(fits["LASER_FWHM"])
        lasers_amplitudes = read(fits["LASER_AMPLITUDE"])
        lasers_dists = read(fits["LASER_DIST"])
        λmap = read(fits["LAMBDA_MAP"]);
        lamp_amplitudes = read(fits["LAMP_AMPLITUDE"])
        (lenslets_models, lasers_fwhms, lasers_amplitudes, lasers_dists, λmap, lamp_amplitudes)
    end
end


function compar(
    (lenslets_modelsA, lasers_fwhmsA, lasers_amplitudesA, lasers_distsA, λmapA, lamp_amplitudesA),
    (lenslets_modelsB, lasers_fwhmsB, lasers_amplitudesB, lasers_distsB, λmapB, lamp_amplitudesB))

    eq = true

    nlens = size(lenslets_modelsA,1)
    nλ = size(lasers_amplitudesA,1)
    nrows_lampAmplitude = size(lamp_amplitudesA,1)
    
    if (nlens,) == size(lenslets_modelsA) == size(lenslets_modelsB)
        for i in 1:nlens
            if isassigned(lenslets_modelsA, i) & isassigned(lenslets_modelsB, i)
                lensA = lenslets_modelsA[i]
                lensB = lenslets_modelsB[i]
                
                bboxA = lensA.bbox
                bboxB = lensB.bbox
                if bboxA != bboxB
                    @warn "different bbox λref lens $i"
                   eq = false
                   
                end
            elseif !isassigned(lenslets_modelsA, i) & !isassigned(lenslets_modelsB, i)
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
        if isassigned(lenslets_modelsA, i) & isassigned(lenslets_modelsB, i)
            dmodelA = lenslets_modelsA[i].disp_model
            dmodelB = lenslets_modelsB[i].disp_model
            if !isapprox(dmodelA.λref, dmodelB.λref)
                @warn "different dmodel λref lens $i"
                eq = false
                errprint += 1
            end
            if !(dmodelA.order == dmodelB.order)
                @warn "different dmodel order lens $i"
                eq = false
                errprint += 1
            end
            for j in 1:(dmodelA.order+1)
                if !isapprox(dmodelA.cxs[j], dmodelB.cxs[j]; rtol=0.05, atol=2)
                    @warn "different dmodel cxs lens $i coeff $j ($(dmodelA.cxs[j]) != $(dmodelB.cxs[j]))"
                    eq = false
                    errprint += 1
                end
                if !isapprox(dmodelA.cys[j], dmodelB.cys[j]; rtol=0.05, atol=2)
                    @warn "different dmodel cys lens $i coeff $j ($(dmodelA.cys[j]) != $(dmodelB.cys[j]))"
                    eq = false
                    errprint += 1
                end
            end
        elseif !isassigned(lenslets_modelsA, i) & !isassigned(lenslets_modelsB, i)
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
        if isassigned(lenslets_modelsA, i) & isassigned(lenslets_modelsB, i)
            profileA = lenslets_modelsA[i].profile_model
            profileB = lenslets_modelsB[i].profile_model
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
        elseif !isassigned(lenslets_modelsA, i) & !isassigned(lenslets_modelsB, i)
            # nothing to do
        else
            continue # already warned in previous tests
        end
    end
    
    if (nλ,nlens) == size(lasers_amplitudesA) == size(lasers_amplitudesB)
        errprint = 0
        for i in 1:nlens
            if isassigned(lenslets_modelsA, i) & isassigned(lenslets_modelsB, i)
                for l in 1:nλ
                    if !isapprox(lasers_amplitudesA[l,i], lasers_amplitudesB[l,i]; atol=2)
                        @warn "lasers_amplitudes lens $i laser $l ($(lasers_amplitudesA[l,i]) != $(lasers_amplitudesB[l,i]))"
                        eq = false
                        errprint += 1
                    end
                    if errprint >= 10
                        break
                    end
                end
            elseif !isassigned(lenslets_modelsA, i) & !isassigned(lenslets_modelsB, i)
                # nothing to do
            else
                @warn "lens $i is assigned in one and unassigned in another"
                eq = false
                errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for lasers_amplitudes, stopping searching them"
                break
            end
        end
    else
        @warn "different lasers_amplitudes sizes"
        eq = false
    end
    
    if (nrows_lampAmplitude,nlens) == size(lamp_amplitudesA) == size(lamp_amplitudesB)
        errprint = 0
        for i in 1:nlens
            if isassigned(lenslets_modelsA, i) & isassigned(lenslets_modelsB, i)
                for r in 1:nrows_lampAmplitude
                    if !isapprox(lamp_amplitudesA[r,i], lamp_amplitudesB[r,i]; atol=1, rtol=0.01, nans=true)
                        @warn "lamp_amplitudes lens $i row $r ($(lamp_amplitudesA[r,i]) != $(lamp_amplitudesB[r,i]))"
                        eq = false
                        errprint += 1
                    end
                    if errprint >= 10
                        break
                    end
                end
            elseif !isassigned(lenslets_modelsA, i) & !isassigned(lenslets_modelsB, i)
                # nothing to do
            else
                @warn "lens $i is assigned in one and unassigned in another"
                eq = false
                errprint += 1
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

    if (nλ,nlens) == size(lasers_fwhmsA) == size(lasers_fwhmsB)
        errprint = 0
        for i in 1:nlens
            if isassigned(lenslets_modelsA, i) & isassigned(lenslets_modelsB, i)
                for l in 1:nλ
                    if !isapprox(lasers_fwhmsA[l,i], lasers_fwhmsB[l,i]; atol=0.05)
                        @warn "lasers_fwhms lens $i laser $l ($(lasers_fwhmsA[l,i]) != $(lasers_fwhmsB[l,i]))"
                        eq = false 
                        errprint += 1
                    end
                end
            elseif !isassigned(lenslets_modelsA, i) & !isassigned(lenslets_modelsB, i)
                # nothing to do
            else
                @warn "lens $i is assigned in one and unassigned in another"
                eq = false
                errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for lasers_fwhms, stopping searching them"
                break
            end
        end
    else
        @warn "different lasers_fwhms sizes"
        eq = false
    end

    if (2048,2048) == size(lasers_distsA) == size(lasers_distsB)
        errprint = 0
        for y in 1:2048, x in 1:2048
            if !isapprox(lasers_distsA[x,y], lasers_distsB[x,y]; atol=0.05, nans=true)
                @warn "lasers_dists x $x y $y ($(lasers_distsA[x,y]) != $(lasers_distsB[x,y]))"
                eq = false 
                errprint +=1
            end
            if errprint >= 10
                @warn "too many errors for lasers_dists, stopping searching them"
                break
            end
        end
    else
        @warn "incorrect lasers_dists sizes"
        eq = false
    end

    if (2048,2048) == size(λmapA) == size(λmapB)
        errprint = 0
        for y in 1:2048, x in 1:2048
            if !isapprox(λmapA[x,y], λmapB[x,y]; atol=0.0001, nans=true)
                @warn "λmap x $x y $y ($(λmapA[x,y]) != $(λmapB[x,y]))"
                eq = false 
                errprint += 1
            end
            if errprint >= 10
                @warn "too many errors for λmap, stopping searching them"
                break
            end
        end
    else
        @warn "incorrect λmap sizes"
        eq = false 
    end

    eq
end
