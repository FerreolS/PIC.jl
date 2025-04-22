function exporte(filepath, A)

    (; lenslets_models, lasers_dists, λmap, lamp_amplitudes) = A
    
    nlens = length(lenslets_models)
    
    assigned_lenslets =
        map(i -> isassigned(lenslets_models, i), eachindex(lenslets_models))
    
    nλs = unique(map(lens -> lens.lasers_model.nλ, lenslets_models[assigned_lenslets]))
    lasers_orders = unique(map(lens -> lens.lasers_model.order, lenslets_models[assigned_lenslets]))
    profile_orders = unique(map(lens -> lens.profile_model.order, lenslets_models[assigned_lenslets]))
    
    length(nλs) == 1 || throw(ArgumentError)
    length(lasers_orders) == 1 || throw(ArgumentError)
    length(profile_orders) == 1 || throw(ArgumentError)
    
    nλ = nλs[1]
    lasers_order = lasers_orders[1]
    profile_order = profile_orders[1]

    bboxs_array = fill(NaN64, 4, nlens)

    lasers_λrefs_array = fill(NaN64, nlens) 
    lasers_cxs_array = fill(NaN64, lasers_order + 1, nlens)
    lasers_cys_array = fill(NaN64, lasers_order + 1, nlens)
    lasers_fwhms_array = fill(NaN64, nλ, nlens)
    lasers_amplitudes_array = fill(NaN64, nλ, nlens)

    profile_λrefs_array = fill(NaN64, nlens) 
    profile_cλs_array = fill(NaN64, profile_order + 1, nlens)
    profile_cxs_array = fill(NaN64, profile_order + 1, nlens)

    for i in eachindex(lenslets_models)
        assigned_lenslets[i] || continue
    
        lens = lenslets_models[i]
        
        bbox = lens.bbox
        bboxs_array[:,i] .= [ bbox.xmin; bbox.xmax; bbox.ymin; bbox.ymax ]
        
        lasers_model = lens.lasers_model
        lasers_λrefs_array[i] = lasers_model.λref
        lasers_cxs_array[:,i] .= lasers_model.cxs
        lasers_cys_array[:,i] .= lasers_model.cys
        lasers_fwhms_array[:,i] .= lasers_model.fwhms
        lasers_amplitudes_array[:,i] .= lasers_model.amplitudes

        profile_model = lens.profile_model
        profile_λrefs_array[i] = profile_model.λref
        profile_cλs_array[:,i] .= profile_model.cλs
        profile_cxs_array[:,i] .= profile_model.cxs
    end
    
    writefits!(filepath,
        FitsHeader("EXTNAME" => "ASSIGNED", "NLENS" => nlens, "NLAMBDA" => nλ,
                   "LASERS_ORDER" => lasers_order, "PROFILE_ORDER" => profile_order),
        reshape(assigned_lenslets, Val(2)),# matrix so DS9 does not crash
        FitsHeader("EXTNAME" => "BBOXS"), bboxs_array,
        FitsHeader("EXTNAME" => "LASERS_LAMBDAREFS"),
        reshape(lasers_λrefs_array, Val(2)), # matrix so DS9 does not crash
        FitsHeader("EXTNAME" => "LASERS_CXS"), lasers_cxs_array,
        FitsHeader("EXTNAME" => "LASERS_CYS"), lasers_cys_array,
        FitsHeader("EXTNAME" => "LASERS_FWHMS"), lasers_fwhms_array,
        FitsHeader("EXTNAME" => "LASERS_AMPLITUDES"), lasers_amplitudes_array,
        FitsHeader("EXTNAME" => "PROFILE_LAMBDAREFS"),
        reshape(profile_λrefs_array, Val(2)), # matrix so DS9 does not crash
        FitsHeader("EXTNAME" => "PROFILE_CLAMBDAS"), profile_cλs_array,
        FitsHeader("EXTNAME" => "PROFILE_CXS"), profile_cxs_array,
        FitsHeader("EXTNAME" => "LASER_DIST"), lasers_dists,
        FitsHeader("EXTNAME" => "LAMBDA_MAP"), λmap,
        FitsHeader("EXTNAME" => "LAMP_AMPLITUDE"), lamp_amplitudes
    )
end

function importe(filepath)
    FitsFile(filepath) do fits
        
        nlens = fits[1]["NLENS"].integer
        nλ = fits[1]["NLAMBDA"].integer
        lasers_order = fits[1]["LASERS_ORDER"].integer
        profile_order = fits[1]["PROFILE_ORDER"].integer
        
        assigned_lenslets = reshape(read(Array{Bool}, fits["ASSIGNED"]), Val(1))
        bboxs_array = read(fits["BBOXS"])
        lasers_λrefs_array = reshape(read(fits["LASERS_LAMBDAREFS"]), Val(1))
        lasers_cxs_array = read(fits["LASERS_CXS"])
        lasers_cys_array = read(fits["LASERS_CYS"])
        lasers_fwhms_array = read(fits["LASERS_FWHMS"])
        lasers_amplitudes_array = read(fits["LASERS_AMPLITUDES"])
        profile_λrefs_array = reshape(read(fits["PROFILE_LAMBDAREFS"]), Val(1))
        profile_cλs_array = read(fits["PROFILE_CLAMBDAS"])
        profile_cxs_array = read(fits["PROFILE_CXS"])
        
        lenslets_models = Array{LensletModel}(undef, nlens)
        for i in 1:nlens
            assigned_lenslets[i] || continue
        
            bbox = BoundingBox(bboxs_array[:,i]...)

            lasers_λref = lasers_λrefs_array[i]
            lasers_cxs = lasers_cxs_array[:,i]
            lasers_cys = lasers_cys_array[:,i]
            lasers_fwhms = lasers_fwhms_array[:,i]
            lasers_amplitudes = lasers_amplitudes_array[:,i]
            lasers_model = LasersModel(
                nλ, lasers_order, lasers_λref, lasers_cxs, lasers_cys, lasers_fwhms, lasers_amplitudes)
            
            profile_λref = profile_λrefs_array[i]
            profile_cλs = profile_cλs_array[:,i]
            profile_cxs = profile_cxs_array[:,i]
            profile_model = ProfileModel(
                profile_λref, profile_order, profile_cλs, profile_cxs)

            lenslets_models[i] = LensletModel(bbox, lasers_model, profile_model)
        end
        
        lasers_dists = read(fits["LASER_DIST"])
        λmap = read(fits["LAMBDA_MAP"]);
        lamp_amplitudes = read(fits["LAMP_AMPLITUDE"])
        (; nlens, nλ, lasers_order, profile_order, assigned_lenslets, lenslets_models,
         lasers_dists, λmap, lamp_amplitudes)
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
    if A.profile_order != B.profile_order
        @warn "different profile order"
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
        dmodelA = A.lenslets_models[i].lasers_model
        dmodelB = B.lenslets_models[i].lasers_model
        if !isapprox(dmodelA.λref, dmodelB.λref)
            @warn "lens $i different lasers λref: $(dmodelA.λref) $(dmodelB.λref)"
            eq = false
            errprint += 1
        end
        if !(dmodelA.order == dmodelB.order)
            @warn "lens $i different lasers order: $(dmodelA.order) $(dmodelB.order)"
            eq = false
            errprint += 1
        end
        for j in 1:(dmodelA.order+1)
            if !isapprox(dmodelA.cxs[j], dmodelB.cxs[j]; rtol=0.05, atol=2)
                @warn "lens $i different lasers cxs[$j]: ($(dmodelA.cxs[j]) != $(dmodelB.cxs[j]))"
                eq = false
                errprint += 1
            end
            if !isapprox(dmodelA.cys[j], dmodelB.cys[j]; rtol=0.05, atol=2)
                @warn "lens $i different lasers cys[$j]: ($(dmodelA.cys[j]) != $(dmodelB.cys[j]))"
                eq = false
                errprint += 1
            end
        end
        for l in 1:nλ
            if !isapprox(dmodelA.fwhms[l], dmodelB.fwhms[l]; atol=0.05)
                @warn "lens $i different lasers fwhms[$l]: ($(dmodelA.fwhms[l]) != $(dmodelB.fwhms[l]))"
                eq = false 
                errprint += 1
            end
        end
        for l in 1:nλ
            if !isapprox(dmodelA.amplitudes[l], dmodelB.amplitudes[l]; atol=2)
                @warn "lens $i different lasers amplitudes[$l]: ($(dmodelA.amplitudes[l]) != $(dmodelB.amplitudes[l]))"
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
            profileA = A.lenslets_models[i].profile_model
            profileB = B.lenslets_models[i].profile_model
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
            if isassigned(A.lenslets_models, i) & isassigned(B.lenslets_models, i)
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
            elseif !isassigned(A.lenslets_models, i) & !isassigned(B.lenslets_models, i)
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

    if (2048,2048) == size(A.lasers_dists) == size(B.lasers_dists)
        errprint = 0
        for y in 1:2048, x in 1:2048
            if !isapprox(A.lasers_dists[x,y], B.lasers_dists[x,y]; atol=0.05, nans=true)
                @warn "lasers_dists x $x y $y ($(A.lasers_dists[x,y]) != $(B.lasers_dists[x,y]))"
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

    if (2048,2048) == size(A.λmap) == size(B.λmap)
        errprint = 0
        for y in 1:2048, x in 1:2048
            if !isapprox(A.λmap[x,y], B.λmap[x,y]; atol=0.0001, nans=true)
                @warn "λmap x $x y $y ($(A.λmap[x,y]) != $(B.λmap[x,y]))"
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
