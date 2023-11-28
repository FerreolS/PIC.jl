using Test,
      TwoDimensional
      
using Base: Fix1, Fix2, splat
using StatsBase: median, mad
using PIC: FitResult

function get_result(filepath::String)

    FitsFile(filepath) do fitsfile
        
        computed_lenses_map = read(Array{Bool}, fitsfile["computed_lenses_map"])

        lenses_boxes = read(fitsfile["lenses_boxes"])
        
        lenses_boxes = map(Base.splat(BoundingBox{Int}), eachcol(lenses_boxes))
        
        wavelamp_order            = fitsfile["wavelamps_fits_cx"]["ORDER"].integer
        λ0                        = fitsfile["wavelamps_fits_cx"]["L0"].float
        wavelamps_fits_cx         = read(fitsfile["wavelamps_fits_cx"])
        wavelamps_fits_cy         = read(fitsfile["wavelamps_fits_cy"])
        wavelamps_fits_fwhm       = read(fitsfile["wavelamps_fits_fwhm"])
        wavelamps_fits_background = read(fitsfile["wavelamps_fits_background"])
        wavelamps_fits_amp        = read(fitsfile["wavelamps_fits_amp"])
        wavelamps_centers_dists   = read(fitsfile["wavelamps_centers_dists"])
        wavelamps_λvals           = read(fitsfile["wavelamps_lvals"])
    
        specpos_order           = fitsfile["specpos_fits_cx"]["ORDER"].integer
        specpos_fits_cx         = read(fitsfile["specpos_fits_cx"])
        specpos_fits_cλ         = read(fitsfile["specpos_fits_cl"])
        specpos_fits_background = read(fitsfile["specpos_fits_background"])
        sepcpos_fits_amps       = read(fitsfile["sepcpos_fits_amps"])
           
        FitResult(lenses_boxes, λ0,
         wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
         wavelamps_fits_fwhm, wavelamps_fits_background, wavelamps_fits_amp,
         wavelamps_centers_dists, wavelamps_λvals,
         specpos_order, specpos_fits_cx, specpos_fits_cλ,
         specpos_fits_background, sepcpos_fits_amps,
         computed_lenses_map)
    end
end

function store_result(filepath, fitresult)

    (lenses_boxes, λ0,
        wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
        wavelamps_fits_fwhm, wavelamps_fits_background, wavelamps_fits_amp,
        wavelamps_centers_dists, wavelamps_λvals,
        specpos_order, specpos_fits_cx, specpos_fits_cλ,
        specpos_fits_background, sepcpos_fits_amps,
        computed_lenses_map) = map(p -> getproperty(fitresult, p), propertynames(fitresult))

    lenses_boxes = reduce(hcat, map(b -> [b...], lenses_boxes))

    fitsfile = FitsFile(filepath, "w!")

    write(fitsfile, FitsHeader("HDUNAME" => "computed_lenses_map"), computed_lenses_map)
    
    write(fitsfile, FitsHeader("EXTNAME" => "lenses_boxes"), lenses_boxes)
    
    write(fitsfile,
        FitsHeader("EXTNAME" => "wavelamps_fits_cx", "ORDER" => wavelamp_order, "L0" => λ0),
        wavelamps_fits_cx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "wavelamps_fits_cy", "ORDER" => wavelamp_order, "L0" => λ0),
        wavelamps_fits_cy)
    write(fitsfile, FitsHeader("EXTNAME" => "wavelamps_fits_fwhm"),     wavelamps_fits_fwhm)
    write(fitsfile, FitsHeader("EXTNAME" => "wavelamps_fits_background"),wavelamps_fits_background)
    write(fitsfile, FitsHeader("EXTNAME" => "wavelamps_fits_amp"),      wavelamps_fits_amp)
    write(fitsfile, FitsHeader("EXTNAME" => "wavelamps_centers_dists"), wavelamps_centers_dists)
    write(fitsfile, FitsHeader("EXTNAME" => "wavelamps_lvals"),         wavelamps_λvals)
    
    write(fitsfile,
        FitsHeader("EXTNAME" => "specpos_fits_cx", "ORDER" => specpos_order, "L0" => λ0),
        specpos_fits_cx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "specpos_fits_cl", "ORDER" => specpos_order, "L0" => λ0),
        specpos_fits_cλ)
        
    write(fitsfile, FitsHeader("EXTNAME" => "specpos_fits_background"), specpos_fits_background)
    write(fitsfile, FitsHeader("EXTNAME" => "sepcpos_fits_amps"),       sepcpos_fits_amps)
    
    close(fitsfile)
end

function keep_numbers(data)
    data |> Fix1(filter, !isnan) |> Fix1(filter, !isinf)
end

function get_summary(fitresult)

    (lenses_boxes, λ0,
     wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
     wavelamps_fits_fwhm, wavelamps_fits_background, wavelamps_fits_amp,
     wavelamps_centers_dists, wavelamps_λvals,
     specpos_order, specpos_fits_cx, specpos_fits_cλ,
     specpos_fits_background, sepcpos_fits_amps,
     computed_lenses_map) = map(p -> getproperty(fitresult, p), propertynames(fitresult))

    goods = findall(computed_lenses_map)
    
    nbgoods = length(goods)
    
    # wavelamps_fits_cx
    quantiles_wavelamps_fits_cx = fill(NaN64, (21, wavelamp_order + 1))
    for a in 1:wavelamp_order+1
        data = wavelamps_fits_cx[a,goods]
        quantiles_wavelamps_fits_cx[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    # wavelamps_fits_cy
    quantiles_wavelamps_fits_cy = fill(NaN64, (21, wavelamp_order + 1))
    for a in 1:wavelamp_order+1
        data = wavelamps_fits_cy[a,goods]
        quantiles_wavelamps_fits_cy[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    # wavelamps_fits_fwhm
    medianMad_wavelamps_fits_fwhm = fill(NaN64, size(wavelamps_fits_fwhm,1), 2)
    for i in 1:size(wavelamps_fits_fwhm,1)
        medianMad_wavelamps_fits_fwhm[i,:] .=
            keep_numbers(wavelamps_fits_fwhm[i,goods]) |> x -> (median(x), mad(x))
    end
    
    # wavelamps_fits_background
    medianMad_wavelamps_fits_background = fill(NaN64, 1, 2)
    m, s = keep_numbers(wavelamps_fits_background[goods]) |> x -> (median(x), mad(x))
    medianMad_wavelamps_fits_background[1,1] = m
    medianMad_wavelamps_fits_background[1,2] = s
    
    # wavelamps_fits_amp
    medianMad_wavelamps_fits_amp = fill(NaN64, size(wavelamps_fits_amp,1), 2)
    for i in 1:size(wavelamps_fits_amp,1)
        medianMad_wavelamps_fits_amp[i,:] .=
            keep_numbers(wavelamps_fits_amp[i,goods]) |> x -> (median(x), mad(x))
    end

    # wavelamps_centers_dists
    medianMad_wavelamps_centers_dists = fill(NaN64, 5, 2)
    for c in 1:5
        cs = [ wavelamps_centers_dists[lenses_boxes[i]][c,:] for i in goods ]
        m, s = keep_numbers(reduce(vcat, cs)) |> x -> (median(x), mad(x))
        medianMad_wavelamps_centers_dists[c,1] = m
        medianMad_wavelamps_centers_dists[c,2] = s
    end
    
    # wavelamps_λvals
    medianMad_wavelamps_λvals = fill(NaN64, 40, 2)
    for l in 1:40
        ls = [ wavelamps_λvals[lenses_boxes[i]][:,l] for i in goods ]
        m, s = keep_numbers(reduce(vcat, ls)) |> x -> (median(x), mad(x))
        medianMad_wavelamps_λvals[l,1] = m
        medianMad_wavelamps_λvals[l,2] = s
    end

    # specpos_fits_cx
    quantiles_specpos_fits_cx = fill(NaN64, (21, specpos_order + 1))
    for a in 1:specpos_order+1
        data = specpos_fits_cx[a,goods]
        quantiles_specpos_fits_cx[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    # specpos_fits_cλ
    quantiles_specpos_fits_cλ = fill(NaN64, (21, specpos_order + 1))
    for a in 1:specpos_order+1
        data = specpos_fits_cλ[a,goods]
        quantiles_specpos_fits_cλ[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end

    # specpos_fits_background
    medianMad_specpos_fits_background = fill(NaN64, 1, 2)
    m, s = keep_numbers(specpos_fits_background[goods]) |> x -> (median(x), mad(x))
    medianMad_specpos_fits_background[1,1] = m
    medianMad_specpos_fits_background[1,2] = s

    # specpos_fits_amps
    medianMad_sepcpos_fits_amps = fill(NaN64, 40, 2)
    for i in 1:40
        m, s = keep_numbers(sepcpos_fits_amps[i, goods]) |> x -> (median(x), mad(x))
        medianMad_sepcpos_fits_amps[i,1] = m
        medianMad_sepcpos_fits_amps[i,2] = s
    end


    (; nbgoods, λ0, wavelamp_order,
       quantiles_wavelamps_fits_cx, quantiles_wavelamps_fits_cy,
       medianMad_wavelamps_fits_fwhm,
       medianMad_wavelamps_fits_background, medianMad_wavelamps_fits_amp,
       medianMad_wavelamps_centers_dists, medianMad_wavelamps_λvals,
       quantiles_specpos_fits_cx, quantiles_specpos_fits_cλ,
       medianMad_specpos_fits_background, medianMad_sepcpos_fits_amps)
end

function display_summary(fitresult, io=stdout)
    summary = get_summary(fitresult)
    foreach(keys(summary)) do k
        print(io, "========== ")
        println(io, k)
        show(io, MIME"text/plain"(), summary[k])
        println(io)
    end
end


function equalOrNans(x, y)
    size(x) == size(y) || return false
    for i in eachindex(x)
        x[i] == y[i] || (isnan(x[i]) && isnan(y[i])) || return false
    end
    true
end

function compare_results(r1, r2)

    if r1.computed_lenses_map != r2.computed_lenses_map
        @warn "different computed_lenses_map"
        @show (sum(r1.computed_lenses_map) , sum(r2.computed_lenses_map))
        @warn "next tests will only be on lenses computed by both sides"
    end
    
    clm = (&).(r1.computed_lenses_map, r2.computed_lenses_map)
     
    if r1.lenses_boxes[clm] != r2.lenses_boxes[clm]
        @warn "different lensbox"
    end
     
    r1.wavelamp_order != r2.wavelamp_order &&
        @warn "different wavelamp_order: $(r1.wavelamp_order) != $(r2.wavelamp_order)"
    r1.λ0 ≈ r2.λ0 || @warn "different λ0: $(r1.λ0) != $(r2.λ0)"
     
    if !equalOrNans(r1.wavelamps_fits_cx[:,clm], r2.wavelamps_fits_cx[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.wavelamps_fits_cx[:,i], r2.wavelamps_fits_cx[:,i])
                bad = i
                break
            end
        end
        @warn "different wavelamps_fits_cx, example lens $bad: $(r1.wavelamps_fits_cx[:,bad]) != $(r2.wavelamps_fits_cx[:,bad])"
    end
     
    if !equalOrNans(r1.wavelamps_fits_cy[:,clm], r2.wavelamps_fits_cy[:,clm])
        @warn "different wavelamps_fits_cy"
    end
    
    if !equalOrNans(r1.wavelamps_fits_fwhm[:,clm], r2.wavelamps_fits_fwhm[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.wavelamps_fits_fwhm[:,i], r2.wavelamps_fits_fwhm[:,i])
                bad = i
                break
            end
        end
        @warn "different wavelamps_fits_fwhm, example lens $bad: "
        show(stdout, MIME"text/plain"(), r1.wavelamps_fits_fwhm[:,bad])
        println("!=")
        show(stdout, MIME"text/plain"(), r2.wavelamps_fits_fwhm[:,bad])
        println()
    end
    
    if !equalOrNans(r1.wavelamps_fits_background[clm], r2.wavelamps_fits_background[clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.wavelamps_fits_background[i], r2.wavelamps_fits_background[i])
                bad = i
                break
            end
        end
        @warn "different wavelamps_fits_background, example lens $bad: $(r1.wavelamps_fits_background[bad]) != $(r2.wavelamps_fits_background[bad])"
    end
    
    if !equalOrNans(r1.wavelamps_fits_amps[:,clm], r2.wavelamps_fits_amps[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.wavelamps_fits_amps[:,i], r2.wavelamps_fits_amps[:,i])
                bad = i
                break
            end
        end
        @warn "different wavelamps_fits_amps, example lens $bad: $(r1.wavelamps_fits_amps[:,bad]) != $(r2.wavelamps_fits_amps[:,bad])"
    end
    
    if !equalOrNans(r1.wavelamps_centers_dists, r2.wavelamps_centers_dists)
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.wavelamps_centers_dists[:,i], r2.wavelamps_centers_dists[:,i])
                bad = i
                break
            end
        end
        @warn string("different wavelamps_centers_dists, example lens $bad: ",
                     r1.wavelamps_centers_dists[:,bad], " != ",  r2.wavelamps_centers_dists[:,bad])
    end
    
    if !equalOrNans(r1.wavelamps_λvals, r2.wavelamps_λvals)
        @warn "different wavelamps_λvals"
    end
    
    r1.specpos_order != r2.specpos_order &&
        @warn "different specpos_order: $r1.specpos_order != $r2.specpos_order"
     
    if !equalOrNans(r1.specpos_fits_cx[:,clm], r2.specpos_fits_cx[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.specpos_fits_cx[:,i], r2.specpos_fits_cx[:,i])
                bad = i
                break
            end
        end
        @warn "different specpos_fits_cx, example lens $bad: $(r1.specpos_fits_cx[:,bad]) != $(r2.specpos_fits_cx[:,bad])"
    end
     
    if !equalOrNans(r1.specpos_fits_cλ[:,clm], r2.specpos_fits_cλ[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.specpos_fits_cλ[:,i], r2.specpos_fits_cλ[:,i])
                bad = i
                break
            end
        end
        @warn "different specpos_fits_cλ, example lens $bad: $(r1.specpos_fits_cλ[:,bad]) != $(r2.specpos_fits_cλ[:,bad])"
    end
    
    if !equalOrNans(r1.specpos_fits_background[clm], r2.specpos_fits_background[clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.specpos_fits_background[i], r2.specpos_fits_background[i])
                bad = i
                break
            end
        end
        @warn "different specpos_fits_background, example lens $bad: $(r1.specpos_fits_background[bad]) != $(r2.specpos_fits_background[bad])"
    end
    
    if !equalOrNans(r1.sepcpos_fits_amps[:,clm], r2.sepcpos_fits_amps[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(r1.sepcpos_fits_amps[:,i], r2.sepcpos_fits_amps[:,i])
                bad = i
                break
            end
        end
        @warn "different sepcpos_fits_amps, example lens $bad: "
        show(stdout, MIME"text/plain"(), r1.sepcpos_fits_amps[:,bad])
        println("!=")
        show(stdout, MIME"text/plain"(), r2.sepcpos_fits_amps[:,bad])
        println()
    end
end


#end


