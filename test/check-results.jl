using Test,
      TwoDimensional
      
using Base: Fix1, Fix2, splat
using StatsBase: median, mad

function get_result(
    (lensboxs, λ0,
     wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
     specpos_order, specpos_fits_cx, specpos_fits_cλ,
     laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
)
    
    lensboxs = reduce(hcat, map(b -> [b...], lensboxs))
    
    (; lensboxs, λ0,
       wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
       specpos_order, specpos_fits_cx, specpos_fits_cλ,
       laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
end

function get_result(filepath::String)

    FitsFile(filepath) do fitsfile
        
        computedlensmap = read(Array{Bool}, fitsfile["computedlensmap"])        

        lensboxs = read(fitsfile["lensboxs"])
        
        wavelamp_order = fitsfile["wavelamps_fits_cx"]["ORDER"].integer
        λ0 = fitsfile["wavelamps_fits_cx"]["L0"].float
        wavelamps_fits_cx = read(fitsfile["wavelamps_fits_cx"])
        wavelamps_fits_cy = read(fitsfile["wavelamps_fits_cy"])
    
        specpos_order = fitsfile["specpos_fits_cx"]["ORDER"].integer
        specpos_fits_cx = read(fitsfile["specpos_fits_cx"])
        specpos_fits_cλ = read(fitsfile["specpos_fits_cl"])
    
        laserAmplitude = read(fitsfile["laserAmp"])
        lampAmplitude = read(fitsfile["lampAmp"])
        laserfwhm = read(fitsfile["laserfwhm"])
        laserdist = read(fitsfile["laserDist"])
        λMap = read(fitsfile["lambdaMap"])

        (; lensboxs, λ0,
           wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
           specpos_order, specpos_fits_cx, specpos_fits_cλ,
           laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
   end
end

function store_result(
    filepath,
    (; lensboxs, λ0,
     wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
     specpos_order, specpos_fits_cx, specpos_fits_cλ,
     laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
)

    fitsfile = FitsFile(filepath, "w!")

    write(fitsfile, FitsHeader("HDUNAME" => "computedlensmap"), computedlensmap)
    
    write(fitsfile, FitsHeader("EXTNAME" => "lensboxs"), lensboxs)
    
    write(fitsfile,
        FitsHeader("EXTNAME" => "wavelamps_fits_cx", "ORDER" => wavelamp_order, "L0" => λ0),
        wavelamps_fits_cx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "wavelamps_fits_cy", "ORDER" => wavelamp_order, "L0" => λ0),
        wavelamps_fits_cy)
    
    write(fitsfile,
        FitsHeader("EXTNAME" => "specpos_fits_cx", "ORDER" => specpos_order, "L0" => λ0),
        specpos_fits_cx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "specpos_fits_cl", "ORDER" => specpos_order, "L0" => λ0),
        specpos_fits_cλ)
        
    write(fitsfile, FitsHeader("EXTNAME" => "laserAmp"), laserAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "lampAmp"), lampAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "laserFWHM"), laserfwhm)
    write(fitsfile, FitsHeader("EXTNAME" => "laserDist"), laserdist)
    write(fitsfile, FitsHeader("EXTNAME" => "lambdaMap"), λMap)
    close(fitsfile)
end

function keep_numbers(data)
    data |> Fix1(filter, !isnan) |> Fix1(filter, !isinf)
end

function get_result_summary(
    (; lensboxs, λ0,
       wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
       specpos_order, specpos_fits_cx, specpos_fits_cλ,
       laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
)

    goods = findall(computedlensmap)
    
    nbgoods = length(goods)
    
    quantiles_wavelamps_fits_cx = fill(NaN, (21, wavelamp_order + 1))
    for a in 1:wavelamp_order+1
        data = wavelamps_fits_cx[a,goods]
        quantiles_wavelamps_fits_cx[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    quantiles_wavelamps_fits_cy = fill(NaN, (21, wavelamp_order + 1))
    for a in 1:wavelamp_order+1
        data = wavelamps_fits_cy[a,goods]
        quantiles_wavelamps_fits_cy[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    quantiles_specpos_fits_cx = fill(NaN, (21, specpos_order + 1))
    for a in 1:specpos_order+1
        data = specpos_fits_cx[a,goods]
        quantiles_specpos_fits_cx[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    quantiles_specpos_fits_cλ = fill(NaN, (21, specpos_order + 1))
    for a in 1:specpos_order+1
        data = specpos_fits_cλ[a,goods]
        quantiles_specpos_fits_cλ[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    # laser amp
    medianLaserAmp1, madLaserAmp1 = keep_numbers(laserAmplitude[1,goods]) |> x -> (median(x), mad(x))
    medianLaserAmp2, madLaserAmp2 = keep_numbers(laserAmplitude[2,goods]) |> x -> (median(x), mad(x))
    medianLaserAmp3, madLaserAmp3 = keep_numbers(laserAmplitude[3,goods]) |> x -> (median(x), mad(x))
    medianMadLasersAmps = [ medianLaserAmp1 madLaserAmp1 ;
                            medianLaserAmp2 madLaserAmp2 ;
                            medianLaserAmp3 madLaserAmp3 ]
    
    # laser FWHM
    medianLaserFWHM1, madLaserFWHM1 = keep_numbers(laserfwhm[1,goods]) |> x -> (median(x), mad(x))
    medianLaserFWHM2, madLaserFWHM2 = keep_numbers(laserfwhm[2,goods]) |> x -> (median(x), mad(x))
    medianLaserFWHM3, madLaserFWHM3 = keep_numbers(laserfwhm[3,goods]) |> x -> (median(x), mad(x))
    medianMadLasersFWHMs = [ medianLaserFWHM1 madLaserFWHM1 ;
                             medianLaserFWHM2 madLaserFWHM2 ;
                             medianLaserFWHM3 madLaserFWHM3 ]

    # lamp Amp
    
    medianMadLampAmp = fill(NaN, 40, 2)
    for i in 1:40
        m, s = keep_numbers(lampAmplitude[1+i, goods]) |> x -> (median(x), mad(x))
        medianMadLampAmp[i,1] = m
        medianMadLampAmp[i,2] = s
    end

    # laser dist
    medianMadDistColumn = fill(NaN, 5, 2)
    for c in 1:5
        cs = [ laserdist[BoundingBox(lensboxs[:,i]...)][c,:] for i in goods ]
        m, s = keep_numbers(reduce(vcat, cs)) |> x -> (median(x), mad(x))
        medianMadDistColumn[c,1] = m
        medianMadDistColumn[c,2] = s
    end
    
    # laser dist
    medianMadλMapLine = fill(NaN, 40, 2)
    for l in 1:40
        ls = [ λMap[BoundingBox(lensboxs[:,i]...)][:,l] for i in goods ]
        m, s = keep_numbers(reduce(vcat, ls)) |> x -> (median(x), mad(x))
        medianMadλMapLine[l,1] = m
        medianMadλMapLine[l,2] = s
    end

    (; nbgoods, λ0,
       wavelamp_order, quantiles_wavelamps_fits_cx, quantiles_wavelamps_fits_cy,
       specpos_order, quantiles_specpos_fits_cx, quantiles_specpos_fits_cλ,
       medianMadLasersAmps,
       medianMadLasersFWHMs,
       medianMadLampAmp,
       medianMadDistColumn, medianMadλMapLine)
end

function display_summary(summary, io=stdout)
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

function compare_results(res_1, res_2)
    
    (lensboxs_1, λ0_1,
     wavelamp_order_1, wavelamps_fits_cx_1, wavelamps_fits_cy_1,
     specpos_order_1, specpos_fits_cx_1, specpos_fits_cλ_1,
     laserAmplitude_1, lampAmplitude_1, laserfwhm_1, laserdist_1, λMap_1,
     computedlensmap_1) = res_1

    (lensboxs_2, λ0_2,
     wavelamp_order_2, wavelamps_fits_cx_2, wavelamps_fits_cy_2,
     specpos_order_2, specpos_fits_cx_2, specpos_fits_cλ_2,
     laserAmplitude_2, lampAmplitude_2, laserfwhm_2, laserdist_2, λMap_2,
     computedlensmap_2) = res_2
     
    if computedlensmap_1 != computedlensmap_2
        @warn "different computedlensmap"
        @warn "next tests will only be on lenses computed by both sides"
    end
    
    clm = (&).(computedlensmap_1, computedlensmap_2)
     
    if lensboxs_1[:,clm] != lensboxs_2[:,clm]
        @warn "different lensbox"
    end
     
    wavelamp_order_1 != wavelamp_order_2 &&
        @warn "different wavelamp_order: $wavelamp_order_1 != $wavelamp_order_2"
    λ0_1 ≈ λ0_2 || @warn "different λ0: $λ0_1 != $λ0_2"
     
    if !equalOrNans(wavelamps_fits_cx_1[:,clm], wavelamps_fits_cx_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(wavelamps_fits_cx_1[:,i], wavelamps_fits_cx_2[:,i])
                bad = i
                break
            end
        end
        @warn "different wavelamps_fits_cx, example lens $bad: $(wavelamps_fits_cx_1[:,bad]) != $(wavelamps_fits_cx_2[:,bad])"
    end
     
    if !equalOrNans(wavelamps_fits_cy_1[:,clm], wavelamps_fits_cy_2[:,clm])
        @warn "different wavelamps_fits_cy"
    end
     
    specpos_order_1 != specpos_order_2 &&
        @warn "different specpos_order: $specpos_order_1 != $specpos_order_2"
     
    if !equalOrNans(specpos_fits_cx_1[:,clm], specpos_fits_cx_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(specpos_fits_cx_1[:,i], specpos_fits_cx_2[:,i])
                bad = i
                break
            end
        end
        @warn "different specpos_fits_cx, example lens $bad: $(specpos_fits_cx_1[:,bad]) != $(specpos_fits_cx_2[:,bad])"
    end
     
    if !equalOrNans(specpos_fits_cλ_1[:,clm], specpos_fits_cλ_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(specpos_fits_cλ_1[:,i], specpos_fits_cλ_2[:,i])
                bad = i
                break
            end
        end
        @warn "different specpos_fits_cλ, example lens $bad: $(specpos_fits_cλ_1[:,bad]) != $(specpos_fits_cλ_2[:,bad])"
    end
     
    if !equalOrNans(laserAmplitude_1[:,clm], laserAmplitude_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(laserAmplitude_1[:,i], laserAmplitude_2[:,i])
                bad = i
                break
            end
        end
        @warn "different laserAmplitude, example lens $bad: $(laserAmplitude_1[:,bad]) != $(laserAmplitude_2[:,bad])"

    end
    
    if !equalOrNans(lampAmplitude_1[:,clm], lampAmplitude_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(lampAmplitude_1[:,i], lampAmplitude_2[:,i])
                bad = i
                break
            end
        end
        @warn "different lampAmplitude, example lens $bad: "
        show(stdout, MIME"text/plain"(), lampAmplitude_1[:,bad])
        println("!=")
        show(stdout, MIME"text/plain"(), lampAmplitude_2[:,bad])
    end
    
    if !equalOrNans(laserfwhm_1[:,clm], laserfwhm_2[:,clm])
        @warn "different laserfwhm"
    end
    
    if !equalOrNans(laserdist_1, laserdist_2)
        bad = 0
        for i in findall(clm)
            if !equalOrNans(laserdist_1[:,i], laserdist_2[:,i])
                bad = i
                break
            end
        end
        @warn string("different laserdist, example lens $bad: ",
                     laserdist_1[:,bad], " != ",  laserdist_2[:,bad])
    end
    
    if !equalOrNans(λMap_1, λMap_2)
        @warn "different λMap"
    end
end


#end


