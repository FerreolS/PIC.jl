using Test,
      TwoDimensional
      
using Base: Fix1, Fix2, splat
using StatsBase: median, mad

function cmpt_goods(lenslettab)
    [ i for i in eachindex(lenslettab)
        if isassigned(lenslettab, i) && all(!isnan, lenslettab[i].dmodel.cx)
                                     && all(!isnan, lenslettab[i].dmodel.cy)
                                     && [1, 0, 0] != lenslettab[i].profile.cx
                                     && [1, 0, 0] != lenslettab[i].profile.cλ   ]
end

function cmpt_clm(lenslettab)
    clm = falses(length(lenslettab))
    goods = cmpt_goods(lenslettab)
    clm[goods] .= true
    clm
end

function store_output(
    filepath,
    (lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, clm))

    fitsfile = FitsFile(filepath, "w!")

    nblenses = length(lenslettab)

    goods = findall(clm)

    bboxs = fill(-1, 4, nblenses)
    dmodelorder = lenslettab[goods[1]].dmodel.order
    dmodelλ0    = lenslettab[goods[1]].dmodel.λ0
    profileorder = lenslettab[goods[1]].profile.order
    profileλ0    = lenslettab[goods[1]].profile.λ0
    dmodelcx = fill(NaN, dmodelorder + 1, nblenses)
    dmodelcy = fill(NaN, dmodelorder + 1, nblenses)
    profilecx = fill(NaN, profileorder + 1, nblenses)
    profilecl = fill(NaN, profileorder + 1, nblenses)

    for good in goods
        lens = lenslettab[good]
        bboxs[:,good] .= [lens.bbox...]
        dmodelcx[:,good] .= lens.dmodel.cx
        dmodelcy[:,good] .= lens.dmodel.cy
        profilecx[:,good] .= lens.profile.cx
        profilecl[:,good] .= lens.profile.cλ
    end

    write(fitsfile, FitsHeader("HDUNAME" => "computedlensmap"), clm)
    write(fitsfile, FitsHeader("EXTNAME" => "bboxs"), bboxs)
    write(fitsfile,
        FitsHeader("EXTNAME" => "dmodelcx", "ORDER" => dmodelorder, "L0" => dmodelλ0), dmodelcx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "dmodelcy", "ORDER" => dmodelorder, "L0" => dmodelλ0), dmodelcy)
    write(fitsfile,
        FitsHeader("EXTNAME" => "profilecx", "ORDER" => profileorder, "L0" => profileλ0), profilecx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "profilecl", "ORDER" => profileorder, "L0" => profileλ0), profilecl)
    write(fitsfile, FitsHeader("EXTNAME" => "laserAmp"), laserAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "lampAmp"), lampAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "laserFWHM"), laserfwhm)
    write(fitsfile, FitsHeader("EXTNAME" => "laserDist"), laserdist)
    write(fitsfile, FitsHeader("EXTNAME" => "lambdaMap"), λMap)
    close(fitsfile)
end

function store_output(
    filepath,
    (λ0, wavelamplkltab, wavelamps_fits_cx, wavelamps_fits_cy, specposlkltab, laserAmplitude,
     lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap))

    fitsfile = FitsFile(filepath, "w!")

    nblenses = length(wavelamplkltab)

    goods = findall(computedlensmap)

    bboxs = fill(-1, 4, nblenses)
    dmodelorder = wavelamplkltab[goods[1]].order
    dmodelλ0    = λ0
    profileorder = specposlkltab[goods[1]].pmodel.order
    profileλ0    = specposlkltab[goods[1]].pmodel.λ0
    dmodelcx = fill(NaN, dmodelorder + 1, nblenses)
    dmodelcy = fill(NaN, dmodelorder + 1, nblenses)
    profilecx = fill(NaN, profileorder + 1, nblenses)
    profilecl = fill(NaN, profileorder + 1, nblenses)

    for good in goods
        wavelamplkl = wavelamplkltab[good]
        profilemodel = specposlkltab[good].pmodel
        bboxs[:,good] .= [wavelamplkl.box...]
        dmodelcx[:,good] .= wavelamps_fits_cx[:,good]
        dmodelcy[:,good] .= wavelamps_fits_cy[:,good]
        profilecx[:,good] .= profilemodel.cx
        profilecl[:,good] .= profilemodel.cλ
    end

    write(fitsfile, FitsHeader("HDUNAME" => "computedlensmap"), computedlensmap)
    write(fitsfile, FitsHeader("EXTNAME" => "bboxs"), bboxs)
    write(fitsfile,
        FitsHeader("EXTNAME" => "dmodelcx", "ORDER" => dmodelorder, "L0" => dmodelλ0), dmodelcx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "dmodelcy", "ORDER" => dmodelorder, "L0" => dmodelλ0), dmodelcy)
    write(fitsfile,
        FitsHeader("EXTNAME" => "profilecx", "ORDER" => profileorder, "L0" => profileλ0), profilecx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "profilecl", "ORDER" => profileorder, "L0" => profileλ0), profilecl)
    write(fitsfile, FitsHeader("EXTNAME" => "laserAmp"), laserAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "lampAmp"), lampAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "laserFWHM"), laserfwhm)
    write(fitsfile, FitsHeader("EXTNAME" => "laserDist"), laserdist)
    write(fitsfile, FitsHeader("EXTNAME" => "lambdaMap"), λMap)
    close(fitsfile)
end

function equalOrNans(x, y)
    size(x) == size(y) || return false
    for i in eachindex(x)
        x[i] == y[i] || (isnan(x[i]) && isnan(y[i])) || return false
    end
    true
end

function approxOrNans(x, y)
    size(x) == size(y) || return false
    for i in eachindex(x)
        x[i] ≈ y[i] || (isnan(x[i]) && isnan(y[i])) || return false
    end
    true
end

function approxOrNans(x, y, atol)
    size(x) == size(y) || return false
    for i in eachindex(x)
        isapprox(x[i], y[i]; atol) || (isnan(x[i]) && isnan(y[i])) || return false
    end
    true
end

function compare_results_strict(res_1, res_2)
    
    (clm_1, bboxs_1,
     dmodelorder_1, dmodelλ0_1, dmodelcx_1, dmodelcy_1,
     profileorder_1, profileλ0_1, profilecx_1, profilecl_1,
     laserAmp_1, lampAmp_1, laserfwhm_1,
     laserdist_1, λMap_1) = res_1
     
    (clm_2, bboxs_2,
     dmodelorder_2, dmodelλ0_2, dmodelcx_2, dmodelcy_2,
     profileorder_2, profileλ0_2, profilecx_2, profilecl_2,
     laserAmp_2, lampAmp_2, laserfwhm_2,
     laserdist_2, λMap_2) = res_2
     
    if clm_1 != clm_2
        @warn "different computedlensmap"
        @warn "next tests will only be on lenses computed by both sides"
    end
    
    clm = (&).(clm_1, clm_2)
     
    if bboxs_1[:,clm] != bboxs_2[:,clm]
        @warn "different bboxs"
    end
     
    dmodelorder_1 != dmodelorder_2 &&
        @warn "different dmodelorder: $dmodelorder_1 != $dmodelorder_2"
    dmodelλ0_1 ≈ dmodelλ0_2 || @warn "different dmodelλ0: $dmodelλ0_1 != $dmodelλ0_2"
     
    if !equalOrNans(dmodelcx_1[:,clm], dmodelcx_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(dmodelcx_1[:,i], dmodelcx_2[:,i])
                bad = i
                break
            end
        end
        @warn "different dmodelcx, example lens $bad: $(dmodelcx_1[:,bad]) != $(dmodelcx_2[:,bad])"
    end
     
    if !equalOrNans(dmodelcy_1[:,clm], dmodelcy_2[:,clm])
        @warn "different dmodelcy"
    end
     
    profileorder_1 != profileorder_2 &&
        @warn "different profileorder: $profileorder_1 != $profileorder_2"
    profileλ0_1 ≈ profileλ0_2 || @warn "different profileλ0: $profileλ0_1 != $profileλ0_2"
     
    if !equalOrNans(profilecx_1[:,clm], profilecx_2[:,clm])
        @warn "different profilecx"
    end
     
    if !equalOrNans(profilecl_1[:,clm], profilecl_2[:,clm])
        @warn "different profilecl"
    end
     
    if !equalOrNans(laserAmp_1[:,clm], laserAmp_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(laserAmp_1[:,i], laserAmp_2[:,i])
                bad = i
                break
            end
        end
        @warn "different laserAmp, example lens $bad: $(laserAmp_1[:,bad]) != $(laserAmp_2[:,bad])"

    end
    
    if !equalOrNans(lampAmp_1[:,clm], lampAmp_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !equalOrNans(lampAmp_1[:,i], lampAmp_2[:,i])
                bad = i
                break
            end
        end
        @warn "different lampAmp, example lens $bad: $(lampAmp_1[:,bad]) != $(lampAmp_2[:,bad])"
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

function compare_results_approx(res_1, res_2)
    
    (clm_1, bboxs_1,
     dmodelorder_1, dmodelλ0_1, dmodelcx_1, dmodelcy_1,
     profileorder_1, profileλ0_1, profilecx_1, profilecl_1,
     laserAmp_1, lampAmp_1, laserfwhm_1,
     laserdist_1, λMap_1) = res_1
     
    (clm_2, bboxs_2,
     dmodelorder_2, dmodelλ0_2, dmodelcx_2, dmodelcy_2,
     profileorder_2, profileλ0_2, profilecx_2, profilecl_2,
     laserAmp_2, lampAmp_2, laserfwhm_2,
     laserdist_2, λMap_2) = res_2
     
    if clm_1 != clm_2
        @warn "different computedlensmap"
        @warn "next tests will only be on lenses computed by both sides"
    end
    
    clm = (&).(clm_1, clm_2)
     
    if bboxs_1[:,clm] != bboxs_2[:,clm]
        @warn "different bboxs"
    end
     
    dmodelorder_1 != dmodelorder_2 &&
        @warn "different dmodelorder: $dmodelorder_1 != $dmodelorder_2"
    dmodelλ0_1 ≈ dmodelλ0_2 || @warn "different dmodelλ0: $dmodelλ0_1 != $dmodelλ0_2"
     
    if !approxOrNans(dmodelcx_1[:,clm], dmodelcx_2[:,clm])
        @warn "different dmodelcx"
    end
     
    if !approxOrNans(dmodelcy_1[:,clm], dmodelcy_2[:,clm])
        @warn "different dmodelcy"
    end
     
    profileorder_1 != profileorder_2 &&
        @warn "different profileorder: $profileorder_1 != $profileorder_2"
    profileλ0_1 ≈ profileλ0_2 || @warn "different profileλ0: $profileλ0_1 != $profileλ0_2"
     
    if !approxOrNans(profilecx_1[:,clm], profilecx_2[:,clm], 0.01)
        bad = 0
        for i in findall(clm)
            if !approxOrNans(profilecx_1[:,i], profilecx_2[:,i], 0.1)
                bad = i
                break
            end
        end
        @warn string("different profilecx, example lens $bad: ",
                     profilecx_1[:,bad], " != ", profilecx_2[:,bad])
    end
     
    if !approxOrNans(profilecl_1[:,clm], profilecl_2[:,clm])
        @warn "different profilecl"
    end
     
    if !approxOrNans(laserAmp_1[:,clm], laserAmp_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !approxOrNans(laserAmp_1[:,i], laserAmp_2[:,i])
                bad = i
                break
            end
        end
        @warn "different laserAmp, example lens $bad: $(laserAmp_1[:,bad]) != $(laserAmp_2[:,bad])"

    end
    
    if !approxOrNans(lampAmp_1[:,clm], lampAmp_2[:,clm])
        bad = 0
        for i in findall(clm)
            if !approxOrNans(lampAmp_1[:,i], lampAmp_2[:,i])
                bad = i
                break
            end
        end
        @warn "different lampAmp, example lens $bad: $(lampAmp_1[:,bad]) != $(lampAmp_2[:,bad])"
    end
    
    if !approxOrNans(laserfwhm_1[:,clm], laserfwhm_2[:,clm])
        @warn "different laserfwhm"
    end
    
    if !approxOrNans(laserdist_1, laserdist_2)
        bad = 0
        for i in findall(clm)
            if !approxOrNans(laserdist_1[:,i], laserdist_2[:,i])
                bad = i
                break
            end
        end
        @warn string("different laserdist, example lens $bad: ",
                     laserdist_1[:,bad], " != ",  laserdist_2[:,bad])
    end
    
    if !approxOrNans(λMap_1, λMap_2)
        @warn "different λMap"
    end
end

function get_result(
    (λ0, wavelamplkltab, wavelamps_fits_cx, wavelamps_fits_cy, specposlkltab, laserAmplitude,
     lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
)

    goods = findall(computedlensmap)
    
    nblenses = length(computedlensmap)
    
    bboxs = fill(-1, 4, nblenses)

    dmodelorder = wavelamplkltab[goods[1]].order
    dmodelλ0    = λ0
    dmodelcx = fill(NaN, dmodelorder + 1, nblenses)
    dmodelcy = fill(NaN, dmodelorder + 1, nblenses)

    profileorder = specposlkltab[goods[1]].pmodel.order
    profileλ0    = specposlkltab[goods[1]].pmodel.λ0
    profilecx = fill(NaN, profileorder + 1, nblenses)
    profilecl = fill(NaN, profileorder + 1, nblenses)

    for good in goods
        wavelamplkl = wavelamplkltab[good]
        bboxs[:,good] .= [wavelamplkl.box...]
        dmodelcx[:,good] .= wavelamps_fits_cx[:,good]
        dmodelcy[:,good] .= wavelamps_fits_cy[:,good]
        specposlkl = specposlkltab[good]
        profilecx[:,good] .= specposlkl.pmodel.cx
        profilecl[:,good] .= specposlkl.pmodel.cλ
    end
    
    (; computedlensmap, bboxs,
       dmodelorder, dmodelλ0, dmodelcx, dmodelcy,
       profileorder, profileλ0, profilecx, profilecl,
       laserAmplitude, lampAmplitude, laserfwhm,
       laserdist, λMap)
end

function get_result(filepath::String)

    FitsFile(filepath) do fitsfile

        hduclm = fitsfile[1]
        hdubboxs = fitsfile["bboxs"]
        hdudmodelcx = fitsfile["dmodelcx"]
        hdudmodelcy = fitsfile["dmodelcy"]
        hduprofilecx = fitsfile["profilecx"]
        hduprofilecl = fitsfile["profilecl"]
        hdulaseramp = fitsfile["laserAmp"]
        hdulampamp = fitsfile["lampAmp"]
        hdulaserfwhm = fitsfile["laserfwhm"]
        hdulaserdist = fitsfile["laserDist"]
        hduλMap = fitsfile["lambdaMap"]
        
        computedlensmap = read(Array{Bool}, hduclm)
        bboxs = read(hdubboxs)
        
        dmodelorder = hdudmodelcx["ORDER"].integer
        dmodelλ0 = hdudmodelcx["L0"].float
        dmodelcx = read(hdudmodelcx)
        dmodelcy = read(hdudmodelcy)
    
        profileorder = hduprofilecx["ORDER"].integer
        profileλ0 = hduprofilecx["L0"].float
        profilecx = read(hduprofilecx)
        profilecl = read(hduprofilecl)
    
        laserAmplitude = read(hdulaseramp)
        laserfwhm = read(hdulaserfwhm)
        lampAmplitude = read(hdulampamp)

        laserdist = read(hdulaserdist)
        λMap = read(hduλMap)

        (; computedlensmap, bboxs,
           dmodelorder, dmodelλ0, dmodelcx, dmodelcy,
           profileorder, profileλ0, profilecx, profilecl,
           laserAmplitude, lampAmplitude, laserfwhm,
           laserdist, λMap)
   end
end

#end

function keep_numbers(data)
    data |> Fix1(filter, !isnan) |> Fix1(filter, !isinf)
end

function keep_90pc(data)
    (lower, upper) = quantile(data, [0.05, 0.95])
    clamp.(data, lower, upper)
end

function get_mean_std(data)
    data = keep_numbers(data)
    data = keep_90pc(data)
    return mean(data), std(data)
end

function get_result_summary(
    (; computedlensmap, bboxs,
       dmodelorder, dmodelλ0, dmodelcx, dmodelcy,
       profileorder, profileλ0, profilecx, profilecl,
       laserAmplitude, lampAmplitude, laserfwhm,
       laserdist, λMap)
)

    goods = findall(computedlensmap)
    
    ratiogoods = length(goods) / length(computedlensmap)
    
    # disp model
    
    quantilesDmodelcx = fill(NaN, (21, dmodelorder + 1))
    for a in 1:dmodelorder+1
        data = dmodelcx[a,goods]
        quantilesDmodelcx[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    quantilesDmodelcy = fill(NaN, (21, dmodelorder+1))
    for a in 1:dmodelorder+1
        data = dmodelcy[a,goods]
        quantilesDmodelcy[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    # profile model
    
    quantilesProfilecx = fill(NaN, (21, profileorder+1))
    for a in 1:profileorder+1
        data = profilecx[a,goods]
        quantilesProfilecx[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
    end
    
    quantilesProfilecl = fill(NaN, (21, profileorder+1))
    for a in 1:profileorder+1
        data = profilecl[a,goods]
        quantilesProfilecl[:,a] .= quantile(keep_numbers(data), (0.0 : 0.05 : 1.0))
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
        cs = [ laserdist[BoundingBox(bboxs[:,i]...)][c,:] for i in goods ]
        m, s = keep_numbers(reduce(vcat, cs)) |> x -> (median(x), mad(x))
        medianMadDistColumn[c,1] = m
        medianMadDistColumn[c,2] = s
    end
    
    # laser dist
    medianMadλMapLine = fill(NaN, 40, 2)
    for l in 1:40
        ls = [ λMap[BoundingBox(bboxs[:,i]...)][:,l] for i in goods ]
        m, s = keep_numbers(reduce(vcat, ls)) |> x -> (median(x), mad(x))
        medianMadλMapLine[l,1] = m
        medianMadλMapLine[l,2] = s
    end

    (; ratiogoods,
       dmodelorder, dmodelλ0, quantilesDmodelcx, quantilesDmodelcy,
       profileorder, profileλ0, quantilesProfilecx, quantilesProfilecl,
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

