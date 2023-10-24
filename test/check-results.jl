using Test

function cmpt_goods(lenslettab)
    [ i for i in eachindex(lenslettab)
        if isassigned(lenslettab, i) && all(!isnan, lenslettab[i].dmodel.cx)
                                     && all(!isnan, lenslettab[i].dmodel.cy)
                                     && [1, 0, 0] != lenslettab[i].profile.cy
                                     && [1, 0, 0] != lenslettab[i].profile.cλ   ]
end

function cmpt_vlm(lenslettab)
    vlm = zeros(length(lenslettab))
    goods = cmpt_goods(lenslettab)
    vlm[goods] .= 1
    vlm
end

function store_result(
    filepath,
    (lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap))

    fitsfile = FitsFile(filepath, "w!")

    nblenses = length(lenslettab)

    goods = cmpt_goods(lenslettab)

    bboxs = fill(-1, 4, nblenses)
    dmodelorder = lenslettab[goods[1]].dmodel.order
    dmodelλ0    = lenslettab[goods[1]].dmodel.λ0
    profileorder = lenslettab[goods[1]].profile.order
    profileλ0    = lenslettab[goods[1]].profile.λ0
    dmodelcx = fill(NaN, dmodelorder + 1, nblenses)
    dmodelcy = fill(NaN, dmodelorder + 1, nblenses)
    profilecy = fill(NaN, profileorder + 1, nblenses)
    profilecl = fill(NaN, profileorder + 1, nblenses)

    for good in goods
        lens = lenslettab[good]
        bboxs[:,good] .= [lens.bbox...]
        dmodelcx[:,good] .= lens.dmodel.cx
        dmodelcy[:,good] .= lens.dmodel.cy
        profilecy[:,good] .= lens.profile.cy
        profilecl[:,good] .= lens.profile.cλ
    end

    vlm = cmpt_vlm(lenslettab)

    write(fitsfile, FitsHeader("HDUNAME" => "validlensesmap"), vlm)
    write(fitsfile, FitsHeader("EXTNAME" => "bboxs"), bboxs)
    write(fitsfile,
        FitsHeader("EXTNAME" => "dmodelcx", "ORDER" => dmodelorder, "L0" => dmodelλ0), dmodelcx)
    write(fitsfile,
        FitsHeader("EXTNAME" => "dmodelcy", "ORDER" => dmodelorder, "L0" => dmodelλ0), dmodelcy)
    write(fitsfile,
        FitsHeader("EXTNAME" => "profilecy", "ORDER" => profileorder, "L0" => profileλ0), profilecy)
    write(fitsfile,
        FitsHeader("EXTNAME" => "profilecl", "ORDER" => profileorder, "L0" => profileλ0), profilecl)
    write(fitsfile, FitsHeader("EXTNAME" => "laserAmp"), laserAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "lampAmp"), lampAmplitude)
    write(fitsfile, FitsHeader("EXTNAME" => "laserFWHM"), laserfwhm)
    write(fitsfile, FitsHeader("EXTNAME" => "laserDist"), laserdist)
    write(fitsfile, FitsHeader("EXTNAME" => "lambdaMap"), λMap)
    close(fitsfile)
end

function test_result_strict(result, fitspath)

    (lenslettab, laserAmp, lampAmp, laserfwhm,laserdist, λMap) = result

    @testset "test_result_strict" begin
    FitsFile(fitspath) do fitsfile

        goods = cmpt_goods(lenslettab)

        @testset "validlensesmap" begin
            vlm = cmpt_vlm(lenslettab)
            goal_vlm = read(fitsfile[1])
            @test vlm == goal_vlm
        end

        @testset "bboxs" begin
            goal_bboxs = read(fitsfile["bboxs"])
            for i in goods
                @test lenslettab[i].bbox.xmin == goal_bboxs[1,i]
                @test lenslettab[i].bbox.xmax == goal_bboxs[2,i]
                @test lenslettab[i].bbox.ymin == goal_bboxs[3,i]
                @test lenslettab[i].bbox.ymax == goal_bboxs[4,i]
            end
        end

        @testset "dmodel" begin
            goal_order = fitsfile["dmodelcx"]["ORDER"].integer
            goal_λ0    = fitsfile["dmodelcx"]["L0"].float
            goal_dmodelcx = read(fitsfile["dmodelcx"])
            goal_dmodelcy = read(fitsfile["dmodelcy"])
            e = 0 
            for i in goods
                @test lenslettab[i].dmodel.order == goal_order || (e += 1 ; false)
                # λ0 is stored as FITS header so it cannot store exactly the value
                # so we need to compare using ≈
                @test lenslettab[i].dmodel.λ0 ≈ goal_λ0 || (e += 1 ; false)
                @test lenslettab[i].dmodel.cx == goal_dmodelcx[:,i] || (e += 1 ; false)
                @test lenslettab[i].dmodel.cy == goal_dmodelcy[:,i] || (e += 1 ; false)
                e > 400 && break
            end
            e > 400 && @warn "skipping tests because too much tests failed"
        end

        @testset "profile" begin
            goal_order = fitsfile["profilecy"]["ORDER"].integer
            goal_λ0    = fitsfile["profilecy"]["L0"].float
            goal_profilecy = read(fitsfile["profilecy"])
            goal_profilecl = read(fitsfile["profilecl"])
            e = 0
            for i in goods
                @test lenslettab[i].profile.order == goal_order || (e += 1 ; false)
                # λ0 is stored as FITS header so it cannot store exactly the value
                # so we need to compare using ≈
                @test lenslettab[i].profile.λ0 ≈ goal_λ0 || (e += 1 ; false)
                @test lenslettab[i].profile.cy == goal_profilecy[:,i] || (e += 1 ; false)
                @test lenslettab[i].profile.cλ == goal_profilecl[:,i] || (e += 1 ; false)
                e > 400 && break
            end
            e > 400 && @warn "skipping tests because too much tests failed"
        end

        @testset "laserAmp" begin
            goal_LaserAmp = read(fitsfile["laserAmp"])
            e = 0
            for i in goods
                @test laserAmp[:,i] == goal_LaserAmp[:,i] || (e += 1 ; false)
                e > 400 && break
            end
            e > 400 && @warn "skipping tests because too much tests failed"
        end

        @testset "laserfwhm" begin
            goal_LaserFWHM = read(fitsfile["laserFWHM"])
            e = 0
            for i in goods
                @test laserfwhm[:,i] == goal_LaserFWHM[:,i] || (e += 1 ; false)
                e > 400 && break
            end
            e > 400 && @warn "skipping tests because too much tests failed"
        end

        @testset "laserdist" begin
            goal_laserdist = read(fitsfile["laserDist"])
            @test size(laserdist) == size(goal_laserdist)
            e = 0
            for y in size(laserdist, 2), x in size(laserdist, 1)
                @test laserdist[x,y] == goal_laserdist[x,y] || (e += 1 ; false)
                e > 400 && break
            end
            e > 400 && @warn "skipping tests because too much tests failed"
        end

        @testset "λMap" begin
            @test λMap == read(fitsfile["lambdaMap"])
        end
    end
    end
end

function test_result_approx(result, fitspath)

    (lenslettab, laserAmp, lampAmp, laserfwhm,laserdist, λMap) = result

    @testset "test_result_strict" begin
    FitsFile(fitspath) do fitsfile

        goods = cmpt_goods(lenslettab)

        @testset "validlensesmap" begin
            vlm = cmpt_vlm(lenslettab)
            goal_vlm = read(fitsfile[1])
            @test vlm == goal_vlm
        end

        @testset "bboxs" begin
            goal_bboxs = read(fitsfile["bboxs"])
            for i in goods
                @test lenslettab[i].bbox.xmin == goal_bboxs[1,i]
                @test lenslettab[i].bbox.xmax == goal_bboxs[2,i]
                @test lenslettab[i].bbox.ymin == goal_bboxs[3,i]
                @test lenslettab[i].bbox.ymax == goal_bboxs[4,i]
            end
        end

        @testset "dmodel" begin
            goal_order = fitsfile["dmodelcx"]["ORDER"].integer
            goal_λ0    = fitsfile["dmodelcx"]["L0"].float
            goal_dmodelcx = read(fitsfile["dmodelcx"])
            goal_dmodelcy = read(fitsfile["dmodelcy"])
            for i in goods
                @test lenslettab[i].dmodel.order == goal_order
                @test lenslettab[i].dmodel.λ0 ≈ goal_λ0
                for o in 1:length(lenslettab[1].dmodel.cx)
                    @test isapprox(lenslettab[i].dmodel.cx[o], goal_dmodelcx[o,i]; atol=0.01)
                    @test isapprox(lenslettab[i].dmodel.cy[o], goal_dmodelcy[o,i]; atol=0.01)
                end
            end
        end

        @testset "profile" begin
            goal_order = fitsfile["profilecy"]["ORDER"].integer
            goal_λ0    = fitsfile["profilecy"]["L0"].float
            goal_profilecy = read(fitsfile["profilecy"])
            goal_profilecl = read(fitsfile["profilecl"])
            for i in goods
                @test lenslettab[i].profile.order == goal_order
                @test lenslettab[i].profile.λ0 ≈ goal_λ0
                for o in 1:length(lenslettab[1].profile.cy)
                    @test isapprox(lenslettab[i].profile.cy[o], goal_profilecy[o,i]; atol=0.01)
                    @test isapprox(lenslettab[i].profile.cλ[o], goal_profilecl[o,i]; atol=0.01)
                end
            end
        end

        @testset "laserAmp" begin
            goal_LaserAmp = read(fitsfile["laserAmp"])
            for i in goods, o in 1:length(laserAmp[:,1])
                @test isapprox(laserAmp[o,i], goal_LaserAmp[o,i]; atol=0.01)
            end
        end

        @testset "laserfwhm" begin
            goal_LaserFWHM = read(fitsfile["laserFWHM"])
            for i in goods, o in 1:length(laserfwhm[:,1])
                @test isapprox(laserfwhm[o,i], goal_LaserFWHM[o,i]; atol=0.05)
            end
        end

        @testset "laserdist" begin
            goal_laserdist = read(fitsfile["laserDist"])
            @test size(laserdist) == size(goal_laserdist)
            for y in size(laserdist, 2), x in size(laserdist, 1)
                @test laserdist[x,y] ≈ goal_laserdist[x,y]
            end
        end

        @testset "λMap" begin
            @test λMap ≈ read(fitsfile["lambdaMap"])
        end
    end
    end
end

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

function get_result_summary(result)

    (lenslettab, laserAmp, lampAmp, laserfwhm,laserdist, λMap) = result

    goods = cmpt_goods(lenslettab)
    
    ratiogoods = length(goods) / length(lenslettab)
    
    # disp model
    
    dmodelorder = lenslettab[goods[1]].dmodel.order
    dmodelλ0 = lenslettab[goods[1]].dmodel.λ0
    
    quantilesDmodelcx = fill(NaN, (9, dmodelorder+1))
    for a in 1:dmodelorder+1
        data = [ lens.dmodel.cx[a] for lens in lenslettab[goods] ]
        quantilesDmodelcx[:,a] .= quantile(filter_numbers(data), 0.1:0.1:0.9)
    end
    
    quantilesDmodelcy = fill(NaN, (9, dmodelorder+1))
    for a in 1:dmodelorder+1
        data = [ lens.dmodel.cy[a] for lens in lenslettab[goods] ]
        quantilesDmodelcy[:,a] .= quantile(filter_numbers(data), 0.1:0.1:0.9)
    end
    
    # profile model
    
    profileorder = lenslettab[goods[1]].profile.order
    
    quantilesProfilecy = fill(NaN, (9, profileorder+1))
    for a in 1:profileorder+1
        data = [ lens.profile.cy[a] for lens in lenslettab[goods] ]
        quantilesProfilecy[:,a] .= quantile(filter_numbers(data), 0.1:0.1:0.9)
    end
    
    quantilesProfilecl = fill(NaN, (9, profileorder+1))
    for a in 1:profileorder+1
        data = [ lens.profile.cλ[a] for lens in lenslettab[goods] ]
        quantilesProfilecl[:,a] .= quantile(filter_numbers(data), 0.1:0.1:0.9)
    end
    
    # laser amp
    meanLaserAmp1, stdLaserAmp1 = get_mean_std(laserAmp[1,goods])
    meanLaserAmp2, stdLaserAmp2 = get_mean_std(laserAmp[2,goods])
    meanLaserAmp3, stdLaserAmp3 = get_mean_std(laserAmp[3,goods])
    
    # laser FWHM
    meanLaserFWHM1, stdLaserFWHM1 = get_mean_std(laserfwhm[1,goods])
    meanLaserFWHM2, stdLaserFWHM2 = get_mean_std(laserfwhm[2,goods])
    meanLaserFWHM3, stdLaserFWHM3 = get_mean_std(laserfwhm[3,goods])

    # lamp Amp
    
    meanLampAmp = fill(NaN, 40)
    stdLampAmp = fill(NaN, 40)
    for i in 1:40
        m, s = get_mean_std(lampAmp[1+i, goods])
        meanLampAmp[i] = m
        stdLampAmp[i] = s
    end

    # laser dist
    meanDistColumn = fill(NaN, 5)
    stdDistColumn = fill(NaN, 5)
    for c in 1:5
        cs = [ laserdist[lens.bbox][c,:] for lens in lenslettab[goods] ]
        m, s = get_mean_std(reduce(vcat, cs))
        meanDistColumn[c] = m
        stdDistColumn[c] = s
    end
    
    # laser dist
    meanλMapLine = fill(NaN, 40)
    stdλMapLine = fill(NaN, 40)
    for l in 1:40
        ls = [ λMap[lens.bbox][:,l] for lens in lenslettab[goods] ]
        m, s = get_mean_std(reduce(vcat, ls))
        meanλMapLine[l] = m
        stdλMapLine[l] = s
    end

    (; ratiogoods, dmodelorder, dmodelλ0, quantilesDmodelcx, quantilesDmodelcy,
       quantilesProfilecy, quantilesProfilecl, meanLaserAmp1, stdLaserAmp1,
       meanLaserAmp2, stdLaserAmp2, meanLaserAmp3, stdLaserAmp3, meanLaserFWHM1, stdLaserFWHM1,
       meanLaserFWHM2, stdLaserFWHM2, meanLaserFWHM3, stdLaserFWHM3,
       meanLampAmp, stdLampAmp, meanDistColumn, stdDistColumn, meanλMapLine, stdλMapLine)
end

function get_result_summary(filepath::String)

    FitsFile(filepath) do fitsfile

        hduvlm = fitsfile[1]
        hdubboxs = fitsfile["bboxs"]
        hdudmodelcx = fitsfile["dmodelcx"]
        hdudmodelcx = fitsfile["dmodelcy"]
        hduprofilecy = fitsfile["profilecy"]
        hduprofilecl = fitsfile["profilecl"]
        hdulaseramp = fitsfile["laserAmp"]
        hdulaserfwhm = fitsfile["laserfwhm"]
        hdulaserdist = fitsfile["laserDist"]
        hduλMap = fitsfile["lambdaMap"]
        
        T = hdudmodelcx.data_eltype

        vlm = read(hduvlm)
        
        ratiogoods = sum(vlm) / length(vlm)
        
        dmodelOrder = hdudmodelcx["ORDER"].integer
        dmodelλ0 = hdudmodelcx["L0"].float
        goal_dmodelcx = read(fitsfile["dmodelcx"])
        goal_dmodelcy = read(fitsfile["dmodelcy"])
        
        # laser amp
        
        laserAmp = read(hdulaseramp)
        
        laserAmp1 = laserAmp[1,:] |> Fix1(filter, !isnan)
        (lower, upper) = quantile(laserAmp1, [0.05, 0.95])
        laserAmp1 = filter(x -> lower ≤ x ≤ upper, laserAmp1)
        meanLaserAmp1 = mean(laserAmp1)
        stdLaserAmp1  = std(laserAmp1)
        
        laserAmp2 = laserAmp[2,:] |> Fix1(filter, !isnan)
        (lower, upper) = quantile(laserAmp2, [0.05, 0.95])
        laserAmp2 = filter(x -> lower ≤ x ≤ upper, laserAmp2)
        meanLaserAmp2 = mean(laserAmp2)
        stdLaserAmp2  = std(laserAmp2)
        
        laserAmp3 = laserAmp[3,:] |> Fix1(filter, !isnan)
        (lower, upper) = quantile(laserAmp3, [0.05, 0.95])
        laserAmp3 = filter(x -> lower ≤ x ≤ upper, laserAmp3)
        meanLaserAmp3 = mean(laserAmp3)
        stdLaserAmp3  = std(laserAmp3)
        
        # laser fwhm
        
        laserFWHM = read(hdulaserfwhm)
        
        laserFWHM1 = laserFWHM[1,:] |> Fix1(filter, !isnan)
        (lower, upper) = quantile(laserFWHM1, [0.05, 0.95])
        laserFWHM1 = filter(x -> lower ≤ x ≤ upper, laserFWHM1)
        meanLaserFWHM1 = mean(laserFWHM1)
        stdlaserFWHM1  = std(laserFWHM1)
        
        laserFWHM2 = laserFWHM[2,:] |> Fix1(filter, !isnan)
        (lower, upper) = quantile(laserFWHM2, [0.05, 0.95])
        laserFWHM2 = filter(x -> lower ≤ x ≤ upper, laserFWHM2)
        meanLaserFWHM2 = mean(laserFWHM2)
        stdLaserFWHM2  = std(laserFWHM2)
        
        laserFWHM3 = laserFWHM[3,:] |> Fix1(filter, !isnan)
        (lower, upper) = quantile(laserFWHM3, [0.05, 0.95])
        laserFWHM3 = filter(x -> lower ≤ x ≤ upper, laserFWHM3)
        meanLaserFWHM3 = mean(laserFWHM3)
        stdLaserFWHM3  = std(laserFWHM3)

        (; ratiogoods,
           meanLaserAmp1, stdLaserAmp1,
           meanLaserAmp2, stdLaserAmp2,
           meanLaserAmp3, stdLaserAmp3,
           meanLaserFWHM1, stdLaserFWHM1,
           meanLaserFWHM2, stdLaserFWHM2,
           meanLaserFWHM3, stdLaserFWHM3)
   end
end

function display_summary(sum)
    foreach(keys(sum)) do k
        print("========== ")
        println(k)
        show(stdout, MIME"text/plain"(),sum[k])
        println()
    end
end

