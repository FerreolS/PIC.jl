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
    write(fitsfile, FitsHeader("EXTNAME" => "laserFWHM"), laserfwhm)
    write(fitsfile, FitsHeader("EXTNAME" => "laserDist"), laserdist)
    write(fitsfile, FitsHeader("EXTNAME" => "lambdaMap"), λMap)
    close(fitsfile)
end

function test_result_strict(result, fitspath)

    (lenslettab, laserAmp, lampAmp, laserfwhm,laserdist, λMap) = result

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
                @test lenslettab[i].dmodel.cx == goal_dmodelcx[:,i]
                @test lenslettab[i].dmodel.cy == goal_dmodelcy[:,i]
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
                @test lenslettab[i].profile.cy == goal_profilecy[:,i]
                @test lenslettab[i].profile.cλ == goal_profilecl[:,i]
            end
        end

        @testset "laserAmp" begin
            goal_LaserAmp = read(fitsfile["laserAmp"])
            for i in goods
                @test laserAmp[:,i] == goal_LaserAmp[:,i]
            end
        end

        @testset "laserfwhm" begin
            goal_LaserFWHM = read(fitsfile["laserFWHM"])
            for i in goods
                @test laserfwhm[:,i] == goal_LaserFWHM[:,i]
            end
        end

        @testset "laserdist" begin
            @test laserdist == read(fitsfile["laserDist"])
        end

        @testset "λMap" begin
            @test λMap == read(fitsfile["lambdaMap"])
        end
    end
end