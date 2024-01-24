module ZeroHalfChvatGom

export ZHOrigSubSys

mutable struct ZeroHalfSystem
    ABar::Array{Int64,2}
    bBar::Array{Int64,1}
    xStar::Array{Float64,1}
    s::Array{Float64,1}
    rowIdx::Vector{BitSet}
    m::Int64
    n::Int64
end
ZeroHalfSystem() = ZeroHalfSystem(zeros(1, 1), zeros(1), zeros(1), zeros(1), [BitSet()], 1, 1)

struct ZeroHalfSepParams
    heurStrength::Int64
    printLevel::Int64
end

struct violatedIneq
    rows::Array{Int64,1}
    violation::Float64
end

function createZeroHalfSystem(A::Array{Int64,2}, b::Array{Int64,1}, xStarIn::Array{Float64,1})
    zhs = ZeroHalfSystem()
  # first construct the matrix ABar and bBar of 0/1 elements
    minElem = minimum([minimum(A) ,minimum(b)])
  # ensure that all numbers in A and b are non-negtive in order for the % operator to give 0/1 values (we don't want -1 entries)
    addElem = 0
    if minElem < 0
        addElem = convert(Int, ceil(-minElem / 2)) * 2
    end
    println("addElem: $addElem")
    zhs.ABar = (A .+ addElem) .% 2
    zhs.bBar = (b .+ addElem) .% 2
    zhs.s = b - A * xStarIn

    (zhs.m, zhs.n) = size(zhs.ABar)
    zhs.xStar = copy(xStarIn)
    zhs.rowIdx = Vector{BitSet}(undef, zhs.m)
    for i = 1:zhs.m
        zhs.rowIdx[i] = BitSet(i)
    end
    return zhs
end

function deleteRow!(zhs::ZeroHalfSystem, row::Int64)
    zhs.ABar[row,:] = zhs.ABar[zhs.m,:]
    zhs.rowIdx[row] = zhs.rowIdx[zhs.m]
    zhs.s[row] = zhs.s[zhs.m]
    zhs.bBar[row] = zhs.bBar[zhs.m]
  ################################################################################################
  # The following four operations are not needed. It's just to make the data structure look nice.
  # later we can remove them to speed up the code a little.
    zhs.ABar[zhs.m,1:zhs.n] = -1 * ones(zhs.n)
    zhs.rowIdx[zhs.m] = BitSet()
    zhs.s[zhs.m] = -1
    zhs.bBar[zhs.m] = -1
  ################################################################################################

    zhs.m -= 1
end

function deleteColumn!(zhs::ZeroHalfSystem, col::Int64)
    zhs.ABar[:,col] = zhs.ABar[:,zhs.n]
    zhs.xStar[col] = zhs.xStar[zhs.n]
  ################################################################################################
  # The following two operations are not needed. It's just to make the data structure look nice.
  # later we can remove them to speed up the code a little.
    zhs.ABar[1:zhs.m,zhs.n] = -1 * ones(zhs.m)
    zhs.xStar[zhs.n] = -1
  ################################################################################################

    zhs.n -= 1
end

# add row i to row ii in zhs
function addrow!(zhs::ZeroHalfSystem, i::Int64, ii::Int64)
    zhs.ABar[ii,1:zhs.n] = (zhs.ABar[ii,1:zhs.n] + zhs.ABar[i,1:zhs.n]) .% 2
    zhs.bBar[ii] = (zhs.bBar[ii] + zhs.bBar[i]) % 2
    zhs.rowIdx[ii] = symdiff(zhs.rowIdx[i], zhs.rowIdx[ii])
    zhs.s[ii] = zhs.s[ii] + zhs.s[i]
end

function reduceMat(zhs::ZeroHalfSystem, EPSVAL::Float64, zhsParams::ZeroHalfSepParams)
    # if zhsParams.printLevel >= 1
        println("Entering reduceMat (m,n)=($(zhs.m),$(zhs.n))")
    # end
    # println("ABar: $(zhs.ABar)")

    j = 1
  # rule 1: delete columns with xStar[j]=0 (we only need to apply this rule once)
    while j <= zhs.n
        if zhs.xStar[j] <= EPSVAL
      # delete column j
            deleteColumn!(zhs, j)
            if zhsParams.printLevel >= 5
                println("Rule 1: Deleting column due to x[$j]==0. (m,n)=($(zhs.m),$(zhs.n))")
            end
        else
            j += 1
        end
    end

  # rule 2: delete rows that are are all zero in ABar and bBar.
    i = 1
    while i <= zhs.m
        allZero = true
        for j = 1:zhs.n
            if zhs.ABar[i,j] == 1
                allZero = false
                break
            end
        end
        if zhs.bBar[i] == 1
            allZero = false
        end
        if allZero
      # delete row i
      # println("row i: ")
      # println(zhs.ABar[i,1:zhs.n])
            deleteRow!(zhs, i)
            if zhsParams.printLevel >= 5
                println("Rule 2: Deleting row $i due to all zero row. (m,n)=($(zhs.m),$(zhs.n))")
            end
        else
            i += 1
        end
    end

  # rule 3: Zero columns in \bar{A} can be removed.
  # rule 4: Identical columns in \bar{A} can be replaced by a single representative
  # with associated variable value as sum of the merged variables.
  # we do rule 5 at the same time:
  # rule 5. Any unit vector column \bar{a}_{i}=e_{j}, 1\leq j\leq m, in \bar{A} can be removed
  # provided that x_{i}^{*} is added to the slack s_{j} of row j .

  # We are checking rules 3, 4 and 5 in the loop below.

  # columnDict is a dictionary with a big integer as key and
  # a column id as value. Each column is converted to a big number
  # (just using binary encoding), if the number already is in the dict then
  # we have found a duplicate column and it is eliminated.
    columnDict = Dict{BigInt,Int64}()
    j = 1
    while j <= zhs.n
        val = BigInt(0)
        onesCount = 0
        lastOne = -1
        for i = 1:zhs.m
            if zhs.ABar[i,j] == 1
                val += BigInt(2)^(i - 1)
                onesCount += 1
                lastOne = i
            end
        end
        if onesCount == 0
      # rule 3:
      # column is all zeros. We can delete it.
            deleteColumn!(zhs, j)
            if zhsParams.printLevel >= 5
                println("Rule 3: Deleting column $j due to all zero column. (m,n)=($(zhs.m),$(zhs.n))")
            end
        elseif onesCount == 1
      # rule 4:
      # column is a unit vector. We can get rid of the column if we add xStar[j] to s[i] where i is the row with the single one.
            zhs.s[lastOne] += zhs.xStar[j]
            deleteColumn!(zhs, j)
            if zhsParams.printLevel >= 5
                println("Rule 4: Deleting column $j and updating slack of row $lastOne as it is a unit vector. s[$lastOne]=$(zhs.s[lastOne]). (m,n)=($(zhs.m),$(zhs.n))")
            end
        elseif haskey(columnDict, val)
      # rule 5:
      # We have a duplicate column. Delete column j and transfer x value to the
      # other "instance" of the column
            othercol = columnDict[val]
            zhs.xStar[othercol] += zhs.xStar[j]
      # delete column j
            deleteColumn!(zhs, j)
            if zhsParams.printLevel >= 5
                println("Rule 5: Deleting column $j due to duplicate column. (m,n)=($(zhs.m),$(zhs.n))")
            end
        else
      # We have not seen this column before. Add it to the dictionary and go on
            columnDict[val] = j
            j += 1
        end
    end

  # Rule 6: Any row 1\leq j\leq m with slack s_{j}\geq1 can be removed.
    i = 1
    while i <= zhs.m
        if zhs.s[i] >= 1
            deleteRow!(zhs, i)
            if zhsParams.printLevel >= 5
                println("Rule 6: Deleting row $i due to slack >= 1. (m,n)=($(zhs.m),$(zhs.n))")
            end
        else
            i += 1
        end
    end

  # Rule 7. Rows identical in \left(\bar{A},\bar{b}\right) can be eliminated except for one with smallest slack value.

  # rowDict is a dictionary with a big integer as key and
  # a column id as value. Each row (including the RHS) is converted to a big number
  # (just using binary encoding), if the number already is in the dict then
  # we have found a duplicate row and the one with largest slack is eleminated
    rowDict = Dict{BigInt,Int64}()

    i = 1
    while i <= zhs.m
        val = BigInt(0)
        for j = 1:zhs.n
            if zhs.ABar[i,j] == 1
                val += BigInt(2)^(j - 1)
            end
        end
        if zhs.bBar[i] == 1
            val += BigInt(2)^zhs.n
        end
        if zhsParams.printLevel >= 5
            println("Testing rule 7. i=$i, val = $val")
        end
        if haskey(rowDict, val)
      # We have a duplicate row.
            otherRow = rowDict[val]
            if (zhs.s[otherRow] < zhs.s[i])
        # otherrow has smallest slack. Then we can just delete row i
        # println("Row $i: ", zhs.ABar[i,1:zhs.n], " RHS: ", zhs.bBar[i], "s: ", zhs.s[i])
        # println("Row $otherRow: ", zhs.ABar[otherRow,1:zhs.n], " RHS: ", zhs.bBar[otherRow], "s: ", zhs.s[otherRow])
                deleteRow!(zhs, i)
                if zhsParams.printLevel >= 5
                    println("Rule 7: Deleting column $i due to duplicate row. (m,n)=($(zhs.m),$(zhs.n))")
                end
            else
        # this row has smallest slack. Copy to other row and delete this one
        # println("Row $i: ", zhs.ABar[i,1:zhs.n], " RHS: ", zhs.bBar[i], "s: ", zhs.s[i])
        # println("Row $otherRow: ", zhs.ABar[otherRow,1:zhs.n], " RHS: ", zhs.bBar[otherRow], "s: ", zhs.s[otherRow])
                zhs.rowIdx[otherRow] = zhs.rowIdx[i]
                zhs.s[otherRow] = zhs.s[i]
                deleteRow!(zhs, i)
                if zhsParams.printLevel >= 5
                    println("Rule 7: Deleting column $i (was $otherRow)due to duplicate row. (m,n)=($(zhs.m),$(zhs.n))")
                end
            end
        else
      # We have not seen this row before. Add it to the dictionary and go on
            rowDict[val] = i
            i += 1
        end
    end
    # if zhsParams.printLevel >= 1
        println("Leaving reduceMat (m,n)=($(zhs.m),$(zhs.n))")
    # end
end

function gaussLikeSimplify(zhs, EPSVAL, zhsParams::ZeroHalfSepParams)
    ineqs = Vector{violatedIneq}()
    # if zhsParams.printLevel >= 1
        println("Entering gaussLikeSimplify (m,n)=($(zhs.m),$(zhs.n))")
    # end
    i = 1
    while i <= zhs.m
        if zhsParams.printLevel >= 5
            println("Gauss: i=$i")
        end
        if (zhs.s[i] <= EPSVAL)
      # row i has zero slack. Three cases:
      # 1) all element in row i of \bar{A} has value zero and \bar{b}_i=1    =>  We have identified a maximally violated zero-half cut
      # 2) all element in row i of \bar{A} has value zero and \bar{b}_i=0    =>  The row can be eleminated according to rule 2 of prop. 3
      # 3) there is at least one "1" in row i of \bar{A} => we use row i to simplify the system based on proposition 5.
      #    Use column j with a[i,j]=1 that has highest value of x^*_j  (suggestion of Koster et al.)
            highX = -99999
            bestJ = -1
            countOnes = 0
            for j = 1:zhs.n
                if zhs.ABar[i,j] == 1
                    countOnes += 1
                    if zhs.xStar[j] > highX
                        highX = zhs.xStar[j]
                        bestJ = j
                    end
                end
            end
            if countOnes == 0
        # AbAr[i,:] is all zeros
                if zhs.bBar[i] == 1
          # maximally violated inequality detected
                    push!(ineqs, violatedIneq(collect(zhs.rowIdx[i]), 0.5))
                    if zhsParams.printLevel >= 2
                        println("Gauss: Maximally violated inequality detected. Rows: $(zhs.rowIdx[i])")
                    end
          # I guess we can delete this row now
                    deleteRow!(zhs, i)
                else
          # all element in row i of \bar{A} has value zero and \bar{b}_i=0    =>  The row can be eleminated according to rule 2 of prop. 3
                    deleteRow!(zhs, i)
                    if zhsParams.printLevel >= 5
                        println("Gauss: deleting row $i according to rule 2 of prop. 3. (m,n)=($(zhs.m),$(zhs.n))")
                    end
                end
            else
        # we are in case 3.
        # We use proposition 5 to get rid of column bestJ
        # Add row i to all rows ii where ABar[ii,bestJ] = 1
                for ii = 1:zhs.m
                    if i != ii && zhs.ABar[ii,bestJ] == 1
                        addrow!(zhs, i, ii)
                    end
                end
                zhs.s[i] += zhs.xStar[bestJ]
                deleteColumn!(zhs, bestJ)
                if zhsParams.printLevel >= 5
                    println("Gauss: deleting column $bestJ according to prop. 5. Slack of row $i updated. (m,n)=($(zhs.m),$(zhs.n))")
                end
                i += 1
            end
        else
            i += 1
        end
    end
    # if zhsParams.printLevel >= 1
        println("Leaving gaussLikeSimplify (m,n)=($(zhs.m),$(zhs.n))")
    # end
    return ineqs
end

function testRowCombination(zhs::ZeroHalfSystem, rows::Vector{Int64}, EPSCUT::Float64)
    # println("in <testRowCombination>")
    sumS = sum(zhs.s[rows])
    # can we disregard this set of rows based on slack alone?
    if (sumS <= 1 - EPSCUT * 2)
        # no, the rows could identify a violated inequality based on the slack
        # make a row vector v with m 0/1 elements. Ones are placed at the position of the used rows
        v = zeros(1, zhs.m)
        # println("rows: $rows")
        for row in rows
            v[row] = 1
        end
        ABar = view(zhs.ABar, 1:zhs.m, 1:zhs.n)
        b = view(zhs.bBar, 1:zhs.m)
        xStar = view(zhs.xStar, 1:zhs.n)
        vb = v * b
        # right hand side need to be odd
        if (vb[1] % 2) == 1
            # Here we compute vs + (v*\bar{A} mod 2) x^* in a compact way. I'm not sure if it is super efficient
            # (is there a lot of overhead with the views?)
            val = sumS .+ ((v * ABar) .% 2 ) * xStar
            # println("val: $val")
            if val[1] < 1 - EPSCUT * 2
                res = zhs.rowIdx[rows[1]]
                for i = 2:length(rows)
                    res = symdiff(res, zhs.rowIdx[rows[i]])
                end
                return (true, violatedIneq(collect(res), 0.5 - val[1] / 2))
            end
        end
    end
    return (false, violatedIneq([], 0))
end

function combineRowsHeur(zhs::ZeroHalfSystem, nRows::Int64, EPSCUT::Float64, zhsParams::ZeroHalfSepParams)
    ineqs = Vector{violatedIneq}()
    counters = zeros(Int64, nRows)
    countersMaxVal = zeros(Int64, nRows)
    m = zhs.m
    if nRows > m
        println("Warning (combineRowsHeur): reducing nRows since zero half system only contains $m rows")
        nRows = m
    end
    for i = 1:nRows
        counters[i] = i
        # Example: if nRows = 3 and m = 10 then countersMaxVal should be [8,9,10]
        countersMaxVal[i] = m - (nRows - i)
    end
    if zhsParams.printLevel >= 10
        println("countersMaxVal: ", countersMaxVal)
    end
    done = false
    while !done
        # println("Counters: $counters")
        # try to combine the rows given by <counters>
        (foundVioaltedIneq, cut) = testRowCombination(zhs, counters, EPSCUT)
        if foundVioaltedIneq
            if zhsParams.printLevel >= 2
                println("Row combination heuristic found violated inequality.")
                println("Rows: $(cut.rows), violation: $(cut.violation)")
            end
            push!(ineqs, cut)
        end
        # here we update counters. This is a mess. Can it be done simpler? If not, then put it into some general utility function maybe?
        cIdx = nRows
        countersNeedUpdate = true
        while countersNeedUpdate
            counters[cIdx] += 1
            # println("Counters: $counters")
            if counters[cIdx] <= countersMaxVal[cIdx]
                countersNeedUpdate = false
            else
                if cIdx <= 1
                    if zhsParams.printLevel >= 10
                        println("cIdx = $cIdx, setting done = true")
                    end
                    countersNeedUpdate = false
                    done = true
                else
                    cIdx -= 1
                end
            end
        end
        # "reset" counters after the last "correctly changed counter"
        for i = (cIdx + 1):nRows
            counters[i] = counters[i - 1] + 1
        end
    end
    # it would be nice to clean up ineqs. Could there be duplicates in the vector?
    return ineqs
end

function ZeroHalfSep(A::Array{Int64,2}, b::Array{Int64,1}, xStarIn::Array{Float64,1}, EPSVAL::Float64, zhsParams::ZeroHalfSepParams)
    println("Entering: ZeroHalfSep. Strength = $(zhsParams.heurStrength)")
    # throw("exit")
    # tic()
    zhs = createZeroHalfSystem(A, b, xStarIn)
    reduceMat(zhs, EPSVAL, zhsParams)
    allCuts = gaussLikeSimplify(zhs, EPSVAL, zhsParams)
    # Can this call ever reduce the system further?
    reduceMat(zhs, EPSVAL, zhsParams)

    println(length(allCuts), " cut(s) found after gauss like elimination")

    if zhs.m >= 1
        if (zhsParams.printLevel >= 10)
            for i = 1:zhs.m
                println("row $i of ZHS is consisting of original rows: ", zhs.rowIdx[i])
            end
        end
        for nRows = 1:minimum([zhsParams.heurStrength,zhs.m])
          #= if length(allCuts) > 0
              break
          end =#
            println("Attempting to combine $nRows rows.")
            newCuts = combineRowsHeur(zhs, nRows, 0.01, zhsParams)
            for cut in newCuts
                push!(allCuts, cut)
            end
        end
        println("$(length(allCuts)) cuts found. Some may be duplicates.")
      # toc()
    else
        println("ZeroHalfSep: Reduction rules completely eliminated system. Skipping heuristic")
    end

    # test one of the cuts
    #= 
    if length(allCuts) > 0
        cutId = rand(1:length(allCuts))
        cut = allCuts[cutId]
        println("Selecting cut: Rows: $(cut.rows), violation: $(cut.violation)")
        (m,n) = size(A)
        u = zeros(1,m)
        #println("rows: $rows")
        for row in cut.rows
            u[row] = 1
        end
        LHSCoeffs = floor.(0.5*u*A)
        RHS = floor.(0.5*u*b)
        println("LHSCoeffs: $LHSCoeffs")
        println("RHS: $RHS")
        println("LHS*xstar = ", LHSCoeffs*xStarIn)
        println("violation: ", LHSCoeffs*xStarIn-RHS )
    end =#
    if (zhsParams.printLevel >= 4)
        for i = 1:length(allCuts)
            println("----- cut $i: rows -----")
            cut = allCuts[i]
            println(cut.rows)
        end
    end

    return allCuts
end
# Since the original MIP may contain real-valued variables A may just contain a subset of the columns of the original constraint matrix
# <vecVarToGlobalVar> maps from a column # in A (a variable in the reduced system) to a column # in the global systems (a variable in the
# global system)
function ZeroHalfSepTranslateVars(A::Array{Int64,2}, b::Array{Int64,1}, xStarIn::Array{Float64,1}, EPSVAL::Float64, zhsParams::ZeroHalfSepParams, vecVarToGlobalVar::Vector{Int64})
    allCuts = ZeroHalfSep(A, b, xStarIn, EPSVAL, zhsParams)
    cutsExpressedInOrigVars = []
    # tic()
    println("Converting cuts to normal form....")
    # @@ This is slower than expected. I guess the matrix multiplications are slow and could be sped up by
    # a more direct approach.
    cutCount = 0
    for cut in allCuts
        (m, n) = size(A)
        u = zeros(1, m)
        # println("rows: $rows")
        for row in cut.rows
            u[row] = 1
        end
        LHSCoeffs = floor.(0.5 * u * A)
        RHS = floor.(0.5 * u * b)
        vars = Vector{Int64}()
        coeffs = Vector{Float64}()
        for i = 1:length(LHSCoeffs)
            if LHSCoeffs[i] != 0
                push!(vars, vecVarToGlobalVar[i])
                push!(coeffs, LHSCoeffs[i])
            end
        end
        push!(cutsExpressedInOrigVars, [vars,coeffs, RHS[1] ])
        cutCount += 1
        if cutCount >= 100
            println("Has 100 cuts now... stopping.")
            break
        end
    end
    println("done converting cuts to normal form. Time spent on conversion: ")
    # toc()
    return cutsExpressedInOrigVars
end

mutable struct ZHOrigSubSys
    ASub::Array{Int64,2}
    bSub::Array{Int64,1}
    varsSorted::Array{Int64,1}
end
ZHOrigSubSys() = ZHOrigSubSys(Array{Int64,2}(undef, 0, 0), Array{Int64,1}(), Array{Int64,1}())


# This is a "near-copy" of a function in AutoDec.jl with the same name. In AutoDec.jl the function takes a
# MIP object as input, here we do not use this object (goal: enable us to use the zero half cut separation outside AutoDec)
function findZeroHalfOrigSubSystem(A, consLB, consUB, varType, addBinVarUBConstraints::Bool, printLevel::Int64, varNames=[])
    (m, n) = size(A)
    println("findZeroHalfOrigSubSystem: orig size: ($m,$n)")
    usedVars = BitSet()
    validRows = Vector{Int64}()
    # bSub : RHS vector
    bSub = Vector{Int64}()
    vecMult = Vector{Int64}()
    for i = 1:m
        # check if row j can be used in the Zero-Half sub-sustem
        # - All coefficients include b_j should be integer (at some point we could lift this assumption, but we would have to
        #   update the separation routine)
        # - All variables that are used in the row should be discrete variables
        rowIsValid = true
        RHS = 0.0
        if consUB[i] == Inf
            RHS = -consLB[i]
            # this is a >= constraint. Convert to <= by multiplying with -1
            push!(vecMult, -1)
        else
            RHS = consUB[i]
            # this is a <= or == constraint. We keep it as a <= constaint
            push!(vecMult, 1)
        end
        if (RHS - round(RHS)) != 0
            rowIsValid = false
        else
            allZeroRow = true
            for j = 1:n
                val = A[i,j]
                if (val != 0)
                    if (val - round(val)) != 0
                        rowIsValid = false
                        break
                    else
                        if (varType[j] == :Int) || (varType[j] == :Bin)
                            push!(usedVars, j)
                            allZeroRow = false
                        else
                            rowIsValid = false
                            break
                        end
                    end
                end
            end
        end
        if rowIsValid && !allZeroRow
            push!(validRows, i)
            push!(bSub, RHS)
        end
    end
    # <validRows> contains the rows we can use for generation of robust zero-half cuts
    # <usedVars> contains the variables that were used in those rows.
    varsSorted = sort(collect(usedVars))
    ASub = convert(Array{Int64,2}, A[validRows, varsSorted])
    for i = 1:length(validRows)
        if vecMult[i] == -1
            ASub[i,:] = -ASub[i,:]
        end
    end
    # generate upper bound rows for binary variables (@@@ we should check if the constraint system already contains these)
    # @@@ We may also want to put in upper bound constraints for general integer variables.
    if addBinVarUBConstraints
        countBin = 0
        for i in usedVars
            if varType[i] == :Bin
                countBin += 1
            end
        end
        if countBin > 0
            (m, n) = size(ASub)
            AUB = zeros(countBin, n)
            bUB = ones(countBin)
            row = 1
            for i = 1:length(varsSorted)
                if varType[varsSorted[i]] == :Bin
                    AUB[row,i] = 1
                    row += 1
                end
            end
            # println(AUB)
            ASub = vcat(ASub, AUB)
            bSub = vcat(bSub, bUB)
        end
    end
    if printLevel >= 10 && length(varNames) > 0
        (m, n) = size(ASub)
        for i = 1:m
            print("row $i: ")
            firstOnLine = true
            for j = 1:n
                if ASub[i,j] != 0
                    if !firstOnLine
                        print(" + ")
                    else
                        firstOnLine = false
                    end
                    print(ASub[i,j], "*$(varNames[varsSorted[j]])")
                end
            end
            println(" <= $(bSub[i])")
        end
    end
    (m, n) = size(ASub)
    println("findZeroHalfOrigSubSystem: subsystem size: ($m,$n)")
    return ZHOrigSubSys(ASub, bSub, varsSorted)
end


end
