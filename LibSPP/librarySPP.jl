# --------------------------------------------------------------------------- #
# collect the un-hidden filenames available in a given folder

function getfname(pathtofolder)

    # recupere tous les fichiers se trouvant dans le repertoire cible
    allfiles = readdir(pathtofolder)

    # vecteur booleen qui marque les noms de fichiers valides
    flag = trues(size(allfiles))

    k=1
    for f in allfiles
        # traite chaque fichier du repertoire
        if f[1] != '.'
            # pas un fichier cache => conserver
            println("fname = ", f)
        else
            # fichier cache => supprimer
            flag[k] = false
        end
        k = k+1
    end

    # extrait les noms valides et retourne le vecteur correspondant
    finstances = allfiles[flag]
    return finstances
end

# --------------------------------------------------------------------------- #
# Loading an instance of SPP (format: OR-library)

function loadSPP(fname)

    f=open(fname)
    # lecture du nbre de contraintes (m) et de variables (n)
    m, n = parse.(Int, split(readline(f)) )
    # lecture des n coefficients de la fonction economique et cree le vecteur d'entiers c
    C = parse.(Int, split(readline(f)) )
    # lecture des m contraintes et reconstruction de la matrice binaire A
    A=zeros(Int, m, n)
    for i=1:m
        # lecture du nombre d'elements non nuls sur la contrainte i (non utilise)
        readline(f)
        # lecture des indices des elements non nuls sur la contrainte i
        for valeur in split(readline(f))
          j = parse(Int, valeur)
          A[i,j]=1
        end
    end
    close(f)
    return C, A
end

# --------------------------------------------------------------------------- #
# Deamon checking the feasibility of a solution given its vector x

function isAdmissible(C, A, x)

    vecSat = zeros(Int, size(A,1))
    vecUnit = ones(Int,size(A,1))
    z::Int64 = 0
    verbose = true
    var1 = findall(isequal(1), x[:])

    for j in var1
        vecSat = vecSat .+ A[:,j]
        z = z + C[j]
    end

    if findfirst(isequal(false), (vecSat .<= vecUnit)) != nothing
        println( "Feasible : no")
        @assert false "no-feasible solution detected"
    end
    println( "Feasible : yes | Σ(x_i) = ", length(var1), " ; z(x) = ", z)
    return true
end
