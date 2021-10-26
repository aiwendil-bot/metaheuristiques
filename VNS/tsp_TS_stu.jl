# -------------------------------------------------------------------------------------------
# X. Gandibleux - Metaheuristiques : TSP TS (version etudiant) - Septembre 2019
# Base de depart « etudiant » pour exercice "TS applique au TSP"

@static if VERSION < v"1.0-"
    error("NOT COMPLIANT WITH JULIA < v1.0.0")
end

using Random
using Printf
using PyPlot # a commenter si le tracage des resultats n'est pas appele

# -------------------------------------------------------------------------------------------
# construit un distancier symetrique aux valeurs de distances uniformément aléatoires dans 1..vSup
function tsp_elaboreDistancierAleatoire(
           n::Int64,    # in : nombre de villes
           vSup::Int64  # in : plus grande distance entre deux villes
         )
    distancier = rand(1:vSup,n,n)
    for i=1:n
        for j=1:n
                distancier[i,j] = distancier[j,i]
        end
        distancier[i,i] = 0
    end
    return distancier # out : matrice carree symetrique de distances entieres
end

# -------------------------------------------------------------------------------------------
# calcule la longueur d’un circuit hamiltonien correspondant a une solution TSP admissible
function tsp_calculeLongueurCircuit(
           d::Matrix{Int64},  # in : distancier symetrique de ville a ville
           x::Array{Int64},   # in : solution admissible représentée par codage de permutation des villes
         )
    n::Int64 = length(x)
    longueur::Int64 = d[ x[n] , x[1] ]
    for i = 1:n-1
        longueur += d[ x[i] , x[i+1] ]
    end
    return longueur # out : longueur du circuit correspondant a la solution admissible
end

# -------------------------------------------------------------------------------------------
# construit une solution aléatoire TSP admissible représentée dans un codage de permutation
function tsp_elaboreSolutionAleatoire(
           n::Int64  # in : nombre de villes
         )
    return shuffle(collect(1:n))  # out : une permutation aleatoire des villes
end

# ------------------------------------------------------------
# Plot a graph with the results obtained by the SA (z for all solutions, best, greedy)

function traceResultat(titre, zAllLst, zBestLst, zBestAllLst)
    figure("Examen d'un run",figsize=(6,6))
    title(titre)
    xlabel("Itérations")
    ylabel("valeurs de z(x)")
    grid()
    plot(zAllLst,label="all solutions")
    plot(zBestAllLst, linewidth=2.0, color="green", label="best solutions")
    legend(loc=1, fontsize ="small")
end

# -------------------------------------------------------------------------------------------
function tsp_solveTS(
           d::Matrix{Int64},      # in : distancier symetrique de ville a ville
           x::Array{Int64},       # in : solution TSP admissible représentée pas un codage de permutation des villes
           z::Int64,              # in : valeur de la solution TSP admissible
           iterTabou::Int64,      # in : nombre d'iterations a realiser par la recherche
           dureeTabou::Int64      # in : valeur de la longueur tabou
         )

    # traces d'activite pour le trace de graphiques
    zAllLst     = [] # all values of f(x) collected for each iteration
    zBestLst    = [] # each best values of f(x) collected
    zBestAllLst = [] # best values of f(x) collected for each iteration
    push!(zAllLst, z)
    push!(zBestLst, z)
    push!(zBestAllLst, z) ; zMin = z

    @printf(" parametres TS : iterTabou= %d, dureeTabou= %d\n\n", iterTabou,  dureeTabou)

    n      = length( x )
    xCur   = copy( x ) # solution courante
    zCur   = z
    xBest  = copy( x ) # meilleure solution
    zBest  = z

    # Initiallisation memoire tabou avec  "tous les mouvements sont autorises"
    listeTabou  = zeros(Int64,n,n)

    # Boucle principale de la TS
    for nr_iteration=1:iterTabou

        trouveAmeliorant = false
        deltaBest = typemax(Int64)
        iDbest = 0  ;  iFbest = 0
        jDbest = 0  ;  jFbest = 0

        # recherche meilleur mouvement autorise en appliquant un 2opt sur un voisinage complet
        for i=0:n-1
            for j=0:n-4

                iD = i%n + 1       ;  iF = (i+1)%n + 1
                jD = (i+j+2)%n + 1 ;  jF = (i+j+3)%n + 1

                d1 = d[ xCur[iD] , xCur[iF] ] + d[ xCur[jD] , xCur[jF] ]
                d2 = d[ xCur[iD] , xCur[jD] ] + d[ xCur[iF] , xCur[jF] ]
                delta = d2 - d1

                voisinTrouve = ( delta  <  deltaBest )
                mouvementAutorise = (listeTabou[ xCur[iD] , xCur[jD] ] <= nr_iteration) && (listeTabou[ xCur[iF] , xCur[jF] ] <= nr_iteration)
                critereAspiration = ( zCur + delta < zBest )

                if ( voisinTrouve  && ( mouvementAutorise || critereAspiration ))
                    trouveAmeliorant = true
                    deltaBest = delta
                    iDbest = iD  ;  iFbest = iF
                    jDbest = jD  ;  jFbest = jF
                end
            end
        end

        if trouveAmeliorant

            # maj memoire tabou
            listeTabou[ xCur[iDbest] , xCur[iFbest] ] = dureeTabou + nr_iteration
            listeTabou[ xCur[jDbest] , xCur[jFbest] ] = dureeTabou + nr_iteration
            listeTabou[ xCur[iFbest] , xCur[iDbest] ] = dureeTabou + nr_iteration
            listeTabou[ xCur[jFbest] , xCur[jDbest] ] = dureeTabou + nr_iteration

            # maj incrementale de la valeur de la solution
            zCur = zCur + deltaBest
            push!(zAllLst,zCur)

            # maj solution : 2opt et inverse les segments du chemin entre les 2 segments traites
            if iDbest < jDbest
                xCur[iDbest+1:jDbest]= reverse(xCur[iDbest+1:jDbest])
            else
                xtwice=vcat(xCur, xCur)
                xtwice[iDbest+1:jDbest+n]= reverse(xtwice[iDbest+1:jDbest+n])
                xCur[1:jDbest] = xtwice[n+1:n+jDbest]
                xCur[iDbest:n] = xtwice[iDbest:n]
            end

            # met a jour la meilleure solution si il y a eu amelioration
            if (zCur < zBest)
                zBest = zCur
                xBest = copy(xCur)
                @printf("  numIter = %3d    z(xBest) = %4d\n", nr_iteration, zBest)
                push!(zBestLst,zBest)
                zMin = zBest
            end

        else
            @assert false "Tous les mouvements sont interdits : liste tabou trop longue"
        end

        push!(zBestAllLst,zMin)
    end # indice sur nr_iteration

    return xBest, zBest, zAllLst, zBestLst, zBestAllLst

end

# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------

function tsp_TS(
           n::Int64            = 20,  # in : nombre de villes
           vSup::Int64         = 25, # in : plus grande distance entre deux villes
           iterTabou::Int64    = 60, # in : nombre d'iterations a realiser par la recherche
           dureeTabou::Int64   = 7   # in : valeur de la longueur tabou
         )

    d = tsp_elaboreDistancierAleatoire(n,vSup) # distancier aleatoire

    xInit = tsp_elaboreSolutionAleatoire(n)
    zInit = tsp_calculeLongueurCircuit(d,xInit)

    println("Init    : z(xInit) = ", zInit)
    println("  ")

    xBest, zBest, zAllLst, zBestLst, zBestAllLst = tsp_solveTS(d, xInit, zInit, iterTabou, dureeTabou)
    println("\nTS      : z(xTS)   = ", zBest)
    println("  ")

    traceResultat("TS : TSP | n=" * string(20) , zAllLst, zBestLst, zBestAllLst)

end

