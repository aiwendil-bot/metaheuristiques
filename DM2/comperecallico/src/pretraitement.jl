function vect_variables(matrix)
    result = [] 
    sizehint!(result, size(matrix, 2))
    for i in 1:size(matrix, 2)
        j = 1
        vect = Vector{Int}(undef, sum(matrix[:,i]))
        for k in 1:size(matrix, 1)
            if matrix[k,i] == 1
                vect[j] = k
                j = j + 1
            end
        end
        push!(result, vect)
    end
    return result
end

# renvoie un vecteur de vecteurs, de taille nb_contraintes
# chaque vecteur contient les indices des variable pour lesquelles la contrainte est liée

function vect_contraintes(matrix)
    result = []
    sizehint!(result, size(matrix, 1))
    for i in 1:size(matrix, 1)
        j = 1
        vect = Vector{Int}(undef, sum(matrix[i,:]))
        for k in 1:size(matrix, 2)
            if matrix[i,k] == 1
                vect[j] = k
                j = j + 1
            end
        end
        push!(result, vect)
    end
    return result
end
