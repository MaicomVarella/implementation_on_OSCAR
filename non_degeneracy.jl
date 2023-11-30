#@doc Markdown.doc"""
#    faces_to_keep(a)
#
#    For a $2\times 3$ matrix $A$ with polynomial entries return "true" if a vector of all_faces(a) contains faces $\sigma_j$ of $\Gamma(a_{i,j})$ such that $\sum_{j=1,2,3} \sigma_j$ belongs to the Newton diagram $\Gamma(a_{1,1}*a_{1,2}*a_{1,3})$ and return "false" otherwise.
#
#    """

function all_faces(a)
        return [[P Q R] for P in Newton_diagram(a[1,1]) for Q in Newton_diagram(a[1,2]) for R in Newton_diagram(a[1,3])]
end

function faces_to_keep(a)
        all = all_faces(a)
        return vcat([all[i][1,1]+all[i][1,2]+all[i][1,3] in Newton_diagram(a[1,1]*a[1,2]*a[1,3]) for i in 1:size(all,1)])
end

#@doc Markdown.doc"""
 #   bounded_faces(a)
#
#    For a $2\times 3$ matrix $A$ with polynomial entries return a subcollection of all_faces(a) containing only the faces $\sigma_j$ of $\Gamma(a_{i,j})$ such that $\sum_{j=1,2,3} \sigma_j$ belongs to the Newton diagram $\Gamma(a_{1,1}*a_{1,2}*a_{1,3})$.

#    """
function bounded_faces(a)
        all = all_faces(a)
        keep = faces_to_keep(a)
        return all(a)[keep .==1]
end

#@doc Markdown.doc"""
 #   matrix_function(a,P,Q,S)

  #  For a $2\times 3$ matrix $A$ with polynomial entries and polyhedra $P$, $Q$ and $R$ return a matrix with entries $(a_{1,1})_P$, $(a_{2,1})_P$, $(a_{1,2})_Q$, $(a_{1,2})_Q$, $(a_{1,3})_R$ and $(a_{2,3})_R$ .

   # """

function matrix_function(a,P,Q,S)
        return R[face_function(a[1,1],P) face_function(a[1,2],Q) face_function(a[1,3],S); face_function(a[2,1],P) face_function(a[2,2],Q) face_function(a[2,3],S)]
end

#@doc Markdown.doc"""
 #   is_non_degenerate(a)
#
 #   For a $2\times 3$ matrix $A$ with polynomial entries return "true" if for all vectors [P Q R] of bounded_faces(a) the determinantal variety defined by matrix_function(a,P,Q,R) is non-singular outside the coordinate hyperplanes $(V(\prod_{i=1}^n x_i),0)$.
  #  """

function is_non_degenerate(a)
        ND1 = Newton_diagram(a[1,1])
        ND2 = Newton_diagram(a[1,2])
        ND3 = Newton_diagram(a[1,3])
        ND_sum = Newton_diagram(a[1,1]*a[1,2]*a[1,3])
        All = [[P Q R] for P in ND1 for Q in ND2 for R in ND3]
        Keep = vcat([All[i][1,1]+All[i][1,2]+All[i][1,3] in ND_sum for i in 1:size(All,1)])
        Bound = All[Keep .==1]
        outside_torus = ideal(R,R([1],[[1 for j in 1:ngens(R)]]))
        Ideal(a) = ideal(R, minors(a,2))
        singular_locus_ideal(a)= modulus(underlying_quotient(OO(singular_locus(Spec(R,ideal(R,minors(a,2)),complement_of_ideal(ideal(R,gens(R))))))))
        saturation_singular =[saturation(singular_locus_ideal(matrix_function(a,Bound[i][1,1],Bound[i][1,2],Bound[i][1,3])),outside_torus)==ideal(R,[1]) for i in 1: size(Bound,1)]
        saturation_dimension = [dim(saturation(ideal(R,minors(matrix_function(a,Bound[i][1,1],Bound[i][1,2],Bound[i][1,3]),2)),outside_torus))==dim(Ideal(a)) for i in 1: size(Bound,1)]
        saturation_empty = [saturation(Ideal(matrix_function(a,Bound[i][1,1],Bound[i][1,2],Bound[i][1,3])),outside_torus)==ideal(R,[1]) for i in 1: size(Bound,1)]
        if [saturation_dimension[i]+saturation_empty[i]+saturation_singular[i] for i in 1:size(Bound,1)] == [2 for i in 1:size(Bound,1)]
                return true
        else
                return false
        end
end
