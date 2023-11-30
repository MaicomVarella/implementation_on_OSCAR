

@doc Markdown.doc"""
    Newton_polytope(f::MPolyElemLoc)

    For a multivariate power series $f = \sum_{a \in NN^n}$ c_a x^a$, the Newton polytope is defined as

    ```math
    \Delta (f) = Conv(\bigcup_{a \in supp(f)} a+{\mathbb R}_{\geq 0}^n),
    ```
    where
    ```math
    supp(f) = \{a \in NN^n : c_a \neq 0 \}.
    ```

    """

function Newton_polytope(f::MPolyElem)
         R = parent(f)
         exponents = reduce(hcat, Oscar.AbstractAlgebra.exponent_vectors(f))'
         positive_orthant = Polyhedron(-identity_matrix(QQ, ngens(R)),vcat([0 for i in 1:ngens(R)]))
         p(i) = positive_orthant + exponents[i,:]
         return foldl(convex_hull,[p(i) for i in 1:size(exponents,1)])
         end

function Newton_polytope(f::MPolyElemLoc)
          g = numerator(f)
          return Newton_polytope(g)
          end

 #@doc Markdown.doc"""
   #  Newton_diagram(f::MPolyElemLoc)
#
#     The Newton diagram of a power series $f$ can be understood as the union of the bounded faces of its Newton polytope.
#
#     """

function Newton_diagram(f::MPolyElem)
         R = parent(f)
         faces_to_keep(x::Polyhedron,d) = vcat([is_bounded(faces(x,d)[i]) for i in 1:length(faces(x,d))])
         bounded_faces(x::Polyhedron,d) = faces(x,d)[faces_to_keep(x,d).==1]
         return foldl(union,[bounded_faces(Newton_polytope(f),d) for d in 0:ngens(R)])
         end

 function Newton_diagram(f::MPolyElemLoc)
          g = numerator(f)
          return Newton_diagram(g)
          end
#@doc Markdown.doc"""
#        is_convenient(f::MPolyElemLoc)
#
#        A power series f is called convenient, if its Newton polytope meets all the coordinate axes.
#
#        """

function touch_positive_axis(x::Polyhedron,i::Int64)
         positive_orthant = Polyhedron(-identity_matrix(QQ, ambient_dim(x)),vcat([0 for i in 1:ambient_dim(x)]))
         positive_axis(i) = faces(positive_orthant,1)[i]
         if dim(intersect(x,positive_axis(i)))<0
         return false
         else
         return true
         end
         end

@attributes fmpq_mpoly

@attr Bool function is_convenient(f::MPolyElem)
         if vcat([touch_positive_axis(Newton_polytope(f),i) for i in 1:ambient_dim(Newton_polytope(f))]) == vcat([1 for i in 1:ambient_dim(Newton_polytope(f))])
         return true
         else
         return false
         end
         end

@attr Bool function is_convenient(f::MPolyElemLoc)
          g = numerator(f)
          return is_convenient(g)
          end

#@doc Markdown.doc"""
#        face_function(f::MPolyElemLoc,P::Polyhedron)
#
#        Let $f = \sum_{a \in NN^n}$ c_a x^a$ be a multivariate, convenient power series and let $\Gamma(f)$ be its Newton-diagram. For each face $\sigma \in \Gamma(f)$ we define a quasihomogeneous polynomial
#```math
#   f_{\sigma}=\sum_{a \in \sigma} c_a x^a
#```
#        """

function face_function(f::MPolyElem,P::Polyhedron)
        g=zero(f)
        for (c,e) in zip(coefficients(f),exponents(f))
        if contains(P,e)
        g = g+R([c],[e])
        end
        end
        return g
        end

function face_function(f::MPolyElemLoc,P::Polyhedron)
        g = numerator(f)
        return face_function(g,P)
        end

#@doc Markdown.doc"""
#        is_Newton_non_degenerate(f::MPolyElemLoc)
#        $f$ is called Newton-non-degenerate, if for all faces $f_{\sigma}$, $\sigma \in \Gamma(f)$, the germ $(V(f_\sigma),0)$ is non-singular outside the coordinate hyperplanes $(V(\prod_{i=1}^n x_i),0)$.
#        """

@attr Bool function is_Newton_non_degenerate(f::MPolyElem)
         R = parent(f)
         singular_locus_ideal(f)= modulus(underlying_quotient(OO(singular_locus(Spec(R,ideal(R,f),complement_of_ideal(ideal(R,gens(R)))))[1])))
         outside_torus = ideal(R,R([1],[[1 for j in 1:ngens(R)]]))
         if [saturation(singular_locus_ideal(face_function(f,P)),outside_torus)==ideal(R,[1]) for P in Newton_diagram(f)] == [1 for P in Newton_diagram(f)]
         return true
         else
         return false
         end
         end

 @attr Bool function is_Newton_non_degenerate(f::MPolyElemLoc)
         g = numerator(f)
         return is_newton_non_degenerate(g)
          end
