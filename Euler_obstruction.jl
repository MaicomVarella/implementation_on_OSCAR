using Combinatorics
using InvertedIndices

function partition(n,k)
  p=with_replacement_combinations(1:n,k)
  return map(A -> [sum(A .== i) for i in 1:n],p)
end

function part(n,k)
           P = partition(n,k)
           element(i,j) = P[i][j] < 2
           vec(i) = sum(element(i,j) for j in 1:(size(P[i],1)))
           keep(i) = vec(i) == size(P[i],1)
           K = [keep(i) for i in 1:size(P,1)]
           return P[K]
           end

#@Markdown.doc """
#           matrix(p,J)
#           For a polyhedron $P$ and a vector $J$ return the matrix with the vertices of the intersection between $P$ and the coordinate hyperplane $\mathbb{C}^J$.
#
#           """

function matrix(p,J)
        V(p) = transpose(reduce(hcat, vertices(p)))
        h(p,J)=V(p)[:, Not(J)]
        l(x)=sum(x)<1
        return V(p)[map(l, eachrow(h(p,J))),J]
end

#@Markdown.doc """
#           NP0(p,J)
#           For a polyhedron $P$ and a vector $J$ return the polyhedron generated by the lines of matrix(p,J) and the origin.
#           """

function NP0(p,J)
        return convex_hull([0; matrix(p,J)])
end

#@Markdown.doc """
#           NP0(p,J)
#           For a polyhedron $P$ and a vector $J$ return the polyhedron generated by the lines of matrix(p,J).
#           """

function NP(p,J)
        return convex_hull(matrix(p,J))
end

#@Markdown.doc """
#           Vol(P,J)
#
#           For a polyhedron $P$ and a vector $J$ return the difference volume(NP0(x,J))-volume(NP(x,J)).
#           """

function Vol(x,J)
        return volume(NP0(x,J))-volume(NP(x,J))
end

function mult(K,x,i)
        if K[i] != 0
        return K[i]*x[i]
        else
        return Polyhedron(-identity_matrix(QQ, ngens(R)),vcat([0 for i in 1:ngens(R)]))
        end
end

#@Markdown.doc """
#           V(J,x...)
#
#           For a set of polyhedra $x_1,\dots,x_{|J|}$ and a vector $J$ return the difference mixed volume MV(x^1^J,\dots,x_|J|^J)
#           """

function MV(J,x...)
        return sum((-1)^(size(J,1)+k)*sum(Vol(sum(mult(K,x,i) for i in 1:size(J,1)),J) for K in part(size(J,1),k)) for k in 1:size(J,1))
end

function MV_3(a,x,b,y,c,z,J)
       V = vcat([x for i in 1:a],[y for i in 1:b],[z for i in 1:c])
       return MV(J,V...)
 end

function MV_4(a,x,b,y,c,z,d,w,J)
        V = vcat([x for i in 1:a],[y for i in 1:b],[z for i in 1:c],[w for i in 1:d])
        return MV(J,V...)
end

function pt(n,k)
        P = partition(n,k)
        element(i,j) = P[i][j] == 0
        vec(i) = sum(element(i,j) for j in 1:(size(P[i],1)))
        keep(i) = vec(i) == 0
        K = [keep(i) for i in 1:size(P,1)]
        return P[K]
end

#@Markdown.doc """
#           Euler_characteristic_Milnor_fiber(f,a)
#
#           For a function $f$ and a $2\times 3$ matrix return the Euler characteristic of the Milnor fiber of $f$ restricted to the determinantal singularity $X_A^2$.
#           """
function Euler_characteristic_Milnor_fiber(f,a)
        x = Newton_polytope(f)
        y = Newton_polytope(a[1,1])
        z = Newton_polytope(a[1,2])
        w = Newton_polytope(a[1,3])
        c(m,p)=collect(powerset([1:m;],p,p))
        return sum(sum(sum((-1)^(p-1)*(MV_3(I[1],x,I[2],y,I[3],z,J) + MV_3(I[1],x,I[2],y,I[3],w,J) + MV_3(I[1],x,I[2],z,I[3],w,J)) for I in pt(3,p)) for J in c(ngens(R),p)) for p in 3:ngens(R)) + sum(sum(sum((-1)^(p+1)*(p-I[1]+1)*MV_4(I[1],x,I[2],y,I[3],z,I[4],w,J) for I in pt(4,p)) for J in c(ngens(R),p)) for p in 4:ngens(R))
       end

       function Euler_characteristic_Milnor_fiber_matrix(f,a)
               x = Newton_polytope(f)
               y = Newton_polytope(a[1,1]+a[1,2]+a[1,3]+a[2,1]+a[2,2]+a[2,3])
               z = y
               w = y
               c(m,p)=collect(powerset([1:m;],p,p))
               return sum(sum(sum((-1)^(p-1)*(MV_3(I[1],x,I[2],y,I[3],z,J) + MV_3(I[1],x,I[2],y,I[3],w,J) + MV_3(I[1],x,I[2],z,I[3],w,J)) for I in pt(3,p)) for J in c(ngens(R),p)) for p in 3:ngens(R)) + sum(sum(sum((-1)^(p+1)*(p-I[1]+1)*MV_4(I[1],x,I[2],y,I[3],z,I[4],w,J) for I in pt(4,p)) for J in c(ngens(R),p)) for p in 4:ngens(R))
              end