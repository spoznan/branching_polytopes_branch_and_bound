from fractions import Fraction
from collections import defaultdict

def get_RNAPolytope( rnapoly_path ):
    return RNAPolytope.construct_from_file( rnapoly_path )

def slice_center(sl,delta=Fraction(1,2)):
    vertex_sum = vector(sl.base_ring(), [0]*sl.ambient_dim())
    #print(vertex_sum)
    for v in sl.vertex_generator():
        vertex_sum += v.vector()
    vertex_sum=vertex_sum / (sl.n_vertices())
    #print("2nd test",vertex_sum)
    # print "Average vertex sum: ", vertex_sum
    Z = IntegerRing()
    delta = Z(1)/Z(2)
    for v in sl.ray_generator():
        # print "Ray:", v.vector()
        #print("delta:", delta)
        for i in range(0,len(vertex_sum)):
            vertex_sum[i] +=delta*v.vector()[i]        
        #vertex_sum += delta*v.vector()
        #print("3rd test", vertex_sum)
        # print "Vertex sum after a ray:", vertex_sum
    vertex_sum.set_immutable()
    return vertex_sum
