from fractions import Fraction
from collections import defaultdict

##def get_bounded_partitions(P,sname,pct=False):
##    #P=RNAPolytope.construct_from_file(sname)
##    bounded = defaultdict(lambda:0)
##    dirs = 'abc'
##    for sl in P.d1_slices():
##        if sl.dimension()!=3:
##            continue
##        rays = sl.rays()
##        unb = ''
##        for i in range(3):
##            unb= unb + (dirs[i] if any([v[i] for v in rays]) else '')
##
##        bounded[unb] += 1
##        #if unb=='a':
##        #    print timestamp("Unbounded only in a: {0}; {1}; {2}".format(sname, sl.original_vertex, slice_center(sl)))
##        #elif unb=='c':
##        #    print timestamp("Unbounded only in c: {0}; {1}; {2}".format(sname, sl.original_vertex, slice_center(sl)))
##        #elif unb=='ab':
##        #    print timestamp("Unbounded only in ab: {0}; {1}; {2}".format(sname, sl.original_vertex, slice_center(sl)))
##        #elif unb == 'bc':
##        #    print timestamp(sname)
##
##    n = len(P.d1_slices())
##
##    total_bounded = bounded['']
##    total_unbounded = n-total_bounded
##
##    if pct:
##        for k in bounded:
##            bounded[k] *= 100.0/n
##
##    bounded['bounded'] = (100.0*total_bounded)/n
##    bounded['unbounded'] = (100.0*total_unbounded)/n
##
##    return bounded
##
##def is_bounded(sl):
##    return len(sl.rays())==0
##def is_unbounded(sl):
##    return len(sl.rays())>0
##def is_wedge(sl):
##    return len(sl.rays())==2
##def is_strip(sl):
##    return len(sl.rays())==1
##
##def get_slice_width_height_depth(sl):
##    from sage.numerical.mip import MIPSolverException
##    I = []
##    for i in range(3):
##        lp,V=sl.to_linear_program(return_variable=True)
##        var = V[i]
##        lp.set_objective(var)
##        try:
##            o1 = lp.solve()
##        except MIPSolverException as e:
##            if str(e).find("unbounded")>-1:
##                o1 = infinity
##            else:
##                raise e
##        lp.set_objective(-var)
##        try:
##            o2 = -lp.solve()
##        except MIPSolverException as e:
##            if str(e).find("unbounded")>-1:
##                o2 = -infinity
##            else:
##                print(e)
##        I.append(o1-o2)
##    return I
##def get_slice_containing_pt(P,pt=(34/10,0,4/10)):
##    slices = []
##    for sl in P.d1_slices():
##        if pt in sl:
##            slices.append(sl)
##    return slices

## USED IN InitializePolytopes
def get_RNAPolytope( rnapoly_path ):
    return RNAPolytope.construct_from_file( rnapoly_path )

## USED IN FinalMerge (data output section) as well as AccuracyOf > CreateAccuracyData
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

##def get_var_intervals(sl,pt=None):
##    if pt is None:
##        pt = slice_center(sl)
##    elif pt=='classic':
##        pt=classic
##    from sage.numerical.mip import MIPSolverException
##    I = []
##    for i in range(3):
##        lp,V=sl.to_linear_program(return_variable=True)
##        var = V[i]
##        lp.set_objective(var)
##        for j in range(3):
##            if j!=i:
##                lp.add_constraint(V[j]==pt[j])
##        try:
##            o1 = lp.solve()
##        except MIPSolverException as e:
##            if str(e).find("unbounded")>-1:
##                o1 = infinity
##            else:
##                print(e)
##        lp.set_objective(-var)
##        try:
##            objective = lp.solve()
##            # objective is the largest value that -1 * var can have. 
##            # So the smallest value that var can have, that satisfies the linear program, is -1 * objective.
##            o2 = -1 * objective
##        except MIPSolverException as e:
##            if str(e).find("unbounded")>-1:
##                o2 = -infinity
##            else:
##                print(e)
##        w = min([abs(o1-pt[i]),abs(o2-pt[i])])
##        I.append(w)
##    return I
