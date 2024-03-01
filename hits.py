import numpy

# the graph
# A->D
# B->C   B->E
# C->A
# D->C
# E->D  E->B   E->F    E->C
# F->C  F->H
# G->A  G->C
# H->A
# where node C is a dangling node
L = numpy.matrix( [ [0., 0., 0., 1., 0., 0., 0., 0.],
                    [0., 0., 1., 0., 1., 0., 0., 0.],
                    [1., 0., 0., 0., 0., 0., 0., 0.],
                    [0., 0., 1., 0., 0., 0., 0., 0.],
                    [0., 1., 1., 1., 0., 1., 0., 0.],
                    [0., 0., 1., 0., 0., 0., 0., 1.],
                    [1., 0., 1., 0., 0., 0., 0., 0.],
                    [1., 0., 0., 0., 0., 0., 0., 0.]
                    ] )



Lt = L.transpose()
print("*******L *************")
print(L)
print("*******Lt*************")
print(Lt)
print("**********************")

no_nodes = L.shape[0]
print("number of nodes in graph L:", no_nodes)

# create two vertical vectors no_nodes x 1, where each element is equal to 1.
a = numpy.matrix( [ [1.]  for i in range(no_nodes) ] )
h = numpy.matrix( [ [1.]  for i in range(no_nodes) ] )
a = a / numpy.sum(a)
h = h / numpy.sum(h)


for i in range(50):
    a_new = Lt * h    # h * L
    h_new = L * a

    # normalize
    a_new = a_new / numpy.sum(a_new)
    h_new = h_new / numpy.sum(h_new)

    # Compute the Euclidean distance between the vectors a and a_new, and h and h_new
    # If the error/distance is very small, terminate the process
    error1 =  numpy.sqrt(numpy.sum(numpy.square((a-a_new))))
    error2 =  numpy.sqrt(numpy.sum(numpy.square((h-h_new))))

    print("ERROR")
    print(error1)
    print(error2)

    a = a_new
    h = h_new
    if error1 < 10**-10 and error2 < 10**-10:
        print("early exit at iteration:", i)
        break

print("final authority scores")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(a)

print("final hub scores")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(h)

