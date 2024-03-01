
import numpy

# the graph
#     A -> B -> C
# where node C is a dangling node
A = numpy.matrix([
    [0., 1., 0.],
    [0., 0., 1.],
    [0., 0., 0.],
])

no_nodes = A.shape[0]
print("number of nodes in graph A:", no_nodes)

# Fix dangling nodes
for row_id in range(no_nodes):
    row_sum = numpy.sum(A[row_id, :])
    if row_sum == 0.0:   # node row_id is dandling
        A[row_id, :] = numpy.matrix([[1 for i in range(no_nodes)]])


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("without dangling nodes")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(A)

# Normalize, making the matrix stochastic
for row_id in range(no_nodes):
    row_sum = numpy.sum(A[row_id, :])  # out degree of node row_id
    A[row_id, :] /= row_sum

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("normalized")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(A)

# Initialize the teleportation matrix T
T = numpy.matrix([[1.0/no_nodes for i in range(no_nodes)] for j in range(no_nodes)])

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("teleportation")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(T)


# Markov matrix
# Final Matrix after removing dangling nodes and adding teleportation
d = 0.85
M = d * A + (1-d) * T

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(" Matrix of the Markov process ")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(M)

# Initial probability distribution
p = numpy.matrix([[1./no_nodes for i in range(no_nodes)]])
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Initial distribution")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(p)

for i in range(1, 101):
    # p*M is equivalent to:    p_new = transpose(M)*p
    p_new = p*M

    if i % 10 == 0:
        print("iter:", i)
        print(p_new)

    # Compute the Euclidean distance between the vectors p and p_new
    # If the error/distance is very small, terminate the Markov process
    error = numpy.sqrt(numpy.sum(numpy.square((p-p_new))))

    p = p_new
    if error < 10**-10:
        print("early exit at iteration:", i)
        break


print("final distribution")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(p)






