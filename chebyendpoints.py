p = 1

n = 5

product = np.zeros((p, n))

for i in range(n):

    for k in range(p-1):

        print((i**2-k**2)/(2*k+1))

        product[k, i] = (i**2-k**2)/(2*k+1)
    
